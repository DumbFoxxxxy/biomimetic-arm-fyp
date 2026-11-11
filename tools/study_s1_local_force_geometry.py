"""Screen GH leverage and endpoint travel using the selected S1-Y2 interfaces."""

import csv
import json

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read


def load_reference(reference):
    path = ROOT / reference["file"]
    if provenance(path) != reference:
        raise ValueError(f"Source fingerprint changed: {path}")
    return read(path)


def span(origin, distal, centre):
    difference = origin - distal
    length = np.linalg.norm(difference)
    direction = difference / length
    return length, direction, np.cross(distal - centre, direction)


def isolated_torque(levers, axis, sign, max_tension):
    perpendicular = null_space(axis.reshape(1, 3))
    result = linprog(-sign * (levers @ axis) / 1000.,
                     A_eq=perpendicular.T @ levers.T / 1000., b_eq=np.zeros(2),
                     bounds=(0., max_tension), method="highs")
    if not result.success:
        raise RuntimeError(result.message)
    return dict(target_torque_Nm=max(0., -result.fun), tensions_N=result.x.tolist(),
                torque_source_Nm=(result.x @ levers / 1000.).tolist())


def main():
    selection_path = ROOT / "design/S1_CURRENT_CANDIDATE.json"
    selection = read(selection_path)
    shoulder_ref = dict(file=selection["assembly_report"], sha256=selection["assembly_sha256"])
    bare_ref = dict(file=selection["bare_motion_report"], sha256=selection["bare_motion_report_sha256"])
    shoulder, bare = load_reference(shoulder_ref), load_reference(bare_ref)
    inventory_path = DEST / "route_inventory.json"
    inventory = read(inventory_path)
    if inventory["shoulder_assembly"] != shoulder_ref:
        raise ValueError("Route inventory and S1 selection refer to different assemblies")
    routes = [row for row in inventory["routes"] if row["region"] == "shoulder"]
    route_ref = routes[0]["local_route_reference"]
    route_document = load_reference(route_ref)
    strand_ref = next(row["current_free_strand"]["source"] for row in routes if "current_free_strand" in row)
    strands = load_reference(strand_ref)
    if strands["assembly_sha256"] != shoulder_ref["sha256"]:
        raise ValueError("Free tendon and S1 selection refer to different assemblies")
    route_ports = {row["head"]: row["ports"] for row in route_document["routes"]}
    centre = np.array(bare["joint_centre_source_mm"])
    axis_names = ["flexion", "abduction", "longitudinal"]
    axes = np.array([bare["axes_source"][name] for name in axis_names])
    if not np.allclose(np.linalg.norm(axes, axis=1), 1.):
        raise ValueError("Expected unit GH axes")

    supports, y_tendons = {}, {}
    for route in routes:
        head = route["id"]
        if "current_free_strand" in route:
            row = next(row for row in strands["heads"] if row["head"] == head)["neutral"]["route"]
            origin, distal = (np.array(row[f"exit_{end}_source_mm"]) for end in ("A", "B"))
            chord, direction, chord_lever = span(origin, distal, centre)
            tangent = np.array(row["exit_tangent_B_toward_free_segment"])
            y_tendons[head] = dict(
                free_terminal_exit_source_mm=origin.tolist(), humeral_exit_source_mm=distal.tolist(),
                nominal_unmeasured_length_mm=row["nominal_length_assumption_mm"], neutral_chord_mm=chord,
                bone_exit_chord_direction=direction.tolist(), bone_exit_tangent_direction=tangent.tolist(),
                tangent_projection_mm=(axes @ np.cross(distal - centre, tangent)).tolist(),
                chord_projection_mm=(axes @ chord_lever).tolist(),
                scope="Prescribed neutral eye orientations. Conditional per-unit-tension projections only."
                      " Actual tension, swivel angle and free TPU terminal equilibrium remain unknown.")
            continue
        distal_key = f"stop:{head}_stop_blank" if route["module"] == "deltoid" else f"cap:{head}/distal"
        if shoulder["intended_body_owners"][distal_key] != "humerus_ideal_attachment":
            raise ValueError(f"{head} has no humerus-owned terminal in this assembly")
        ports = route_ports[head]
        points = [np.array(ports[role]["outer_face_source_mm"]) for role in ("origin", "distal")]
        supports[head] = (route, points, distal_key)

    neutral, samples, summary = {}, [], {}
    for head, (route, (origin, distal), distal_key) in supports.items():
        length, direction, lever = span(origin, distal, centre)
        neutral[head] = dict(origin_body=route["anatomical_origin_body"], distal_body="humerus",
                             proximal_port_source_mm=origin.tolist(), distal_port_source_mm=distal.tolist(),
                             distal_geometry_key=distal_key, chord_mm=length, distal_force_per_unit_tension=direction.tolist(),
                             proximal_reaction_per_unit_tension=(-direction).tolist(),
                             moment_arms_mm=(axes @ lever).tolist())

    for path in bare["results"]:
        name, sign = path["axis"], int(path["sign"])
        axis_index = axis_names.index(name)
        axis = axes[axis_index]
        limit = abs(path["last_contiguous_sample_clear_deg"])
        angles = np.unique(np.r_[np.arange(0., limit, 1.), limit]) * sign
        direction_summary = dict(last_prior_bare_sample_deg=float(angles[-1]), channels={})
        all_levers = []
        for angle in angles:
            rotation = Rotation.from_rotvec(np.deg2rad(angle) * axis)
            pose_levers = []
            for head, (_, (origin, distal), _) in supports.items():
                offset = rotation.apply(distal - centre)
                length, direction, lever = span(origin, centre + offset, centre)
                projected = axes @ lever
                epsilon = 1e-6
                distances = [np.linalg.norm(origin - centre - Rotation.from_rotvec(delta * axis).apply(offset))
                             for delta in (-epsilon, epsilon)]
                numerical_lever = (distances[0] - distances[1]) / (2. * epsilon)
                if not np.isclose(numerical_lever, projected[axis_index], atol=1e-6, rtol=0.):
                    raise RuntimeError(f"Moment arm sign or units disagree for {head}")
                pose_levers.append(lever)
                samples.append(dict(axis=name, direction=sign, angle_deg=angle, channel=head, chord_mm=length,
                                    delta_from_neutral_mm=length-neutral[head]["chord_mm"],
                                    flexion_arm_mm=projected[0], abduction_arm_mm=projected[1],
                                    longitudinal_arm_mm=projected[2]))
            all_levers.append(pose_levers)
        for head in supports:
            rows = [row for row in samples if row["channel"] == head and row["axis"] == name
                    and row["direction"] == sign]
            arms = np.array([sign * row[f"{name}_arm_mm"] for row in rows])
            travel = [row["delta_from_neutral_mm"] for row in rows]
            sign_changes = [(float(angles[i]), float(angles[i + 1])) for i in range(len(angles) - 1)
                            if (arms[i] > 0.) != (arms[i + 1] > 0.)]
            direction_summary["channels"][head] = dict(
                signed_lever_range_in_motion_direction_mm=[float(arms.min()), float(arms.max())],
                chord_change_range_mm=[min(travel), max(travel)],
                chord_change_at_last_sample_mm=travel[-1], sign_change_brackets_deg=sign_changes)
        probe_indices = {0, len(angles) // 2, len(angles) - 1}
        if sign == 1 and name in ("flexion", "abduction"):
            probe_indices.update(np.flatnonzero(angles == 96.))
        direction_summary["isolated_active_torque_probes"] = {
            str(float(angles[index])): isolated_torque(np.array(all_levers[index]), axis, sign, 25.)
            for index in sorted(probe_indices)}
        summary[f"{name}_{sign:+d}"] = direction_summary

    output = DEST.parent / "s1_local_muscle_load"
    output.mkdir(exist_ok=True)
    with (output / "gh_endpoint_geometry.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
        writer.writeheader()
        writer.writerows(samples)
    report = dict(generator=provenance(ROOT / "tools/study_s1_local_force_geometry.py"), selection=provenance(selection_path),
                  shoulder_assembly=shoulder_ref, bare_motion=bare_ref, port_geometry=route_ref,
                  free_strand_geometry=strand_ref, axes=axis_names, axes_source=axes.tolist(),
                  gh_centre_source_mm=centre.tolist(), step_max_deg=1.,
                  signed_longitudinal_anatomical_mapping_verified=bare["anatomical_axial_sign_mapping_verified"],
                  neutral=neutral, Y_neutral_conditional=y_tendons, paths=summary, sample_rows=len(samples),
                  isolated_torque_probe_bound_N_per_channel=25., tension_bound_is_assumed=True,
                  scope="S1-Y2 fixed shoulder girdle, single GH axes, prescribed diagnostic poses."
                        " Chords are local endpoint force surrogates; clearance and TPU strain are not predicted."
                        " Four-channel torque fit excludes gravity, passive load and the two free-terminal Y muscles."
                        " Prior bare-motion limits are probe bounds, not local actuator or complete assembly limits.")
    (output / "gh_endpoint_geometry.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(samples)} local span samples and {len(y_tendons)} conditional Y projections")


if __name__ == "__main__":
    main()
