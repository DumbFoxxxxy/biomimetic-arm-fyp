"""Sample the six-millimetre FJ1500 candidate along GH flexion and return."""

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def summary(rows, loads):
    summaries = {}
    for load in loads:
        poses = [pose for pose in rows if pose["input_N"] == load]
        fail_index = next((index for index, pose in enumerate(poses) if pose["hits"]), len(poses))
        summaries[str(int(load))] = dict(
            last_clear=poses[fail_index - 1] if fail_index else None,
            first_collision=poses[fail_index] if fail_index < len(poses) else None,
            failed_poses=sum(bool(pose["hits"]) for pose in poses))
    return summaries


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    rest_file = folder / "lateral_rest_tpu_geometry.json"
    rest = read(rest_file)
    screen = load_reference(rest["bias_screen"])
    sleeve = load_reference(screen["compact_sleeve"])
    previous = load_reference(screen["reference_equilibrium"])
    geometry = load_reference(previous["local_geometry"])
    assembly = load_reference(previous["shoulder_assembly"])
    model = next(row for row in load_reference(previous["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    if provenance(ROOT / rest["candidate"]["file"])["sha256"] != rest["candidate"]["sha256"]:
        raise ValueError("TPU candidate changed after readback")
    sleeve_key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    sources = dict(assembly["sources"])
    if sources[sleeve_key] != sleeve["original"]:
        raise ValueError("Sleeve candidate uses a different S1 assembly")
    sources[sleeve_key] = sleeve["selected"]["source"]
    solids = {key: load_component(source) for key, source in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    pin_a, pin_b = (np.array(model[f"neutral_{end}_pin_source_mm"]) for end in ("cap", "bone"))
    axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
    offset = np.array(rest["rest_offset_source_mm"])
    adjusted_model, ports = rest_geometry(model, route["ports"], offset)

    def evaluate(angle, levels):
        gh = Rotation.from_rotvec(flexion * np.deg2rad(angle))
        rows = []
        for state in solve_branch(adjusted_model, ports, centre, gh, screen["stiffness"], levels, retain_leads=True):
            pin, hits = collision_state(solids, assembly["intended_body_owners"], state, offset, gh,
                                        centre, pin_a, pin_b, axis_a, axis_b)
            rows.append(dict(angle_deg=angle, input_N=state["commanded_muscle_tension_N"],
                             cap_pin_source_mm=pin.tolist(), cap_rotation_vector_deg=state["cap_rotation_vector_deg"],
                             swivel_angles_deg=state["swivel_angles_deg"], tension_N=state["short_tendon_tension_N"],
                             equilibrium_residual_N=state["stationarity_residual_N"], hits=hits))
        return rows

    loads = [0., 5., 10.]
    outward = []
    for angle in range(97):
        outward.extend(evaluate(angle, loads))
        if angle % 20 == 0:
            print("Outward through", angle, "degrees", flush=True)
    fine = {}
    for key, limit in summary(outward, loads).items():
        clear, hit = limit["last_clear"], limit["first_collision"]
        if not (clear and hit):
            continue
        for index in range(1, 10):
            angle = round(clear["angle_deg"] + index * .1, 1)
            fine.setdefault(angle, None)
    for angle in fine:
        fine[angle] = evaluate(angle, loads)
    out_samples = sorted(outward + [row for rows in fine.values() for row in rows], key=lambda row: row["angle_deg"])
    limits = summary(out_samples, loads)
    ceiling = int(np.ceil(max((row["last_clear"]["angle_deg"] if row["last_clear"] else 0.)
                              for row in limits.values())))
    angles = sorted({row["angle_deg"] for row in out_samples if row["angle_deg"] <= ceiling}, reverse=True)
    returning = []
    for angle in angles:
        # Each return pose is re-solved from rest with the reverse load sequence.
        returning.extend(evaluate(angle, loads[::-1]))
    outward_lookup = {(row["angle_deg"], row["input_N"]): row for row in out_samples}
    direction_agreement = {}
    for load in loads:
        compared = [row for row in returning if row["input_N"] == load]
        deltas = []
        disagreements = []
        for row in compared:
            peer = outward_lookup[row["angle_deg"], load]
            deltas.append(np.linalg.norm(np.array(row["cap_pin_source_mm"]) - peer["cap_pin_source_mm"]))
            if bool(row["hits"]) != bool(peer["hits"]):
                disagreements.append(row["angle_deg"])
        contiguous = [row for row in out_samples if row["input_N"] == load and
                      row["angle_deg"] <= (limits[str(int(load))]["last_clear"] or {"angle_deg": -1.})["angle_deg"]]
        max_step = max((np.linalg.norm(np.array(a["cap_pin_source_mm"]) - b["cap_pin_source_mm"])
                        for a, b in zip(contiguous[:-1], contiguous[1:])), default=0.)
        direction_agreement[str(int(load))] = dict(
            same_angle_comparisons=len(compared), mismatched_collision_angles_deg=disagreements,
            max_same_angle_pin_delta_mm=max(deltas, default=0.),
            max_adjacent_outward_clear_pin_step_mm=max_step,
            return_collision_angles_inside_outward_clear_range_deg=[
                row["angle_deg"] for row in compared if row["hits"]
                and (limits[str(int(load))]["last_clear"] or {"angle_deg": -1.})["angle_deg"] >= row["angle_deg"]])
    report = dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_flexion_path.py"),
        rest_geometry=provenance(rest_file), sleeve=screen["compact_sleeve"], rest_offset_mm=rest["rest_offset_mm"],
        loads_N=loads, fixed_GH_search_limit_deg=96., coarse_step_deg=1., boundary_step_deg=.1,
        rigid_interference_threshold_mm3=.001, limits=limits, return_agreement=direction_agreement,
        outward=out_samples, return_samples=returning,
        scope="Prescribed GH angles. Local unloaded/loaded equilibrium sampled in both directions with independent"
              " starts and compared at matching angles. Excludes dynamic motion, compliant contact, steel solids"
              " and bend stiffness. Collision thresholds and input force do not establish physical reach.")
    save(folder / "FJ1500_flexion_path.json", report)
    for load, limit in limits.items():
        print(load, "N last", (limit["last_clear"] or {}).get("angle_deg"), "first",
              (limit["first_collision"] or {}).get("angle_deg"), "return", direction_agreement[load], flush=True)


if __name__ == "__main__":
    main()
