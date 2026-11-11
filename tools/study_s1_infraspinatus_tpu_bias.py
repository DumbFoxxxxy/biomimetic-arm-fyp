"""Compare shared FJ1500 TPU rest offsets with free-cap equilibrium and rigid readback."""

from copy import deepcopy
from itertools import combinations

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save, load_component
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_terminal_reposition import collisions, with_bounds
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def rest_geometry(model, ports, offset):
    adjusted_model, adjusted_ports = deepcopy(model), deepcopy(ports)
    for key in ("neutral_cap_pin_source_mm", "neutral_crimp_outer_exit_A_source_mm"):
        adjusted_model[key] = (np.array(model[key]) + offset).tolist()
    path = np.array(model["active_strand_centreline_source_mm"])
    path[:2] += offset
    adjusted_model["active_strand_centreline_source_mm"] = path.tolist()
    key = "outer_face_source_mm"
    adjusted_ports["distal"][key] = (np.array(ports["distal"][key]) + offset).tolist()
    return adjusted_model, adjusted_ports


def collision_state(solids, owners, state, offset, gh, centre, pin_a, pin_b, axis_a, axis_b):
    cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
    phi_a, phi_b = np.deg2rad(state["swivel_angles_deg"])
    cap_pin = pin_a + offset + state["cap_translation_source_mm"]
    cap_pose = rigid_transform(cap, pin_a, cap_pin)
    eye_pose = rigid_transform(cap * Rotation.from_rotvec(axis_a * phi_a), pin_a, cap_pin)
    bone_eye_pose = rigid_transform(gh * Rotation.from_rotvec(axis_b * phi_b), pin_b, centre + gh.apply(pin_b - centre))
    arm_pose = rigid_transform(gh, centre, centre)
    transforms = dict(cap=cap_pose, eye=eye_pose, bone_eye=bone_eye_pose, humerus=arm_pose)
    groups = {name: {} for name in (*transforms, "fixed")}
    for key, part in solids.items():
        group = "fixed"
        if owners[key] == f"free_terminal:{HEAD}":
            group = "cap"
        elif key in (f"link:{HEAD}/muscle_groove_ring", f"link:{HEAD}/muscle_uncrimped_sleeve"):
            group = "eye"
        elif key in (f"link:{HEAD}/bone_groove_ring", f"link:{HEAD}/bone_uncrimped_sleeve"):
            group = "bone_eye"
        elif owners[key] == "humerus_ideal_attachment":
            group = "humerus"
        groups[group][key] = part if group == "fixed" else part.transform(transforms[group])
    bounds = {name: with_bounds(parts) for name, parts in groups.items()}
    hits = []
    moving_groups = {"cap", "eye", "bone_eye"}
    for left, right in combinations(bounds, 2):
        if {left, right} & moving_groups:
            hits.extend(collisions(bounds[left], bounds[right], np.zeros(3)))
    return cap_pin, sorted(hits, key=lambda row: -row["overlap_mm3"])


def main():
    folder = DEST.parent / "s1_local_muscle_load"
    output = folder / "FJ1500_cap_relief"
    previous_path = folder / "y_terminal_retained_lead_equilibrium.json"
    previous = read(previous_path)
    geometry = load_reference(previous["local_geometry"])
    assembly = load_reference(previous["shoulder_assembly"])
    model = next(row for row in load_reference(previous["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    sleeve_file = output / "compact_sleeve_screen.json"
    sleeve_report = read(sleeve_file)
    sources = dict(assembly["sources"])
    sleeve_key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    if sources[sleeve_key] != sleeve_report["original"] or sleeve_report["selected"]["hits"]:
        raise ValueError("Compact sleeve does not match this source model")
    sources[sleeve_key] = sleeve_report["selected"]["source"]
    solids = {key: load_component(row) for key, row in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    owners = assembly["intended_body_owners"]
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    pin_a, pin_b = (np.array(model[f"neutral_{end}_pin_source_mm"]) for end in ("cap", "bone"))
    axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
    origin = np.array(route["ports"]["origin"]["inner_face_source_mm"])
    stop = np.array(route["ports"]["distal"]["outer_face_source_mm"])
    axial = (stop - origin) / np.linalg.norm(stop - origin)
    outward = np.array([-1., 0., 0.])
    outward -= (outward @ axial) * axial
    outward /= np.linalg.norm(outward)
    stiffness = [.2, .4, 50., 1.]
    levels = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 8., 6., 4., 2., 0.]
    results = []
    for distance in (0., 3., 6., 9., 12.):
        offset = outward * distance
        adjusted_model, adjusted_ports = rest_geometry(model, route["ports"], offset)
        for angle in (0., 30., 60., 96.):
            gh = Rotation.from_rotvec(flexion * np.deg2rad(angle))
            states = solve_branch(adjusted_model, adjusted_ports, centre, gh, stiffness, levels, retain_leads=True)
            for state in states:
                pin, hits = collision_state(solids, owners, state, offset, gh, centre, pin_a, pin_b, axis_a, axis_b)
                state["resolved_cap_pin_source_mm"] = pin.tolist()
                state["hits"] = hits
            failing = [index for index, state in enumerate(states) if state["hits"]]
            results.append(dict(rest_offset_mm=distance, rest_offset_source_mm=offset.tolist(), GH_flexion_deg=angle,
                                rows=states, first_failure_index=failing[0] if failing else None,
                                rigid_conflict_states=len(failing)))
            print(f"Rest offset {distance:g} mm, GH {angle:g} deg: {len(failing)}/{len(states)} conflicts", flush=True)
    save(output / "tpu_rest_offset_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_tpu_bias.py"),
        equilibrium_solver=provenance(ROOT / "tools/study_s1_y_terminal_equilibrium.py"),
        reference_equilibrium=provenance(previous_path), compact_sleeve=provenance(sleeve_file),
        trial_direction_source=outward.tolist(), stiffness=stiffness, results=results,
        scope="A common lateral rest offset at the proximal-anchored TPU end, not a pose-by-pose imposed shift."
              " Each cap position, orientation and two swivels are solved under cable load and TPU restoring forces."
              " The original folded tendon lead length hypothesis and the same point-bend lower bound remain."
              " Rigid overlaps are checked after the solve; there is no contact force, flexible mesh or support CAD."))


if __name__ == "__main__":
    main()
