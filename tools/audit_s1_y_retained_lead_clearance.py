"""Screen original S1-Y2 hardware at representative lead-constrained states."""

from itertools import combinations

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_local_force_geometry import load_reference


def rigid_transform(rotation, old_pin, new_pin):
    matrix = rotation.as_matrix()
    return np.c_[matrix, new_pin - matrix @ old_pin]


def screen(components, moving):
    bounds = {key: np.array(part.bounding_box()) for key, part in components.items()}
    hits, nearest = [], None
    gap_cap = 1.
    for left, right in combinations(components, 2):
        if not ({left, right} & moving):
            continue
        lb, rb = bounds[left], bounds[right]
        distance = np.linalg.norm(np.maximum(np.maximum(lb[:3] - rb[3:], rb[:3] - lb[3:]), 0.))
        if distance > gap_cap:
            continue
        first, second = components[left], components[right]
        volume = float(A.I.intersection_volume(first, second)) if distance == 0. else 0.
        if volume > .001:
            hits.append(dict(left=left, right=right, overlap_mm3=volume))
        elif distance <= gap_cap:
            gap = A.I.gap(first, second, gap_cap)
            if gap < gap_cap:
                gap_cap, nearest = float(gap), [left, right]
    return dict(overlaps=sorted(hits, key=lambda row: -row["overlap_mm3"]),
                min_gap_capped_1mm=0. if hits else gap_cap, nearest_pair=nearest)


def main():
    folder = DEST.parent / "s1_local_muscle_load"
    input_file = folder / "y_terminal_retained_lead_equilibrium.json"
    equilibrium = read(input_file)
    geometry = load_reference(equilibrium["local_geometry"])
    assembly = load_reference(equilibrium["shoulder_assembly"])
    models = {row["head"]: row for row in load_reference(equilibrium["kinematics"])["models"]}
    axes = dict(zip(geometry["axes"], np.array(geometry["axes_source"])))
    centre = np.array(geometry["gh_centre_source_mm"])
    owners = assembly["intended_body_owners"]
    sources = {key: row for key, row in assembly["sources"].items()
               if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
               and not key.endswith("/tendon_wire")}
    solids = {key: load_component(row) for key, row in sources.items()}
    records = []
    for case in equilibrium["cases"]:
        head = case["head"]
        model = models[head]
        pin_a = np.array(model["neutral_cap_pin_source_mm"])
        pin_b = np.array(model["neutral_bone_pin_source_mm"])
        axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
        gh = Rotation.from_rotvec(axes.get(case["prescribed_axis"], axes["flexion"])
                                  * np.deg2rad(case["prescribed_GH_deg"]))
        shoulder = rigid_transform(gh, centre, centre)
        humeral = {key: part.transform(shoulder) if owners[key] == "humerus_ideal_attachment" else part
                   for key, part in solids.items()}
        sample_steps = (0, 5, 10, len(case["rows"]) - 1)
        if head == "FJ1500" and case["prescribed_axis"] == "flexion" and case["prescribed_GH_deg"] >= 60.:
            sample_steps = range(len(case["rows"]))
        for step in sample_steps:
            state = case["rows"][step]
            cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
            swivel_a, swivel_b = np.deg2rad(state["swivel_angles_deg"])
            moving_pin = pin_a + state["cap_translation_source_mm"]
            cap_transform = rigid_transform(cap, pin_a, moving_pin)
            eye_a = rigid_transform(cap * Rotation.from_rotvec(axis_a * swivel_a), pin_a, moving_pin)
            eye_b = rigid_transform(gh * Rotation.from_rotvec(axis_b * swivel_b), pin_b,
                                    centre + gh.apply(pin_b - centre))
            posed = dict(humeral)
            moving = set()
            for key, part in solids.items():
                if owners[key] == f"free_terminal:{head}":
                    posed[key] = part.transform(cap_transform)
                elif key in {f"link:{head}/muscle_{end}" for end in ("groove_ring", "uncrimped_sleeve")}:
                    posed[key] = part.transform(eye_a)
                elif key in {f"link:{head}/bone_{end}" for end in ("groove_ring", "uncrimped_sleeve")}:
                    posed[key] = part.transform(eye_b)
                else:
                    continue
                moving.add(key)
            hardware = screen(posed, moving)
            tip_turn_angles = None
            tip_reactions = None
            if state["short_tendon_tension_N"] > 1e-5:
                exit_a, exit_b = np.array(state["tendon_exits_source_mm"])
                tip_a, tip_b = np.array(state["lead_tips_source_mm"])
                vectors = np.diff([exit_a, tip_a, tip_b, exit_b], axis=0)
                directions = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
                tip_turn_angles = np.rad2deg(np.arccos(np.clip(np.sum(directions[:-1] * directions[1:], axis=1),
                                                               -1., 1.))).tolist()
                tip_reactions = (state["short_tendon_tension_N"]
                                 * np.linalg.norm(np.diff(directions, axis=0), axis=1)).tolist()
            records.append(dict(head=head, prescribed_axis=case["prescribed_axis"],
                                prescribed_GH_deg=case["prescribed_GH_deg"], step=step,
                                assumed_length_offset_mm=case["assumed_length_offset_mm"],
                                commanded_muscle_tension_N=state["commanded_muscle_tension_N"],
                                short_tendon_tension_N=state["short_tendon_tension_N"],
                                tip_point_turn_angles_deg=tip_turn_angles,
                                hypothetical_tip_support_reactions_N=tip_reactions, **hardware))
        print(head, case["prescribed_axis"], case["prescribed_GH_deg"], case["assumed_length_offset_mm"],
              "audited", flush=True)
    save(folder / "y_retained_lead_clearance.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_y_retained_lead_clearance.py"),
        equilibrium=provenance(input_file), shoulder_assembly=equilibrium["shoulder_assembly"],
        interference_threshold_mm3=.001, rigid_sources_count=len(sources), records=records,
        scope="Selected prescribed GH poses and load steps. Compare the moving cap, its hardware and both eye"
              " packages against original S1-Y2 rigid hardware. Opposite free cap and swivels stay at their"
              " neutral positions, but bone-owned parts follow GH. TPU, wire solids and bench fixtures are"
              " excluded. Adjacent supports can legitimately have near-zero gaps. Tip turns represent point"
              " bends in the lower-bound model; there is no designed support at those tips."))
    print(f"{sum(bool(row['overlaps']) for row in records)} of {len(records)} samples have rigid overlaps")


if __name__ == "__main__":
    main()
