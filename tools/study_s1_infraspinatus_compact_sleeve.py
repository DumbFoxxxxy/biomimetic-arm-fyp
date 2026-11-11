"""Screen compact FJ1500 crimp sleeves while keeping their outlet and bores."""

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import CAP, HEAD
from study_s1_infraspinatus_terminal_reposition import collisions, with_bounds
from study_s1_local_force_geometry import load_reference


SLEEVE = f"link:{HEAD}/muscle_uncrimped_sleeve"


def sleeve(centre, u, v, length, radius):
    tips = [centre + u * 8.9 + v * sign * .45 for sign in (-1, 1)]
    barrels = [rod(tip - u * length, tip, radius) for tip in tips]
    bores = [rod(tip - u * 4., tip + u, .4) for tip in tips]
    return barrels[0] + barrels[1] - bores[0] - bores[1]


def main():
    folder = DEST.parent / "s1_local_muscle_load"
    output = folder / "FJ1500_cap_relief"
    audit_file = folder / "y_retained_lead_clearance.json"
    audit = read(audit_file)
    equilibrium = load_reference(audit["equilibrium"])
    geometry = load_reference(equilibrium["local_geometry"])
    assembly = load_reference(equilibrium["shoulder_assembly"])
    model = next(row for row in load_reference(equilibrium["kinematics"])["models"] if row["head"] == HEAD)
    axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
    pin_a, pin_b = (np.array(model[f"neutral_{end}_pin_source_mm"]) for end in ("cap", "bone"))
    path = np.array(model["active_strand_centreline_source_mm"])
    u = path[1] - path[0]
    u /= np.linalg.norm(u)
    v = np.cross(axis_a, u)
    centre = pin_a + axis_a * 2.2
    if not np.allclose(centre + u * 8.9 + v * .45, path[0], atol=1e-8, rtol=0.):
        raise ValueError("Reconstructed sleeve no longer matches the actual source outlet")
    axes = dict(zip(geometry["axes"], np.array(geometry["axes_source"])))
    gh_centre = np.array(geometry["gh_centre_source_mm"])
    owners = assembly["intended_body_owners"]
    solids = {key: load_component(row) for key, row in assembly["sources"].items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    original = solids.pop(SLEEVE)
    reference = sleeve(centre, u, v, 2.5, 1.)
    source_difference = float((original - reference).volume() + (reference - original).volume())
    if source_difference > .001:
        raise ValueError(f"Sleeve source mismatch {source_difference:.6g} mm³")
    cases = {(row["prescribed_axis"], row["prescribed_GH_deg"], row["assumed_length_offset_mm"]): row["rows"]
             for row in equilibrium["cases"] if row["head"] == HEAD}
    probes = []
    for row in audit["records"]:
        if row["head"] != HEAD:
            continue
        name, angle, offset = (row[key] for key in ("prescribed_axis", "prescribed_GH_deg", "assumed_length_offset_mm"))
        state = cases[name, angle, offset][row["step"]]
        gh = Rotation.from_rotvec(axes.get(name, axes["flexion"]) * np.deg2rad(angle))
        cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
        phi_a, phi_b = np.deg2rad(state["swivel_angles_deg"])
        free_pin = pin_a + state["cap_translation_source_mm"]
        cap_pose = rigid_transform(cap, pin_a, free_pin)
        eye_pose = rigid_transform(cap * Rotation.from_rotvec(axis_a * phi_a), pin_a, free_pin)
        bone_eye_pose = rigid_transform(gh * Rotation.from_rotvec(axis_b * phi_b), pin_b,
                                        gh_centre + gh.apply(pin_b - gh_centre))
        arm_pose = rigid_transform(gh, gh_centre, gh_centre)
        context = {}
        for key, part in solids.items():
            if owners[key] == f"free_terminal:{HEAD}":
                part = part.transform(cap_pose)
            elif key == f"link:{HEAD}/muscle_groove_ring":
                part = part.transform(eye_pose)
            elif key in {f"link:{HEAD}/bone_{end}" for end in ("groove_ring", "uncrimped_sleeve")}:
                part = part.transform(bone_eye_pose)
            elif owners[key] == "humerus_ideal_attachment":
                part = part.transform(arm_pose)
            context[key] = part
        info = dict(axis=name, GH_deg=angle, offset_mm=offset, step=row["step"],
                    input_N=row["commanded_muscle_tension_N"])
        probes.append((info, eye_pose, with_bounds(context)))

    def screen(part):
        rows = []
        for info, eye_pose, context in probes:
            hits = collisions(with_bounds({SLEEVE: part.transform(eye_pose)}), context, np.zeros(3))
            if hits:
                rows.append(dict(**info, hits=hits))
        return rows

    def minimum_cap_gap(part):
        nearest = []
        for info, eye_pose, context in probes:
            cap = next(body for key, body, _ in context if key == CAP)
            gap = float(A.I.gap(part.transform(eye_pose), cap, 1.))
            nearest.append(dict(**info, gap_capped_1mm=gap))
        return min(nearest, key=lambda row: row["gap_capped_1mm"])

    baseline = screen(original)
    trials = []
    for length in (2.5, 2.2, 1.9, 1.6, 1.3):
        for radius in (1., .9, .8, .7):
            part = sleeve(centre, u, v, length, radius)
            if str(part.status()).split(".")[-1] != "NoError" or len(part.decompose()) != 1:
                raise ValueError("Candidate sleeve lost connected material")
            hits = screen(part)
            trials.append(dict(length_mm=length, outer_radius_mm=radius, outer_wall_nominal_mm=radius - .4,
                               volume_mm3=float(part.volume()), interfering_states=len(hits), hits=hits,
                               minimum_cap_gap=minimum_cap_gap(part) if not hits else None))
    selected = next((row for row in trials if not row["hits"]
                     and row["minimum_cap_gap"]["gap_capped_1mm"] > 0.), None)
    if selected:
        part = sleeve(centre, u, v, selected["length_mm"], selected["outer_radius_mm"])
        actual, source = export(part, output / "FJ1500_compact_muscle_crimp_sleeve.obj",
                                np.array(assembly["sources"][SLEEVE]["translation_source_mm"]))
        if screen(actual):
            raise RuntimeError("Sleeve export changed the clearance result")
        selected = {**selected, "source": source, "minimum_cap_gap": minimum_cap_gap(actual)}
    save(output / "compact_sleeve_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_compact_sleeve.py"),
        pose_audit=provenance(audit_file), original=assembly["sources"][SLEEVE],
        source_reconstruction_difference_mm3=source_difference, sampled_poses=len(probes),
        baseline_hits=baseline, trials=trials, selected=selected,
        sleeve_outlet_source_mm=path[0].tolist(), bore_radius_mm=.4, wire_radius_mm=.3,
        nominal_ligament_between_bores_mm=.1,
        scope="Only sleeve-to-rigid comparisons are checked at the prior diagnostic states."
              " Pin, outlet, dual bore centers and assumed tendon length remain fixed."
              " Grip, manufacturing and other rigid or compliant conflicts are not checked."))
    print("Poses", len(probes), "baseline hit states", len(baseline), "selected", selected, flush=True)
    print("Trials (length, radius, hit count):",
          [(row["length_mm"], row["outer_radius_mm"], row["interfering_states"]) for row in trials], flush=True)
    print("Minimum cap gaps:", [(row["length_mm"], row["outer_radius_mm"],
                                row["minimum_cap_gap"]["gap_capped_1mm"])
                               for row in trials if row["minimum_cap_gap"]], flush=True)


if __name__ == "__main__":
    main()
