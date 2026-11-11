"""Screen geometric endpoint alternatives at fixed short-tendon length."""

import argparse

import numpy as np
from scipy.linalg import null_space
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_local_force_geometry import load_reference


def bounds_overlap(left, right):
    return np.all(left[:3] <= right[3:]) and np.all(right[:3] <= left[3:])


def collisions(mobile, context, translation, stop_at_first=False):
    hits = []
    for mobile_key, part, bounds in mobile:
        shifted = bounds + np.tile(translation, 2)
        posed = None
        for other_key, other, other_bounds in context:
            if not bounds_overlap(shifted, other_bounds):
                continue
            if posed is None:
                posed = part.translate(translation)
            volume = float(A.I.intersection_volume(posed, other))
            if volume > .001:
                hits.append(dict(moving=mobile_key, obstruction=other_key, overlap_mm3=volume))
                if stop_at_first:
                    return hits
    return hits


def with_bounds(parts):
    return [(key, part, np.array(part.bounding_box())) for key, part in parts.items()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact-sleeve", action="store_true")
    args = parser.parse_args()
    output = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    audit_file = output.parent / "y_retained_lead_clearance.json"
    audit = read(audit_file)
    equilibrium = load_reference(audit["equilibrium"])
    geometry = load_reference(equilibrium["local_geometry"])
    assembly = load_reference(equilibrium["shoulder_assembly"])
    model = next(row for row in load_reference(equilibrium["kinematics"])["models"] if row["head"] == HEAD)
    owners = assembly["intended_body_owners"]
    sources = dict(assembly["sources"])
    sleeve_report_file = output / "compact_sleeve_screen.json"
    sleeve_key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    if args.compact_sleeve:
        sleeve_report = read(sleeve_report_file)
        if (sleeve_report["pose_audit"] != provenance(audit_file)
                or sleeve_report["original"] != sources[sleeve_key]
                or sleeve_report["selected"]["hits"]):
            raise ValueError("Sleeve clearance is not verified against these source poses")
        sources[sleeve_key] = sleeve_report["selected"]["source"]
    solids = {key: load_component(row) for key, row in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    axis = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    centre = np.array(geometry["gh_centre_source_mm"])
    pin_a = np.array(model["neutral_cap_pin_source_mm"])
    pin_b = np.array(model["neutral_bone_pin_source_mm"])
    axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
    outputs = []
    for case in equilibrium["cases"]:
        if (case["head"] != HEAD or case["prescribed_axis"] != "flexion"
                or case["prescribed_GH_deg"] not in (60., 96.) or case["assumed_length_offset_mm"] != 0.):
            continue
        angle = case["prescribed_GH_deg"]
        gh = Rotation.from_rotvec(axis * np.deg2rad(angle))
        arm = rigid_transform(gh, centre, centre)
        for step in range(11):
            state = case["rows"][step]
            cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
            phi_a, phi_b = np.deg2rad(state["swivel_angles_deg"])
            free_pin = pin_a + state["cap_translation_source_mm"]
            cap_frame = rigid_transform(cap, pin_a, free_pin)
            eye_a = rigid_transform(cap * Rotation.from_rotvec(axis_a * phi_a), pin_a, free_pin)
            eye_b = rigid_transform(gh * Rotation.from_rotvec(axis_b * phi_b), pin_b,
                                    centre + gh.apply(pin_b - centre))
            moving, context = {}, {}
            for key, part in solids.items():
                if owners[key] == f"free_terminal:{HEAD}":
                    moving[key] = part.transform(cap_frame)
                elif key in {f"link:{HEAD}/muscle_{end}" for end in ("groove_ring", "uncrimped_sleeve")}:
                    moving[key] = part.transform(eye_a)
                elif key in {f"link:{HEAD}/bone_{end}" for end in ("groove_ring", "uncrimped_sleeve")}:
                    context[key] = part.transform(eye_b)
                else:
                    context[key] = part.transform(arm) if owners[key] == "humerus_ideal_attachment" else part
            baseline = next(row for row in audit["records"] if row["head"] == HEAD
                            and row["prescribed_GH_deg"] == angle and row["step"] == step
                            and row["assumed_length_offset_mm"] == 0.)
            baseline_hits = [hit for hit in baseline["overlaps"]
                             if not (args.compact_sleeve and sleeve_key in (hit["left"], hit["right"]))]
            invariant = [hit for hit in baseline_hits
                         if (hit["left"] in moving) == (hit["right"] in moving)]
            row = dict(GH_deg=angle, input_N=state["commanded_muscle_tension_N"], invariant_overlaps=invariant,
                       baseline_overlaps=baseline_hits, candidate=None, tested_count=0)
            if not invariant:
                moving_bounds, context_bounds = with_bounds(moving), with_bounds(context)
                p, q = np.array(state["lead_tips_source_mm"])
                radius = np.linalg.norm(p - q)
                radial = (p - q) / radius
                tangent_plane = null_space(radial.reshape(1, 3))
                azimuths = np.linspace(0., 2. * np.pi, 32, endpoint=False)
                for distance in np.arange(0., 21., 1.):
                    tilt = 2. * np.arcsin(distance / (2. * radius))
                    for azimuth in azimuths[:1] if distance == 0. else azimuths:
                        tangent = tangent_plane @ np.array([np.cos(azimuth), np.sin(azimuth)])
                        shift = radius * ((np.cos(tilt) - 1.) * radial + np.sin(tilt) * tangent)
                        row["tested_count"] += 1
                        if collisions(moving_bounds, context_bounds, shift, stop_at_first=True):
                            continue
                        error = abs(np.linalg.norm(p + shift - q) - radius)
                        if error > 1e-8:
                            raise RuntimeError("Length-preserving displacement failed")
                        row["candidate"] = dict(displacement_mm=float(np.linalg.norm(shift)),
                                                shift_source_mm=shift.tolist(), chord_length_error_mm=error)
                        break
                    if row["candidate"] is not None:
                        break
            outputs.append(row)
            print(angle, row["input_N"], row["candidate"], len(invariant), row["tested_count"], flush=True)
    filename = ("terminal_reposition_with_compact_sleeve.json" if args.compact_sleeve
                else "terminal_reposition_screen.json")
    save(output / filename, dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_terminal_reposition.py"),
        pose_audit=provenance(audit_file), rows=outputs, shift_step_mm=1., azimuth_count=32,
        compact_sleeve=provenance(sleeve_report_file) if args.compact_sleeve else None,
        shift_limit_mm=20.,
        scope="Independently translate the original moving terminal on a sphere about the distal lead tip."
              " Keeps the span length, cap/eye orientations and prescribed GH unchanged."
              " Checks rigid terminal component intersections only. No path interpolation, equilibrium,"
              " stiffness, clearance margin, muscle deformation or actual wire bends are verified."))


if __name__ == "__main__":
    main()
