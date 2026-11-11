"""Map S1-Y2 free cap collisions into its neutral geometry and screen cuts."""

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_s1_local_force_geometry import load_reference


HEAD = "FJ1500"
CAP = f"free_terminal:{HEAD}"


def union(solids):
    manifold = A.I._manifold()
    return manifold.Manifold.batch_boolean(list(solids), manifold.OpType.Add)


def inverse(transform):
    rotation = transform[:, :3].T
    return np.c_[rotation, -rotation @ transform[:, 3]]


def guard_solids(model, route):
    axis = np.array(model["cap_pin_axis_source"])
    pin = np.array(model["neutral_cap_pin_source_mm"])
    port = route["ports"]["distal"]
    inner, outer, normal = (np.array(port[key]) for key in
                            ("inner_face_source_mm", "outer_face_source_mm", "normal_source"))
    tangent = np.array(route["distal_terminal_override"]["tangent_source"])
    tangent -= normal * (tangent @ normal)
    tangent /= np.linalg.norm(tangent)
    frame = np.c_[tangent, np.cross(normal, tangent), normal, inner]
    manifold = A.I._manifold()
    hook_seats = union(manifold.Manifold.cube([2., 2.4, 3.]).translate([-1., y, 0.]).transform(frame)
                       for y in (2.6, -5.))
    return dict(pin_bearing=rod(pin-axis*1.5, pin+axis*1.5, 3.2),
                stop_housing=rod(inner, outer, 2.4), TPU_hook_seats=hook_seats)


def main():
    base = DEST.parent / "s1_local_muscle_load"
    output = base / "FJ1500_cap_relief"
    output.mkdir(exist_ok=True)
    audit_file = base / "y_retained_lead_clearance.json"
    audit = read(audit_file)
    equilibrium = load_reference(audit["equilibrium"])
    assembly = load_reference(audit["shoulder_assembly"])
    geometry = load_reference(equilibrium["local_geometry"])
    model = next(row for row in load_reference(equilibrium["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    axes = dict(zip(geometry["axes"], np.array(geometry["axes_source"])))
    centre = np.array(geometry["gh_centre_source_mm"])
    pin = np.array(model["neutral_cap_pin_source_mm"])
    cap_axis = np.array(model["cap_pin_axis_source"])
    original = load_component(assembly["sources"][CAP])
    guards = {key: original ^ solid for key, solid in guard_solids(model, route).items()}
    obstacles = {}
    cases = {(row["prescribed_axis"], row["prescribed_GH_deg"], row["assumed_length_offset_mm"]): row["rows"]
             for row in equilibrium["cases"] if row["head"] == HEAD}
    patches = []
    rows = []
    for record in audit["records"]:
        if record["head"] != HEAD or not record["overlaps"]:
            continue
        name, angle, offset = (record[key] for key in
                               ("prescribed_axis", "prescribed_GH_deg", "assumed_length_offset_mm"))
        state = cases[name, angle, offset][record["step"]]
        cap_rotation = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
        gh = Rotation.from_rotvec(axes[name] * np.deg2rad(angle))
        free_pin = pin + state["cap_translation_source_mm"]
        cap_frame = rigid_transform(cap_rotation, pin, free_pin)
        cap_eye = cap_rotation * Rotation.from_rotvec(cap_axis * np.deg2rad(state["swivel_angles_deg"][0]))
        other_overlaps = []
        cap_overlaps = []
        for hit in record["overlaps"]:
            if CAP not in (hit["left"], hit["right"]):
                other_overlaps.append(hit)
                continue
            key = hit["right"] if hit["left"] == CAP else hit["left"]
            if key not in obstacles:
                obstacles[key] = load_component(assembly["sources"][key])
            body = obstacles[key]
            owner = assembly["intended_body_owners"][key]
            if owner == "humerus_ideal_attachment":
                body = body.transform(rigid_transform(gh, centre, centre))
            elif key in (f"link:{HEAD}/muscle_groove_ring", f"link:{HEAD}/muscle_uncrimped_sleeve"):
                body = body.transform(rigid_transform(cap_eye, pin, free_pin))
            elif owner not in ("fixed_context", "scapula_ideal_fastened"):
                raise ValueError(f"Unmapped obstruction: {key} ({owner})")
            patch = original ^ body.transform(inverse(cap_frame))
            if abs(patch.volume() - hit["overlap_mm3"]) > .0001:
                raise ValueError(f"Mismatch in inverse-mapped overlap: {key}")
            patches.append((angle, key, patch))
            cap_overlaps.append(dict(obstacle=key, overlap_mm3=float(patch.volume())))
        rows.append(dict(GH_deg=angle, offset_mm=offset, step=record["step"],
                         input_N=record["commanded_muscle_tension_N"], cap_overlaps=cap_overlaps,
                         other_overlaps=other_overlaps))

    results = []
    for label, selected in (("flexion_60", [p for p in patches if p[0] == 60.]),
                            ("flexion_96", [p for p in patches if p[0] == 96.]),
                            ("combined", patches)):
        required = union(patch for _, _, patch in selected)
        guards_hit = {key: float((guard ^ required).volume()) for key, guard in guards.items()}
        safe_cut = required - union(guards.values())
        candidate = original - safe_cut
        positive = [float(part.volume()) for part in candidate.decompose() if part.volume() > 1e-8]
        residuals = [(obstacle, float((candidate ^ patch).volume())) for _, obstacle, patch in selected]
        residuals = [dict(obstacle=obstacle, overlap_mm3=volume) for obstacle, volume in residuals if volume > .001]
        mesh = None
        if len(positive) == 1 and not residuals:
            actual, mesh = export(candidate, output / f"FJ1500_{label}_local_cut.obj",
                                  np.array(assembly["sources"][CAP]["translation_source_mm"]))
            if any((actual ^ patch).volume() > .001 for _, _, patch in selected):
                raise RuntimeError("Readback reintroduced a conflict")
        results.append(dict(label=label, patches=len(selected), required_removal_mm3=float(required.volume()),
                            protected_material_in_conflict_mm3=guards_hit,
                            removal_outside_protected_mm3=float(safe_cut.volume()),
                            candidate_positive_components_mm3=positive,
                            protected_cut_residuals=residuals, candidate_mesh=mesh))
        print(label, results[-1]["required_removal_mm3"], guards_hit, len(residuals), flush=True)
    save(output / "cap_relief_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_cap_relief.py"),
        pose_audit=provenance(audit_file), original_cap=assembly["sources"][CAP],
        original_cap_volume_mm3=float(original.volume()),
        protected_material_mm3={key: float(part.volume()) for key, part in guards.items()},
        results=results, pose_conflicts=rows,
        scope="Cap removal evaluated against prior sampled diagnostic poses with zero extra clearance."
              " Original bearing envelope, a 2.4 mm outer stop housing guard and hook landing footprints"
              " are protected for design comparison. Their dimensions are not strength tolerances."
              " Other collided hardware remains unchanged; a cap-only cut cannot resolve those contacts."))


if __name__ == "__main__":
    main()
