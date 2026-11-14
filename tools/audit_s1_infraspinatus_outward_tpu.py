"""Check whether the solved FJ1500 cap path fits the existing rest-shaped lattice."""

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from study_elbow_series_package import bounded, overlaps
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1500_outward_path_screen.json"
    screen = read(path_file)
    rest_path = folder / "lateral_rest_tpu_geometry.json"
    rest = read(rest_path)
    reference = load_reference(load_reference(rest["bias_screen"])["reference_equilibrium"])
    geometry = load_reference(reference["local_geometry"])
    assembly = load_reference(screen["shoulder_assembly"])
    if rest["original"] != assembly["sources"]["tpu:infraspinatus"] or screen["selected"]["rest_mm"] != rest["rest_offset_mm"]:
        raise ValueError("Existing mesh does not have the selected rest shape")
    route = next(row for row in load_reference(screen["input_route"])["routes"] if row["head"] == HEAD)
    model = next(row for row in load_reference(reference["kinematics"])["models"] if row["head"] == HEAD)
    origin, stop = (np.array(route["ports"][name][face]) for name, face in
                    (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm")))
    length = np.linalg.norm(stop-origin)
    axis = (stop-origin) / length
    pin = np.array(model["neutral_cap_pin_source_mm"]) + rest["rest_offset_source_mm"]
    mesh = load_component(rest["candidate"])
    points, faces = mesh_geom.manifold_mesh_to_arrays(mesh.to_mesh64())
    mf = A.I._manifold()
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    start, end = rest["original_mesh_deflection_fraction_limits"]
    s = np.clip(((points-origin)@axis/length-start)/(end-start), 0., 1.)
    w = s**2 * (3.-2.*s)
    dw = (6.*s*(1.-s)/((end-start)*length))[:, None]*axis
    transition = (s > 0.) & (s < 1.)
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    targets = {key: load_component(source) for key, source in assembly["sources"].items()
               if key.startswith(("bone:", "GH_cap_"))}
    rows = []
    for state in screen["outward"]:
        angle = state["GH_deg"]
        if angle not in (0., 30., 60., 68., 75., 85., 96.):
            continue
        gh = Rotation.from_rotvec(flexion*np.deg2rad(angle))
        cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
        delta = cap.as_matrix()-np.eye(3)
        shift = (points-pin) @ delta.T + state["cap_translation_source_mm"]
        jacobian = np.eye(3) + w[transition, None, None]*delta + shift[transition, :, None]*dw[transition, None, :]
        determinant = np.linalg.det(jacobian)
        singular = np.linalg.svd(jacobian, compute_uv=False)
        record = dict(GH_deg=angle, input_N=state["commanded_muscle_tension_N"],
                      minimum_sampled_deformation_determinant=float(determinant.min()),
                      maximum_local_stretch_ratio=float(singular.max()), minimum_local_stretch_ratio=float(singular.min()),
                      nonpositive_sample_count=int((determinant <= 0.).sum()), checked_transition_vertices=len(determinant))
        if record["nonpositive_sample_count"] == 0:
            moved = mf.Manifold(mf.Mesh64(np.ascontiguousarray(points+w[:, None]*shift), faces))
            record["mesh_status"] = str(moved.status())
            if moved.status() == mf.Error.NoError:
                transform = rigid_transform(gh, centre, centre)
                posed = {key: solid.transform(transform) if assembly["intended_body_owners"][key] == "humerus_ideal_attachment" else solid
                         for key, solid in targets.items()}
                record["bone_and_cover_hits_mm3"] = overlaps(moved, bounded(posed))
        rows.append(record)
        print(angle, record["input_N"], record["minimum_sampled_deformation_determinant"], record.get("bone_and_cover_hits_mm3"), flush=True)
    save(folder / "FJ1500_outward_TPU_compatibility.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_infraspinatus_outward_tpu.py"),
        trajectory=provenance(path_file), rest_geometry=provenance(rest_path), cases=rows,
        scope="Imposed smooth axial blend from fixed proximal mesh to actual solved rigid cap motion."
              " Negative Jacobian rejects this particular interpolation, not the physical muscle."
              " Stretch and solid overlaps are geometric indicators; no material law or contact equilibrium."))


if __name__ == "__main__":
    main()
