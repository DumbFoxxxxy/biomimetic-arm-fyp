"""Compare a tension-only series take-up element at the FJ1500 short tendon."""

from itertools import product

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from check_shoulder_front_liner_keepers import is_compliant
from study_elbow_series_package import bounded, overlaps
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1500_outward_path_screen.json"
    outward = read(path_file)
    tpu_file = folder / "FJ1500_outward_TPU_compatibility.json"
    original_check = read(tpu_file)
    if original_check["trajectory"] != provenance(path_file):
        raise ValueError("Baseline mesh check uses a different cap trajectory")
    smooth = load_reference(outward["parent"])
    coupling = load_reference(smooth["coupling"])
    flexion_path = load_reference(coupling["prior_path"])
    rest = load_reference(original_check["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    prior = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(prior["local_geometry"])
    assembly = load_reference(outward["shoulder_assembly"])
    model = next(row for row in load_reference(prior["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(outward["input_route"])["routes"] if row["head"] == HEAD)
    sleeve_key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    sleeve = load_reference(flexion_path["sleeve"])
    if sleeve["original"] != assembly["sources"][sleeve_key]:
        raise ValueError("Shoulder sleeve does not match the reference assembly")
    sources = {**assembly["sources"], sleeve_key: sleeve["selected"]["source"]}
    parts = {key: load_component(source) for key, source in sources.items()
             if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
             and not key.endswith("/tendon_wire")}
    owners = assembly["intended_body_owners"]
    bone_keys = [key for key in parts if key.startswith(("bone:", "GH_cap_"))]
    shoulder_centre = np.array(geometry["gh_centre_source_mm"])
    flexion_axis = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    pins = [np.array(model[f"neutral_{k}_pin_source_mm"]) for k in ("cap", "bone")]
    pin_axes = [np.array(model[f"{k}_pin_axis_source"]) for k in ("cap", "bone")]
    origin, stop = (np.array(route["ports"][k][face]) for k, face in
                    (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm")))
    span = np.linalg.norm(stop-origin)
    axis = (stop-origin) / span
    choice = outward["selected"]
    offset = np.array(rest["rest_offset_source_mm"])
    direction = Rotation.from_rotvec(-axis*np.deg2rad(choice["coupling_roll_deg"])).apply(outward["direction_source"])
    adjusted, ports = rest_geometry(model, route["ports"], offset)
    mesh = load_component(rest["candidate"])
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(mesh.to_mesh64())
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    mf = A.I._manifold()
    lo, hi = rest["original_mesh_deflection_fraction_limits"]
    fraction = np.clip(((vertices-origin) @ axis/span-lo)/(hi-lo), 0., 1.)
    w = fraction**2*(3.-2.*fraction)
    dw = (6.*fraction*(1.-fraction)/(span*(hi-lo)))[:, None]*axis
    in_blend = (fraction > 0.) & (fraction < 1.)

    def evaluate(spring, multiplier, extra, angles, loads):
        stiffness = list(bias["stiffness"])
        stiffness[:2] = [multiplier*k for k in stiffness[:2]]
        cases = []
        for angle in angles:
            gh = Rotation.from_rotvec(flexion_axis*np.deg2rad(angle))
            arm_pose = rigid_transform(gh, shoulder_centre, shoulder_centre)
            bones = bounded({key: parts[key].transform(arm_pose) if owners[key] == "humerus_ideal_attachment" else parts[key]
                             for key in bone_keys})
            states = solve_branch(adjusted, ports, shoulder_centre, gh, stiffness, loads,
                                  retain_leads=True, installed_length_offset=extra,
                                  lateral_coupling=(direction, choice["gain"]), strand_bending_Nmm2=10.,
                                  series_stiffness_N_per_mm=spring)
            for state in states:
                pin, rigid = collision_state(parts, owners, state, offset, gh, shoulder_centre, *pins, *pin_axes)
                rotation = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
                delta = rotation.as_matrix()-np.eye(3)
                shift = (vertices-pins[0]-offset) @ delta.T + state["cap_translation_source_mm"]
                jacobian = np.eye(3) + w[in_blend, None, None]*delta + shift[in_blend, :, None]*dw[in_blend, None, :]
                minimum = float(np.linalg.det(jacobian).min())
                soft = None
                if minimum > 0.:
                    posed = mf.Manifold(mf.Mesh64(np.ascontiguousarray(vertices+w[:, None]*shift), faces))
                    if posed.status() == mf.Error.NoError:
                        soft = overlaps(posed, bones)
                cases.append(dict(GH_deg=angle, cap_pin_source_mm=pin.tolist(), rigid_hits=rigid,
                                  minimum_TPU_map_determinant=minimum, interpolated_TPU_bone_hits_mm3=soft, **state))
        return cases

    def metrics(rows):
        return dict(rigid_count=sum(bool(row["rigid_hits"]) for row in rows),
                    folded_count=sum(row["interpolated_TPU_bone_hits_mm3"] is None for row in rows),
                    TPU_hit_count=sum(bool(row["interpolated_TPU_bone_hits_mm3"]) for row in rows),
                    peak_TPU_overlap_mm3=max((max((row["interpolated_TPU_bone_hits_mm3"] or {}).values(), default=0.) for row in rows)),
                    max_cap_translation_mm=max(float(np.linalg.norm(row["cap_translation_source_mm"])) for row in rows),
                    max_series_extension_mm=max(row.get("series_extension_mm", 0.) for row in rows),
                    min_tendon_tension_at_10N=min(row["short_tendon_tension_N"] for row in rows if row["commanded_muscle_tension_N"] == 10.),
                    peak_tendon_tension_N=max(row["short_tendon_tension_N"] for row in rows))

    angles = (0., 30., 60., 68., 75., 85., 96.)
    loads = (0., 5., 10.)
    reference = evaluate(None, 1., 0., angles, loads)
    for old, new in zip(original_check["cases"], reference, strict=True):
        np.testing.assert_allclose(new["minimum_TPU_map_determinant"], old["minimum_sampled_deformation_determinant"], atol=2e-6, rtol=0.)
        if new["interpolated_TPU_bone_hits_mm3"] is not None:
            matches = old["bone_and_cover_hits_mm3"]
            if set(matches) != set(new["interpolated_TPU_bone_hits_mm3"]):
                raise ValueError("Baseline TPU collision pairs changed")
    variants = []
    for spring, scale, extra in product((.25, .5, 1.), (1., 2., 4.), (0., 6.)):
        rows = evaluate(spring, scale, extra, angles, loads)
        result = dict(series_stiffness_N_per_mm=spring, TPU_translation_stiffness_multiplier=scale,
                      extra_free_length_mm=extra, metrics=metrics(rows), cases=rows)
        variants.append(result)
        print(spring, scale, extra, result["metrics"], flush=True)
    # Keep the complete comparisons. A slack tendon cannot be treated as
    # a functioning replacement merely because its geometry clears.
    candidates = [trial for trial in variants if trial["metrics"]["min_tendon_tension_at_10N"] > 0.]
    def order(trial):
        m = trial["metrics"]
        return (m["rigid_count"], m["folded_count"], m["peak_TPU_overlap_mm3"],
                m["max_cap_translation_mm"], -m["min_tendon_tension_at_10N"])
    best = min(candidates, key=order) if candidates else None
    validation = []
    if best:
        for branch, samples, force in (("outward", tuple(range(0, 97, 4)), loads),
                                       ("return", tuple(range(96, -1, -4)), loads[::-1])):
            rows = evaluate(best["series_stiffness_N_per_mm"], best["TPU_translation_stiffness_multiplier"],
                            best["extra_free_length_mm"], samples, force)
            validation.append(dict(branch=branch, metrics=metrics(rows), cases=rows))
    save(folder / "FJ1500_series_compliance_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_series_compliance.py"),
        solver=provenance(ROOT / "tools/study_s1_y_terminal_equilibrium.py"),
        baseline_trajectory=provenance(path_file), baseline_mesh_check=provenance(tpu_file), baseline_metrics=metrics(reference),
        shoulder_probe_angles_deg=list(angles), comparisons=variants,
        comparison_for_validation={k: v for k, v in best.items() if k != "cases"} if best else None,
        validation=validation, scope="Hypothetical tension-only series take-up. Effective extension is required mechanism travel, not steel strain."
        " Exact sampled bone/cover overlap uses one imposed TPU interpolation; no TPU constitutive/contact model."
        " Prescribed shoulder poses and hypothetical stiffness, no installed spring, guide or free drive."))
    print("Validation", [(row["branch"], row["metrics"]) for row in validation], flush=True)


if __name__ == "__main__":
    main()
