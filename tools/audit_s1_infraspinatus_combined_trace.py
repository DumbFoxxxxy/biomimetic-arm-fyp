"""Compare two FJ1500 terminals on recorded shoulder poses from the combined drive."""

from copy import deepcopy

import mujoco
import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_series_elbow_spring_integration import load_candidate
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from audit_s1_y_retained_lead_clearance import rigid_transform
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from check_shoulder_front_liner_keepers import is_compliant
from study_elbow_series_package import bounded, overlaps
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    exit_path = folder / "FJ1500_humeral_exit_screen.json"
    revision = read(exit_path)
    selected = revision["selected"]
    baseline = load_reference(revision["baseline"])
    coupling = load_reference(load_reference(baseline["parent"])["coupling"])
    previous = load_reference(coupling["prior_path"])
    rest = load_reference(previous["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    reference = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(reference["local_geometry"])
    assembly = load_reference(baseline["shoulder_assembly"])
    head = "FJ1500"
    specification = next(r for r in load_reference(reference["kinematics"])["models"] if r["head"] == head)
    route = next(r for r in load_reference(baseline["input_route"])["routes"] if r["head"] == head)
    sleeve_key = f"link:{head}/muscle_uncrimped_sleeve"
    sleeve = load_reference(previous["sleeve"])
    if sleeve["original"] != assembly["sources"][sleeve_key]:
        raise ValueError("Unexpected sleeve source")
    sources = {**assembly["sources"], sleeve_key: sleeve["selected"]["source"]}
    parts = {k: load_component(src) for k, src in sources.items()
             if not is_compliant(k) and not k.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
             and not k.endswith("/tendon_wire")}
    shifted_parts = {**parts, **{k: load_component(src) for k, src in selected["meshes"].items()}}
    owners = assembly["intended_body_owners"]
    centre = np.array(geometry["gh_centre_source_mm"])
    _, _, longitudinal = np.array(geometry["axes_source"])
    cap_bias = np.array(rest["rest_offset_source_mm"])
    adjusted, ports = rest_geometry(specification, route["ports"], cap_bias)
    origin, stop = (np.array(route["ports"][k][face]) for k, face in
                    (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm")))
    direction = Rotation.from_rotvec(A.D.unit(origin-stop)*np.deg2rad(baseline["selected"]["coupling_roll_deg"])).apply(baseline["direction_source"])
    pin_a, pin_b = (np.array(specification[f"neutral_{k}_pin_source_mm"]) for k in ("cap", "bone"))
    axis_a, axis_b = (np.array(specification[f"{k}_pin_axis_source"]) for k in ("cap", "bone"))
    points, faces = mesh_geom.manifold_mesh_to_arrays(load_component(rest["candidate"]).to_mesh64())
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    span = np.linalg.norm(stop-origin)
    axial = (stop-origin)/span
    lower, upper = rest["original_mesh_deflection_fraction_limits"]
    u = np.clip(((points-origin)@axial/span-lower)/(upper-lower), 0., 1.)
    weights = u**2*(3.-2.*u)
    derivative = (6.*u*(1.-u)/((upper-lower)*span))[:,None]*axial
    transition = (u>0.) & (u<1.)
    mf = A.I._manifold()

    integration = DEST.parent / "integration_review/s1_series_elbow"
    _, extended, _, _, _, _, _, retention = load_candidate(integration)
    model_path = ROOT / extended["model"]["file"]
    if provenance(model_path) != extended["model"]:
        raise ValueError("Combined dynamics model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    pivot = model.joint("GH")
    np.testing.assert_allclose(model.body("humerus").pos*1000., centre, atol=1e-5, rtol=0.)
    trace = checked_trace(extended["trajectory"])
    address = pivot.qposadr.item()
    rotations = Rotation.from_quat(np.array([[float(row[f"qpos_{i}"]) for i in range(address,address+4)]
                                             for row in trace]), scalar_first=True)
    frames = sorted(set(range(0, len(trace), 20)) | {len(trace)-1}
                    | {int(np.argmax([float(row[key]) for row in trace]))
                       for key in ("GH_delta_deg", "flexion_projection_deg", "elbow_deg", "slide_mm")})
    scenarios = []
    for label, shift, correction, solids in (("original", np.zeros(3), 0., parts),
                                            ("offset_16mm", np.array(selected["parameters"]["shift_source_mm"]),
                                             selected["length_adjustment_mm"], shifted_parts)):
        test_model = deepcopy(adjusted)
        for key in ("neutral_bone_pin_source_mm", "neutral_crimp_outer_exit_B_source_mm"):
            test_model[key] = (np.array(test_model[key])+shift).tolist()
        strand = np.array(test_model["active_strand_centreline_source_mm"])
        strand[-2:] += shift
        test_model["active_strand_centreline_source_mm"] = strand.tolist()
        a = np.array(test_model["neutral_crimp_outer_exit_A_source_mm"])
        b = centre+rotations.apply(np.array(test_model["neutral_crimp_outer_exit_B_source_mm"])-centre)
        chord = np.linalg.norm(a-b, axis=1)
        arms = np.einsum("ij,ij->i", np.cross(b-centre, (a-b)/chord[:,None]), rotations.apply(longitudinal))
        rows = []
        for i in frames:
            gh = rotations[i]
            arm_pose = rigid_transform(gh, centre, centre)
            obstacles = bounded({k: solid.transform(arm_pose) if owners[k] == "humerus_ideal_attachment" else solid
                                 for k, solid in solids.items() if k.startswith(("bone:", "GH_cap_"))})
            states = solve_branch(test_model, ports, centre, gh, bias["stiffness"], [0.,5.,10.], retain_leads=True,
                                  installed_length_offset=correction, lateral_coupling=(direction, baseline["selected"]["gain"]),
                                  strand_bending_Nmm2=10.)
            for state in states:
                pin, hits = collision_state(solids, owners, state, cap_bias, gh, centre, pin_a, pin_b+shift, axis_a, axis_b)
                torque = float(gh.apply(longitudinal)@state["GH_torque_source_Nm"])
                delta = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"])).as_matrix()-np.eye(3)
                displacement = (points-pin_a-cap_bias)@delta.T+state["cap_translation_source_mm"]
                gradient = np.eye(3)+weights[transition,None,None]*delta+displacement[transition,:,None]*derivative[transition,None,:]
                determinant = float(np.linalg.det(gradient).min())
                soft_hits = None
                if determinant>0.:
                    mesh = mf.Manifold(mf.Mesh64(np.ascontiguousarray(points+weights[:,None]*displacement), faces))
                    if mesh.status() == mf.Error.NoError:
                        soft_hits = overlaps(mesh, obstacles)
                rows.append(dict(frame=i, time_s=float(trace[i]["time_s"]), cap_pin_source_mm=pin.tolist(),
                                 terminal_hits=hits, longitudinal_torque_Nm=torque, min_TPU_map_determinant=determinant,
                                 interpolated_TPU_hits_mm3=soft_hits, **state))
        high = [row for row in rows if row["commanded_muscle_tension_N"] == 10.]
        summary = dict(chord_range_mm=[float(chord.min()),float(chord.max())], chord_span_mm=float(np.ptp(chord)),
                       moment_arm_range_mm=[float(arms.min()),float(arms.max())],
                       ten_newton_torque_range_Nm=[min(q["longitudinal_torque_Nm"] for q in high),max(q["longitudinal_torque_Nm"] for q in high)],
                       ten_newton_tension_range_N=[min(q["short_tendon_tension_N"] for q in high),max(q["short_tendon_tension_N"] for q in high)],
                       max_cap_shift_mm=max(float(np.linalg.norm(q["cap_translation_source_mm"])) for q in rows),
                       failing_terminal_states=sum(bool(q["terminal_hits"]) for q in rows),
                       invalid_TPU_states=sum(q["interpolated_TPU_hits_mm3"] is None for q in rows),
                       TPU_bone_cover_contact_states=sum(bool(q["interpolated_TPU_hits_mm3"]) for q in rows),
                       peak_TPU_overlap_mm3=max(max((q["interpolated_TPU_hits_mm3"] or {}).values(), default=0.) for q in rows))
        scenarios.append(dict(label=label, summary=summary, equilibrium=rows))
        print(label, summary, flush=True)
    save(folder / "FJ1500_combined_trace_screen.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_infraspinatus_combined_trace.py"),
        exit_revision=provenance(exit_path), model=extended["model"], trace=extended["trajectory"], retention=retention,
        chord_sample_count=len(trace), equilibrium_frames=frames, scenarios=scenarios,
        scope="Conditional local FJ1500 equilibrium on the recorded shoulder poses."
              " Moment arm projected on the rotating humeral longitudinal axis."
              " The 0/5/10 N loads are probes; the recorded combined trajectory was not driven by FJ1500."
              " Local terminal checks use shoulder geometry only. TPU is checked against bones and GH covers"
              " with a prescribed mesh interpolation. Wire contact and live dynamics are outside scope."))


if __name__ == "__main__":
    main()
