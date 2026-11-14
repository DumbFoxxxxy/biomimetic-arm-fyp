"""Test distributed FJ1500 deformation against the recorded loaded cap poses."""

import mujoco
import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def blend(stations, lower, upper, ramp):
    x = np.clip((stations-lower)/(upper-lower), 0., 1.)
    if ramp == 0.:
        return x*x*(3.-2.*x), 6.*x*(1.-x)/(upper-lower)
    # Constant slope through the belly, with quadratic end ramps.
    start = x < ramp
    end = x > 1.-ramp
    integral = np.where(start, x*x/(2.*ramp), np.where(end, 1.-ramp-(1.-x)**2/(2.*ramp), x-ramp/2.))
    slope = np.minimum(np.minimum(x/ramp, (1.-x)/ramp), 1.)
    return integral/(1.-ramp), slope/((upper-lower)*(1.-ramp))


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    parent_file = folder / "FJ1500_combined_trace_screen.json"
    parent = read(parent_file)
    revision = load_reference(parent["exit_revision"])
    outward = load_reference(revision["baseline"])
    coupling = load_reference(load_reference(outward["parent"])["coupling"])
    old = load_reference(coupling["prior_path"])
    rest = load_reference(old["rest_geometry"])
    prior = load_reference(load_reference(rest["bias_screen"])["reference_equilibrium"])
    geometry = load_reference(prior["local_geometry"])
    assembly = load_reference(outward["shoulder_assembly"])
    spec = next(r for r in load_reference(prior["kinematics"])["models"] if r["head"] == "FJ1500")
    route = next(r for r in load_reference(outward["input_route"])["routes"] if r["head"] == "FJ1500")
    origin, stop = (np.array(route["ports"][k][face]) for k, face in
                    (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm")))
    axis = A.D.unit(stop-origin)
    span = float(np.linalg.norm(stop-origin))
    pin = np.array(spec["neutral_cap_pin_source_mm"])+rest["rest_offset_source_mm"]
    centre = np.array(geometry["gh_centre_source_mm"])
    wire_path = np.array(rest["rest_path_source_mm"])
    solid = load_component(rest["candidate"])
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    coords = (vertices-origin)@axis
    mf = A.I._manifold()
    # Protect full axial extent of the actual end pieces where their padded
    # envelopes touch TPU. Other hard pieces are checked in the motion screen.
    sources = assembly["sources"]
    distal_keys = [k for k, owner in assembly["intended_body_owners"].items() if owner == "free_terminal:FJ1500"]
    proximal_keys = [k for k in sources if k == "cap:FJ1500/origin" or k.startswith("origin:FJ1500/")]
    guards = {}
    for end, keys in (("proximal", proximal_keys), ("distal", distal_keys)):
        values = []
        for key in keys:
            part = load_component(sources[key])
            if end == "distal":
                part = part.translate(rest["rest_offset_source_mm"])
            intersection = part.minkowski_sum(mf.Manifold.cube([1.]*3, True)) ^ solid
            if intersection.volume() > .001:
                points, _ = mesh_geom.manifold_mesh_to_arrays(intersection.to_mesh64())
                stations = (points-origin)@axis
                values.append(dict(part=key, interval_mm=[float(stations.min()), float(stations.max())]))
        if not values:
            raise ValueError(f"No {end} pad interface found")
        guards[end] = values
    lower = max(row["interval_mm"][1] for row in guards["proximal"])
    upper = min(row["interval_mm"][0] for row in guards["distal"])
    if upper <= lower:
        raise ValueError("Protected interfaces leave no deformation span")
    rows = next(case["equilibrium"] for case in parent["scenarios"] if case["label"] == "offset_16mm")
    trace = checked_trace(parent["trace"])
    model_path = ROOT/parent["model"]["file"]
    if provenance(model_path) != parent["model"]:
        raise ValueError("Recorded dynamics model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    address = model.joint("GH").qposadr.item()
    bones = {k: load_component(v) for k,v in sources.items() if k.startswith(("bone:", "GH_cap_"))}
    obstacles = {}
    for frame in parent["equilibrium_frames"]:
        gh = Rotation.from_quat([float(trace[frame][f"qpos_{i}"]) for i in range(address,address+4)], scalar_first=True)
        arm = rigid_transform(gh, centre, centre)
        obstacles[frame] = bounded({k: part.transform(arm) if assembly["intended_body_owners"][k] == "humerus_ideal_attachment"
                                    else part for k, part in bones.items()})
    nominal_lo, nominal_hi = np.array(rest["original_mesh_deflection_fraction_limits"])*span
    scenarios = [("previous", nominal_lo, nominal_hi, 0.), ("extended_cubic", lower, upper, 0.),
                 ("extended_ramps_20pct", lower, upper, .20), ("extended_ramps_10pct", lower, upper, .10)]
    results = []
    for name, lo, hi, ramp in scenarios:
        w, derivative = blend(coords, lo, hi, ramp)
        dw = derivative[:,None]*axis
        selected = derivative>0.
        wire_weight, _ = blend((wire_path-origin)@axis, lo, hi, ramp)
        checks = []
        for row in rows:
            delta = Rotation.from_rotvec(np.deg2rad(row["cap_rotation_vector_deg"])).as_matrix()-np.eye(3)
            shift = np.array(row["cap_translation_source_mm"])
            displacement = (vertices-pin)@delta.T+shift
            jacobian = np.eye(3)+w[selected,None,None]*delta+displacement[selected,:,None]*dw[selected,None,:]
            determinant = float(np.linalg.det(jacobian).min())
            if name == "previous":
                np.testing.assert_allclose(determinant, row["min_TPU_map_determinant"], rtol=0., atol=1e-9)
            record = dict(frame=row["frame"], input_N=row["commanded_muscle_tension_N"], minimum_sampled_determinant=determinant)
            if determinant>0.:
                posed = mf.Manifold(mf.Mesh64(np.ascontiguousarray(vertices+w[:,None]*displacement), faces))
                record["mesh_status"] = str(posed.status())
                if posed.status() == mf.Error.NoError:
                    wire_points = wire_path+wire_weight[:,None]*((wire_path-pin)@delta.T+shift)
                    wire = A.D.outer_envelope(wire_points, .30)
                    record.update(bone_and_cover_hits_mm3=overlaps(posed, obstacles[row["frame"]]),
                                  nominal_internal_wire_overlap_mm3=float((posed^wire).volume()))
            checks.append(record)
        results.append(dict(name=name, blend_range_mm=[lo,hi], ramp_fraction=ramp, cases=checks))
        print(name, 'range', lo, hi, 'invalid', sum('bone_and_cover_hits_mm3' not in r for r in checks),
              'bone hits',sum(bool(r.get('bone_and_cover_hits_mm3')) for r in checks),
              'wire hits',sum(r.get('nominal_internal_wire_overlap_mm3',0.)>.001 for r in checks),flush=True)
    save(folder/"FJ1500_distributed_deformation_screen.json", dict(
        generator=provenance(ROOT/"tools/audit_s1_infraspinatus_distributed_deformation.py"), combined=provenance(parent_file),
        rest_geometry=old["rest_geometry"], pad_guard_intervals_mm=guards, muscle_axis_span_mm=span, results=results,
        scope="Prescribed interpolation of actual rest-shaped lattice and nominal inner route on sampled local equilibria."
              " Protected pad bounds derived from actual geometry with a 0.5 mm cube dilation."
              " No constitutive equilibrium, stiffness identification, actual inner cable sliding or contact forces."))


if __name__ == "__main__":
    main()
