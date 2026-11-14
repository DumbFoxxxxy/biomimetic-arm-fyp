"""Measure the brief near-zero extensor–humerus contact on the actual drive trace."""

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, rod, save
from channel_first_bone_contact import BoneContactSurface
from refine_elbow_short_guide_relief import indexed_solid, inverse_map
from study_s1_local_force_geometry import load_reference
from study_s1_series_elbow_release_timing import run


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    source_report = folder / "interface_relief_elbow_muscle_geometry.json"
    report = read(source_report)
    extended = load_reference(report["extended_drive"])
    budget = load_reference(extended["budget"])
    timing = load_reference(budget["timing"])
    baseline = load_reference(timing["baseline"])
    assembly = load_reference(budget["assembly"])
    active = load_reference(timing["input_cycles"])
    start = load_reference(active["start"])
    shoulder = load_reference(baseline["shoulder_endpoints"])
    release = timing["cases"][timing["selected"]]
    model_file = ROOT / extended["model"]["file"]
    if provenance(model_file) != extended["model"]:
        raise ValueError("Drive model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    dense = run(model, np.array(start["candidate"]["quaternion_wxyz"]),
                [model.actuator(name).id for name in active["shoulder_channels"]],
                np.array(active["preload_N"]), .6 * np.array(active["full_flexion_increment_N"]),
                release["shoulder_release_start_s"], release["elbow_release_start_s"],
                np.array(shoulder["axes_source"][0]), sample_period_s=.0002)
    parent = checked_trace(report["trace"])
    period = dense[1, 0]
    parent_indices = np.rint(np.array([float(row["time_s"]) for row in parent]) / period).astype(int)
    np.testing.assert_allclose(dense[parent_indices, -model.nq:],
                               [[float(row[f"qpos_{i}"]) for i in range(model.nq)] for row in parent], atol=1e-10, rtol=0.)

    elbow = load_reference(load_reference(report["passive"])["assembly"])["muscles"]["extensor"]
    definition = load_reference(elbow["source_definition"])
    path = np.array(elbow["wire_path_source_mm"])[1:-1]
    axis = np.array(definition["basis_columns_display"])[:, 2]
    gradient = axis / ((path[-1] - path[0]) @ axis)
    mf = A.I._manifold()

    def mesh_data(solid):
        points, triangles = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
        return points, weigh(points, path[0], gradient), np.ascontiguousarray(triangles, dtype=np.uint64)

    source = report["mesh_sources"]["extensor"]
    _, vertices, faces = indexed_solid(source)
    lattice = vertices, weigh(vertices, path[0], gradient), faces
    bone = load_component(assembly["sources"]["bone:humerus"])
    short = read(DEST.parent / "elbow_local_muscle_load/elbow_slider_short_guides.json")
    guide_design = load_reference(short["short_guide_design"])
    guides = guide_design["actual_short_tube_guides"]
    np.testing.assert_allclose(np.array([row["centre_display_mm"] for row in guides]) + [-140., -56., 1100.],
                               path, atol=1e-7, rtol=0.)
    guide_meshes = []
    for centre, spec in zip(path, guides, strict=True):
        half = np.array(spec["axis_display"]) * spec["length_mm"] / 2.
        guide_meshes.append(mesh_data(rod(centre - half, centre + half, spec["outer_radius_mm"])))

    def pose(index):
        data.qpos[:] = dense[index, -model.nq:]
        mujoco.mj_kinematics(model, data)
        placement = transforms(model, data, neutral)
        h = placement["humerus"]
        d = placement["FJ1477_distal_pivot"]
        carriage = placement["FJ1477_carriage"]
        r, t = h[:, :3].T, h[:, 3]
        return r @ d[:, :3], r @ (d[:, 3] - t), r @ (carriage[:, 3] - t)

    def placed(mesh, affine):
        points, weights, triangles = mesh
        current = displacement_map(points, weights, *affine)
        solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(current), triangles))
        if solid.status() != mf.Error.NoError:
            raise ValueError("Invalid posed extensor geometry")
        return solid

    seeds = {int(parent_indices[i]): float(row["muscles"]["extensor"]["rigid_TPU_contacts"]["bone:humerus"])
             for i, row in enumerate(report["records"])
             if "bone:humerus" in row["muscles"]["extensor"].get("rigid_TPU_contacts", {})}
    if not seeds:
        raise ValueError("Expected reported return contact")
    half_window = round(.020 / period)
    begin, end = min(seeds) - half_window, max(seeds) + half_window
    records = {}
    peak_contact = None

    def check(indices):
        nonlocal peak_contact
        for index in sorted(set(indices) - records.keys()):
            contact = placed(lattice, pose(index)) ^ bone
            volume = float(contact.volume())
            if volume > max(records.values(), default=0.):
                peak_contact = contact
            records[index] = volume

    check(set(range(begin, end + 1, 5)) | seeds.keys() | {end})
    while records[begin] > .001 and begin > 0:
        left = max(0, begin - half_window)
        check(set(range(left, begin, 5)) | {left})
        begin = left
    while records[end] > .001 and end < len(dense) - 1:
        right = min(len(dense) - 1, end + half_window)
        check(set(range(end + 5, right + 1, 5)) | {right})
        end = right
    initial_hits = [i for i, volume in records.items() if volume > .001]
    print("Coarse contact", len(initial_hits), "window", dense[begin, 0], dense[end, 0],
          "peak", max(records.values()), flush=True)
    check(i for hit in initial_hits for i in range(max(0, hit - 5), min(len(dense), hit + 6)))
    for index, value in seeds.items():
        np.testing.assert_allclose(records[index], value, atol=1e-7, rtol=0.)
    if records[begin] > .001 or records[end] > .001:
        raise ValueError("Contact reaches the scan boundary")
    index = max(records, key=records.get)
    affine = pose(index)
    contact = peak_contact
    surface = BoneContactSurface(bone)
    # Distance to the humeral surface from the overlap boundary. This is a
    # sampled geometric depth and has no constitutive material interpretation.
    depths = []
    for step in (.20, .10):
        refined = contact.refine_to_length(step)
        vertices, triangles = mesh_geom.manifold_mesh_to_arrays(refined.to_mesh64())
        points = np.vstack([vertices, vertices[triangles].mean(axis=1)])
        distance, _ = surface.query(points)
        deepest = int(distance.argmin())
        depths.append(dict(subdivision_target_edge_mm=step, point_count=len(points),
                           sampled_depth_mm=float(-distance[deepest]),
                           deepest_humerus_frame_mm=points[deepest].tolist()))
    raw_points, _ = mesh_geom.manifold_mesh_to_arrays(contact.to_mesh64())
    source_points = inverse_map(raw_points, *affine, path[0], gradient)
    bounds = np.array([source_points.min(axis=0), source_points.max(axis=0)])
    station = (source_points - path[0]) @ axis
    wire = A.D.outer_envelope(displacement_map(path, weigh(path, path[0], gradient), *affine), .30)
    guide_gaps = [float(contact.min_gap(placed(mesh, affine), 40.)) for mesh in guide_meshes]
    guide_index = int(np.argmin(guide_gaps))
    results = dict(
        generator=provenance(ROOT / "tools/measure_s1_series_elbow_extensor_contact.py"),
        run_source=provenance(ROOT / "tools/study_s1_series_elbow_release_timing.py"),
        assembly_check=provenance(source_report), extensor=source, humerus=assembly["sources"]["bone:humerus"],
        guide_definition=short["short_guide_design"], integrator_timestep_s=float(model.opt.timestep),
        scan_window_s=[float(dense[begin, 0]), float(dense[end, 0])], coarse_period_s=period * 5,
        refined_period_s=period, recorded_pose_match_count=len(parent), overlap_threshold_mm3=.001,
        samples=[dict(time_s=float(dense[i, 0]), elbow_deg=float(dense[i, 3]), volume_mm3=records[i]) for i in sorted(records)],
        peak=dict(time_s=float(dense[index, 0]), elbow_deg=float(dense[index, 3]),
                  volume_mm3=records[index], depths=depths, source_bounds_mm=bounds.tolist(),
                  axial_range_from_first_guide_mm=[float(station.min()), float(station.max())],
                  nominal_strut_diameter_mm=guide_design["parameters"]["strut_diameter_mm"],
                  wire_gap_mm_capped_40=float(contact.min_gap(wire, 40.)),
                  nearest_guide_index_zero_based=guide_index, guide_outer_envelope_gap_mm_capped_40=guide_gaps[guide_index]),
        scope="Actual dynamic states at smaller output intervals; TPU geometry still uses cap displacement interpolation."
              " Depth is sampled on the overlap boundary; no material compression, contact force or physical clearance claim.")
    save(folder / "FJ1477_return_humerus_contact.json", results)
    print(results["peak"], "frames", len(records), flush=True)


if __name__ == "__main__":
    main()
