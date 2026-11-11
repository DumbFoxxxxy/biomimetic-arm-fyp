"""Evaluate FJ1477 end displacement, channel strain demand and series travel."""

import csv
from pathlib import Path

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from simulate_elbow_antagonistic_load import JOINTS


def polyline_lengths(points):
    return np.linalg.norm(np.diff(points, axis=0), axis=1)


def turn_radius(points):
    left, right = np.diff(points, axis=0)[:-1], np.diff(points, axis=0)[1:]
    product = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1) * np.linalg.norm(left + right, axis=1)
    twice_cross = 2. * np.linalg.norm(np.cross(left, right), axis=1)
    curved = twice_cross > 1e-8
    return float((product[curved] / twice_cross[curved]).min()) if curved.any() else None


def map_points(points, rotation, translation, anchor, station_gradient):
    weights = np.clip((points - anchor) @ station_gradient, 0., 1.)
    displacement = points @ (rotation - np.eye(3)).T + translation
    return points + weights[:, None] * displacement


def channel_record(path, posed, inner_radius, wire_radius):
    neutral = polyline_lengths(path)
    lengths = polyline_lengths(posed)
    ratios = lengths / neutral
    min_radius = inner_radius / np.sqrt(ratios.max())
    chord = np.linalg.norm(posed[-1] - posed[0])
    return dict(
        bore_span_mm=float(chord),
        straight_lower_bound_extension_mm=float(chord - neutral.sum()),
        interpolated_path_length_mm=float(lengths.sum()),
        interpolated_path_change_mm=float(lengths.sum() - neutral.sum()),
        max_segment_stretch_ratio=float(ratios.max()),
        min_segment_stretch_ratio=float(ratios.min()),
        max_segment_stretch_index=int(ratios.argmax()),
        interpolated_polyline_min_three_point_radius_mm=turn_radius(posed),
        ideal_round_lumen_min_radius_mm=float(min_radius),
        ideal_round_lumen_min_radial_clearance_mm=float(min_radius - wire_radius),
    )


def mesh_record(vertices, faces, edges, edge_lengths, rotation, translation, anchor,
                gradient, rigid_bones, forearm_affine):
    posed = map_points(vertices, rotation, translation, anchor, gradient)
    ratios = np.linalg.norm(posed[edges[:, 1]] - posed[edges[:, 0]], axis=1) / edge_lengths
    station = (vertices - anchor) @ gradient
    interior = (station > 0.) & (station < 1.)
    displacement = vertices[interior] @ (rotation - np.eye(3)).T + translation
    jacobians = (np.eye(3) + station[interior, None, None] * (rotation - np.eye(3))
                 + displacement[:, :, None] * gradient[None, None, :])
    determinants = np.linalg.det(jacobians)
    mf = A.I._manifold()
    solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), faces))
    valid = solid.status() == mf.Error.NoError
    return dict(
        surface_edge_length_ratio_quantiles=np.quantile(ratios, [0., .05, .5, .95, 1.]).tolist(),
        min_interpolation_jacobian=float(determinants.min()),
        nonpositive_jacobians=int(np.count_nonzero(determinants <= 0.)),
        manifold_status=str(solid.status()),
        interpolated_TPU_overlap_mm3=(
            {"humerus": float(A.I.intersection_volume(solid, rigid_bones["humerus"])),
             "forearm": float(A.I.intersection_volume(solid, rigid_bones["forearm"].transform(forearm_affine)))}
            if valid and determinants.min() > 0. else None),
    )


def short_guide_record(vertices, faces, centres, axes, half_lengths, posed_path,
                       rotation, translation, anchor, gradient, wire_radius):
    entrances = map_points(centres - axes * half_lengths, rotation, translation, anchor, gradient)
    exits = map_points(centres + axes * half_lengths, rotation, translation, anchor, gradient)
    gaps = np.linalg.norm(entrances[1:] - exits[:-1], axis=1)
    mf = A.I._manifold()
    posed = map_points(vertices, rotation, translation, anchor, gradient)
    solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), faces))
    valid = solid.status() == mf.Error.NoError
    return dict(
        source_manifold_status_after_mapping=str(solid.status()),
        max_gap_between_guide_mouths_mm=float(gaps.max()),
        nominal_wire_overlap_TPU_mm3=float(A.I.intersection_volume(
            solid, A.D.outer_envelope(posed_path, wire_radius))) if valid else None,
    )


def main():
    output = DEST.parent / "elbow_local_muscle_load"
    response_file = output / "elbow_attachment_response.json"
    response = read(response_file)
    model_path = ROOT / response["local_model"]["file"]
    passive_path = DEST / "elbow_passive_model/model_check.json"
    passive = read(passive_path)
    assembly_path = ROOT / passive["assembly"]["file"]
    assembly = read(assembly_path)
    case = next(row for row in response["cases"]
                if row["arrangement"] == "original" and row["stiffness_N_per_m"] == [200., 50.])
    trace_path = ROOT / case["trace"]["file"]
    for path, expected in ((model_path, response["local_model"]), (trace_path, case["trace"]),
                           (assembly_path, passive["assembly"])):
        if provenance(path) != expected:
            raise RuntimeError(f"Input fingerprint changed: {path.name}")

    tpu = assembly["display_only_TPU_sources"]["tpu:elbow_extensor"]
    design_path = (ROOT / tpu["file"]).with_name("design.json")
    design = read(design_path)
    parent_path = DEST.parent / "anatomical_concept/direct_muscle_interfaces/extensor_reanchor_extended/muscle/FJ1477"
    short_mesh_path = parent_path / "TPU_muscle.obj"
    short_design_path = parent_path / "design.json"
    if provenance(short_mesh_path)["sha256"] != design["parent_OBJ_sha256"]:
        raise RuntimeError("Short-guide source differs from the continuous-channel parent")
    short_design = read(short_design_path)
    short_vertices, short_faces = mesh_geom.load_obj(str(short_mesh_path))
    short_vertices += np.asarray(tpu["translation_source_mm"])
    short_faces = np.ascontiguousarray(short_faces, dtype=np.uint64)
    repair_path = DEST.parent / "integration_review/s1_elbow_shared_humerus/elbow_tpu_readback_repair.json"
    repair = next(row for row in read(repair_path)["rows"] if row["key"] == "tpu:elbow_extensor")
    if repair["original_source"]["sha256"] != tpu["sha256"] or provenance(ROOT / tpu["file"])["sha256"] != tpu["sha256"]:
        raise RuntimeError("Extensor repair does not match original muscle geometry")
    mesh = load_component(repair["candidate_source"])
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(mesh.to_mesh64())
    faces = faces.astype(np.uint64)
    edges = np.unique(np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1), axis=0)
    edge_lengths = np.linalg.norm(vertices[edges[:, 1]] - vertices[edges[:, 0]], axis=1)
    edge_mask = edge_lengths >= .2
    edges, edge_lengths = edges[edge_mask], edge_lengths[edge_mask]
    muscle = assembly["muscles"]["extensor"]
    pads = np.asarray(muscle["pad_centres_source_mm"])
    path = np.asarray(muscle["wire_path_source_mm"][1:-1])
    guides = short_design["actual_short_tube_guides"]
    guide_centres = np.asarray([row["centre_display_mm"] for row in guides]) + tpu["translation_source_mm"]
    guide_axes = np.asarray([row["axis_display"] for row in guides])
    guide_half_lengths = .5 * np.asarray([row["length_mm"] for row in guides])[:, None]
    np.testing.assert_allclose(guide_centres, path, atol=1e-8)
    axis = np.asarray(design["basis_columns_display"])[:, 2]
    pad_half_thickness = 2.
    active_span = (pads[1] - pads[0]) @ axis - 2. * pad_half_thickness
    # The bore seats move with their caps. Hold the entire end regions outside
    # the two bore planes rigid in this imposed displacement map.
    anchor = path[0]
    blend_span = (path[-1] - path[0]) @ axis
    if not 0. < blend_span <= active_span:
        raise RuntimeError("Unexpected extensor bore planes relative to the two pad faces")
    gradient = axis / blend_span

    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    cap = data.body("FJ1477_distal_pivot")
    forearm = data.body("forearm")
    neutral_cap, neutral_forearm = cap.xpos.copy() * 1000., forearm.xpos.copy() * 1000.
    with trace_path.open(newline="", encoding="utf-8") as stream:
        trace = list(csv.DictReader(stream))
    peak_index = max(range(len(trace)), key=lambda index: float(trace[index]["elbow_deg"]))
    samples = {0, peak_index, len(trace) - 1}
    for angle in (30., 60., 90.):
        for window in (range(peak_index + 1), range(peak_index, len(trace))):
            samples.add(min(window, key=lambda index: abs(float(trace[index]["elbow_deg"]) - angle)))
    stop_geometry_path = DEST / "elbow_extension_stop/assembly.json"
    stop_geometry = read(stop_geometry_path)
    bones = {name: load_component(stop_geometry["sources"][f"bone:{name}"]) for name in ("humerus", "forearm")}
    lumen_radius = float(design["tube_inner_radius_mm"])
    wall_radius = float(design["tube_outer_radius_mm"])
    wire_radius = .3
    channel_rows, mesh_rows = [], []
    for index, row in enumerate(trace):
        pose = [float(row[f"{name}_deg"]) for name in JOINTS]
        data.qpos[addresses] = np.deg2rad(pose)
        mujoco.mj_kinematics(model, data)
        rotation = cap.xmat.reshape(3, 3).copy()
        translation = cap.xpos * 1000. - rotation @ neutral_cap
        posed = map_points(path, rotation, translation, anchor, gradient)
        record = dict(time_s=float(row["time_s"]), elbow_deg=pose[0], cap_deg=pose[2],
                      **channel_record(path, posed, lumen_radius, wire_radius))
        np.testing.assert_allclose(record["bore_span_mm"] - np.linalg.norm(path[-1] - path[0]),
                                   float(row["FJ1477_length_change_mm"]), atol=1e-6)
        channel_rows.append(record)
        if index in samples:
            forearm_r = forearm.xmat.reshape(3, 3)
            forearm_affine = np.c_[forearm_r, forearm.xpos * 1000. - forearm_r @ neutral_forearm]
            mesh_rows.append(dict(
                **record, **mesh_record(vertices, faces, edges, edge_lengths, rotation, translation,
                                        anchor, gradient, bones, forearm_affine),
                short_guide_screen=short_guide_record(short_vertices, short_faces, guide_centres,
                                                     guide_axes, guide_half_lengths, posed, rotation,
                                                     translation, anchor, gradient, wire_radius),
            ))

    peak = channel_rows[peak_index]
    total = float(trace[peak_index]["FJ1477_length_change_mm"])
    force = 50. * total * .001
    tube_length = polyline_lengths(path).sum()
    area = np.pi * (wall_radius**2 - lumen_radius**2)
    # Uniform straight-tube small-strain axial stiffness before lattice contributions.
    tube_stiffness = {str(e): float(1000. * e * area / tube_length) for e in (5., 15., 30.)}
    splits = []
    for series_travel in (0., 25., 50., 60.):
        lattice_travel = total - series_travel
        splits.append(dict(lattice_travel_mm=lattice_travel, series_travel_mm=series_travel,
                           reference_lattice_span_ratio=lattice_travel / active_span,
                           lattice_secant_stiffness_N_per_m=1000. * force / lattice_travel,
                           series_secant_stiffness_N_per_m=1000. * force / series_travel if series_travel else None))
    report = dict(
        generator=provenance(Path(__file__).resolve()), response=provenance(response_file),
        trace=provenance(trace_path), extensor_mesh=provenance(ROOT / tpu["file"]),
        channel_definition=provenance(design_path), rigid_bones=provenance(stop_geometry_path),
        mesh_repair=provenance(repair_path), extensor_collision_mesh=repair["candidate_source"],
        short_guide_parent_mesh=provenance(short_mesh_path), short_guide_definition=provenance(short_design_path),
        nominal_bore_span_mm=float(np.linalg.norm(path[-1] - path[0])),
        nominal_internal_polyline_length_mm=float(tube_length),
        nominal_active_span_mm=float(active_span), bore_plane_blend_span_mm=float(blend_span),
        lumen_radius_mm=lumen_radius,
        wall_outer_radius_mm=wall_radius, assumed_wire_radius_mm=wire_radius, peak_pose=peak,
        tube_tangent_stiffness_N_per_m_by_assumed_modulus_MPa=tube_stiffness,
        tube_only_modulus_for_50_N_per_m_MPa=float(50. * tube_length / (1000. * area)),
        series_design_budget_at_peak=dict(total_travel_mm=total, assumed_passive_spring_force_N=force, splits=splits),
        sampled_mesh_poses=mesh_rows, channel_poses=channel_rows,
        scope="Kinematic demand from original-anchor 25 N and 50 N/m extensor response. Mesh and "
        "channel are prescribed by interpolating the installed end-cap transforms. Edge ratios and "
        "overlap measure that deformation map, not actual TPU strain or contact. Lumen estimate assumes "
        "uniform incompressible round tube stretch and ignores flattening, bending and wire motion. "
        "Tube stiffness is an axial linear estimate with assumed moduli; series budget is not hardware.")
    save(output / "elbow_extensor_long_stroke.json", report)
    print("peak", peak, flush=True)
    print("tube stiffness", tube_stiffness, "series force", force, flush=True)
    for row in mesh_rows:
        print(round(row["elbow_deg"], 2), row["min_interpolation_jacobian"],
              row["interpolated_TPU_overlap_mm3"], row["short_guide_screen"], flush=True)


if __name__ == "__main__":
    main()
