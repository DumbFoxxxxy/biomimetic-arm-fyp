"""Check FJ1477 short guides and the local wire against the free-slider trajectories."""

import csv
from pathlib import Path

import mujoco
import numpy as np

from audit_elbow_cap_travel import load_geometry
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_spring_mount import spring_envelope
from build_s1_teres_origin_mount import A, export, load_component, mesh_geom, rod, save
from simulate_elbow_antagonistic_load import JOINTS
from study_elbow_extensor_long_stroke import polyline_lengths, turn_radius
from study_elbow_series_package import ADAPTER, PREFIX, bounded, overlaps, pose_solids, union


def displacement_map(points, weight, rotation, distal_translation, proximal_shift):
    return (points + (1. - weight[:, None]) * proximal_shift
            + weight[:, None] * (points @ (rotation - np.eye(3)).T + distal_translation))


def weigh(points, first, gradient):
    return np.clip((points - first) @ gradient, 0., 1.)


def mesh_edges(faces):
    return np.unique(np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]],
                                             faces[:, [2, 0]])), axis=1), axis=0)


def sample_indices(rows):
    angles = [float(row["elbow_deg"]) for row in rows]
    peak = int(np.argmax(angles))
    growth = [float(row["FJ1477_span_change_mm"]) for row in rows]
    sample = {0, peak, len(rows) - 1, int(np.argmax([float(row["carriage_mm"]) for row in rows])),
              int(np.argmin(growth)), int(np.argmax(growth))}
    for angle in (5., 10., 15., 20., 25., 30., 60., 90.):
        for window in (range(peak + 1), range(peak, len(rows))):
            sample.add(min(window, key=lambda index: abs(angles[index] - angle)))
    return sample


def relief_meshes(vertices, faces, origin, centres, half_axes, stations, rail_base, rail_axis, travel):
    mf = A.I._manifold()
    initial = mf.Manifold(mf.Mesh64(np.ascontiguousarray(vertices), faces))
    guide_keepout = union(rod(centre - half_axis, centre + half_axis, row["outer_radius_mm"] + .25)
                         for centre, half_axis, row in zip(centres, half_axes, stations, strict=True))
    protected = initial ^ guide_keepout
    options = []
    for radius in (2.5, 2.8, 3.1, 3.5, 4.5, 5.5):
        cutter = rod(rail_base - 28. * rail_axis, rail_base + (travel - 5.) * rail_axis, radius)
        mesh = (initial - cutter).as_original().simplify(1e-6)
        removed = float((initial - mesh).volume())
        protected_loss = float((protected - mesh).volume())
        components = sum(piece.volume() > 1e-5 for piece in mesh.decompose())
        row = dict(cutter_radius_mm=radius, removed_mm3=removed, protected_guide_removed_mm3=protected_loss,
                   positive_components=components, source_origin_mm=origin.tolist())
        options.append((row, mesh))
    return options


def export_lattice(mesh, path, origin):
    try:
        return export(mesh, path, origin)
    except AssertionError:
        # Retain the mesh's native indexed edges when welding coincident points
        # would join separate lattice faces, as with the parent TPU source.
        vertices, faces = mesh_geom.load_obj(str(path))
        _, boundary, nonmanifold = mesh_geom.edge_census(faces)
        mf = A.I._manifold()
        result = mf.Manifold(mf.Mesh64(np.ascontiguousarray(vertices + origin),
                                       np.ascontiguousarray(faces, dtype=np.uint64)))
        error = abs(result.volume() - mesh.volume()) / mesh.volume()
        if (boundary or nonmanifold or result.status() != mf.Error.NoError
                or sum(piece.volume() > 1e-5 for piece in result.decompose()) != 1 or error > 1e-5):
            raise ValueError("Native lattice mesh failed indexed readback") from None
        return result, dict(**provenance(path), translation_source_mm=origin.tolist(),
                            mesh_import_mode="native_indexed_without_welding", roundtrip_relative_volume_error=error)


def main():
    out = DEST.parent / "elbow_local_muscle_load"
    response_file = out / "elbow_series_free_slider.json"
    long_file = out / "elbow_extensor_long_stroke.json"
    mount_file = out / "elbow_series_spring_mount.json"
    guide_file = out / "elbow_series_guide.json"
    response, prior, mount, guide = map(read, (response_file, long_file, mount_file, guide_file))
    if response["guide"] != provenance(guide_file) or mount["guide"] != provenance(guide_file):
        raise RuntimeError("Guide changed between the spring mount and force simulation")
    sources = [response["local_model"], prior["short_guide_parent_mesh"], prior["short_guide_definition"],
               guide["stop_assembly"], guide["package"]]
    for source in sources:
        if provenance(ROOT / source["file"]) != source:
            raise RuntimeError(f"Source changed: {source['file']}")
    model = mujoco.MjModel.from_xml_string((ROOT / response["local_model"]["file"]).read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    slider = model.joint("FJ1477_series_travel").qposadr.item()
    extensor = data.body("FJ1477_distal_pivot")
    neutral_distal = extensor.xpos.copy() * 1000.
    neutral_slide = data.body("FJ1477_carriage").xpos.copy() * 1000.
    slide_axis = np.asarray(guide["axis_source"])
    passive = read(DEST / "elbow_passive_model/model_check.json")
    assembly_file = ROOT / passive["assembly"]["file"]
    if provenance(assembly_file) != passive["assembly"]:
        raise RuntimeError("Original mount assembly changed")
    assembly = read(assembly_file)
    original = assembly["muscles"]["extensor"]
    internal_path = np.asarray(original["wire_path_source_mm"][1:-1])
    tpu_translation = assembly["display_only_TPU_sources"]["tpu:elbow_extensor"]["translation_source_mm"]
    short_design = read(ROOT / prior["short_guide_definition"]["file"])
    stations = short_design["actual_short_tube_guides"]
    centres = np.asarray([row["centre_display_mm"] for row in stations]) + tpu_translation
    np.testing.assert_allclose(centres, internal_path, atol=1e-8)
    half_axes = .5 * np.asarray([row["length_mm"] for row in stations])[:, None] * np.asarray(
        [row["axis_display"] for row in stations])
    muscle_axis = np.asarray(short_design["basis_columns_display"])[:, 2]
    gradient = muscle_axis / ((internal_path[-1] - internal_path[0]) @ muscle_axis)
    vertices, faces = mesh_geom.load_obj(str(ROOT / prior["short_guide_parent_mesh"]["file"]))
    vertices = vertices + tpu_translation
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    edges = mesh_edges(faces)
    edge_lengths = np.linalg.norm(vertices[edges[:, 1]] - vertices[edges[:, 0]], axis=1)
    edges, edge_lengths = edges[edge_lengths >= .2], edge_lengths[edge_lengths >= .2]
    neutral_length = np.linalg.norm(internal_path[-1] - internal_path[0])
    nominal_segments = polyline_lengths(internal_path)
    samples = {key: (points, weigh(points, internal_path[0], gradient))
               for key, points in (("path", internal_path), ("vertices", vertices),
                                   ("entry", centres - half_axes), ("exit", centres + half_axes))}
    rigid, origins = load_geometry(passive, read(ROOT / guide["stop_assembly"]["file"]))
    package = read(ROOT / guide["package"]["file"])
    moving_names = package["proximal_moving_parts"]
    fixed_rigid = {key: part for key, (owner, part) in rigid.items()
                   if owner not in origins and key not in moving_names and key != ADAPTER}
    moving_rigid = {key: row for key, row in rigid.items() if row[0] in origins}
    spring = mount["spring_mount"]
    installed = {key: load_component(row) for key, row in spring["sources"].items()}
    fixed_rigid.update({"series_adapter": installed["adapter"], "spring_fixed_pin": installed["fixed_pin"]})
    fixed_rigid.update({f"rail_{side}_{name}": load_component(source)
                        for side, lane in mount["retention"].items()
                        for name, source in lane["sources"].items()})
    carried_rigid = {key: rigid[key][1] for key in moving_names if key != PREFIX + "proximal_rigid_cap"}
    carried_rigid.update({"series_cap": installed["cap"], "spring_moving_pin": installed["moving_pin"]})
    fixed_bounds = bounded(fixed_rigid)
    spring_start = np.asarray(spring["fixed_eye_source_mm"])
    spring_end = np.asarray(spring["moving_eye_neutral_source_mm"])
    mf = A.I._manifold()
    neutral_muscle = mf.Manifold(mf.Mesh64(np.ascontiguousarray(vertices), faces))
    old_cap = rigid[PREFIX + "proximal_rigid_cap"][1]
    cap_overlap = neutral_muscle ^ installed["cap"]
    cap_cutter = cap_overlap.hull().minkowski_sum(mf.Manifold.cube([.3] * 3, True))
    trimmed_cap = (installed["cap"] - (cap_cutter - old_cap)).as_original().simplify(1e-6)
    cap_trim = installed["cap"] - trimmed_cap
    protected_sleeve_losses = [float((cap_trim ^ rod(np.asarray(base) - 17. * slide_axis,
                                                    np.asarray(base) - 10. * slide_axis, 4.)).volume())
                               for base in guide["selected"]["rail_origins_source_mm"]]
    if max(protected_sleeve_losses) > .001 or (old_cap - trimmed_cap).volume() > .001:
        raise RuntimeError("Cap relief cuts through original bearing or cap material")
    trimmed_cap, trimmed_cap_source = export(trimmed_cap, out / "series_cap_short_guide_relief.obj", internal_path[0])
    rail_base = np.asarray(guide["selected"]["rail_origins_source_mm"][0])
    relief = relief_meshes(vertices, faces, internal_path[0], centres, half_axes, stations,
                           rail_base, slide_axis, guide["travel_mm"])
    possible_reliefs = {}
    for row, mesh in relief:
        if row["positive_components"] != 1 or row["protected_guide_removed_mm3"] > .001:
            continue
        radius = row["cutter_radius_mm"]
        for tolerance in (1e-6, 1e-5, 1e-4):
            try:
                exported, source = export_lattice(mesh.as_original().simplify(tolerance),
                                                  out / f"elbow_short_guides_rail_relief_r{radius:g}.obj",
                                                  internal_path[0])
            except (AssertionError, ValueError) as error:
                row["export_error"] = str(error) or type(error).__name__
                continue
            row.pop("export_error", None)
            row["candidate_source"] = source
            points, triangles = mesh_geom.manifold_mesh_to_arrays(exported.to_mesh64())
            possible_reliefs[radius] = (points, weigh(points, internal_path[0], gradient), triangles.astype(np.uint64))
            break
    results = {}
    for index in (2, 4):
        case = response["cases"][index]
        trace_path = ROOT / case["trace"]["file"]
        if provenance(trace_path) != case["trace"]:
            raise RuntimeError("Slider trace fingerprint changed")
        with trace_path.open(newline="", encoding="utf-8") as stream:
            trace = list(csv.DictReader(stream))
        selected = sample_indices(trace)
        motion, mesh_poses = [], []
        for step, row in enumerate(trace):
            angles = [float(row[f"{name}_deg"]) for name in JOINTS]
            travel = float(row["carriage_mm"])
            data.qpos[addresses] = np.deg2rad(angles)
            data.qpos[slider] = travel * .001
            mujoco.mj_kinematics(model, data)
            rotation = extensor.xmat.reshape(3, 3)
            distal_translation = extensor.xpos * 1000. - rotation @ neutral_distal
            proximal_shift = data.body("FJ1477_carriage").xpos * 1000. - neutral_slide

            def pose(name):
                return displacement_map(*samples[name], rotation, distal_translation, proximal_shift)

            path, entry, exit_ = pose("path"), pose("entry"), pose("exit")
            np.testing.assert_allclose(path[[0, -1]], data.site_xpos[
                [model.site("FJ1477_origin").id, model.site("FJ1477_distal").id]] * 1000., atol=1e-7)
            chord = np.linalg.norm(path[-1] - path[0])
            np.testing.assert_allclose(chord - neutral_length, float(row["FJ1477_span_change_mm"]), atol=1e-7)
            segments = polyline_lengths(path)
            gap = np.linalg.norm(entry[1:] - exit_[:-1], axis=1)
            record = dict(time_s=float(row["time_s"]), elbow_deg=angles[0], carriage_mm=travel,
                          bore_growth_mm=float(chord - neutral_length),
                          wire_span_change_mm=float(segments.sum() - nominal_segments.sum()),
                          max_guide_spacing_ratio=float((segments / nominal_segments).max()),
                          max_between_guide_gap_mm=float(gap.max()),
                          wire_polyline_min_three_point_radius_mm=turn_radius(path))
            motion.append(record)
            if step not in selected:
                continue
            posed = pose("vertices")
            interior = (samples["vertices"][1] > 0.) & (samples["vertices"][1] < 1.)
            displacement = vertices[interior] @ (rotation - np.eye(3)).T + distal_translation - proximal_shift
            jacobian = (np.eye(3) + samples["vertices"][1][interior, None, None] * (rotation - np.eye(3))
                        + displacement[:, :, None] * gradient[None, None, :])
            determinant = float(np.linalg.det(jacobian).min())
            solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), faces))
            if solid.status() != mf.Error.NoError or determinant <= 0.:
                raise RuntimeError(f"Imposed geometry is invalid at t={row['time_s']}")
            obstacles = {**fixed_bounds, **bounded(pose_solids(data, moving_rigid, origins)),
                         **bounded({key: part.translate(proximal_shift) for key, part in carried_rigid.items()}),
                         **bounded({"spring_envelope": spring_envelope(spring_start, spring_end + proximal_shift, slide_axis)})}
            wire = A.D.outer_envelope(path, .3)
            wire_extra = A.D.outer_envelope(path, .45)
            ratios = np.linalg.norm(posed[edges[:, 1]] - posed[edges[:, 0]], axis=1) / edge_lengths
            proposed = {}
            revised_obstacles = {**obstacles, **bounded({"series_cap": trimmed_cap.translate(proximal_shift)})}
            for radius, (points, weights, triangles) in possible_reliefs.items():
                new_points = displacement_map(points, weights, rotation, distal_translation, proximal_shift)
                new_solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(new_points), triangles))
                proposed[str(radius)] = dict(
                    geometry_status=str(new_solid.status()),
                    tpu_rigid_overlap_mm3=overlaps(new_solid, revised_obstacles)
                    if new_solid.status() == mf.Error.NoError else None,
                    local_wire_overlap_mm3=float(A.I.intersection_volume(wire, new_solid)),
                )
            mesh_poses.append(dict(
                **record, source_pose_index=step, min_map_jacobian=determinant,
                surface_edge_ratio_quantiles=np.quantile(ratios, [.05, .5, .95]).tolist(),
                tpu_rigid_overlap_mm3=overlaps(solid, obstacles),
                local_wire_rigid_overlap_mm3=overlaps(wire, obstacles),
                local_wire_tpu_overlap_mm3=float(A.I.intersection_volume(wire, solid)),
                enlarged_wire_tpu_overlap_mm3=float(A.I.intersection_volume(wire_extra, solid)),
                relief_probes=proposed,
            ))
        results[str(index)] = dict(trace=provenance(trace_path), pose_count=len(motion),
                                   mesh_samples=mesh_poses, motion=motion)
        print(index, "peak bore growth", max(row["bore_growth_mm"] for row in motion),
              "max gap", max(row["max_between_guide_gap_mm"] for row in motion), flush=True)
        for sample in mesh_poses:
            print(round(sample["time_s"], 4), round(sample["elbow_deg"], 3),
                  sample["local_wire_tpu_overlap_mm3"], sample["local_wire_rigid_overlap_mm3"],
                  sample["tpu_rigid_overlap_mm3"], flush=True)
    selected_relief = None
    samples_for_selection = [pose for case in results.values() for pose in case["mesh_samples"]]
    for row, _ in relief:
        radius = row["cutter_radius_mm"]
        if radius not in possible_reliefs:
            continue
        tests = [pose["relief_probes"][str(radius)] for pose in samples_for_selection]
        row["checked_poses"] = len(tests)
        row["possible_contacts"] = sorted({
            key for test in tests for key in (test["tpu_rigid_overlap_mm3"] or {})})
        if (any(name not in ("bone:humerus", "bone:forearm") for name in row["possible_contacts"])
                or any(test["geometry_status"] != str(mf.Error.NoError) or test["local_wire_overlap_mm3"] > .001
                       for test in tests)):
            continue
        selected_relief = radius
        break
    save(out / "elbow_slider_short_guides.json", dict(
        generator=provenance(Path(__file__).resolve()), response=provenance(response_file),
        short_guide_mesh=prior["short_guide_parent_mesh"], short_guide_design=prior["short_guide_definition"],
        spring_mount=provenance(mount_file), guide_count=len(stations), local_wire_diameter_mm=.6,
        enlarged_test_wire_diameter_mm=.9, rigid_overlap_threshold_mm3=.001, cases=results,
        shaft_relief_options=[row for row, _ in relief], selected_relief_radius_mm=selected_relief,
        trimmed_cap=dict(source=trimmed_cap_source, neutral_original_overlap_mm3=float(cap_overlap.volume()),
                         removed_mm3=float(cap_trim.volume()), protected_sleeve_removed_mm3=protected_sleeve_losses),
        scope="Only the local guide path between existing proximal and distal bores. Finite poses "
        "use prescribed interpolation between the actual sliding and rotating end frames. Straight "
        "spans connect physical guide centres; friction, wire departure within holes, connector "
        "tails, restoring force and flexible contact are not solved. Interference with spring "
        "envelope means possible space conflict, not real TPU damage or rigid obstruction. "
        "Relief candidates cut original short-guide stock only around shaft A and protect all guide "
        "stations. A separate cap trim removes overlap from new bracket stock outside the sleeve regions."))
    print("relief options", [row for row, _ in relief], flush=True)


if __name__ == "__main__":
    main()
