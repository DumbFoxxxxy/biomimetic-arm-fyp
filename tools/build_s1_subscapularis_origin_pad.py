"""Refit the FJ1504 TPU origin transition to its recessed cap."""

from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_elbow_slider_short_guides import export_lattice
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, rod, save
from study_elbow_extensor_long_stroke import turn_radius
from study_elbow_series_package import bounded, overlaps, union


def nearest_station(points, path):
    delta = np.diff(path, axis=0)
    length = np.linalg.norm(delta, axis=1)
    station = np.r_[0., np.cumsum(length)]
    best_dist = np.full(len(points), np.inf)
    result = np.zeros(len(points))
    for index, (p, d, size) in enumerate(zip(path[:-1], delta, length, strict=True)):
        if size < 1e-7:
            continue
        w = np.clip((points - p) @ d / size**2, 0., 1.)
        d2 = ((points - p - w[:, None] * d)**2).sum(axis=1)
        better = d2 < best_dist
        best_dist[better] = d2[better]
        result[better] = station[index] + size * w[better]
    return result


def fade(station, span):
    f = np.clip(station / span, 0., 1.)
    return 1. - f * f * (3. - 2. * f)


def reconnect(body, rigid, lumen):
    mf = A.I._manifold()
    solids = sorted((part for part in body.decompose() if part.volume() > 1e-8),
                    key=lambda part: -part.volume())
    if len(solids) <= 1:
        return body, []
    main = solids[0]
    links = []
    for island in solids[1:]:
        bulk = mesh_geom.manifold_mesh_to_arrays(main.to_mesh64())[0]
        pieces = mesh_geom.manifold_mesh_to_arrays(island.to_mesh64())[0]
        distance, index = cKDTree(bulk).query(pieces)
        for i in np.argsort(distance)[:16]:
            if distance[i] > 5.:
                break
            tip, root = pieces[i], bulk[index[i]]
            for radius in (.8, 1.1, 1.4):
                bridge = A.D.outer_envelope(np.array([tip, root]), radius) - lumen
                if overlaps(bridge, rigid):
                    continue
                if min(float((bridge ^ part).volume()) for part in (main, island)) <= 1e-5:
                    continue
                links.append(dict(radius_mm=radius, points_mm=[tip.tolist(), root.tolist()],
                                  added_mm3=float((bridge - main - island).volume())))
                main = union([main, island, bridge])
                break
            else:
                continue
            break
        else:
            return body, []
        if distance[i] > 5.:
            return body, []
    return union([body, main]).as_original().simplify(1e-6), links


def remove_terminal_scraps(body, terminal_region, protected_wall):
    parts = sorted((part for part in body.decompose() if part.volume() > 1e-8),
                   key=lambda part: -part.volume())
    scraps = parts[1:]
    if (not scraps or sum(part.volume() for part in scraps) >= 1.
            or any((part - terminal_region).volume() > 1e-6
                   or (part ^ protected_wall).volume() > 1e-6 for part in scraps)):
        return body, []
    return (body - union(scraps)).as_original().simplify(1e-6), [
        dict(volume_mm3=float(part.volume()), bounds_source_mm=list(part.bounding_box()))
        for part in scraps
    ]


def smooth_neck(path, normal, world):
    proposals = []
    u = np.linspace(0., 1., 24)[:, None]
    for index in (5, 7, 9, 11, 13, 15):
        direction = path[index + 1] - path[index]
        direction /= np.linalg.norm(direction)
        for first in (2., 3., 4., 5., 6.):
            for last in (2., 3., 4., 5., 6.):
                p0, p3 = path[0], path[index]
                p1, p2 = p0 - first * normal, p3 - last * direction
                arc = (1. - u)**3 * p0 + 3. * u * (1. - u)**2 * p1 + 3. * u**2 * (1. - u) * p2 + u**3 * p3
                curve = np.vstack([arc, path[index + 1:]])
                radius = turn_radius(curve[:len(arc) + 1])
                if radius is not None:
                    proposals.append((radius, curve, len(arc)))
    for radius, curve, arc_count in sorted(proposals, key=lambda row: -row[0]):
        if not overlaps(A.D.outer_envelope(curve[:arc_count], .3), world):
            return curve, radius
    return path, turn_radius(path[:16])


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    assembly_file = folder / "assembly_scapular_relief.json"
    revision_file = folder / "scapular_thorax_relief.json"
    assembly, revision = map(read, (assembly_file, revision_file))
    if assembly["scapular_mount_revision"] != provenance(revision_file):
        raise RuntimeError("Scapular revision changed")
    route_file = ROOT / revision["route_reference"]["file"]
    if provenance(route_file) != revision["route_reference"]:
        raise RuntimeError("Origin route reference changed")
    route = next(row for row in read(route_file)["routes"] if row["head"] == "FJ1504")
    port = route["ports"]["origin"]
    normal = np.asarray(port["normal_source"])
    original_inner = np.asarray(port["inner_face_source_mm"])
    shift = next(row["selected"]["shift_mm"] for row in revision["reports"] if row["head"] == "FJ1504") * normal
    inner = original_inner + shift
    old_path = np.asarray(route["full_wire_path_source_mm"])
    np.testing.assert_allclose(old_path[2], original_inner, atol=1e-8)
    old_path = old_path[2:]
    muscle_source = assembly["flexible_geometry_sources"]["tpu:subscapularis"]
    original = load_component(muscle_source)
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(original.to_mesh64())
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    stations = nearest_station(vertices, old_path)
    wire_stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(old_path, axis=0), axis=1))]
    fixed = {key: load_component(row) for key, row in assembly["sources"].items()
             if key.startswith(("thorax_reference:", "origin:FJ1504/"))
             or key in ("bone:scapula", "cap:FJ1504/origin")}
    world = bounded(fixed)
    cap = fixed["cap:FJ1504/origin"]
    support = fixed["origin:FJ1504/front_support"]
    mf = A.I._manifold()
    cutter_margin = mf.Manifold.cube([.2] * 3, True)
    trials = []
    selected = None
    for span in (26., 27., 28., 29., 30., 24., 32., 36.):
        moved = vertices + fade(stations, span)[:, None] * shift
        core = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
        row = dict(blend_span_mm=span, warped_mesh_status=str(core.status()))
        trials.append(row)
        if core.status() != mf.Error.NoError:
            continue
        path, bend_radius = smooth_neck(old_path + fade(wire_stations, span)[:, None] * shift, normal, world)
        # At the unchanged axial cap bore the transition begins on its centre.
        # The lumen drill extends outside the pad to preserve an open end.
        lumen = A.D.outer_envelope(np.vstack([inner + 1. * normal, path]), .86)
        new_stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
        active_transition = path[:np.searchsorted(new_stations, span + 8.) + 1]
        protected_wall = A.D.outer_envelope(active_transition, 1.65) - lumen
        proximal_vertices = moved[stations <= span]
        lower, upper = proximal_vertices.min(axis=0) - 1., proximal_vertices.max(axis=0) + 1.
        terminal_region = mf.Manifold.cube(upper - lower).translate(lower)
        wire = A.D.outer_envelope(active_transition, .3)
        wire_hits = overlaps(wire, world)
        row["local_wire_rigid_overlap_mm3"] = wire_hits
        for radius in (10.8, 9.6, 8.4):
            pad = rod(inner - 4. * normal, inner, radius)
            body = (core + pad) - lumen
            initial_hits = overlaps(body, world)
            removed_contact = union(part.hull().minkowski_sum(cutter_margin)
                                    for key in initial_hits for part in (body ^ fixed[key]).decompose()
                                    if part.volume() > 1e-7)
            candidate = (body - removed_contact).as_original().simplify(1e-6)
            candidate, links = reconnect(candidate, world, lumen)
            candidate, scraps = remove_terminal_scraps(candidate, terminal_region, protected_wall)
            components = sorted([float(part.volume()) for part in candidate.decompose() if part.volume() > 1e-6],
                                reverse=True)
            result = dict(pad_radius_mm=radius, input_overlap_mm3=initial_hits,
                          local_clearance_cut_mm3=float((body - candidate).volume()),
                          positive_component_volumes_mm3=components, retained_lattice_bridges=links,
                          removed_detached_terminal_scraps=scraps)
            row.setdefault("pad_trials", []).append(result)
            if len(components) != 1:
                continue
            # Estimate available load bearing from a small normal displacement
            # into each rigid face, outside the candidate's free contact gap.
            compression = .35
            result["cap_bearing_probe_mm3"] = float((candidate.translate(compression * normal) ^ cap).volume())
            result["support_bearing_probe_mm3"] = float((candidate.translate(-compression * normal) ^ support).volume())
            result["residual_overlap_mm3"] = overlaps(candidate, world)
            if (result["residual_overlap_mm3"] or wire_hits
                    or min(result["cap_bearing_probe_mm3"], result["support_bearing_probe_mm3"]) <= .001):
                continue
            try:
                actual, source = export_lattice(candidate, folder / "FJ1504_TPU_recessed_origin.obj", inner)
            except (AssertionError, ValueError) as error:
                result["export_error"] = str(error) or type(error).__name__
                continue
            roundtrip = overlaps(actual, world)
            wire_overlap = float((actual ^ wire).volume())
            if roundtrip or wire_overlap > .001:
                result["readback_hits_mm3"] = roundtrip
                result["readback_wire_overlap_mm3"] = wire_overlap
                continue
            selected = dict(**result, blend_span_mm=span, source=source, path_source_mm=path.tolist(),
                            proximal_wire_sampled_bend_radius_mm=bend_radius,
                            readback_wire_overlap_mm3=wire_overlap)
            break
        if selected:
            break
    save(folder / "FJ1504_recessed_TPU_origin.json", dict(
        generator=provenance(Path(__file__).resolve()), assembly=provenance(assembly_file),
        original_TPU=muscle_source, cap_inner_source_mm=inner.tolist(), proximal_shift_mm=shift.tolist(),
        selected=selected, trials=trials, scope="Neutral design mesh with nominal axial interpolation "
        "and a reconstructed end pad. This is a rest-shape design study, not a stress prediction. "
        "Cap and saddle bearing probes show potential material engagement; friction, compression, "
        "force, actual TPU movement and assembly remain unverified."))
    print("selected TPU", selected, flush=True)
    for trial in trials:
        print(trial["blend_span_mm"], trial["local_wire_rigid_overlap_mm3"], [
            (row["pad_radius_mm"], len(row["positive_component_volumes_mm3"]), row.get("residual_overlap_mm3"))
            for row in trial["pad_trials"]], flush=True)


if __name__ == "__main__":
    main()
