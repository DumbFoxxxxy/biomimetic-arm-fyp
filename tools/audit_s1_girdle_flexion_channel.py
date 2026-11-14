"""Replay the FJ1468 channel on the checked SC/AC and GH flexion paths."""

import argparse

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_ac_blind_coupled_motion import load_assembly
from audit_elbow_slider_short_guides import export_lattice
from audit_s1_infraspinatus_distributed_deformation import blend
from audit_s1_series_elbow_integration import MOBILE, compose
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, export, load_component, save
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_series_package import bounded, overlaps, union
from study_s1_coupled_scapula import ShoulderAssembly, affine_at
from study_s1_infraspinatus_wire_travel import route_metrics
from study_s1_local_force_geometry import load_reference


def neutral_overlap(intersection, clavicle, humerus, origin, axis, lower, upper):
    """Invert the same axial blend used to pose the lattice."""
    target = np.asarray(intersection.to_mesh64().vert_properties)[:, :3]
    a, b = clavicle[:, :3], humerus[:, :3]
    ta, tb = clavicle[:, 3], humerus[:, 3]
    neutral = (target - ta) @ a
    for _ in range(24):
        station = (neutral - origin) @ axis
        weight, slope = blend(station, lower, upper, .1)
        delta = neutral @ (b - a).T + tb - ta
        error = neutral @ a.T + ta + weight[:, None] * delta - target
        residual = float(np.linalg.norm(error, axis=1).max())
        if residual < 1e-7:
            break
        jac = a + weight[:, None, None] * (b - a) + slope[:, None, None] * delta[:, :, None] * axis
        neutral -= np.linalg.solve(jac, error[:, :, None])[:, :, 0]
    else:
        raise ValueError(f"Overlap inverse did not converge: {residual:g} mm")
    return neutral, residual


def screen_detour(points, stations, pose_targets):
    """Compare local smooth offsets against the two obstructing fasteners."""
    # Enclose segments 57–59 with a smooth taper and retain both port axes.
    start, stop = 55, 61
    chord = points[stop] - points[start]
    chord /= np.linalg.norm(chord)
    u = np.cross(chord, np.eye(3)[np.argmin(np.abs(chord))])
    u /= np.linalg.norm(u)
    v = np.cross(chord, u)
    local_stations = np.unique(np.r_[stations[start:stop + 1], np.linspace(stations[start], stations[stop], 25)])
    nominal = np.column_stack([np.interp(local_stations, stations, points[:, i]) for i in range(3)])
    fraction = (local_stations - local_stations[0]) / (local_stations[-1] - local_stations[0])
    taper = np.sin(np.pi * fraction) ** 2
    trials = []
    for amplitude in (0., 1., 2., 3., 4., 6.):
        for azimuth in (range(0, 360, 45) if amplitude else (0,)):
            theta = np.deg2rad(azimuth)
            direction = u * np.cos(theta) + v * np.sin(theta)
            proposed = nominal + amplitude * taper[:, None] * direction
            candidate = np.vstack((points[:start], proposed, points[stop + 1:]))
            score = 0.
            blocked = 0
            for deform, targets in pose_targets:
                envelope = A.D.outer_envelope(deform(proposed), .40)
                hits = overlaps(envelope, targets)
                blocked += bool(hits)
                score = max(score, sum(hits.values()))
            length, _, _ = route_metrics(candidate)
            trials.append((dict(amplitude_mm=amplitude, azimuth_deg=azimuth, blocked_poses=blocked,
                                worst_overlap_mm3=score, neutral_arc_mm=length), candidate))
    trials.sort(key=lambda trial: (trial[0]["blocked_poses"], trial[0]["worst_overlap_mm3"],
                                   trial[0]["amplitude_mm"]))
    return trials


def locate_wire_hits(route, targets, hits, neutral_stations, weights):
    locations = {}
    for part in hits:
        solid, box = targets[part]
        lower, upper = np.asarray(box).reshape(2, 3)
        starts, ends = route[:-1], route[1:]
        nearby = np.all(np.minimum(starts, ends) - .30 <= upper, axis=1)
        nearby &= np.all(np.maximum(starts, ends) + .30 >= lower, axis=1)
        segments = []
        for index in np.flatnonzero(nearby):
            wire = A.D.outer_envelope(route[index:index + 2], .30)
            intersection = wire ^ solid
            volume = float(intersection.volume())
            if volume > 1e-6:
                segments.append(dict(segment=int(index), neutral_station_mm=neutral_stations[index:index + 2].tolist(),
                                     distal_blend_weights=weights[index:index + 2].tolist(),
                                     overlap_mm3=volume, intersection_bounds_mm=list(intersection.bounding_box())))
        locations[part] = segments
    return locations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--locate", action="store_true", help="Localize rigid wire hits without replaying TPU contact")
    mode.add_argument("--detour", action="store_true", help="Screen local centreline offsets and replay the best route")
    mode.add_argument("--adapter-relief", action="store_true", help="Screen a local relief in the elbow flexor adapter")
    mode.add_argument("--lined-detour", action="store_true", help="Build and replay the detour with a TPU liner")
    mode.add_argument("--liner-contact", action="store_true", help="Localize contacts around the exported liner")
    mode.add_argument("--liner-assembly", action="store_true", help="Replay the liner with the 1.10 mm adapter relief")
    mode.add_argument("--protected-belly", action="store_true", help="Replay the protected 0.25 mm inset lattice")
    mode.add_argument("--neighbor-wires", action="store_true", help="Compare neighboring wire passage before and after TPU relief")
    mode.add_argument("--posterior-exits", action="store_true", help="Locate posterior wire rigid overlaps along its route")
    mode.add_argument("--map-belly-overlap", action="store_true", help="Map the worst belly overlaps to neutral TPU")
    mode.add_argument("--liner-relief", action="store_true", help="Screen adapter cuts around the exported liner")
    mode.add_argument("--liner-hull-relief", action="store_true", help="Screen conservative convex relief cuts")
    args = parser.parse_args()
    replay_protected = args.protected_belly or args.neighbor_wires or args.posterior_exits
    local = DEST.parent / "s1_local_muscle_load"
    folder = local / "FJ1500_cap_relief"
    rigid_ref = provenance(local / "FJ1468_girdle_flexion_clearance.json")
    guide_ref = provenance(folder / "FJ1468_revised_channel_motion.json")
    rigid, guide = load_reference(rigid_ref), load_reference(guide_ref)
    path = load_reference(rigid["girdle_path"])
    span = load_reference(rigid["span_study"])
    endpoints = load_reference(span["local_endpoints"])
    if path["closure_audit"] != provenance(folder / "FJ1504_AC_blind_closure_audit.json"):
        raise ValueError("Rigid assembly provenance changed")
    origins = load_reference(path["pose_reference"])["neutral_origins_m"]
    assembly, _, sources = load_assembly(folder)
    model = ShoulderAssembly({}, {}, origins)
    axis = np.asarray(endpoints["axes_source"][0])

    wall = load_reference(guide["wall_audit"])
    transition = load_reference(wall["source"])
    points = np.array(next(row["centreline_source_mm"] for row in transition["recuts"] if "mesh" in row))
    _, vertices, faces = indexed_solid(guide["merged_TPU"])
    length, _, _ = route_metrics(points)
    np.testing.assert_allclose(length, guide["neutral_arc_mm"], atol=1e-8, rtol=0.)
    toward = (points[-1] - points[0]) / np.linalg.norm(points[-1] - points[0])
    lower, upper = guide["deformation_window_mm"]
    weights, derivative = blend((vertices - points[0]) @ toward, lower, upper, .1)
    path_weights, _ = blend((points - points[0]) @ toward, lower, upper, .1)
    stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    solids = {key: load_component(meta) for key, meta in sources.items()}
    owners = {key: assembly["owners"].get(key, "scapula") for key in solids}
    mf = A.I._manifold()
    liner = {}
    check_neighborhood = args.liner_contact or args.liner_assembly or replay_protected
    if args.lined_detour or check_neighborhood or args.liner_relief or args.liner_hull_relief or args.map_belly_overlap:
        relief_ref = provenance(local / "FJ1486_adapter_FJ1468_relief.json")
        relief = load_reference(relief_ref)
        detour = load_reference(relief["detour"])
        if relief["rigid_audit"] != rigid_ref or detour["guide_motion"] != guide_ref:
            raise ValueError("Liner inputs do not match the active assembly")
        choice = relief["cases"]["detour"]
        if choice["readback_hit_poses"] or choice["bolt_hit_poses"]:
            raise ValueError("Adapter clearance candidate has unresolved wire hits")
        adapter_key = "elbow:FJ1486_bone_adapter"
        solids[adapter_key] = load_component(choice["mesh"])
        candidate = np.asarray(detour["neutral_path_source_mm"])
        # The exact shared prefix and suffix keep the modified lumen connected
        # to the unmodified channel. Include one common segment at each end.
        prefix = next(i for i, (a, b) in enumerate(zip(points, candidate)) if not np.array_equal(a, b))
        suffix = next(i for i, (a, b) in enumerate(zip(points[::-1], candidate[::-1])) if not np.array_equal(a, b))
        liner_slice = slice(prefix - 2, len(candidate) - suffix + 2)
        if check_neighborhood or args.liner_relief or args.liner_hull_relief or args.map_belly_overlap:
            lined_ref = provenance(local / "FJ1468_girdle_flexion_channel_lined_detour.json")
            saved = load_reference(lined_ref)
            if saved["adapter_relief"] != relief_ref or saved["rigid_audit"] != rigid_ref:
                raise ValueError("Exported liner does not match the active assembly")
            np.testing.assert_array_equal(candidate, saved["neutral_path_source_mm"])
            mesh = saved["lined_tpu"]
            liner = dict(lined_detour=lined_ref, lined_tpu=mesh)
        else:
            body, _, _ = indexed_solid(guide["merged_TPU"])
            outer = A.D.outer_envelope(candidate[liner_slice], .95)
            bore = A.D.outer_envelope(candidate, .50)
            wall = outer - bore
            bond = float((wall ^ body).volume())
            lined = ((body + outer) - bore).as_original().simplify(1e-7)
            if bond <= .001 or len([part for part in lined.decompose() if part.volume() > 1e-5]) != 1:
                raise ValueError("The added liner is not integrated into the lattice")
            _, mesh = export_lattice(lined, local / "FJ1468_lined_detour_TPU.obj", points[0])
            liner = dict(adapter_relief=relief_ref, lined_tpu=mesh, bore_radius_mm=.50, outer_radius_mm=.95,
                         shared_wall_volume_mm3=bond, liner_net_TPU_volume_change_mm3=float(lined.volume() - body.volume()))
        _, vertices, faces = indexed_solid(mesh)
        weights, derivative = blend((vertices - points[0]) @ toward, lower, upper, .1)
        path_weights, _ = blend((candidate - points[0]) @ toward, lower, upper, .1)
        points = candidate
        stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    if args.liner_assembly or args.map_belly_overlap or replay_protected:
        hull_ref = provenance(local / "FJ1486_adapter_FJ1468_liner_relief_hull.json")
        hull = load_reference(hull_ref)
        prior_contact = load_reference(hull["liner_contact"])
        if prior_contact["lined_detour"] != lined_ref or not hull["convex_cutter"]:
            raise ValueError("Adapter relief does not match the exported liner")
        choice = next(row for row in hull["cases"] if row["radius_mm"] == 1.10)
        if choice["components"] != 1 or choice["guard_hits_mm3"] or choice["readback_hit_poses"]:
            raise ValueError("Expected a checked adapter relief")
        solids[adapter_key] = load_component(choice["mesh"])
        liner.update(adapter_hull_relief=hull_ref, adapter_mesh=choice["mesh"])
    if replay_protected:
        trial_ref = provenance(local / "FJ1468_belly_overlap_neutral.json")
        trial_report = load_reference(trial_ref)
        baseline = load_reference(trial_report["assembly_audit"])
        if baseline["adapter_mesh"] != choice["mesh"] or baseline["lined_tpu"] != mesh or baseline["rigid_audit"] != rigid_ref:
            raise ValueError("Trial uses different assembly geometry")
        trial = next(row for row in trial_report["protected_combined_trials"] if row["cutter_inset_mm"] == .25)
        if trial["guard_hits_mm3"] or "mesh" not in trial:
            raise ValueError("Protected trial did not qualify for replay")
        _, vertices, faces = indexed_solid(trial["mesh"])
        weights, derivative = blend((vertices - points[0]) @ toward, lower, upper, .1)
        liner.update(protected_belly_trial=trial_ref, replayed_tpu=trial["mesh"])
    if args.posterior_exits:
        inventory_ref = trial_report["neighbor_channel_reference"]
        inventory = load_reference(inventory_ref)
        posterior = next(row for row in inventory["routes"] if row["id"] == "FJ1467")
        if (posterior["anatomical_origin_body"], posterior["distal_body"]) != ("scapula", "humerus"):
            raise ValueError("Unexpected FJ1467 route owners")
        cable = np.asarray(posterior["local_route"]["points_source_mm"])
        weight, _ = blend((cable - points[0]) @ toward, lower, upper, .1)
        arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(cable, axis=0), axis=1))]
        previous_ref = provenance(local / "FJ1468_girdle_flexion_channel_neighbor_wires.json")
        previous = load_reference(previous_ref)
        if (previous["protected_belly_trial"] != trial_ref or previous["rigid_audit"] != rigid_ref
                or previous["adapter_mesh"] != choice["mesh"]):
            raise ValueError("Neighbor replay uses a different assembly")
        rows = {}
        for motion in previous["cases"].values():
            for sample in motion:
                rows[tuple(sample["pose_deg"])] = sample["neighbor_wires"]["FJ1467"]["rigid_hits_mm3"]
        results = []
        scapula_cuts = []
        scapula_wires = []
        scapular_stock = solids["bone:scapula"]
        for pose, expected in rows.items():
            angle, *sc, flex = pose
            clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
            gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
            humerus = compose(scapula, gh)
            a, b = scapula[:, :3], humerus[:, :3]
            ta, tb = scapula[:, 3], humerus[:, 3]
            posed = cable @ a.T + ta + weight[:, None] * (cable @ (b - a).T + tb - ta)
            frames = {"scapula": scapula, "clavicle": clavicle, **dict.fromkeys(MOBILE, humerus)}
            targets = bounded({key: solids[key].transform(frames[owners[key]]) for key in expected})
            wire = A.D.outer_envelope(posed, .30)
            hits = overlaps(wire, targets)
            if hits.keys() != expected.keys():
                raise ValueError(f"Posterior wire replay changed hit pairs at {pose}")
            np.testing.assert_allclose(list(hits.values()), list(expected.values()), atol=1e-6, rtol=0.)
            results.append(dict(pose_deg=list(pose), overlaps_mm3=hits,
                                locations=locate_wire_hits(posed, targets, hits, arc, weight)))
            proximal = (posed[:9] - ta) @ a
            envelope = A.D.outer_envelope(proximal, .50)
            scapula_wires.append(envelope)
            cut = envelope ^ scapular_stock
            if cut.volume() > 1e-8:
                scapula_cuts.append(cut)
        guard_keys = [key for key in solids if key.startswith(("origin:FJ1467/", "cap:FJ1467/"))]
        guard_stock = bounded({key: solids[key].minkowski_sum(mf.Manifold.cube([1.] * 3, True))
                               for key in guard_keys})
        cutter = mf.Manifold.hull_points(np.concatenate(
            [np.asarray(part.to_mesh64().vert_properties)[:, :3] for part in scapula_cuts]))
        removed = scapular_stock ^ cutter
        relieved = (scapular_stock - cutter).as_original().simplify(1e-7)
        proximal_trial = dict(radius_mm=.50, removed_mm3=float(removed.volume()),
                              remaining_components=len(relieved.decompose()),
                              guard_keys=guard_keys, guard_hits_mm3=overlaps(removed, guard_stock))
        if proximal_trial["remaining_components"] == 1 and not proximal_trial["guard_hits_mm3"]:
            actual, info = export(relieved, local / "FJ1467_scapula_exit_trial.obj", model.ac)
            proximal_trial.update(mesh=info, max_readback_overlap_mm3=max(
                float((envelope ^ actual).volume()) for envelope in scapula_wires))
        save(local / "FJ1467_posterior_rigid_locations.json", dict(
            generator=provenance(ROOT / "tools/audit_s1_girdle_flexion_channel.py"),
            neighbor_replay=previous_ref, route_inventory=inventory_ref, source_ports=posterior["endpoint_ports_reference"],
            neutral_source_route_mm=cable.tolist(), neutral_arc_mm=float(arc[-1]), cases=results,
            proximal_clearance_trial=proximal_trial,
            scope="Localizes previously observed rigid overlaps at sampled poses. Segment spheres overlap and "
                  "their volumes are not additive. Source ports are historical nominal coordinates. "
                  "Overlap with a bone is reported as metal-envelope intersection with actual solid stock.",
        ))
        print(f"Localized posterior wire interference in {len(results)} poses", flush=True)
        return
    if args.neighbor_wires:
        inventory_ref = trial_report["neighbor_channel_reference"]
        inventory = load_reference(inventory_ref)
        original_mesh, original_vertices, original_faces = indexed_solid(mesh)
        original_weights, _ = blend((original_vertices - points[0]) @ toward, lower, upper, .1)
        neighbors = {}
        for item in inventory["routes"]:
            if item["id"] in ("FJ1513", "FJ1467"):
                if item["anatomical_origin_body"] != "scapula" or item["distal_body"] != "humerus":
                    raise ValueError("Expected scapula-to-humerus neighbor route")
                cable = np.asarray(item["local_route"]["points_source_mm"])
                w, _ = blend((cable - points[0]) @ toward, lower, upper, .1)
                neighbors[item["id"]] = (cable, w)
        if len(neighbors) != 2 or original_mesh.status() != mf.Error.NoError:
            raise ValueError("Missing neighbor routes or invalid original TPU")
        liner["neighbor_channel_reference"] = inventory_ref
    if args.map_belly_overlap:
        reference = provenance(local / "FJ1468_girdle_flexion_channel_liner_assembly.json")
        prior = load_reference(reference)
        if prior["adapter_mesh"] != choice["mesh"] or prior["lined_tpu"] != mesh or prior["rigid_audit"] != rigid_ref:
            raise ValueError("Overlap report does not match the restored assembly")
        zones = {}
        inventory_ref = provenance(DEST / "route_inventory.json")
        inventory = load_reference(inventory_ref)
        neutral_tpu, _, _ = indexed_solid(mesh)
        screws_and_caps = {key: solids[key] for key in solids if key.startswith(
            ("origin:FJ1468/", "distal:FJ1468/", "cap:FJ1468/"))}
        guard_paths = {row["id"]: np.asarray(row["local_route"]["points_source_mm"])
                       for row in inventory["routes"] if row["id"] in {"FJ1513", "FJ1467"}}
        if len(guard_paths) != 2:
            raise ValueError("Missing neighboring deltoid reference channels")
        guards = bounded({**screws_and_caps, "channel:FJ1468": A.D.outer_envelope(points, 1.10),
                          **{f"channel:{name}": A.D.outer_envelope(route, 1.10)
                             for name, route in guard_paths.items()}})
        protected = union([solid for solid, _ in guards.values()])
        cuts = []
        peak_cases = []
        rows = [row for case in prior["cases"].values() for row in case]
        for key in ("bone:scapula", "elbow:FJ1477_bone_adapter"):
            worst = max(rows, key=lambda row: row["TPU_rigid_hits_mm3"].get(key, 0.))
            angle, *sc, flex = worst["pose_deg"]
            clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
            gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
            humerus = compose(scapula, gh)
            a, b = clavicle[:, :3], humerus[:, :3]
            ta, tb = clavicle[:, 3], humerus[:, 3]
            moved = vertices @ a.T + ta + weights[:, None] * (vertices @ (b - a).T + tb - ta)
            shaped = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
            frames = {"scapula": scapula, **dict.fromkeys(MOBILE, humerus)}
            hit = shaped ^ solids[key].transform(frames[owners[key]])
            np.testing.assert_allclose(hit.volume(), worst["TPU_rigid_hits_mm3"][key], atol=.001, rtol=0.)
            samples, residual = neutral_overlap(hit, clavicle, humerus, points[0], toward, lower, upper)
            stations_hit = (samples - points[0]) @ toward
            # Use the closest point on every nominal route segment, not only its nodes.
            offsets = samples[:, None, :] - points[:-1]
            segments = np.diff(points, axis=0)
            fraction = np.clip(np.einsum("ijk,jk->ij", offsets, segments) / np.sum(segments**2, axis=1), 0., 1.)
            distances = np.linalg.norm(offsets - fraction[:, :, None] * segments, axis=2).min(axis=1)
            cutter = mf.Manifold.hull_points(np.ascontiguousarray(samples))
            # Screen the actual removed material against channel and attachment
            # guards. A hull can bridge gaps between separate lattice contacts.
            removed = neutral_tpu ^ cutter
            remaining = neutral_tpu - cutter
            pieces = [part for part in hit.decompose() if part.volume() > 1e-6]
            local_hulls = []
            for piece in pieces:
                local_samples, _ = neutral_overlap(piece, clavicle, humerus, points[0], toward, lower, upper)
                local_hulls.append(mf.Manifold.hull_points(np.ascontiguousarray(local_samples)))
            local_cut = union(local_hulls)
            cuts.append(local_cut)
            peak_cases.append((key, clavicle, humerus, frames, worst))
            local_removed = neutral_tpu ^ local_cut
            local_remaining = neutral_tpu - local_cut
            zones[key] = dict(pose_deg=worst["pose_deg"], overlap_mm3=float(hit.volume()),
                              neutral_bounds_mm=[samples.min(axis=0).tolist(), samples.max(axis=0).tolist()],
                              axial_station_bounds_mm=[float(stations_hit.min()), float(stations_hit.max())],
                              minimum_sampled_distance_to_channel_mm=float(distances.min()),
                              vertices_within_1p1_mm=int(np.count_nonzero(distances <= 1.10)),
                              vertices_before_blend=int(np.count_nonzero(stations_hit <= lower)),
                              vertices_after_blend=int(np.count_nonzero(stations_hit >= upper)),
                              sampled_vertices=len(samples), maximum_inverse_residual_mm=residual,
                              hull_trial_removed_mm3=float(removed.volume()),
                              hull_trial_guard_hits_mm3=overlaps(removed, guards),
                              hull_trial_remaining_components=sum(piece.volume() > 1e-5
                                                                  for piece in remaining.decompose()),
                              component_hull_count=len(local_hulls),
                              component_hull_removed_mm3=float(local_removed.volume()),
                              component_hull_guard_hits_mm3=overlaps(local_removed, guards),
                              component_hull_remaining_components=sum(piece.volume() > 1e-5
                                                                      for piece in local_remaining.decompose()))
        protected_trials = []
        combined = union(cuts)
        for inset in (0., .25, .50, 1., 2.):
            core = (combined if inset == 0. else
                    combined.minkowski_difference(mf.Manifold.cube([2 * inset] * 3, True)))
            cut = core - protected
            remaining = (neutral_tpu - cut).as_original().simplify(1e-7)
            removed = neutral_tpu - remaining
            volumes = sorted((float(piece.volume()) for piece in remaining.decompose()), reverse=True)
            guard_hits = overlaps(removed, guards)
            row = dict(cutter_inset_mm=inset, removed_mm3=float(removed.volume()),
                       guard_hits_mm3=guard_hits, remaining_component_volumes_mm3=volumes)
            if row["removed_mm3"] > .001 and sum(v > 1e-5 for v in volumes) == 1 and not guard_hits:
                _, metadata = export_lattice(remaining, local / f"FJ1468_protected_belly_{inset:g}.obj", points[0])
                _, verts, triangles = indexed_solid(metadata)
                weight, _ = blend((verts - points[0]) @ toward, lower, upper, .1)
                row["mesh"] = metadata
                row["peak_pose_overlap_mm3"] = {}
                for key, clavicle, humerus, frames, worst in peak_cases:
                    a, b = clavicle[:, :3], humerus[:, :3]
                    ta, tb = clavicle[:, 3], humerus[:, 3]
                    moved = verts @ a.T + ta + weight[:, None] * (verts @ (b - a).T + tb - ta)
                    replay = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), triangles))
                    if replay.status() != mf.Error.NoError:
                        raise ValueError(f"Invalid protected TPU after posing {worst['pose_deg']}")
                    row["peak_pose_overlap_mm3"][key] = float(
                        (replay ^ solids[key].transform(frames[owners[key]])).volume())
            protected_trials.append(row)
            print("Protected trial", inset, row["removed_mm3"], len(volumes), row.get("peak_pose_overlap_mm3"), flush=True)
        save(local / "FJ1468_belly_overlap_neutral.json", dict(
            generator=provenance(ROOT / "tools/audit_s1_girdle_flexion_channel.py"),
            assembly_audit=reference, neighbor_channel_reference=inventory_ref,
            blend_window_mm=[lower, upper], zones=zones, protected_combined_trials=protected_trials,
            scope="Only the maximum sampled pose for each named pair. Inverse-mapped intersection vertices "
                  "give sampled distances, not guaranteed cutting clearance. End-blend zones are screening "
                  "regions and do not identify the exact attachment surface. Hulls are independent "
                  "single-pose trials. Neighbor cable guards are baseline route candidates. "
                  "Combined protected trials preserve the guarded stock. Eligible meshes are replayed "
                  "at the two source peak poses for the corresponding pair.",
        ))
        print(zones, flush=True)
        return
    if args.liner_relief or args.liner_hull_relief:
        contact_ref = provenance(local / "FJ1468_girdle_flexion_channel_liner_contact.json")
        contact = load_reference(contact_ref)
        if contact["lined_detour"] != lined_ref:
            raise ValueError("Liner contact study uses different geometry")
        result = screen_liner_relief(model, rigid, solids, owners, points[liner_slice],
                                     path_weights[liner_slice], axis, local, convex=args.liner_hull_relief)
        suffix = "_hull" if args.liner_hull_relief else ""
        save(local / f"FJ1486_adapter_FJ1468_liner_relief{suffix}.json", dict(
            generator=provenance(ROOT / "tools/audit_s1_girdle_flexion_channel.py"),
            liner_contact=contact_ref, convex_cutter=args.liner_hull_relief, cases=result,
            scope="Pose-dependent radial corridors in humerus coordinates. Clearance is verified "
                  "against swept centreline envelopes rather than solving TPU deformation. "
                  "Mechanical strength of remaining stock is not assessed.",
        ))
        return
    if args.adapter_relief:
        detour_ref = provenance(local / "FJ1468_girdle_flexion_channel_detour.json")
        detour_report = load_reference(detour_ref)
        if detour_report["rigid_audit"] != rigid_ref or detour_report["guide_motion"] != guide_ref:
            raise ValueError("Adapter relief inputs changed")
        results = screen_adapter_relief(model, rigid, solids, owners, points, toward, lower, upper, axis,
                                        np.asarray(detour_report["neutral_path_source_mm"]))
        for name, row in results.items():
            candidate = row.pop("candidate")
            wires = row.pop("wires")
            if candidate is not None:
                actual, mesh = export(candidate, local / f"FJ1486_adapter_FJ1468_{name}.obj", model.gh)
                residuals = [float((wire ^ actual).volume()) for wire in wires]
                row.update(mesh=mesh, readback_maximum_overlap_mm3=max(residuals),
                           readback_hit_poses=sum(volume > .001 for volume in residuals),
                           readback_added_mm3=float((actual - solids["elbow:FJ1486_bone_adapter"]).volume()))
            print(name, row, flush=True)
        save(local / "FJ1486_adapter_FJ1468_relief.json", dict(
            generator=provenance(ROOT / "tools/audit_s1_girdle_flexion_channel.py"),
            rigid_audit=rigid_ref, detour=detour_ref, cases=results,
            scope="Sampled cable clearance against the humerus-side adapter and installed extensor bolt. "
                  "Protected volumes are 0.5 mm padded envelopes of nearby mounting hardware and bone. "
                  "Connectivity and guard tests apply to the removed adapter material. "
                  "The bolt and TPU remain fixed; the report is a geometric candidate screen.",
        ))
        return
    if args.detour:
        obstacles = ("elbow:FJ1486_bone_adapter", "elbow:FJ1477_bone_0_bolt")
        targets = bounded({key: solids[key] for key in obstacles})
        if any(owners[key] not in MOBILE for key in obstacles):
            raise ValueError("Local detour requires the baseline humerus-side obstacles")
        seeds = {}
        for motion in rigid["cases"].values():
            for pose in motion["clear_path_deg"][::5] + motion["clear_path_deg"][-1:]:
                seeds[tuple(pose)] = pose
        probes = []
        for angle, *sc, flex in seeds.values():
            clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
            gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
            humerus = compose(scapula, gh)

            def deform(route, clavicle=clavicle, humerus=humerus):
                # Apply the same blend as the full replay, in the obstacle frame.
                weight, _ = blend((route - points[0]) @ toward, lower, upper, .1)
                start = route @ clavicle[:, :3].T + clavicle[:, 3]
                end = route @ humerus[:, :3].T + humerus[:, 3]
                return (start + weight[:, None] * (end - start) - humerus[:, 3]) @ humerus[:, :3]

            probes.append((deform, targets))
        trials = screen_detour(points, stations, probes)
        selected, candidate = trials[0]
        print("Detour screen", len(seeds), "poses", selected, flush=True)
        points = candidate
        path_weights, _ = blend((points - points[0]) @ toward, lower, upper, .1)
        stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    cases = {}
    cache = {}
    for name, motion in rigid["cases"].items():
        if motion["reverse_replay_hits"]:
            raise ValueError(f"Rigid return check failed: {name}")
        rows = []
        for pose in motion["clear_path_deg"]:
            key = tuple(pose)
            if key not in cache:
                angle, *sc, flex = pose
                clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
                gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
                humerus = compose(scapula, gh)
                a, b = clavicle[:, :3], humerus[:, :3]
                ta, tb = clavicle[:, 3], humerus[:, 3]
                delta = vertices @ (b - a).T + tb - ta
                moved = vertices @ a.T + ta + weights[:, None] * delta
                jacobian = a + weights[:, None, None] * (b - a) + derivative[:, None, None] * delta[:, :, None] * toward
                determinant = float(np.linalg.det(jacobian).min())
                route = points @ a.T + ta + path_weights[:, None] * (points @ (b - a).T + tb - ta)
                arc, turns, radii = route_metrics(route)
                wire = A.D.outer_envelope(route, .30)
                frames = {"clavicle": clavicle, "scapula": scapula, **dict.fromkeys(MOBILE, humerus)}
                targets = bounded({part: solid.transform(frames[owners[part]]) if owners[part] in frames else solid
                                   for part, solid in solids.items()})
                hits = overlaps(wire, targets)
                row = dict(pose_deg=pose, arc_mm=arc, arc_change_mm=arc - length, total_turn_rad=float(turns.sum()),
                           minimum_polyline_radius_mm=float(radii.min()), minimum_vertex_jacobian=determinant,
                           wire_rigid_hits_mm3=hits,
                           wire_rigid_locations=locate_wire_hits(route, targets, hits, stations, path_weights))
                if args.neighbor_wires:
                    original_moved = (original_vertices @ a.T + ta + original_weights[:, None]
                                      * (original_vertices @ (b - a).T + tb - ta))
                    before = mf.Manifold(mf.Mesh64(np.ascontiguousarray(original_moved), original_faces))
                    after = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
                    if any(part.status() != mf.Error.NoError for part in (before, after)):
                        raise ValueError(f"Invalid neighbor replay at {pose}")
                    row["neighbor_wires"] = {}
                    w0, t0 = scapula[:, :3], scapula[:, 3]
                    for head, (cable, weight) in neighbors.items():
                        cable_world = cable @ w0.T + t0 + weight[:, None] * (cable @ (b - w0).T + tb - t0)
                        envelope = A.D.outer_envelope(cable_world, .30)
                        row["neighbor_wires"][head] = dict(
                            baseline_TPU_mm3=float((envelope ^ before).volume()),
                            protected_TPU_mm3=float((envelope ^ after).volume()),
                            rigid_hits_mm3=overlaps(envelope, targets),
                            neutral_origin_body="scapula", distal_body="humerus")
                elif not args.locate:
                    shaped = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
                    valid = shaped.status() == mf.Error.NoError and determinant > 0.
                    row.update(valid_interpolation=bool(valid),
                               wire_TPU_overlap_mm3=float((wire ^ shaped).volume()) if valid else None,
                               TPU_rigid_hits_mm3=overlaps(shaped, targets) if valid else None)
                    if check_neighborhood:
                        if not valid:
                            raise ValueError(f"Invalid liner pose: {pose}")
                        corridor = A.D.outer_envelope(route[liner_slice], 1.10)
                        inner = A.D.outer_envelope(route[liner_slice], .65)
                        for label, zone in (("liner_neighborhood", corridor), ("inner_wall_neighborhood", inner)):
                            neighborhood = shaped ^ zone
                            row[f"{label}_hits_mm3"] = overlaps(
                                neighborhood, {key: targets[key] for key in row["TPU_rigid_hits_mm3"]})
                cache[key] = row
            rows.append(cache[key])
        cases[name] = rows
        print(name, len(rows), "arc range", min(row["arc_change_mm"] for row in rows),
              max(row["arc_change_mm"] for row in rows), flush=True)
    suffix = ("_neighbor_wires" if args.neighbor_wires else "_protected_belly" if args.protected_belly else "_liner_assembly" if args.liner_assembly else "_liner_contact" if args.liner_contact else "_lined_detour" if args.lined_detour
              else "_detour" if args.detour else "_locations" if args.locate else "")
    detour = (dict(screened_seed_poses=list(seeds.values()), selected=selected,
                   trials=[row for row, _ in trials]) if args.detour else {})
    save(local / f"FJ1468_girdle_flexion_channel{suffix}.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_girdle_flexion_channel.py"), rigid_audit=rigid_ref,
        guide_motion=guide_ref, cases=cases, unique_poses=len(cache), neutral_arc_mm=length,
        neutral_path_source_mm=points.tolist(), tpu_contact_checked=not args.locate, **detour, **liner,
        scope="FJ1468 local centreline and merged lattice moved by the same nominal axial blend. "
              "Cable radius 0.30 mm. All rigid mesh intersections recorded including baseline mounting contacts. "
              "The return revisits identical interpolated geometries. Contact deformation, tension, "
              "external routing and actuator motion are outside this calculation."
              + (" Detour trials use a 0.40 mm radius on a subset of poses and two obstacles; "
                 "full replay uses 0.30 mm against all rigid parts and unchanged TPU. Arc changes "
                 "are relative to the original neutral route." if args.detour else "")
              + (" Local neighborhoods are current centrelines swept at radii 1.10 and 0.65 mm. "
                 "They classify posed lattice material and are not an exact partition by manufacturing feature."
                 if check_neighborhood and not args.neighbor_wires else "")
              + (" Neighbor wires use stock route centreline geometry with two-end blend weights "
                 "and their own scapular origins. The shared TPU uses the front bundle blend "
                 "from clavicle to humerus. This compares pre-/post-cut geometry; actual multi-origin "
                 "deformation and tension are not solved." if args.neighbor_wires else ""),
    ))


def adapter_guards(solids, owners, adapter_key):
    mf = A.I._manifold()
    return bounded({key: part.minkowski_sum(mf.Manifold.cube([1.] * 3, True))
                    for key, part in solids.items()
                    if (key == "bone:humerus" or key.startswith(("elbow:FJ1486_bone_", "elbow:FJ1486_proximal_")))
                    and key != adapter_key and owners[key] in MOBILE})


def screen_liner_relief(model, rigid, solids, owners, route, weights, axis, output, *, convex=False):
    key = "elbow:FJ1486_bone_adapter"
    if owners[key] not in MOBILE:
        raise ValueError("Adapter must follow the humerus")
    adapter = solids[key]
    guards = adapter_guards(solids, owners, key)
    poses = sorted({tuple(pose) for case in rigid["cases"].values() for pose in case["clear_path_deg"]})
    paths = []
    for angle, *sc, flex in poses:
        clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
        gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
        humerus = compose(scapula, gh)
        start = (route @ clavicle[:, :3].T + clavicle[:, 3] - humerus[:, 3]) @ humerus[:, :3]
        paths.append(start + weights[:, None] * (route - start))
    results = []
    for radius in (.70, .85, 1.00, 1.10):
        envelopes = [A.D.outer_envelope(path, radius) for path in paths]
        cuts = [envelope ^ adapter for envelope in envelopes]
        cuts = [cut for cut in cuts if cut.volume() > 1e-8]
        if convex:
            vertices = np.concatenate([np.asarray(cut.to_mesh64().vert_properties)[:, :3] for cut in cuts])
            cutter = A.I._manifold().Manifold.hull_points(vertices)
        else:
            cutter = union(cuts)
        revised = (adapter - cutter).as_original().simplify(1e-7)
        removed = cutter ^ adapter
        count = len([piece for piece in revised.decompose() if piece.volume() > 1e-5])
        hits = overlaps(removed, guards)
        row = dict(radius_mm=radius, checked_poses=len(poses), components=count,
                   incremental_removed_mm3=float(removed.volume()), guard_hits_mm3=hits)
        if count == 1 and not hits:
            try:
                suffix = "_hull" if convex else ""
                actual, mesh = export(revised, output / f"FJ1486_adapter_liner_relief_{radius:.2f}{suffix}.obj", model.gh)
            except (AssertionError, ValueError) as error:
                row["export_error"] = str(error) or type(error).__name__
            else:
                residuals = [float((envelope ^ actual).volume()) for envelope in envelopes]
                row.update(mesh=mesh, maximum_readback_overlap_mm3=max(residuals),
                           readback_hit_poses=sum(volume > .001 for volume in residuals))
        results.append(row)
        print(radius, count, row["incremental_removed_mm3"], hits,
              row.get("readback_hit_poses", "no readback"), flush=True)
    return results


def screen_adapter_relief(model, rigid, solids, owners, points, toward, lower, upper, axis, detour):
    adapter_key = "elbow:FJ1486_bone_adapter"
    bolt_key = "elbow:FJ1477_bone_0_bolt"
    if any(owners[key] not in MOBILE for key in (adapter_key, bolt_key)):
        raise ValueError("Relief requires both obstacles in the humerus frame")
    adapter, bolt = solids[adapter_key], solids[bolt_key]
    guards = adapter_guards(solids, owners, adapter_key)
    poses = {tuple(pose) for case in rigid["cases"].values() for pose in case["clear_path_deg"]}
    results = {}
    for name, route in (("original", points), ("detour", detour)):
        weights, _ = blend((route - points[0]) @ toward, lower, upper, .1)
        cutters = []
        wires = []
        bolt_hits = 0
        for angle, *sc, flex in sorted(poses):
            clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
            gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
            humerus = compose(scapula, gh)
            start = route @ clavicle[:, :3].T + clavicle[:, 3]
            end = route @ humerus[:, :3].T + humerus[:, 3]
            posed = (start + weights[:, None] * (end - start) - humerus[:, 3]) @ humerus[:, :3]
            # Keep the surveyed distal transition and both adjacent segments.
            local = posed[(posed - points[0]) @ toward > 90.]
            wire = A.D.outer_envelope(local, .40)
            wires.append(wire)
            bolt_hits += float((wire ^ bolt).volume()) > .001
            cut = wire ^ adapter
            if cut.volume() > 1e-8:
                cutters.append(cut)
        removed = union(cutters)
        revised = (adapter - removed).as_original().simplify(1e-7)
        guard_hits = overlaps(removed, guards)
        components = len(revised.decompose())
        results[name] = dict(checked_poses=len(poses), guard_keys=list(guards), radius_mm=.40,
                             removed_mm3=float(removed.volume()), components=components,
                             bolt_hit_poses=int(bolt_hits), guard_hits_mm3=guard_hits,
                             candidate=revised if components == 1 and not guard_hits else None, wires=wires)
    return results


if __name__ == "__main__":
    main()
