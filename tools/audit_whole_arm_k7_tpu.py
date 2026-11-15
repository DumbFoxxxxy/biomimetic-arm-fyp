"""Screen all eight TPU bellies on the checked K7 ten-line trajectory."""

import argparse
from itertools import combinations

import mujoco
import numpy as np
from scipy.spatial import cKDTree

from audit_elbow_slider_short_guides import export_lattice, weigh
from audit_s1_infraspinatus_distributed_deformation import blend
from audit_s1_series_elbow_integration import transforms
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, mesh_geom, rod, save
from build_whole_arm_girdle_tpu import wire_frame
from project_runtime import load_mujoco_model
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def affine(points, frame):
    return points @ frame[:, :3].T + frame[:, 3]


def mapped_contact_cutters(body, obstacles, neutral_vertices, posed_vertices, minimum_volume, padding):
    """Map contact components back with the nearest local vertex displacements."""
    mf = A.I._manifold()
    tree = cKDTree(posed_vertices)
    displacement = posed_vertices - neutral_vertices
    cutters = []
    for obstacle in obstacles.values():
        for component in (body ^ obstacle).decompose():
            if component.volume() <= minimum_volume:
                continue
            points = np.asarray(component.to_mesh64().vert_properties)[:, :3]
            distance, nearby = tree.query(points, k=6)
            weights = 1. / np.maximum(distance, .03)
            estimate = points - np.einsum("ij,ijk->ik", weights, displacement[nearby]) / weights.sum(axis=1)[:, None]
            hull = mf.Manifold.hull_points(np.ascontiguousarray(estimate))
            if hull.volume() > .001:
                cutters.append(hull.minkowski_sum(mf.Manifold.cube([padding] * 3, True)))
    return cutters


def local_route_field(vertices, route):
    """Sample the cable's motion at the nearest station of its neutral path."""
    points = np.asarray(route["points_source_mm"])
    weight = np.asarray(route["distal_weights"])
    distance2 = np.full(len(vertices), np.inf)
    station_weight = np.zeros(len(vertices))
    for index, (start, end) in enumerate(zip(points[:-1], points[1:])):
        span = end - start
        offset = vertices - start
        fraction = np.clip(np.einsum("ij,j->i", offset, span) / (span @ span), 0., 1.)
        delta = offset - fraction[:, None] * span
        squared = np.einsum("ij,ij->i", delta, delta)
        closer = squared < distance2
        distance2[closer] = squared[closer]
        station_weight[closer] = weight[index] + fraction[closer] * (weight[index + 1] - weight[index])
    falloff = 1. - blend(np.sqrt(distance2), 1.25, 8., 0.)[0]
    return station_weight[:, None], falloff


def reconnect_lattice(candidate, targets, bores):
    """Reconnect cut lattice struts with short braces clear of the line bores."""
    mf = A.I._manifold()
    targets = bounded(targets)
    repairs = []
    for _ in range(12):
        pieces = sorted((p for p in candidate.decompose() if p.volume() > 1e-5),
                        key=lambda p: -p.volume())
        if len(pieces) < 2:
            break
        major = pieces[0]
        main_v = np.asarray(major.to_mesh64().vert_properties)[:, :3]
        tree = cKDTree(main_v)
        options = []
        for island in pieces[1:]:
            vertices = np.asarray(island.to_mesh64().vert_properties)[:, :3]
            distances, nearest = tree.query(vertices, k=12)
            used = set()
            for rank in np.argsort(distances.ravel()):
                i, j = np.unravel_index(rank, distances.shape)
                distance = distances[i, j]
                if distance > 14. or len(used) == 72:
                    break
                if distance <= .05:
                    continue
                a, b = vertices[i], main_v[nearest[i, j]]
                # Avoid using the whole trial budget on near-duplicate points
                # along one face of a cut strut.
                voxel = tuple(np.floor(np.r_[a, b] / 1.5).astype(int))
                if voxel not in used:
                    used.add(voxel)
                    options.append((distance, a, b, island))
        for _, a, b, island in sorted(options, key=lambda option: option[0]):
            for radius in (.65, .45):
                brace = union([rod(a, b, radius, 16),
                               mf.Manifold.sphere(radius, 12).translate(a),
                               mf.Manifold.sphere(radius, 12).translate(b)])
                if overlaps(brace, targets) or A.I.intersection_volume(brace, bores) > .001:
                    continue
                if min((brace ^ major).volume(), (brace ^ island).volume()) < .01:
                    continue
                candidate = (candidate + brace).as_original().simplify(1e-6)
                repairs.append(dict(start_source_mm=a.tolist(), end_source_mm=b.tolist(), radius_mm=radius))
                break
            else:
                continue
            break
        else:
            # A pair of straight segments can turn around an obstacle when
            # every direct connection crosses the elbow package.
            for _, a, b, island in sorted(options, key=lambda option: option[0])[:24]:
                direction = (b - a) / np.linalg.norm(b - a)
                side = np.cross(direction, np.eye(3)[np.argmin(abs(direction))])
                side /= np.linalg.norm(side)
                up = np.cross(direction, side)
                radius = .55
                found = False
                for offset in (2., 4., 6., 9.):
                    for theta in np.arange(8) * np.pi / 4.:
                        mid = (a + b) / 2 + offset * (np.cos(theta) * side + np.sin(theta) * up)
                        brace = A.D.outer_envelope(np.array([a, mid, b]), radius)
                        if overlaps(brace, targets) or A.I.intersection_volume(brace, bores) > .001:
                            continue
                        if min((brace ^ major).volume(), (brace ^ island).volume()) < .01:
                            continue
                        candidate = (candidate + brace).as_original().simplify(1e-6)
                        repairs.append(dict(start_source_mm=a.tolist(), bend_source_mm=mid.tolist(),
                                            end_source_mm=b.tolist(), radius_mm=radius))
                        found = True
                        break
                    if found:
                        break
                if found:
                    break
            else:
                break
    return candidate, repairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    variants = parser.add_mutually_exclusive_group()
    variants.add_argument("--deltoid-relief", action="store_true",
                        help="Test a guarded neutral relief against the shared elbow")
    variants.add_argument("--reconnect-deltoid", action="store_true")
    variants.add_argument("--track-deltoid-channels", action="store_true")
    variants.add_argument("--repair-deltoid-bore", action="store_true")
    variants.add_argument("--sweep-deltoid-relief", action="store_true",
                          help="Test a local swept relief against the elbow over K7")
    variants.add_argument("--verify-swept-deltoid", action="store_true")
    variants.add_argument("--track-teres-channel", action="store_true",
                          help="Test local alignment with the FJ1507 route on the swept-relief assembly")
    variants.add_argument("--relieve-teres-ends", action="store_true")
    variants.add_argument("--check-moving-cuffs", action="store_true")
    parser.add_argument("--render-replay", action="store_true",
                        help="Render the existing moving-cuff audit with its exact interpolation")
    parser.add_argument("--k1-spring-relief", action="store_true",
                        help="Check the refined K1 spring clearance parts in K7")
    parser.add_argument("--sweep-elbow-min-volume", type=float, default=.005,
                        help="Minimum per-component contact used to build a swept-relief cutter")
    args = parser.parse_args()
    if args.render_replay and not args.check_moving_cuffs:
        parser.error("--render-replay requires --check-moving-cuffs")
    if args.k1_spring_relief and not args.check_moving_cuffs:
        parser.error("--k1-spring-relief requires --check-moving-cuffs")
    folder = DEST.parent / "whole_arm_first_pass"
    combo_ref = provenance(folder / "whole_arm_combined_cuff_k7.json")
    combo = load_reference(combo_ref)
    motion = load_reference(combo["motion"])
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    routed = load_reference(combo["routed"])["routes"]
    support = load_reference(combo["support"])
    revision = load_reference(load_reference(support["audit"])["humeral_revision"])["selected"]["meshes"]
    liner_ref = provenance(folder / "whole_arm_posterior_liner.json")
    liner = load_reference(liner_ref)
    if liner["audit"] != combo["routed"]:
        raise ValueError("Deltoid liner belongs to another wire path")
    girdle_ref = provenance(folder / "girdle_tpu_geometry.json")
    girdle = load_reference(girdle_ref)
    elbow_ref = provenance(folder / "whole_arm_elbow_flexible_audit.json")
    elbow = load_reference(elbow_ref)
    base_ref = provenance(folder / "tabletop_base_check.json")
    base = load_reference(base_ref)
    if girdle["motion"] != combo["motion"] or base["motion"] != combo["motion"]:
        raise ValueError("Girdle or base evidence does not match K7")
    scene_sources = load_reference(scene["neutral_tpu_reference"])["flexible_geometry_sources"]
    tpu_sources = {**scene_sources, "tpu:deltoid": liner["source"],
                   **{f"tpu:elbow_{key}": source for key, source in elbow["source_meshes"].items()},
                   **{key: row["source"] for key, row in girdle["muscles"].items()}}
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
               "elbow:FJ1477_bone_adapter": support["source"], "base:T12": base["source"]}
    sources.update({key: meta for key, meta in revision.items()
                    if key.startswith("link:") and key in sources})
    spring_trial_ref = None
    if args.k1_spring_relief:
        spring_trial_ref = provenance(folder / "whole_arm_k1_spring_relief_refined.json")
        spring_trial = load_reference(spring_trial_ref)
        k1_path = load_reference(spring_trial["path"])
        if k1_path["combined"] != combo_ref or not all(row["source"] for row in spring_trial["trial"].values()):
            raise ValueError("The spring clearance inputs do not match the K7 assembly")
        spring_definition = load_reference(spring_trial["spring_definition"])
        spring_coil = load_component(spring_definition["exports"]["neutral"]["mesh"])
        for key in ("elbow:FJ1477_bone_adapter", "elbow:FJ1477_proximal_rigid_cap"):
            sources[key] = spring_trial["trial"][key]["source"]
    rigid = {} if args.render_replay else {key: load_component(source) for key, source in sources.items()}
    owners = {**transfer["owner_by_source"], "base:T12": "thorax"}
    origins = {key: np.asarray(scene["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    mf = A.I._manifold()
    sections = {}
    neutral_design = None
    tracked_ref = None
    track_teres = args.track_teres_channel or args.relieve_teres_ends or args.check_moving_cuffs
    use_swept = args.verify_swept_deltoid or track_teres
    track_channels = (args.track_deltoid_channels or args.repair_deltoid_bore
                      or args.sweep_deltoid_relief or use_swept)
    if track_channels:
        tracked_ref = provenance(folder / "whole_arm_k7_deltoid_reconnected_audit.json")
        tracked = load_reference(tracked_ref)
        if tracked["combined"] != combo_ref:
            raise ValueError("Reconnected mesh belongs to another K7 trajectory")
        tpu_sources["tpu:deltoid"] = tracked["tpu_sources"]["tpu:deltoid"]
    if args.sweep_deltoid_relief or use_swept:
        bore_ref = provenance(folder / "whole_arm_k7_deltoid_bore_audit.json")
        bore_audit = load_reference(bore_ref)
        if bore_audit["combined"] != combo_ref:
            raise ValueError("The repaired bore belongs to another assembled motion")
        tpu_sources["tpu:deltoid"] = bore_audit["tpu_sources"]["tpu:deltoid"]
    swept_ref = None
    if use_swept:
        swept_ref = provenance(folder / "whole_arm_k7_deltoid_swept_relief.json")
        swept = load_reference(swept_ref)
        if swept["baseline"] != bore_ref or not swept["source"]:
            raise ValueError("No connected swept-relief mesh for the current assembly")
        tpu_sources["tpu:deltoid"] = swept["source"]
    teres_ref = None
    teres_design = None
    cuff_ref = None
    if args.check_moving_cuffs:
        teres_ref = provenance(folder / "whole_arm_k7_teres_end_relief.json")
        end_relief = load_reference(teres_ref)
        cuff_ref = provenance(folder / "whole_arm_k7_moving_cuff_geometry.json")
        cuff_motion = load_reference(cuff_ref)
        if end_relief["combined"] != combo_ref or cuff_motion["combined"] != combo_ref:
            raise ValueError("TPU and cuff hardware must use the same checked trajectory")
        tpu_sources = end_relief["tpu_sources"]
        if spring_trial_ref:
            tpu_sources["tpu:deltoid"] = spring_trial["trial"]["tpu:deltoid"]["source"]
        if not args.render_replay:
            moving_sliders = {key: load_component(source) for key, source in cuff_motion["moving_sources"].items()}
    if args.render_replay:
        from render_whole_arm_k7_replay import K7Replay

        replay = K7Replay(folder, combo_ref, tpu_sources, cuff_ref, cuff_motion, spring_trial_ref)
    if args.relieve_teres_ends:
        teres_ref = provenance(folder / "whole_arm_k7_teres_channel_tracking.json")
        baseline = load_reference(teres_ref)
        if baseline["combined"] != combo_ref or baseline["tpu_sources"] != tpu_sources:
            raise ValueError("Teres relief input does not match the assembly")
    if args.repair_deltoid_bore:
        key = "tpu:deltoid"
        mesh, _, _ = indexed_solid(tpu_sources[key])
        points = np.asarray(routed["FJ1467"]["points_source_mm"])
        wire = A.D.outer_envelope(points, .30)
        intrusion = mesh ^ wire
        patch = intrusion.minkowski_sum(mf.Manifold.cube([.65] * 3, True))
        cutter = A.D.outer_envelope(points, .52) ^ patch
        candidate = (mesh - cutter).as_original().simplify(1e-7)
        removed = mesh - candidate
        guards = {name: A.D.outer_envelope(np.asarray(routed[name]["points_source_mm"]), .95)
                  for name in ("FJ1468", "FJ1513")}
        affected = {name: A.I.intersection_volume(removed, guard) for name, guard in guards.items()}
        residual = A.I.intersection_volume(candidate, wire)
        if residual > .001 or max(affected.values()) > .001:
            raise ValueError(f"Bore repair leaves wire overlap {residual} or alters neighbor guards {affected}")
        origin = np.asarray(tpu_sources[key]["translation_source_mm"])
        actual, tpu_sources[key] = export_lattice(candidate, folder / "K7_deltoid_reconnected_bore.obj", origin)
        neutral_design = dict(original_overlap_mm3=float(intrusion.volume()),
                              reloaded_overlap_mm3=A.I.intersection_volume(actual, wire),
                              removed_mm3=float(removed.volume()), neighbor_material_loss_mm3=affected,
                              cutter_radius_mm=.52,
                              affected_bounds_source_mm=list(intrusion.bounding_box()))
    for key, source in tpu_sources.items():
        if source["file"].lower().endswith(".obj"):
            solid, vertices, faces = indexed_solid(source)
        else:
            solid = load_component(source)
            vertices, faces = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
        if solid.status() != mf.Error.NoError:
            raise ValueError(f"Invalid neutral TPU geometry: {key}")
        sections[key] = (vertices, np.array(faces, dtype=np.uint64, order="C", copy=True))
    if args.deltoid_relief or args.reconnect_deltoid:
        key = "tpu:deltoid"
        original, _, _ = indexed_solid(tpu_sources[key])
        extensor, _, _ = indexed_solid(tpu_sources["tpu:elbow_extensor"])
        targets = {name: shape for name, shape in rigid.items() if name.startswith(("series:", "elbow:"))}
        targets["tpu:elbow_extensor"] = extensor
        exact = overlaps(original, bounded(targets))
        cutters = union((original ^ targets[name]).minkowski_sum(mf.Manifold.cube([.5]*3, True))
                        for name in exact)
        protected = union(A.D.outer_envelope(np.asarray(routed[name]["points_source_mm"]), .95)
                          for name in ("FJ1468", "FJ1467", "FJ1513"))
        candidate = (original - (cutters - protected)).as_original().simplify(1e-6)
        repairs = []
        raw_components = sum(part.volume() > 1e-5 for part in candidate.decompose())
        if args.reconnect_deltoid:
            bores = union(A.D.outer_envelope(np.asarray(routed[name]["points_source_mm"]), .52)
                          for name in ("FJ1468", "FJ1467", "FJ1513"))
            candidate, repairs = reconnect_lattice(candidate, targets, bores)
        discarded = []
        if args.reconnect_deltoid:
            pieces = sorted((piece for piece in candidate.decompose() if piece.volume() > 1e-5),
                            key=lambda piece: -piece.volume())
            if len(pieces) > 1:
                scraps = pieces[1:]
                if (sum(piece.volume() for piece in scraps) < .0002 * pieces[0].volume()
                        and all(A.I.intersection_volume(piece, protected) <= .001 for piece in scraps)):
                    # Intersect with the dominant stock to preserve any internal void shells.
                    candidate = (candidate ^ pieces[0]).as_original().simplify(1e-6)
                    discarded = [float(piece.volume()) for piece in scraps]
            # A protected outer annulus can protrude slightly. Keep the bore
            # and a 0.2 mm radial wall at the local trimming site.
            inner_guard = union(A.D.outer_envelope(np.asarray(routed[name]["points_source_mm"]), .72)
                                for name in ("FJ1468", "FJ1467", "FJ1513"))
            for name in overlaps(candidate, bounded(targets)):
                patch = (candidate ^ targets[name]).minkowski_sum(mf.Manifold.cube([.12]*3, True))
                candidate = (candidate - (patch - inner_guard)).as_original().simplify(1e-6)
        connected = sum(part.volume() > 1e-5 for part in candidate.decompose())
        residual = overlaps(candidate, bounded(targets))
        positive_parts = sorted((p for p in candidate.decompose() if p.volume() > 1e-5),
                                key=lambda p: -p.volume())
        neutral_design = dict(overlap_before_mm3=exact, overlap_after_mm3=residual,
                              positive_components=connected,
                              components_before_bracing=raw_components, repairs=repairs,
                              discarded_off_channel_fragments_mm3=discarded,
                              component_volumes_mm3=[float(p.volume()) for p in positive_parts],
                              detached_parts=[dict(volume_mm3=float(p.volume()),
                                                   bounds_source_mm=list(p.bounding_box()),
                                                   wire_guard_overlap_mm3=A.I.intersection_volume(p, protected))
                                              for p in positive_parts[1:]],
                              residual_inside_wire_guards_mm3={
                                  name: A.I.intersection_volume(candidate ^ targets[name], protected)
                                  for name in residual},
                              removed_mm3=float((original - candidate).volume()),
                              removed_in_wire_guards_mm3=A.I.intersection_volume(original - candidate, protected))
        guarded_loss = (A.I.intersection_volume(original - candidate, inner_guard) if args.reconnect_deltoid
                        else neutral_design["removed_in_wire_guards_mm3"])
        neutral_design["removed_in_minimum_wall_guard_mm3"] = guarded_loss
        if connected != 1 or residual or guarded_loss > .001:
            filename = ("whole_arm_k7_deltoid_reconnection_trial.json" if args.reconnect_deltoid else
                        "whole_arm_k7_deltoid_relief_trial.json")
            save(folder / filename, dict(
                baseline=provenance(folder / "whole_arm_k7_tpu_audit.json"), design=neutral_design,
                adopted=False, scope="Neutral cut trial; preserve three channel guard volumes."))
            print("Neutral relief rejected", neutral_design, flush=True)
            return
        source_origin = np.asarray(tpu_sources[key]["translation_source_mm"])
        name = "K7_deltoid_reconnected.obj" if args.reconnect_deltoid else "K7_deltoid_elbow_relief.obj"
        _, tpu_sources[key] = export_lattice(candidate, folder / name, source_origin)
        _, vertices, faces = indexed_solid(tpu_sources[key])
        sections[key] = vertices, faces

    front = np.asarray(routed["FJ1468"]["points_source_mm"])
    axis = front[-1] - front[0]
    axis /= np.linalg.norm(axis)
    window = load_reference(routed["FJ1468"]["guide_motion_reference"])["deformation_window_mm"]
    deltoid = sections["tpu:deltoid"][0]
    deltoid_weights = blend((deltoid - front[0]) @ axis, *window, .1)[0]
    if track_channels:
        fields = {key: local_route_field(deltoid, routed[key])
                  for key in ("FJ1468", "FJ1467", "FJ1513")}
        total = sum(influence for _, influence in fields.values())
        # Mix overlapping local neighborhoods without amplifying displacement.
        fields = {key: (weight, (influence / np.maximum(total, 1.))[:, None])
                  for key, (weight, influence) in fields.items()}
    if track_teres:
        teres_weights, teres_influence = local_route_field(sections["tpu:teres_major"][0], routed["FJ1507"])

    model_file = ROOT / scene["model"]["file"]
    if provenance(model_file) != scene["model"]:
        raise ValueError("Shared model changed")
    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    poses = motion["paths"]["K7_SC1_pos"]["accepted"]
    rows = []
    relief_hulls = []
    if args.relieve_teres_ends:
        key = "tpu:teres_major"
        first, last = np.asarray(routed["FJ1507"]["points_source_mm"])[[0, -1]]
        span = last - first
        vertices, faces = sections[key]
        weight = weigh(vertices, first, span / (span @ span))[:, None]
        weight += teres_influence[:, None] * (teres_weights - weight)
        guard = A.D.outer_envelope(np.asarray(routed["FJ1507"]["points_source_mm"]), .85)
        training = (0, 4, 8, len(combo["combined"]) - 1)
        cuts = []
        for index in training:
            data.qpos[:] = poses[index]["qpos"]
            mujoco.mj_forward(model, data)
            frames = transforms(model, data, origins)
            moved = (1. - weight) * affine(vertices, frames["scapula"]) + weight * affine(vertices, frames["humerus"])
            body = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
            if body.status() != mf.Error.NoError:
                raise ValueError("Invalid teres major training shape")
            names = ("bone:scapula", "cap:FJ1507/origin", "cap:FJ1507/distal")
            targets = {name: rigid[name].transform(frames[owners[name]]) for name in names}
            cuts.extend(mapped_contact_cutters(body, targets, vertices, moved, .001, .32))
        original = load_component(tpu_sources[key])
        candidate = (original - (union(cuts) - guard)).as_original().simplify(1e-6)
        lost = original - candidate
        protected_loss = A.I.intersection_volume(lost, guard)
        if protected_loss > .001:
            raise ValueError("Teres relief alters the channel protection volume")
        origin = np.asarray(tpu_sources[key]["translation_source_mm"])
        _, tpu_sources[key] = export_lattice(candidate, folder / "K7_teres_end_relief.obj", origin)
        _, vertices, faces = indexed_solid(tpu_sources[key])
        sections[key] = vertices, faces
        teres_weights, teres_influence = local_route_field(vertices, routed["FJ1507"])
        teres_design = dict(training_pose_indices=training, cutters=len(cuts), removed_mm3=float(lost.volume()),
                            guard_radius_mm=.85, protected_loss_mm3=protected_loss)
    for index, case in enumerate(combo["combined"]):
        pose = poses[index]
        if pose["input_deg"] != case["input_deg"]:
            raise ValueError("Mismatched K7 pose")
        data.qpos[:] = pose["qpos"]
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        obstacles = bounded({key: part.transform(frames[owners[key]]) for key, part in rigid.items()})
        wires, paths = {}, {}
        for key, route in routed.items():
            points = np.asarray(route["points_source_mm"])
            a, b = (frames[owner] for owner in route["owners"])
            w = np.asarray(route["distal_weights"])[:, None]
            moved = (1. - w) * affine(points, a) + w * affine(points, b)
            paths[key] = moved
            if not args.render_replay:
                wires[key] = rod(*moved, .30) if len(moved) == 2 else A.D.outer_envelope(moved, .30)
        high = next(step for step in case["steps"] if step["load_N"] == 10.)
        for key, points in high["central_tendon_tips_world_mm"].items():
            paths[key] = np.asarray(points)
            if not args.render_replay:
                wires[key] = A.D.outer_envelope(paths[key], .30)
        tpu = {}
        invalid = []
        render_meshes = {}
        posed_deltoid_vertices = None
        for key, (vertices, faces) in sections.items():
            if key == "tpu:deltoid":
                w = deltoid_weights[:, None]
                moved = (1. - w) * affine(vertices, frames["clavicle"]) + w * affine(vertices, frames["humerus"])
                if track_channels:
                    base_motion = moved.copy()
                    for head, (weight, influence) in fields.items():
                        a, b = (frames[owner] for owner in routed[head]["owners"])
                        guide_motion = (1. - weight) * affine(vertices, a) + weight * affine(vertices, b)
                        moved += influence * (guide_motion - base_motion)
                posed_deltoid_vertices = moved
            elif key in ("tpu:infraspinatus", "tpu:subscapularis"):
                head = "FJ1500" if key == "tpu:infraspinatus" else "FJ1504"
                state = next(row for row in combo["cases"][head]["K7_SC1_pos"][index]["states"]
                             if row["command_N"] == 10.)
                start = np.asarray(state["origin_port_source_mm"])
                delta = np.asarray(state["cap_pin_neutral_source_mm"]) - start
                w = weigh(vertices, start, delta / (delta @ delta))[:, None]
                moved = affine((1. - w) * vertices + w * affine(vertices, np.asarray(state["cap_affine_source_mm"])),
                               frames["scapula"])
            elif key in girdle["muscles"]:
                first, last = np.asarray(routed[key]["points_source_mm"])
                rest, length = wire_frame(first, last, np.array([0., 1., 0.]))
                posed, span = wire_frame(*paths[key], rest[:, 0])
                stretch = (span - 4.) / (length - 4.)
                local = (vertices - first) @ rest
                local[:, 2] = 2. + (local[:, 2] - 2.) * stretch
                moved = local @ posed.T + paths[key][0]
            else:
                head = {"tpu:elbow_flexor": "FJ1486", "tpu:elbow_extensor": "FJ1477",
                        "tpu:teres_major": "FJ1507"}[key]
                route = routed[head]
                first, last = np.asarray(route["points_source_mm"])[[0, -1]]
                span = last - first
                w = weigh(vertices, first, span / (span @ span))[:, None]
                if track_teres and head == "FJ1507":
                    w += teres_influence[:, None] * (teres_weights - w)
                a, b = (frames[owner] for owner in route["owners"])
                moved = (1. - w) * affine(vertices, a) + w * affine(vertices, b)
            if args.render_replay:
                render_meshes[key] = moved, faces
                continue
            solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
            if solid.status() != mf.Error.NoError or solid.volume() <= 0.:
                invalid.append(key)
            else:
                tpu[key] = solid
        if args.render_replay:
            replay.render(index, pose, render_meshes, paths)
            continue
        wire_targets = bounded(wires)
        findings = {key: dict(rigid_hits_mm3=overlaps(part, obstacles), wire_hits_mm3=overlaps(part, wire_targets))
                    for key, part in tpu.items()}
        if args.check_moving_cuffs:
            state = cuff_motion["frames"][index]
            if state["input_deg"] != pose["input_deg"] or state["matrices"].keys() != moving_sliders.keys():
                raise ValueError("Moving cuff assemblies and TPU pose are out of sync")
            moving = bounded({key: part.transform(np.asarray(state["matrices"][key]))
                              for key, part in moving_sliders.items()})
            for key, part in tpu.items():
                findings[key]["moving_cuff_hits_mm3"] = overlaps(part, moving)
        pair_hits = [dict(left=a, right=b, overlap_mm3=volume)
                     for (a, x), (b, y) in combinations(tpu.items(), 2)
                     if (volume := A.I.intersection_volume(x, y)) > .001]
        row = dict(input_deg=pose["input_deg"], humerus_deg=case["humerus_deg"], muscles=findings,
                   tpu_pair_hits=pair_hits, invalid_meshes=invalid)
        if spring_trial_ref:
            slide_index = model.joint("FJ1477_series_travel").qposadr.item()
            if abs(data.qpos[slide_index]) > 1e-10:
                raise ValueError("K7 spring check requires the carriage at its neutral position")
            coil = spring_coil.transform(frames["humerus"])
            row["spring_contacts_mm3"] = dict(rigid=overlaps(coil, obstacles), wire=overlaps(coil, wire_targets),
                                             tpu=overlaps(coil, bounded(tpu)), moving_cuff=overlaps(coil, moving))
        rows.append(row)
        if args.sweep_deltoid_relief and index in (3, 6, 9, 13):
            torso_mesh = tpu["tpu:deltoid"]
            neighbors = {"tpu:elbow_extensor": tpu["tpu:elbow_extensor"],
                         **{name: obstacles[name][0] for name in findings["tpu:deltoid"]["rigid_hits_mm3"]
                            if name.startswith(("series:", "elbow:"))}}
            relief_hulls.extend(mapped_contact_cutters(torso_mesh, neighbors, deltoid, posed_deltoid_vertices,
                                                       args.sweep_elbow_min_volume, .36))
        print(pose["input_deg"], "muscles", len(tpu), "TPU pairs", len(pair_hits), flush=True)
    if args.render_replay:
        replay.finish()
        return
    if args.sweep_deltoid_relief:
        if not relief_hulls:
            raise ValueError("No swept contact samples to map to the neutral mesh")
        original, _, _ = indexed_solid(tpu_sources["tpu:deltoid"])
        guard = union(A.D.outer_envelope(np.asarray(routed[key]["points_source_mm"]), .72)
                      for key in ("FJ1468", "FJ1467", "FJ1513"))
        sweep = union(relief_hulls)
        candidate = (original - (sweep - guard)).as_original().simplify(1e-6)
        parts = sorted((p for p in candidate.decompose() if p.volume() > 1e-5), key=lambda p: -p.volume())
        discarded = []
        # These fragments are disconnected by the swept cut. Keep their
        # removal explicit and reject any loss within the channel walls.
        if len(parts) > 1 and (sum(p.volume() for p in parts[1:]) < .002 * parts[0].volume()
                              and all(A.I.intersection_volume(p, guard) <= .001 for p in parts[1:])):
            discarded = [float(p.volume()) for p in parts[1:]]
            candidate = (candidate ^ parts[0]).as_original().simplify(1e-6)
        lost = original - candidate
        protected_loss = A.I.intersection_volume(lost, guard)
        components = [float(p.volume()) for p in candidate.decompose() if p.volume() > 1e-5]
        source = None
        if len(components) == 1 and protected_loss <= .001:
            _, source = export_lattice(candidate, folder / "K7_deltoid_swept_relief.obj",
                                       np.asarray(tpu_sources["tpu:deltoid"]["translation_source_mm"]))
        save(folder / "whole_arm_k7_deltoid_swept_relief.json", dict(
            baseline=bore_ref, combined=combo_ref, source=source, hull_count=len(relief_hulls),
            training_pose_indices=[3, 6, 9, 13], removed_mm3=float(lost.volume()),
            minimum_contact_volume_mm3=args.sweep_elbow_min_volume,
            protected_radius_mm=.72, protected_loss_mm3=protected_loss,
            component_volumes_mm3=components, discarded_off_channel_fragments_mm3=discarded,
            scope="Approximate inverse map using six neighboring mesh vertices. Candidate generation; verify full K7 replay."))
        print("Swept relief:", len(components), "components;", lost.volume(), "mm3 removed", flush=True)
        return
    output = ("whole_arm_k7_spring_relief.json" if spring_trial_ref else
              "whole_arm_k7_tpu_moving_cuff_audit.json" if args.check_moving_cuffs else
              "whole_arm_k7_teres_end_relief.json" if args.relieve_teres_ends else
              "whole_arm_k7_teres_channel_tracking.json" if args.track_teres_channel else
              "whole_arm_k7_deltoid_swept_audit.json" if args.verify_swept_deltoid else
              "whole_arm_k7_deltoid_bore_audit.json" if args.repair_deltoid_bore else
              "whole_arm_k7_deltoid_channel_tracking.json" if args.track_deltoid_channels else
              "whole_arm_k7_deltoid_reconnected_audit.json" if args.reconnect_deltoid else
              "whole_arm_k7_deltoid_relief_audit.json" if args.deltoid_relief else "whole_arm_k7_tpu_audit.json")
    save(folder / output, dict(
        combined=combo_ref, liner=liner_ref, girdle=girdle_ref, elbow=elbow_ref, base=base_ref,
        tpu_sources=tpu_sources, samples=rows, neutral_relief=neutral_design, tracking_reference=tracked_ref,
        swept_reference=swept_ref,
        teres_tracking_reference=teres_ref if not args.check_moving_cuffs else None,
        teres_end_relief=teres_design, moving_cuff_reference=cuff_ref,
        latest_mesh_reference=teres_ref if args.check_moving_cuffs else None,
        spring_relief_reference=spring_trial_ref,
        generator=provenance(ROOT / "tools/audit_whole_arm_k7_tpu.py"),
        scope="Prescribed affine endpoint interpolation of existing and first-pass TPU geometry, "
              "on the same K7 poses at the 10 N cuff probe state. Contacts are geometric diagnostics, "
              "not a solution of deformation or force. "
              + ("Moving cuff hardware contacts recorded separately." if args.check_moving_cuffs
                 else "The rigid list has fixed cuff hardware only.")
              + (" The neutral spring is included with the K1 clearance parts." if spring_trial_ref
                 else " Spring geometry is outside this comparison.")))


if __name__ == "__main__":
    main()
