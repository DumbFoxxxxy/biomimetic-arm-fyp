"""Screen FJ1467 proximal paths and matching support clearance."""

import argparse

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial.transform import Rotation

from audit_elbow_slider_short_guides import export_lattice
from audit_s1_ac_blind_coupled_motion import load_assembly
from audit_s1_infraspinatus_distributed_deformation import blend
from audit_s1_series_elbow_integration import MOBILE, compose
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, export, load_component, save, shoulder_pads
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_series_package import bounded, overlaps
from study_s1_coupled_scapula import ShoulderAssembly, affine_at
from study_s1_infraspinatus_wire_travel import route_metrics
from study_s1_local_force_geometry import load_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--smooth-port", action="store_true", help="Fit a smooth transition and test a matching bore")
    mode.add_argument("--distal-detour", action="store_true", help="Screen a smooth detour around the elbow cap")
    mode.add_argument("--distal-guide", action="store_true", help="Build and replay a TPU liner for the checked distal detour")
    mode.add_argument("--distal-guide-motion", action="store_true",
                      help="Replay the distal liner with local scapular support")
    mode.add_argument("--distal-guide-relief", action="store_true", help="Trim the liner where it meets the rail")
    args = parser.parse_args()
    local = DEST.parent / "s1_local_muscle_load"
    mount_ref = provenance(local / "FJ1467_stepped_mount.json")
    mount = load_reference(mount_ref)
    stock = load_reference(mount["stock_study"])
    locations = load_reference(stock["route"])
    neighbor = load_reference(locations["neighbor_replay"])
    rigid = load_reference(neighbor["rigid_audit"])
    girdle = load_reference(rigid["girdle_path"])
    endpoints = load_reference(load_reference(rigid["span_study"])["local_endpoints"])
    guide = load_reference(neighbor["guide_motion"])
    assembly, _, sources = load_assembly(local / "FJ1500_cap_relief")
    sources["elbow:FJ1486_bone_adapter"] = neighbor["adapter_mesh"]
    owners = {key: assembly["owners"].get(key, "scapula") for key in sources}
    owners.update(dict.fromkeys(mount["meshes"], "scapula"))
    solids = {key: load_component(meta) for key, meta in {**sources, **mount["meshes"]}.items()}
    if "support" not in mount["meshes"]:
        raise ValueError("Stepped mount is missing the readback support")
    model = ShoulderAssembly({}, {}, load_reference(girdle["pose_reference"])["neutral_origins_m"])
    axis = np.asarray(endpoints["axes_source"][0])
    path = np.asarray(locations["neutral_source_route_mm"])
    front = np.asarray(neighbor["neutral_path_source_mm"])
    toward = front[-1] - front[0]
    toward /= np.linalg.norm(toward)
    window = guide["deformation_window_mm"]
    poses = []
    for sample in locations["cases"]:
        angle, *sc, flex = sample["pose_deg"]
        clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
        gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
        poses.append((sample["pose_deg"], clavicle, scapula, compose(scapula, gh)))
    if args.distal_guide or args.distal_guide_motion or args.distal_guide_relief:
        audit_distal_guide(local, mount_ref, neighbor, solids, owners, poses, front, toward, window,
                           local_support=args.distal_guide_motion or args.distal_guide_relief,
                           relief=args.distal_guide_relief)
        return
    start, stop = 1, 25
    proximal_ref = None
    if args.distal_detour:
        proximal_ref = provenance(local / "FJ1467_smooth_port_exit_route.json")
        proximal = load_reference(proximal_ref)
        if proximal["mount"] != mount_ref or proximal["matching_port"]["envelope_hit_poses"]:
            raise ValueError("Expected checked proximal support")
        path = np.asarray(proximal["neutral_candidate_mm"])
        solids["support"] = load_component(proximal["matching_port"]["mesh"])
        start, stop = 66, 74
    local_base = path[start:stop + 1]
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(local_base, axis=0), axis=1))]
    fraction = arc / arc[-1]
    taper = np.sin(np.pi * fraction)**2
    taper[[0, -1]] = 0.
    normal = np.asarray(stock["normal_source"])
    basis = shoulder_pads.pad_basis(normal)
    local_targets = bounded({key: solids[key] for key in mount["meshes"]})
    azimuth_step = 45
    sizes = (0., .5, 1., 1.5, 2., 3., 4.)
    if args.distal_detour:
        # Screen in humeral coordinates because the interfering cap follows it.
        local_targets = bounded({key: solid for key, solid in solids.items() if owners[key] in MOBILE})
        tangent = path[stop] - path[start]
        tangent /= np.linalg.norm(tangent)
        normal = np.cross(tangent, np.eye(3)[np.argmin(abs(tangent))])
        normal /= np.linalg.norm(normal)
        basis = np.column_stack((np.cross(tangent, normal), tangent, normal))
        azimuth_step = 45
        sizes = (0., 1., 2., 3., 4., 5., 6.)
    if args.smooth_port or args.distal_detour:
        if args.smooth_port:
            del local_targets["support"]
        tangent = np.array([path[start] - path[start - 1], path[stop + 1] - path[stop]])
        tangent /= np.linalg.norm(tangent, axis=1)[:, None]
        curve = CubicHermiteSpline([0., 1.], local_base[[0, -1]], tangent * arc[-1])
        local_base = curve(fraction)
        taper = 16. * fraction**2 * (1. - fraction)**2
    seeds = [poses[i] for i in sorted({0, len(poses) - 1, *range(0, len(poses), 15)})]

    def pose_route(route, scapula, humerus, weight):
        start = route @ scapula[:, :3].T + scapula[:, 3]
        finish = route @ humerus[:, :3].T + humerus[:, 3]
        return start + weight[:, None] * (finish - start)

    trials = []
    for size in sizes:
        for azimuth in (range(0, 360, azimuth_step) if size else (0,)):
            theta = np.deg2rad(azimuth)
            direction = normal * np.cos(theta) + basis[:, 0] * np.sin(theta)
            candidate = path.copy()
            candidate[start:stop + 1] = local_base + size * taper[:, None] * direction
            weight, _ = blend((candidate - front[0]) @ toward, *guide["deformation_window_mm"], .1)
            short, w = candidate[start - 1:stop + 2], weight[start - 1:stop + 2]
            length, _, radius = route_metrics(candidate)
            blocked = []
            for pose, _, scapula, humerus in seeds:
                frame = humerus if args.distal_detour else scapula
                posed = (pose_route(short, scapula, humerus, w) - frame[:, 3]) @ frame[:, :3]
                hits = overlaps(A.D.outer_envelope(posed, .40), local_targets)
                if hits:
                    blocked.append(dict(pose_deg=pose, hits_mm3=hits))
            row = dict(offset_mm=size, azimuth_deg=azimuth, neutral_arc_mm=length,
                       neutral_minimum_edited_radius_mm=float(radius[start - 1:stop + 1].min()), blocked=blocked)
            trials.append((row, candidate, weight))
    ranked = sorted(trials, key=lambda entry: (len(entry[0]["blocked"]),
                    max((sum(row["hits_mm3"].values()) for row in entry[0]["blocked"]), default=0.),
                    -entry[0]["neutral_minimum_edited_radius_mm"], entry[0]["offset_mm"]))
    best, candidate, weight = ranked[0]
    print("Best direction:", best["offset_mm"], best["azimuth_deg"], "blocked seeds:", len(best["blocked"]), flush=True)
    port = {}
    if args.smooth_port:
        support = solids["support"]
        proximal_envelopes = []
        cuts = []
        for _, _, scapula, humerus in poses:
            posed = (pose_route(candidate[:stop + 2], scapula, humerus, weight[:stop + 2])
                     - scapula[:, 3]) @ scapula[:, :3]
            envelope = A.D.outer_envelope(posed, .50)
            proximal_envelopes.append(envelope)
            overlap = envelope ^ support
            if overlap.volume() > 1e-8:
                cuts.append(np.asarray(overlap.to_mesh64().vert_properties)[:, :3])
        if cuts:
            mf = A.I._manifold()
            cutter = mf.Manifold.hull_points(np.ascontiguousarray(np.concatenate(cuts)))
            removed = cutter ^ support
            revised = (support - cutter).as_original().simplify(1e-7)
            guards = bounded({key: solids[key].minkowski_sum(mf.Manifold.cube([1.] * 3, True))
                              for key in mount["meshes"] if key.startswith("fastener_")})
            guard_hits = overlaps(removed, guards)
            components = len(revised.decompose())
            port.update(removed_mm3=float(removed.volume()), guard_hits_mm3=guard_hits, components=components)
            if components == 1 and not guard_hits:
                actual, source = export(revised, local / "FJ1467_smooth_port_support.obj", candidate[1])
                residuals = [float((envelope ^ actual).volume()) for envelope in proximal_envelopes]
                port.update(mesh=source, envelope_hit_poses=sum(value > .001 for value in residuals),
                            maximum_envelope_overlap_mm3=max(residuals))
                solids["support"] = actual
    rows = []
    for pose, clavicle, scapula, humerus in poses:
        frames = {"clavicle": clavicle, "scapula": scapula, **dict.fromkeys(MOBILE, humerus)}
        targets = bounded({key: solid.transform(frames[owners[key]]) if owners[key] in frames else solid
                           for key, solid in solids.items()})
        posed = pose_route(candidate, scapula, humerus, weight)
        length, _, radius = route_metrics(posed)
        rows.append(dict(pose_deg=pose, arc_mm=length,
                         edited_minimum_radius_mm=float(radius[start - 1:stop + 1].min()),
                         proximal_hits_mm3=overlaps(A.D.outer_envelope(posed[:27], .30), targets),
                         full_hits_mm3=overlaps(A.D.outer_envelope(posed, .30), targets)))
    name = ("FJ1467_distal_detour_route.json" if args.distal_detour else
            "FJ1467_smooth_port_exit_route.json" if args.smooth_port else "FJ1467_stepped_exit_route.json")
    save(local / name, dict(
        generator=provenance(ROOT / "tools/study_s1_posterior_exit_route.py"),
        mount=mount_ref, seed_poses=len(seeds), trials=[entry[0] for entry in trials],
        selected_for_full_replay=best, neutral_candidate_mm=candidate.tolist(), cases=rows,
        matching_port=port, proximal_exit=proximal_ref, edited_nodes=[start, stop],
        scope="Screen uses a 0.40 mm envelope on a pose subset. Smooth-port mode excludes the support "
              "from initial screening, then tests a 0.50 mm swept clearance with fastener guards. "
              "Full replay uses a 0.30 mm envelope against the rigid assembly. Offset ends retain "
              "nominal anchor points and the incoming lead. Distal mode screens humeral obstacles. "
              "TPU support and tension are not solved.",
    ))
    print(f"Replayed {len(rows)} poses, {sum(bool(row['proximal_hits_mm3']) for row in rows)} proximal hits", flush=True)


def guide_local_weight(vertices, path, fade_radius):
    """Return a smooth scapular influence and its spatial gradient."""
    start = path[:-1]
    edge = np.diff(path, axis=0)
    relative = vertices[:, None, :] - start
    fraction = np.clip(np.einsum("ijk,jk->ij", relative, edge) / np.sum(edge**2, axis=1), 0., 1.)
    offsets = relative - fraction[:, :, None] * edge
    closest = offsets[np.arange(len(vertices)), np.argmin(np.sum(offsets**2, axis=2), axis=1)]
    distance = np.linalg.norm(closest, axis=1)
    weight, slope = blend(distance, 1.25, fade_radius, 0.)
    direction = closest / np.maximum(distance[:, None], 1e-12)
    return 1. - weight, -slope[:, None] * direction


def audit_distal_guide(local, mount_ref, neighbor, solids, owners, poses, front, toward, window, *, local_support, relief):
    route_ref = provenance(local / "FJ1467_distal_detour_route.json")
    route = load_reference(route_ref)
    proximal = load_reference(route["proximal_exit"])
    if route["mount"] != mount_ref or proximal["mount"] != mount_ref or any(
            row["full_hits_mm3"] for row in route["cases"]):
        raise ValueError("Expected the rigid-checked posterior route and mount")
    solids["support"] = load_component(proximal["matching_port"]["mesh"])
    path = np.asarray(route["neutral_candidate_mm"])
    first, last = route["edited_nodes"]
    region = slice(first - 1, last + 2)
    original, _, _ = indexed_solid(neighbor["replayed_tpu"])
    mf = A.I._manifold()
    outer = A.D.outer_envelope(path[region], .95)
    lumen = A.D.outer_envelope(path[region], .50)
    wall = outer - lumen
    part = ((original + outer) - lumen).as_original().simplify(1e-7)
    bond = float((wall ^ original).volume())
    baseline_ref = None
    if local_support:
        baseline_ref = provenance(local / "FJ1467_distal_guide_TPU.json")
        baseline = load_reference(baseline_ref)
        if baseline["route"] != route_ref or baseline["parent_TPU"] != neighbor["replayed_tpu"]:
            raise ValueError("Previous distal guide does not match the selected route and lattice")
        mesh = baseline["mesh"]
    else:
        if bond <= .001 or sum(piece.volume() > 1e-5 for piece in part.decompose()) != 1:
            raise ValueError("Distal guide did not join the lattice")
        _, mesh = export_lattice(part, local / "FJ1467_distal_guide_TPU.obj", path[first])
    actual, vertices, faces = indexed_solid(mesh)
    inventory = load_reference(neighbor["neighbor_channel_reference"])
    middle = np.asarray(next(row["local_route"]["points_source_mm"] for row in inventory["routes"] if row["id"] == "FJ1513"))
    neighbors = {"FJ1468": (front, "clavicle"), "FJ1513": (middle, "scapula")}
    removed = original - actual
    preserved = {name: float((removed ^ A.D.outer_envelope(points, 1.10)).volume())
                 for name, (points, _) in neighbors.items()}

    def weights(points):
        return blend((points - front[0]) @ toward, *window, .1)[0]

    def deform(points, weight, start, end):
        beginning = points @ start[:, :3].T + start[:, 3]
        finish = points @ end[:, :3].T + end[:, 3]
        return beginning + weight[:, None] * (finish - beginning)

    relief_study = {}
    if relief:
        supported_ref = provenance(local / "FJ1467_distal_guide_supported_motion.json")
        supported = load_reference(supported_ref)
        if supported["mesh"] != mesh or supported["route"] != route_ref:
            raise ValueError("Supported replay does not match the guide")
        # Map the nearby shaft into the undeformed TPU at the neutral,
        # worst and terminal poses. Check every saved pose after the cut.
        from scipy.optimize import root
        shaft = "series:a_shaft"
        peak = max(supported["cases"], key=lambda row: row["neighborhood_rigid_hits_mm3"].get(shaft, 0.))
        pose_keys = {tuple(poses[0][0]), tuple(peak["pose_deg"]), tuple(poses[-1][0])}
        samples = []
        for pose, clavicle, scapula, humerus in poses:
            if tuple(pose) not in pose_keys:
                continue
            w = weights(vertices)
            influence, _ = guide_local_weight(vertices, path[region], 24.)
            moved = deform(vertices, w, clavicle, humerus)
            displacement = (vertices @ (scapula[:, :3] - clavicle[:, :3]).T
                            + scapula[:, 3] - clavicle[:, 3])
            moved += (1. - w[:, None]) * influence[:, None] * displacement
            posed = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
            cable = deform(path[region], weights(path[region]), scapula, humerus)
            contact = (posed ^ A.D.outer_envelope(cable, 1.20)) ^ solids[shaft].transform(humerus)
            if contact.volume() <= .001:
                continue
            for point in np.asarray(contact.to_mesh64().vert_properties)[:, :3]:
                def residual(value):
                    x = value[None, :]
                    axial = weights(x)
                    share, _ = guide_local_weight(x, path[region], 24.)
                    delta = x @ (scapula[:, :3] - clavicle[:, :3]).T + scapula[:, 3] - clavicle[:, 3]
                    return (deform(x, axial, clavicle, humerus)
                            + (1. - axial[:, None]) * share[:, None] * delta)[0] - point

                estimate = (point - humerus[:, 3]) @ humerus[:, :3]
                solution = root(residual, estimate, tol=1e-9)
                if np.linalg.norm(residual(solution.x)) > 1e-6:
                    raise ValueError("Could not map shaft overlap back to neutral")
                samples.append(solution.x)
        if not samples:
            raise ValueError("No shaft overlap samples for relief")
        cutter = mf.Manifold.hull_points(np.ascontiguousarray(samples))
        cutter = cutter.minkowski_sum(mf.Manifold.cube([.30] * 3, True))
        protected = A.D.outer_envelope(path[region], .72)
        revised = (actual - (cutter - protected)).as_original().simplify(1e-7)
        if sum(piece.volume() > 1e-5 for piece in revised.decompose()) != 1:
            raise ValueError("Relief disconnected the TPU lattice")
        _, repaired = export_lattice(revised, local / "FJ1467_distal_guide_rail_relief.obj", path[first])
        repaired_solid, vertices, faces = indexed_solid(repaired)
        lost = actual - repaired_solid
        guard_hits = {name: float((lost ^ A.D.outer_envelope(points, 1.10)).volume())
                      for name, (points, _) in neighbors.items()}
        protected_loss = max(0., float((lost ^ protected).volume()))
        if max([protected_loss, *guard_hits.values()]) > .001:
            raise ValueError("Relief removed protected channel material")
        relief_study = dict(supported_motion=supported_ref, mesh_before_relief=mesh,
                            removed_mm3=float(lost.volume()), guarded_radius_mm=.72,
                            guard_material_loss_mm3=protected_loss, neighbor_guard_loss_mm3=guard_hits,
                            cutter_pose_keys=sorted(pose_keys))
        mesh = repaired
    vertex_weights = weights(vertices)
    if local_support:
        fade_radius = 24.
        influence, spatial_gradient = guide_local_weight(vertices, path[region], fade_radius)
        _, axial_slope = blend((vertices - front[0]) @ toward, *window, .1)
    path_weights = weights(path)
    neighbor_weights = {name: weights(points) for name, (points, _) in neighbors.items()}
    probes = []
    for pose, clavicle, scapula, humerus in poses:
        frames = {"clavicle": clavicle, "scapula": scapula, **dict.fromkeys(MOBILE, humerus)}
        moved = deform(vertices, vertex_weights, clavicle, humerus)
        minimum_jacobian = None
        if local_support:
            a, b, c = (frame[:, :3] for frame in (clavicle, scapula, humerus))
            ta, tb, tc = (frame[:, 3] for frame in (clavicle, scapula, humerus))
            difference = vertices @ (b - a).T + tb - ta
            base_rotation = a + influence[:, None, None] * (b - a)
            origin_world = vertices @ a.T + ta + influence[:, None] * difference
            start_jacobian = base_rotation + difference[:, :, None] * spatial_gradient[:, None, :]
            moved += (1. - vertex_weights[:, None]) * influence[:, None] * difference
            jacobian = (start_jacobian + vertex_weights[:, None, None] * (c - start_jacobian)
                        + axial_slope[:, None, None] * (vertices @ c.T + tc - origin_world)[:, :, None] * toward)
            minimum_jacobian = float(np.linalg.det(jacobian).min())
        shaped = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
        if shaped.status() != mf.Error.NoError or (local_support and minimum_jacobian <= 0.):
            if not local_support:
                raise ValueError(f"Invalid interpolated TPU at {pose}")
            probes.append(dict(pose_deg=pose, minimum_vertex_jacobian=minimum_jacobian,
                               invalid_mesh_reason=str(shaped.status())))
            continue
        posterior = deform(path[region], path_weights[region], scapula, humerus)
        envelope = A.D.outer_envelope(posterior, .30)
        zone = shaped ^ A.D.outer_envelope(posterior, 1.10)
        targets = bounded({key: solid.transform(frames[owners[key]]) if owners[key] in frames else solid
                           for key, solid in solids.items()})
        adjacent = {}
        for name, (points, owner) in neighbors.items():
            cable = deform(points, neighbor_weights[name], frames[owner], humerus)
            adjacent[name] = float((A.D.outer_envelope(cable, .30) ^ zone).volume())
        probes.append(dict(pose_deg=pose, posterior_wire_TPU_mm3=float((envelope ^ shaped).volume()),
                           neighborhood_rigid_hits_mm3=overlaps(zone, targets), neighbor_wire_zone_mm3=adjacent,
                           minimum_vertex_jacobian=minimum_jacobian))
    filename = ("FJ1467_distal_guide_rail_relief.json" if relief else
                "FJ1467_distal_guide_supported_motion.json" if local_support else "FJ1467_distal_guide_TPU.json")
    save(local / filename, dict(
        generator=provenance(ROOT / "tools/study_s1_posterior_exit_route.py"),
        route=route_ref, parent_TPU=neighbor["replayed_tpu"], mesh=mesh, prior_guide_replay=baseline_ref,
        bore_radius_mm=.50, outer_radius_mm=.95, nominal_wall_mm=.45,
        neutral_wall_bond_volume_mm3=bond, neighboring_guard_removed_mm3=preserved, cases=probes,
        relief=relief_study,
        scope=("Shared TPU follows scapular support within 1.25 mm of the neutral distal route, "
               f"fading into the clavicular support field by {fade_radius:g} mm. " if local_support else
               "Candidate liner with nominal blend of shared TPU from clavicle to humerus. ") + "Posterior and "
              "middle steel follow scapula-to-humerus interpolation. Check only the edited posterior "
              "segment and 1.10 mm neighborhood. Neighbor intersections are local, while neutral "
              "stock guards follow complete neighbor routes. Contact and load capacity are unverified.",
    ))
    print(f"Checked distal guide in {len(probes)} poses", flush=True)


if __name__ == "__main__":
    main()
