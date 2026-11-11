"""Quantify imposed boundary and lattice geometry demand for the original elbow."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, mesh_geom, read, save, sha


def screen_deltoid_relief(folder, deltoid, extensor_poses, obstacles,
                         filename='deltoid_anatomical_elbow_motion_relief.obj',
                         capture_threshold_mm3=.001, near_capture_mm=0.):
    mf = A.I._manifold()
    contact_patches = []
    directions = [np.zeros(3)]
    if near_capture_mm:
        directions.extend([*np.eye(3), *(-np.eye(3))])
    for body in extensor_poses:
        if near_capture_mm and A.I.gap(deltoid, body, near_capture_mm) >= near_capture_mm:
            continue
        for direction in directions:
            shifted = body.translate(near_capture_mm*direction)
            if A.I.intersection_volume(deltoid, shifted) > capture_threshold_mm3:
                contact_patches.append(deltoid ^ shifted)
    if not contact_patches:
        return dict(trials=[], scope='No overlap exceeded the screening threshold.')
    hull = mf.Manifold.batch_hull(contact_patches)
    origin = np.asarray(hull.bounding_box()).reshape(2, 3).mean(axis=0)
    route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    heads = {'FJ1468', 'FJ1513', 'FJ1467'}
    routes = {row['head']: row for row in read(route_file)['routes'] if row['head'] in heads}
    assert routes.keys() == heads
    guards = {head: A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), 2.)
              for head, row in routes.items()}
    guards.update({key: part for key, part in obstacles.items()
                   if key.startswith('cap:') and key.split('/')[0][4:] in heads})
    trials = []
    for margin in (.20, .10, .05):
        cutter = hull.minkowski_sum(mf.Manifold.cube([2*margin]*3, True))
        candidate = (deltoid-cutter).as_original().simplify(1e-6)
        removed = deltoid-candidate
        losses = {name: A.I.intersection_volume(removed, guard) for name, guard in guards.items()}
        row = dict(padding_mm=margin, removed_mm3=float(removed.volume()), protected_stock_loss_mm3=losses)
        trials.append(row)
        if max(losses.values()) > .001:
            continue
        try:
            actual, source = export(candidate, folder/filename, origin)
        except (AssertionError, ValueError) as error:
            row['export_error'] = str(error) or type(error).__name__
            continue
        row['peak_readback_overlap_mm3'] = max(A.I.intersection_volume(actual, body)
                                               for body in extensor_poses)
        row['min_readback_clearance_mm'] = min(A.I.gap(actual, body, 1.) for body in extensor_poses)
        if row['peak_readback_overlap_mm3'] <= .001 and row['min_readback_clearance_mm'] > 0:
            row['source'] = source
            break
    return dict(route_sha256=sha(route_file), tested_pose_count=len(extensor_poses), trials=trials,
                cutter_capture_threshold_mm3=capture_threshold_mm3,
                near_capture_axis_translation_mm=near_capture_mm,
                scope='Shoulder-side subtractive pocket against prescribed elbow TPU interpolations. '
                      'No coupled shoulder motion, contact response or load verification.')


def map_humerus_relief(obstacles, extensor_poses, motion):
    inventory_file = A.DESIGN.parent/'anatomical_concept/inventory.json'
    reference = next(row for row in read(inventory_file)['models'] if row['id'] == 'FJ3368')
    source_path = inventory_file.parents[2]/'输入'/reference['source']
    anatomy, _, _ = A.load_checked(source_path, expected=reference['sha256'])
    humerus = obstacles['bone:humerus']
    mf = A.I._manifold()
    contacts = []
    patches = []
    for angle, body in extensor_poses.items():
        if A.I.intersection_volume(humerus, body) <= .001:
            continue
        patch = humerus ^ body
        native_volume = A.I.intersection_volume(patch, anatomy)
        contacts.append(dict(elbow_deg=angle, overlap_mm3=float(patch.volume()),
                             inside_source_anatomy_mm3=native_volume,
                             outside_source_anatomy_mm3=max(0., float(patch.volume())-native_volume)))
        patches.append(patch)
    if not patches:
        return dict(native_humerus_sha256=reference['sha256'], contact_poses=contacts)
    envelope = mf.Manifold.batch_hull(patches)
    padding = .2
    stock = humerus ^ envelope.minkowski_sum(mf.Manifold.cube([2*padding]*3, True))
    native_stock = A.I.intersection_volume(stock, anatomy)
    joint_keys = ('elbow:elbow_shoulder_pin', 'elbow:elbow_spacer_left', 'elbow:elbow_spacer_right',
                  'elbow:elbow_bearing_bush', 'bone:forearm')
    gaps = {}
    for key in joint_keys:
        part = obstacles[key]
        if key in motion['moving_keys']:
            parts = (part.transform(A.I.rotation_about(motion['centre_source_mm'],
                                                       motion['axis_source'], angle))
                     for angle in range(121))
        else:
            parts = (part,)
        gaps[key] = min(A.I.gap(stock, posed, 5.) for posed in parts)
    return dict(native_humerus_sha256=reference['sha256'], contact_poses=contacts,
                pocket_padding_mm=padding, pocket_stock_mm3=float(stock.volume()),
                pocket_inside_source_anatomy_mm3=native_stock, pocket_to_joint_gap_mm_capped5=gaps,
                scope='Stock mapping only. Humerus retained because this pocket includes original anatomy. '
                      'TPU contact and strength remain unverified.')


def screen_extensor_bypass(vertices, faces, stations, obstacles, before):
    """Find a low-displacement geometric bypass with both pad slabs fixed."""
    mf = A.I._manifold()
    humerus = obstacles['bone:humerus']
    threshold = .001
    initial = before['bone:humerus']
    best = None
    checked = 0
    # Search around the distal shaft on the posterior side of the arm.
    for amplitude in (1., 2., 3., 4., 5., 6., 8., 10., 12.):
        for centre, width in ((.85, .15), (.90, .10), (.80, .20), (.75, .25)):
            phase = np.clip((stations-centre)/width, -1., 1.)
            bump = .5*(1.+np.cos(np.pi*phase))
            for azimuth in (90., 60., 120., 30., 150.):
                azimuth_rad = np.deg2rad(azimuth)
                direction = np.array([np.cos(azimuth_rad), np.sin(azimuth_rad), 0.])
                posed = vertices+amplitude*bump[:, None]*direction
                body = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), faces))
                if body.status() != mf.Error.NoError or body.is_empty():
                    continue
                checked += 1
                residual = A.I.intersection_volume(body, humerus)
                if best is not None and residual >= best['humerus_overlap_mm3'] and residual > threshold:
                    continue
                hits = {key: volume for key, part in obstacles.items() if key != 'bone:humerus'
                        and (volume := A.I.intersection_volume(body, part)) > threshold}
                increased = {key: volume-before.get(key, 0.) for key, volume in hits.items()
                             if volume > before.get(key, 0.)+threshold}
                if increased:
                    continue
                gap = A.I.gap(body, humerus, 1.) if residual <= threshold else 0.
                row = dict(amplitude_mm=amplitude, centre_station=centre, half_width_station=width,
                           direction_source=direction.tolist(), humerus_overlap_mm3=residual,
                           humerus_clearance_mm_capped1=gap, other_overlaps_mm3=hits)
                if best is None or residual < best['humerus_overlap_mm3']:
                    best = row
                if residual <= threshold and gap > 0.:
                    return dict(initial_humerus_overlap_mm3=initial, tested_candidates=checked,
                                selected=row, positive_clearance=True)
    return dict(initial_humerus_overlap_mm3=initial, tested_candidates=checked,
                selected=best, positive_clearance=False)


def bypass_displacement(angle, stations, bypass_poses):
    knots = [row['elbow_deg'] for row in bypass_poses]
    selected = [row['selected'] for row in bypass_poses]
    assert all(row['positive_clearance'] for row in bypass_poses)
    assert len(knots) == 2 and all(np.allclose(row['direction_source'], [0, 1, 0]) for row in selected)
    start = knots[0]-10
    knot_angles = [start, *knots]
    amplitude = np.interp(angle, knot_angles, [0., *[row['amplitude_mm'] for row in selected]])
    centre = np.interp(angle, knots, [row['centre_station'] for row in selected])
    half_width = np.interp(angle, knots, [row['half_width_station'] for row in selected])
    phase = np.clip((stations-centre)/half_width, -1., 1.)
    return amplitude*.5*(1.+np.cos(np.pi*phase)), float(amplitude)


def verify_combined_bypass(extensor_meshes, bypass_poses, deltoid, humerus):
    mf = A.I._manifold()
    rows = []
    for angle, (vertices, faces, stations) in extensor_meshes.items():
        offset, amplitude = bypass_displacement(angle, stations, bypass_poses)
        posed = vertices.copy()
        posed[:, 1] += offset
        body = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), faces))
        assert body.status() == mf.Error.NoError and not body.is_empty()
        rows.append(dict(elbow_deg=angle, imposed_bypass_mm=float(amplitude),
                         humerus_overlap_mm3=A.I.intersection_volume(body, humerus),
                         deltoid_overlap_mm3=A.I.intersection_volume(body, deltoid),
                         humerus_clearance_mm_capped1=A.I.gap(body, humerus, 1.),
                         deltoid_clearance_mm_capped1=A.I.gap(body, deltoid, 1.)))
    return rows


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    assembly_file = folder/'assembly_anatomical_elbow_angled_candidate.json'
    motion_file = folder/'anatomical_elbow_rigid_motion.json'
    wire_file = folder/'anatomical_elbow_constrained_wire_motion.json'
    baseline_file = folder/'anatomical_elbow_baseline.json'
    motion, wire = read(motion_file), read(wire_file)
    assert motion['assembly_sha256'] == sha(assembly_file)
    assert wire['rigid_motion_sha256'] == sha(motion_file)
    assert wire['anatomical_baseline_sha256'] == sha(baseline_file)
    sources = read(assembly_file)['sources']
    identity = np.eye(3)
    results = {}
    mf = A.I._manifold()
    posed_muscles = {}
    extensor_meshes = {}
    swivel_transforms = {}
    # The original generator makes each TPU pad from a 4 mm thick disk.
    pad_half_thickness = 2.
    for label, muscle in read(baseline_file)['muscles'].items():
        source = sources[f"elbow:{muscle['anatomical_id']}_cap_rotor"]
        screen_file = (ROOT/source['file']).parent/'screen.json'
        assert sha(screen_file) == wire['swivel_screen_sha256'][label]
        screen = read(screen_file)
        pivot = np.asarray(screen['pivot_zero_display_mm'])+np.asarray(source['translation_source_mm'])
        hinge = screen['trunnion_axis_zero_display']
        definition = ROOT/muscle['definition']
        assert sha(definition) == muscle['definition_sha256']
        axis = np.asarray(read(definition)['basis_columns_display'])[:, 2]
        proximal, distal = np.asarray(muscle['pad_centres_source_mm'])
        neutral_span = np.linalg.norm(distal-proximal)
        free_span_origin = proximal+axis*pad_half_thickness
        station_grad = axis/((distal-proximal)@axis-2*pad_half_thickness)
        proxy = load_component(muscle['collision_proxy'])
        vertices, faces = mesh_geom.manifold_mesh_to_arrays(proxy.to_mesh64())
        edges = np.unique(np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]],
                                                   faces[:, [2, 0]])), axis=1), axis=0)
        lengths = np.linalg.norm(vertices[edges[:, 1]]-vertices[edges[:, 0]], axis=1)
        mask = lengths >= .2
        edges, lengths = edges[mask], lengths[mask]
        stations = (vertices-free_span_origin)@station_grad
        weights = np.clip(stations, 0., 1.)
        sample_vertices = vertices[(stations > 0.) & (stations < 1.)]
        sample_weights = weights[(stations > 0.) & (stations < 1.)]
        poses = []
        for row in wire['poses'][::5] if label != 'extensor' else wire['poses']:
            angle = row['elbow_deg']
            wire_pose = row['muscles'][label]
            swivel = A.I.rotation_about(pivot, hinge, wire_pose['swivel_deg'])
            bend = A.I.rotation_about(motion['centre_source_mm'], motion['axis_source'], angle)
            rotation = bend[:, :3]@swivel[:, :3]
            translation = bend[:, :3]@swivel[:, 3]+bend[:, 3]
            distal_posed = rotation@distal+translation
            span = np.linalg.norm(distal_posed-proximal)
            delta_rotation = rotation-identity
            displacement = vertices@delta_rotation.T+translation
            posed = vertices+weights[:, None]*displacement
            if label == 'extensor':
                extensor_meshes[angle] = (posed, faces.astype(np.uint32), stations)
            if angle % 5:
                continue
            posed_muscles[angle, label] = mf.Manifold(
                mf.Mesh64(np.ascontiguousarray(posed), faces.astype(np.uint32)))
            swivel_transforms[angle, muscle['anatomical_id']] = swivel
            ratios = np.linalg.norm(posed[edges[:, 1]]-posed[edges[:, 0]], axis=1)/lengths
            # Gradient of x + w(x) * ((R-I)x + t), on the interpolated span.
            gradient = (identity+sample_weights[:, None, None]*delta_rotation
                        +(sample_vertices@delta_rotation.T+translation)[:, :, None]*station_grad)
            det = np.linalg.det(gradient)
            poses.append(dict(elbow_deg=angle, distal_pad_displacement_mm=(distal_posed-distal).tolist(),
                              pad_distance_mm=float(span), pad_distance_change_mm=float(span-neutral_span),
                              pad_distance_change_percent=float(100.*(span/neutral_span-1.)),
                              wire_length_change_mm=wire_pose['length_change_mm'],
                              surface_edge_length_ratios=np.quantile(ratios, [0., .05, .5, .95, 1.]).tolist(),
                              min_interpolation_jacobian=float(det.min()),
                              nonpositive_jacobian_vertex_count=int(np.count_nonzero(det <= 0))))
        results[label] = dict(neutral_pad_distance_mm=float(neutral_span), sampled_surface_edges=len(edges),
                              sampled_interior_span_vertices=len(sample_vertices), poses=poses)
        print(label, poses[-1])
    demand_file = folder/'anatomical_elbow_tpu_clamped_boundary_demand.json'
    save(demand_file, dict(
        assembly_sha256=sha(assembly_file), wire_motion_sha256=sha(wire_file),
        baseline_sha256=sha(baseline_file), pose_step_deg=5, min_surface_edge_length_mm=.2,
        surface_edge_ratio_quantiles=[0., .05, .5, .95, 1.], muscles=results,
        rigid_pad_half_thickness_mm=pad_half_thickness,
        scope='Imposed end-cap boundary motion at neutral shoulder. Surface edge ratios and local map '
              'Jacobian characterize the interpolation, not strain, allowable stroke or contact of real TPU.'))
    moving = set(motion['moving_keys'])
    obstacles = {key: load_component(source) for key, source in sources.items()
                 if key not in motion['excluded_unposed_elbow_tpu']}
    rotating = {}
    for head in {head for _, head in swivel_transforms}:
        rotating[head] = {key for key in moving if key.startswith(f'elbow:{head}_distal_')
                          or key in {f'elbow:{head}_cap_rotor', f'elbow:{head}_active_wire_crimp_barrel'}}
    contact_poses = []
    bypass_poses = []
    threshold = .001
    for angle in range(0, 121, 5):
        bend = A.I.rotation_about(motion['centre_source_mm'], motion['axis_source'], angle)
        posed_obstacles = {**obstacles, **{key: obstacles[key].transform(bend) for key in moving}}
        for head, pack in rotating.items():
            swivel = swivel_transforms[angle, head]
            posed_obstacles.update({key: obstacles[key].transform(swivel).transform(bend) for key in pack})
        bodies = {label: posed_muscles[angle, label] for label in results}
        status = {label: str(body.status()) for label, body in bodies.items()}
        valid = {label: body for label, body in bodies.items() if body.status() == mf.Error.NoError
                 and not body.is_empty()}
        contacts = {label: {key: volume for key, part in posed_obstacles.items()
                            if (volume := A.I.intersection_volume(body, part)) > threshold}
                    for label, body in valid.items()}
        mutual = A.I.intersection_volume(*valid.values()) if len(valid) == 2 else None
        extensor_regions = {
            key: np.asarray((valid['extensor'] ^ posed_obstacles[key]).bounding_box()).reshape(2, 3).tolist()
            for key in ('bone:humerus', 'tpu:deltoid') if key in contacts.get('extensor', {})}
        contact_poses.append(dict(elbow_deg=angle, mesh_status=status, muscle_to_parts_mm3=contacts,
                                  between_elbow_muscles_mm3=mutual,
                                  extensor_overlap_bounds_source_mm=extensor_regions))
        if 'bone:humerus' in extensor_regions:
            bypass_poses.append(dict(elbow_deg=angle, **screen_extensor_bypass(
                *extensor_meshes[angle], {**posed_obstacles, 'tpu:elbow_flexor': valid['flexor']},
                contacts['extensor'])))
    contact_file = folder/'anatomical_elbow_tpu_clamped_contact_screen.json'
    save(contact_file, dict(
        boundary_demand_sha256=sha(demand_file), rigid_pad_half_thickness_mm=pad_half_thickness,
        fixed_and_moving_obstacles=len(obstacles), pose_step_deg=5, threshold_mm3=threshold,
        poses=contact_poses, scope='Boolean overlap of imposed TPU interpolations versus posed assembly. '
                                  'Mating interfaces included. No contact response, self-contact, or load model.'))
    print('contact poses:', sum(any(row['muscle_to_parts_mm3'].values()) for row in contact_poses),
          'invalid mesh poses:', sum(len(row['muscle_to_parts_mm3']) < 2 for row in contact_poses))
    assert all(len(row['muscle_to_parts_mm3']) == 2 for row in contact_poses)
    extensor_poses = {angle: body for (angle, label), body in posed_muscles.items() if label == 'extensor'}
    relief = screen_deltoid_relief(folder, obstacles['tpu:deltoid'], extensor_poses.values(), obstacles)
    relief['contact_screen_sha256'] = sha(contact_file)
    save(folder/'anatomical_elbow_deltoid_motion_relief.json', relief)
    print(relief['trials'])
    bone_map = map_humerus_relief(obstacles, extensor_poses, motion)
    bone_map['contact_screen_sha256'] = sha(contact_file)
    save(folder/'anatomical_elbow_humerus_contact_stock.json', bone_map)
    print(bone_map)
    bypass_file = folder/'anatomical_elbow_extensor_bypass_screen.json'
    save(bypass_file, dict(
        contact_screen_sha256=sha(contact_file), poses=bypass_poses,
        scope='Cosine-shaped lateral displacement at flagged poses, with original rigid pad motion. '
              'No force equilibrium, self-contact, cable fit, or continuous bypass trajectory.'))
    print(bypass_poses)
    deltoid_source = next(row['source'] for row in relief['trials'] if 'source' in row)
    combined = verify_combined_bypass(extensor_meshes, bypass_poses, load_component(deltoid_source),
                                     obstacles['bone:humerus'])
    save(folder/'anatomical_elbow_combined_bypass.json', dict(
        bypass_screen_sha256=sha(bypass_file), deltoid_source=deltoid_source, sample_step_deg=1,
        path_deg=[0, 120, 0], poses=combined,
        scope='Continuous prescribed bypass law, sampled at neutral shoulder against humerus and '
              'relieved deltoid only. Return retraces the same geometry. No cable or equilibrium check.'))
    for name in ('humerus', 'deltoid'):
        print(name, 'peak', max(row[f'{name}_overlap_mm3'] for row in combined),
              'min gap', min(row[f'{name}_clearance_mm_capped1'] for row in combined))


if __name__ == '__main__':
    main()
