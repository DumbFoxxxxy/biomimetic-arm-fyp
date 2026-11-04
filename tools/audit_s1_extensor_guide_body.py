"""Screen a finite guide body at the selected fixed extensor waypoint."""
import numpy as np
from itertools import combinations

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha
from study_s1_relocated_elbow_bridges import check_fasteners
from study_s1_elbow_proximal_stations import proximal_groups


def screen_mount(bridge, pair, guide, incoming, relief, folder):
    reaction = next(row for row in bridge['proximal_reaction_seats'] if row['muscle'] == 'extensor')
    cap = load_component(next(flat['source'] for trial in reaction['trials']
                              for flat in trial.get('flat_reaction_cap', []) if 'source' in flat))
    insert = load_component(pair['insert'])
    tpu = load_component(pair['tpu'])
    adapters = {key: load_component(row['source']) for key, row in bridge['selected_adapters'].items()}
    u = A.D.unit(np.cross(incoming, np.eye(3)[np.argmin(abs(incoming))]))
    v = np.cross(incoming, u)
    rows = []
    for angle in range(0, 360, 15):
        theta = np.radians(angle)
        radial = u*np.cos(theta)+v*np.sin(theta)
        for reach in (2., 3., 4., 5.):
            outer = guide+radial*reach
            ray = cap.ray_cast(outer, outer-incoming*14.)
            if not ray:
                continue
            root = np.array(ray[0].position)-incoming*.8
            bridge_stock = rod(guide+radial*1.15, outer, .65)+rod(outer+incoming*.35, root, .65)
            mount = (cap+insert+(bridge_stock-relief)).as_original().simplify(1e-6)
            added = mount-cap
            row = dict(azimuth_deg=angle, radial_reach_mm=reach, components=len(mount.decompose()),
                       wire_overlap_mm3=float(A.I.intersection_volume(added, relief)),
                       tpu_overlap_mm3=float(A.I.intersection_volume(added, tpu)),
                       adapter_hits={key: vol for key, part in adapters.items()
                                     if (vol := float(A.I.intersection_volume(added, part))) > .001})
            rows.append(row)
            if row['components'] != 1 or row['wire_overlap_mm3'] > .001 or row['tpu_overlap_mm3'] > .001 or row['adapter_hits']:
                continue
            try:
                _, row['source'] = export(mount, folder/'extensor_cap_with_fixed_guide.obj', guide)
            except (AssertionError, ValueError) as error:
                row['export_error'] = str(error)
                continue
            return rows
    return rows


def screen_cap_after_tpu(cap, tpu, bridge, incoming, extra_obstacles=None):
    """Test straight entry of a rigid component into the local assembly."""
    obstacles = {'extensor_tpu': tpu}
    obstacles.update({key: load_component(row['source'])
                      for key, row in bridge['selected_adapters'].items()})
    obstacles.update(extra_obstacles or {})
    u = A.D.unit(np.cross(incoming, np.eye(3)[np.argmin(abs(incoming))]))
    v = np.cross(incoming, u)
    mf = A.I._manifold()
    results = []
    for pitch in (-90, -45, 0, 45, 90):
        for yaw in ([0] if abs(pitch) == 90 else range(0, 360, 30)):
            theta, phi = np.radians([pitch, yaw])
            direction = incoming*np.sin(theta)+np.cos(theta)*(u*np.cos(phi)+v*np.sin(phi))
            delta = direction*12.
            hull = mf.Manifold.batch_hull([cap, cap.translate(delta)])
            nearby = {key: part for key, part in obstacles.items()
                      if A.I.intersection_volume(hull, part) > .001}
            first_hit = None
            for distance in np.linspace(0., 12., 49):
                pose = cap.translate(direction*distance)
                hits = {key: volume for key, part in nearby.items()
                        if (volume := float(A.I.intersection_volume(pose, part))) > .001}
                if hits:
                    first_hit = dict(distance_from_seat_mm=float(distance), hits=hits)
                    break
            results.append(dict(offset_mm=delta.tolist(), first_hit=first_hit))
    return dict(travel_mm=12., step_mm=.25, checked_obstacles=list(obstacles),
                passing=sum(row['first_hit'] is None for row in results), trials=results)


def screen_side_opening(removable, tpu, cap, bridge, paths, folder, guide):
    mf = A.I._manifold()
    rigid = {'cap': cap, **{key: load_component(row['source'])
                          for key, row in bridge['selected_adapters'].items()}}
    trials = []
    for path in paths:
        delta = np.array(path['offset_mm'])
        envelope = mf.Manifold.batch_hull([removable, removable.translate(delta)])
        hits = {key: volume for key, part in rigid.items()
                if (volume := float(A.I.intersection_volume(envelope, part))) > .001}
        row = dict(offset_mm=delta.tolist(), rigid_envelope_hits=hits)
        trials.append(row)
        if hits:
            blocked = next((dict(fraction=float(t), part=key, overlap_mm3=volume)
                            for t in np.linspace(0., 1., 49) for key in hits
                            if (volume := float(A.I.intersection_volume(
                                removable.translate(delta*t), rigid[key]))) > .001), None)
            row['rigid_sampled_first_hit'] = blocked
            if blocked is not None:
                continue
        opening = envelope.minkowski_sum(mf.Manifold.cube([.12]*3, True))
        revised = (tpu-opening).as_original().simplify(1e-6)
        row.update(components=len(revised.decompose()),
                   removed_TPU_mm3=float((tpu-revised).volume()),
                   residual_sweep_overlap_mm3=float(A.I.intersection_volume(revised, envelope)))
        if row['components'] != 1 or row['residual_sweep_overlap_mm3'] > .001:
            continue
        try:
            _, row['source'] = export(revised, folder/'extensor_lattice_guide_side_opening.obj', guide)
        except (AssertionError, ValueError) as error:
            row['export_error'] = str(error)
            continue
        break
    return trials


def load_relocated_assembly(base, bridge, channels, tpu_source, extensor_cap_source=None):
    manifest = base/'assembly_tpu_readback_candidate.json'
    assert sha(manifest) == bridge['assembly_sha256']
    groups = proximal_groups(read(base/'elbow_proximal_rotated_station_screen.json')['moving_groups'])
    moving = set(sum(groups.values(), []))
    replaced = moving | {'bone:humerus', 'tpu:elbow_flexor', 'tpu:elbow_extensor'}
    sources = read(manifest)['sources']
    solids = {key: load_component(row) for key, row in sources.items() if key not in replaced}
    frame = read(A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json')['frame_source']
    axial = np.array(frame['basis_columns'])[:, 2]
    pair = channels['choice']['pair']
    for muscle, head in [('flexor', 'FJ1486'), ('extensor', 'FJ1477')]:
        pose = A.I.rotation_about(frame['origin_mm'], axial, pair[f'{muscle}_angle_deg'])
        pose[:, 3] += axial*pair[f'{muscle}_shift_mm']
        for key in groups[head]:
            if '_bone_' not in key and not key.endswith('_proximal_rigid_cap'):
                solids[key] = load_component(sources[key]).transform(pose)
    replacements = {'bone:humerus': read(base/'elbow_relocated_bone_seats/bone_seat_check.json')['source'],
                    'tpu:extensor_relocated': tpu_source}
    replacements.update({key: row['source'] for key, row in bridge['selected_adapters'].items()})
    replacements['tpu:flexor_relocated'] = next(
        terminal['source'] for trial in bridge['flexor_lattice_midspan']
        for terminal in trial.get('terminal_trials', []) if 'source' in terminal)
    for reaction, head in zip(bridge['proximal_reaction_seats'], ('FJ1486', 'FJ1477'), strict=True):
        replacements[f'elbow:{head}_proximal_rigid_cap'] = (
            extensor_cap_source if head == 'FJ1477' and extensor_cap_source else next(
                flat['source'] for trial in reaction['trials']
                for flat in trial.get('flat_reaction_cap', []) if 'source' in flat))
    replacements.update({f'extensor_plate:{key}': row
                         for key, row in bridge['extensor_thread_plate']['selected']['sources'].items()})
    solids.update({key: load_component(row) for key, row in replacements.items()})
    seats = read(base/'elbow_relocated_bone_seats/bone_seat_check.json')['seats']
    axes = read(base/'elbow_relocated_seat_audit.json')['anchors']
    flexor_hardware = {}
    check_fasteners(seats, axes, {}, bridge['selected_adapters'], installed_flexor=flexor_hardware)
    solids.update(flexor_hardware)
    extensor_bolts = bridge['extensor_thread_plate']['selected']['bolts']
    for seat, axis_row in zip(seats, axes, strict=True):
        if seat['muscle'] != 'extensor':
            continue
        index = seat['anchor']
        n = np.array(axis_row['axis'])
        front = np.array(seat['front_face_mm'])
        head = front+n*3.5
        length = extensor_bolts[index]['nominal_length_mm']
        solids[f'extensor_plate:bolt_{index}'] = rod(head, head-n*length, 1.5)+rod(head, head+n*3, 2.75)
        solids[f'extensor_plate:washer_{index}'] = rod(front+n*3, head, 2.75)-rod(front, head+n, 1.7)
    return solids


def verify_side_opening_assembly(base, bridge, channels, removable, side, wire_paths):
    solids = load_relocated_assembly(base, bridge, channels, side['source'])
    delta = np.array(side['offset_mm'])
    hull = A.I._manifold().Manifold.batch_hull([removable, removable.translate(delta)])
    nearby = {key: part for key, part in solids.items() if A.I.intersection_volume(hull, part) > .001}
    failures = []
    for fraction in np.linspace(1., 0., 49):
        pose = removable.translate(delta*fraction)
        hits = {key: volume for key, part in nearby.items()
                if (volume := float(A.I.intersection_volume(pose, part))) > .001}
        if hits:
            failures.append(dict(distance_mm=float(np.linalg.norm(delta)*fraction), hits=hits))
    prefix = 'elbow:FJ1477_proximal_'
    washer_key, nut_key = (prefix+suffix for suffix in ('nut_washer_right', 'nut_right'))
    delayed = {washer_key, nut_key}
    remaining = [row for row in failures if set(row['hits'])-delayed]
    centre = lambda part: np.array(part.bounding_box()).reshape(2, 3).mean(axis=0)
    washer, nut = solids[washer_key], solids[nut_key]
    normal = A.D.unit(centre(nut)-centre(washer))
    installed = {key: part for key, part in solids.items() if key not in delayed}
    installed['detachable_guide'] = removable
    sequence = []
    for key, part in ((washer_key, washer), (nut_key, nut)):
        sweep = A.I._manifold().Manifold.batch_hull([part, part.translate(normal*12.)])
        nearby = {name: obstacle for name, obstacle in installed.items()
                  if name != prefix+'bolt_right' and A.I.intersection_volume(sweep, obstacle) > .001}
        hits = []
        for distance in np.linspace(12., 0., 49):
            pose = part.translate(normal*distance)
            collisions = {name: vol for name, obstacle in nearby.items()
                          if (vol := float(A.I.intersection_volume(pose, obstacle))) > .001}
            if collisions:
                hits.append(dict(distance_mm=float(distance), hits=collisions))
        sequence.append(dict(key=key, failures=hits))
        installed[key] = part
    nut_centre = centre(nut)
    socket = rod(nut_centre+normal*1.5, nut_centre+normal*22., 4.2)
    socket_hits = {name: volume for name, part in installed.items() if name not in delayed
                   and name != prefix+'bolt_right'
                   and (volume := float(A.I.intersection_volume(socket, part))) > .001}
    # A keyed washer can carry a stop tongue while the original screw
    # continues to clamp the cap. Screen the tongue and washer as one part.
    direction = A.D.unit(delta)
    washer_centre = centre(washer)
    vertices = np.asarray(removable.to_mesh64().vert_properties)[:, :3]
    surface = vertices[np.argmax(vertices@direction)]
    tangent = A.D.unit(np.cross(normal, direction))
    retention = []
    nut_clearance = nut.minkowski_sum(A.I._manifold().Manifold.cube([.16]*3, True))
    for gap in (.15, .3, .5):
        for lateral in (-5., -4., 4., 5.):
            radius = .65
            tip = surface+direction*(gap+radius)
            root = washer_centre+tangent*np.sign(lateral)*2.7
            bend = washer_centre+tangent*lateral
            bridge_stock = (rod(root, bend, radius)+rod(bend, tip, radius)
                            +A.I._manifold().Manifold.sphere(radius, 24).translate(bend))
            bolt_keepout = rod(washer_centre-normal*3., washer_centre+normal*3., 1.7)
            keeper = (washer+(bridge_stock-nut_clearance)-bolt_keepout).as_original().simplify(1e-6)
            added = keeper-washer
            collisions = {name: vol for name, part in installed.items() if name != washer_key
                          and (vol := float(A.I.intersection_volume(added, part))) > .001}
            row = dict(gap_parameter_mm=gap, tongue_radius_mm=radius, lateral_offset_mm=lateral,
                       components=len(keeper.decompose()), static_hits=collisions)
            retention.append(row)
            if collisions or row['components'] != 1:
                continue
            first = next((float(distance) for distance in np.arange(.05, 1.01, .05)
                          if A.I.intersection_volume(removable.translate(direction*distance), added) > .001), None)
            row['first_reverse_stop_mm'] = first
            if first is None:
                continue
            row['axial_insertion_hits'] = {
                name: vol for name, part in installed.items() if name not in (washer_key, nut_key)
                and (vol := max(float(A.I.intersection_volume(added.translate(normal*distance), part))
                                for distance in np.linspace(0., 12., 49))) > .001}
            row['nut_insertion_overlap_mm3'] = max(
                float(A.I.intersection_volume(added, nut.translate(normal*distance)))
                for distance in np.linspace(0., 12., 49))
            row['socket_overlap_mm3'] = float(A.I.intersection_volume(added, socket))
            if (not row['axial_insertion_hits'] and row['nut_insertion_overlap_mm3'] <= .001
                    and row['socket_overlap_mm3'] <= .001 and not any('source' in r for r in retention)):
                try:
                    _, row['source'] = export(
                        keeper, base/'elbow_relocated_bone_seats/extensor_guide_stop_washer.obj', washer_centre)
                except (AssertionError, ValueError) as error:
                    row['export_error'] = str(error)
    localization = None
    selected = next((row['source'] for row in retention if 'source' in row), None)
    if selected is not None:
        fixed = {**solids, washer_key: load_component(selected)}
        origin = centre(removable)
        axes = (direction, tangent, A.D.unit(np.cross(direction, tangent)))
        limits = []
        for index, axis in enumerate(axes):
            for sign in (-1, 1):
                for motion, step, count in (('translation_mm', .05, 20), ('rotation_deg', .25, 20)):
                    first = None
                    for number in range(1, count+1):
                        value = sign*step*number
                        posed = (removable.translate(axis*value) if motion == 'translation_mm'
                                 else removable.transform(A.I.rotation_about(origin, axis, value)))
                        hits = {name: vol for name, part in fixed.items()
                                if (vol := float(A.I.intersection_volume(posed, part))) > .001}
                        if hits:
                            first = dict(value=value, hits=hits)
                            break
                    limits.append(dict(axis=index, sign=sign, motion=motion, step=step,
                                       tested_limit=step*count, first_contact=first))
        localization = dict(pivot_source_mm=origin.tolist(), axes=[axis.tolist() for axis in axes],
                            independent_axis_probes=limits,
                            scope='Independent perturbations only; compliant TPU contacts are recorded separately by part name.')
        cap_key = prefix+'rigid_cap'
        cap = fixed[cap_key]
        mf = A.I._manifold()
        sweep = mf.Manifold.batch_hull([removable, removable.translate(delta)])
        keepout = sweep.minkowski_sum(mf.Manifold.cube([.08]*3, True))
        nearby_stock = mf.Manifold.cube([12.]*3, True).translate(origin)
        trials = []
        for reach in (.5, 1., 1.5, 2.):
            expanded = (cap^nearby_stock).minkowski_sum(mf.Manifold.cube([2*reach]*3, True))
            added = (expanded-keepout)-cap
            candidate = (cap+added).as_original().simplify(1e-6)
            hits = {key: vol for key, part in fixed.items() if key != cap_key
                    and (vol := float(A.I.intersection_volume(added, part))) > .001}
            stops = []
            for axis_id, motion, sign in ((1, 'translation_mm', 1), (1, 'rotation_deg', -1),
                                          (2, 'rotation_deg', 1)):
                axis = axes[axis_id]
                step = .05 if motion == 'translation_mm' else .25
                first = None
                for i in range(1, 21):
                    value = sign*step*i
                    posed = (removable.translate(axis*value) if motion == 'translation_mm'
                             else removable.transform(A.I.rotation_about(origin, axis, value)))
                    if A.I.intersection_volume(posed, added) > .001:
                        first = value
                        break
                stops.append(dict(axis=axis_id, motion=motion, sign=sign, first_stop=first))
            trials.append(dict(expansion_mm=reach, added_mm3=float(added.volume()),
                               components=len(candidate.decompose()), static_hits=hits, stops=stops))
        localization['rigid_locator_trials'] = trials
        # Register a rectangular tongue at a cap-facing vertex. Its matching
        # groove runs parallel to the already tested side-entry direction.
        frame_axes = np.stack(axes, axis=1)
        key_trials = []
        viable_keys = []
        contacts = []
        for vertex in vertices:
            hits = cap.ray_cast(vertex+direction*.05, vertex-direction*.6)
            if hits and np.linalg.norm(np.array(hits[0].position)-vertex) < .2:
                contacts.append(vertex)
        contacts.sort(key=lambda vertex: float(np.linalg.norm(vertex-origin)))
        for point in contacts[::max(1, len(contacts)//10)][:10]:
            for width in (.8, 1.2):
                frame = np.c_[frame_axes, point]
                tongue = mf.Manifold.cube([1.4, width, width]).translate(
                    [-1., -width/2, -width/2]).transform(frame)
                keyed = (removable+tongue).as_original().simplify(1e-6)
                side_clearance = .015
                groove = mf.Manifold.cube([13.52, width+2*side_clearance, width+2*side_clearance]).translate(
                    [-1.06, -width/2-side_clearance, -width/2-side_clearance]).transform(frame)
                cut_cap = (cap-groove).as_original().simplify(1e-6)
                key_stock = keyed-removable
                swept_key = mf.Manifold.batch_hull([key_stock, key_stock.translate(delta)])
                hits = {name: vol for name, part in fixed.items()
                        if name not in (cap_key, nut_key, washer_key)
                        and (vol := float(A.I.intersection_volume(swept_key, part))) > .001}
                groove_volume = float((cap-cut_cap).volume())
                row = dict(root_mm=point.tolist(), width_mm=width, cap_removed_mm3=groove_volume,
                           guide_components=len(keyed.decompose()), cap_components=len(cut_cap.decompose()),
                           key_entry_hits=hits, static_cap_overlap_mm3=float(A.I.intersection_volume(keyed, cut_cap)),
                           delayed_fastener_static_hits={
                               name: vol for name in (nut_key, washer_key)
                               if (vol := float(A.I.intersection_volume(key_stock, fixed[name]))) > .001})
                key_trials.append(row)
                if (hits or groove_volume < .1 or row['guide_components'] != 1 or row['cap_components'] != 1
                        or row['static_cap_overlap_mm3'] > .001 or row['delayed_fastener_static_hits']):
                    continue
                row['rigid_contact_probes'] = []
                viable_keys.append((point, tongue, groove))
                for axis_id, value, rotate in ((1, .1, False), (1, -1., True), (2, 1., True)):
                    axis = axes[axis_id]
                    posed = (keyed.transform(A.I.rotation_about(origin, axis, value)) if rotate
                             else keyed.translate(axis*value))
                    row['rigid_contact_probes'].append(float(A.I.intersection_volume(posed, cut_cap)))
        localization['sliding_key_trials'] = key_trials
        dual = []
        ranked = []
        for left, right in combinations(viable_keys, 2):
            separation = float(np.linalg.norm(left[0]-right[0]))
            if separation < .7:
                continue
            keyed = (removable+left[1]+right[1]).as_original().simplify(1e-6)
            cut_cap = (cap-left[2]-right[2]).as_original().simplify(1e-6)
            if len(keyed.decompose()) != 1 or len(cut_cap.decompose()) != 1:
                continue
            stops = []
            for axis_id, sign in ((1, -1), (2, 1)):
                first = next((angle for angle in np.arange(.25, 5.01, .25)
                              if A.I.intersection_volume(keyed.transform(
                                  A.I.rotation_about(origin, axes[axis_id], sign*angle)), cut_cap) > .001), None)
                stops.append(first)
            row = dict(separation_mm=separation, nominal_side_clearance_mm=.015, first_rotation_stop_deg=stops,
                       roots_mm=[left[0].tolist(), right[0].tolist()])
            dual.append(row)
            if all(stop is not None for stop in stops):
                ranked.append((row, keyed, cut_cap))
        localization['dual_key_trials'] = dual
        entry_obstacles = {key: part for key, part in fixed.items() if key not in delayed}
        for row, keyed, cut_cap in sorted(ranked, key=lambda item: max(item[0]['first_rotation_stop_deg'])):
            entry_obstacles[cap_key] = cut_cap
            hull = mf.Manifold.batch_hull([keyed, keyed.translate(delta)])
            nearby = {key: part for key, part in entry_obstacles.items()
                      if A.I.intersection_volume(hull, part) > .001}
            first = None
            for fraction in np.linspace(1., 0., 49):
                pose = keyed.translate(delta*fraction)
                hits = {key: vol for key, part in nearby.items()
                        if (vol := float(A.I.intersection_volume(pose, part))) > .001}
                if hits:
                    first = dict(distance_mm=float(np.linalg.norm(delta)*fraction), hits=hits)
                    break
            row['entry_readback'] = dict(step_mm=.25, first_hit=first)
            if first is not None:
                continue
            folder = base/'elbow_relocated_bone_seats'
            try:
                guide_actual, guide_source = export(keyed, folder/'extensor_dual_key_guide.obj', origin)
                _, cap_source = export(cut_cap, folder/'extensor_dual_key_cap.obj', origin)
            except (ValueError, AssertionError) as error:
                row['export_error'] = str(error)
                continue
            row['sources'] = dict(guide=guide_source, cap=cap_source)
            wire = mf.Manifold.batch_boolean(
                [A.D.outer_envelope(path, .3) for path in wire_paths], mf.OpType.Add)
            checks = []
            for axis_id, sign, limit in zip((1, 2), (-1, 1), row['first_rotation_stop_deg'], strict=True):
                for angle in np.arange(0., limit+.125, .25):
                    pose = guide_actual.transform(A.I.rotation_about(origin, axes[axis_id], sign*angle))
                    checks.append(dict(axis=axis_id, rotation_deg=float(sign*angle),
                                       wire_overlap_mm3=float(A.I.intersection_volume(pose, wire))))
            row['wire_motion_clearance'] = dict(
                checks=checks, wire_radius_mm=.3, elbow_step_deg=1., guide_step_deg=.25,
                scope='Fixed nominal cable path envelope; individual guide rotations about nominal pivot.')
            elbow_manifest = (ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces'
                              /'mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
                              /'continuous_channel_candidate/assembly.json')
            moving = {'bone:forearm'}
            moving.update(f"elbow:{part['key']}" for part in read(elbow_manifest)['instances']
                          if part['owner'] == 'elbow' or part['owner'].endswith('_distal_pivot'))
            moving_solids = {key: part for key, part in fixed.items() if key in moving}
            for head in ('FJ1486', 'FJ1477'):
                source = next((item['source'] for cable in bridge['combined_terminal_wire_audit']['rows']
                               for item in (cable['distal_feed_trial'], *cable['distal_feed_trial'].get(
                                   'aligned_sleeves', [])) if 'source' in item and cable['muscle'] ==
                               ('flexor' if head == 'FJ1486' else 'extensor')), None)
                key = f'elbow:{head}_active_wire_crimp_barrel'
                if head == 'FJ1477':
                    source = next(item['source'] for cable in bridge['combined_terminal_wire_audit']['rows']
                                  for item in cable['distal_feed_trial'].get('aligned_sleeves', [])
                                  if 'source' in item)
                moving_solids[key] = load_component(source) if source else fixed[key]
            elbow_centre, elbow_axis = A.frames.elbow_axis()
            motion_hits = []
            shoulder = {key: part for key, part in fixed.items()
                        if key not in moving_solids and not key.startswith(
                            ('elbow:', 'flexor_', 'extensor_plate:', 'tpu:', 'channel:'))
                        and key != 'bone:humerus'}
            shoulder_hits = []
            for angle in range(121):
                pose = A.I.rotation_about(elbow_centre, elbow_axis, angle)
                posed_parts = {key: part.transform(pose) for key, part in moving_solids.items()}
                hits = {key: vol for key, part in posed_parts.items()
                        if (vol := float(A.I.intersection_volume(guide_actual, part))) > .001}
                if hits:
                    motion_hits.append(dict(angle_deg=angle, hits=hits))
                cross_hits = [dict(moving=key, shoulder=target, overlap_mm3=volume)
                              for key, part in posed_parts.items() for target, obstacle in shoulder.items()
                              if (volume := float(A.I.intersection_volume(part, obstacle))) > .001]
                if cross_hits:
                    shoulder_hits.append(dict(angle_deg=angle, pairs=cross_hits))
            row['elbow_guide_sweep'] = dict(step_deg=1., samples=121, moving_solids=sorted(moving_solids),
                                           hits=motion_hits,
                                           scope='Shoulder neutral; forearm rigid movement; distal swivel angles held neutral.')
            row['elbow_to_shoulder_rigid_sweep'] = dict(
                step_deg=1., samples=121, shoulder_solids=sorted(shoulder), hits=shoulder_hits,
                scope='Rigid cross-module pairs only. Shared humerus and TPU are separate checks.')
            reference = A.DESIGN/'s1_y_muscle_pair/bare_motion_review/bare_joint_motion_checks.json'
            bare = read(reference)
            owner_file = A.DESIGN/'s1_integrated_humerus_anchors/assembly_check.json'
            owners = read(owner_file)['intended_body_owners']
            stationary = {key: part for key, part in shoulder.items()
                          if owners.get(key) in ('fixed_context', 'clavicle_ideal_fastened',
                                                 'scapula_ideal_fastened', 'nominal_scapula_clamp',
                                                 'nominal_clavicle_bench_clamp', 'nominal_bench_fastened',
                                                 'assumed_external_fixed_bench_region', 'front_liner_nominal_retainer')}
            coupled = []
            wire_routes = {}
            guide_point = np.array(bridge['extensor_crimp_elbow_sweep']
                                   ['extensor_fixed_guide_screen']['selected']['guide_source_mm'])
            for channel in channels['channels']:
                label = channel['label']
                start, end = np.array(channels['choice'][f'{label}_endpoints_source_mm'])
                neutral = (np.stack([start, guide_point, end]) if label == 'extensor'
                           else np.vstack([start, channel['path_source_mm'], end]))
                length = np.r_[0., np.cumsum(np.linalg.norm(np.diff(neutral, axis=0), axis=1))]
                fraction = np.array([0., 0., 1.]) if label == 'extensor' else length/length[-1]
                wire_routes[label] = neutral, fraction
            moving_wire_hits = []
            for angle in range(121):
                elbow = A.I.rotation_about(elbow_centre, elbow_axis, angle)
                posed_parts = {key: part.transform(elbow) for key, part in moving_solids.items()}
                for muscle, (neutral, fraction) in wire_routes.items():
                    end = elbow[:, :3]@neutral[-1]+elbow[:, 3]
                    route = neutral+np.outer(fraction, end-neutral[-1])
                    cable = A.D.outer_envelope(route, .3)
                    hits = {key: vol for key, part in posed_parts.items()
                            if (vol := float(A.I.intersection_volume(cable, part))) > .001}
                    if hits:
                        moving_wire_hits.append(dict(elbow_deg=angle, muscle=muscle, hits=hits))
            row['wire_to_moving_elbow_sweep'] = dict(
                samples=121, step_deg=1., obstacle_keys=sorted(moving_solids), hits=moving_wire_hits,
                scope='Nominal cable from anchor to anchor. Distal swivels held neutral. Shared GH transform cancels in these relative poses.')
            for label in ('flexion', 'abduction'):
                results = []
                wire_hits = []
                for angle in range(121):
                    elbow = A.I.rotation_about(elbow_centre, elbow_axis, angle)
                    gh = A.I.rotation_about(bare['joint_centre_source_mm'], bare['axes_source'][label], angle/4)
                    combined = gh[:, :3]@elbow
                    combined[:, 3] += gh[:, 3]
                    hits = []
                    for key, part in moving_solids.items():
                        posed = part.transform(combined)
                        hits.extend(dict(moving=key, shoulder=target, overlap_mm3=volume)
                                    for target, obstacle in stationary.items()
                                    if (volume := float(A.I.intersection_volume(posed, obstacle))) > .001)
                    if hits:
                        results.append(dict(gh_deg=angle/4, elbow_deg=angle, pairs=hits))
                    for muscle, (neutral, fraction) in wire_routes.items():
                        end = elbow[:, :3]@neutral[-1]+elbow[:, 3]
                        route = neutral+np.outer(fraction, end-neutral[-1])
                        route = route@gh[:, :3].T+gh[:, 3]
                        wire = A.D.outer_envelope(route, .3)
                        interference = {key: vol for key, part in stationary.items()
                                        if (vol := float(A.I.intersection_volume(wire, part))) > .001}
                        if interference:
                            wire_hits.append(dict(gh_deg=angle/4, elbow_deg=angle,
                                                  muscle=muscle, hits=interference))
                coupled.append(dict(gh_axis=label, poses=121, hits=results, wire_hits=wire_hits))
            row['coupled_rigid_sweep'] = dict(
                axis_reference_sha256=sha(reference), owner_reference_sha256=sha(owner_file),
                shoulder_keys=sorted(stationary), paths=coupled, gh_step_deg=.25, elbow_step_deg=1.,
                scope='Prescribed GH 0–30 degrees with elbow 0–120 degrees. Forearm solids and nominal elbow wires versus identified stationary shoulder solids. Flexor uses interpolated route, extensor uses fixed guide.')
            break
    return dict(checked_solids=len(solids), samples=49, step_mm=.25, failures=failures,
                first_reverse_path_block=failures[-1] if failures else None,
                delayed_fastener_trial=dict(guide_failures_remaining=remaining, fastener_entry=sequence,
                                            socket_proxy_hits=socket_hits, travel_mm=12.),
                washer_tongue_trials=retention,
                guide_localization=localization,
                scope='Includes both revised lattices, caps, adapters and both nominal bone-seat fastener sets. Dedicated guide retention unverified.')


def main():
    base = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    bridge_file = base / 'elbow_relocated_bridge_screen.json'
    bridge = read(bridge_file)
    selected = bridge['extensor_crimp_elbow_sweep']['extensor_fixed_guide_screen']['selected']
    channel_file = base / 'elbow_relocated_channel_cores/channel_core_check.json'
    channels = read(channel_file)
    proximal, distal = np.array(channels['choice']['extensor_endpoints_source_mm'])
    guide = np.array(selected['guide_source_mm'])
    incoming = A.D.unit(guide-proximal)
    centre, axis = A.frames.elbow_axis()
    mf = A.I._manifold()
    paths = []
    for angle in range(121):
        affine = A.I.rotation_about(centre, axis, angle)
        far = affine[:, :3]@distal+affine[:, 3]
        outgoing = A.D.unit(far-guide)
        paths.append(np.stack([guide-incoming*5., guide, guide+outgoing*5.]))
    relief = mf.Manifold.batch_boolean(
        [A.D.outer_envelope(path, .5) for path in paths], mf.OpType.Add)
    extensor = next(trial['source'] for row in bridge['extensor_lattice_midspan']['trials']
                    for terminal in row.get('terminal_trials', [])
                    for trial in terminal.get('entry_relief_trials', []) if 'source' in trial)
    muscle = load_component(extensor)
    channel = next(load_component(row['source']) for row in channels['channels']
                   if row['label'] == 'extensor')
    trials = []
    pockets = []
    folder = base/'elbow_relocated_bone_seats'
    for radius in (1.4, 1.6, 1.8, 2.):
        blank = rod(guide-incoming*1.5, guide+incoming*1.5, radius)
        body = (blank-relief).as_original().simplify(1e-6)
        trials.append(dict(outer_radius_mm=radius, length_mm=3., components=len(body.decompose()),
                           body_volume_mm3=float(body.volume()),
                           lattice_overlap_mm3=float(A.I.intersection_volume(body, muscle)),
                           continuous_channel_overlap_mm3=float(A.I.intersection_volume(body, channel)),
                           path_overlap_mm3=float(A.I.intersection_volume(body, relief))))
        pocket = rod(guide-incoming*1.55, guide+incoming*1.55, radius+.08)
        # Also clear the angular wire fan where it exits the short insert.
        local = rod(guide-incoming*2.5, guide+incoming*2.5, radius+1.)
        revised = (muscle-pocket-(relief^local)).as_original().simplify(1e-6)
        row = dict(outer_radius_mm=radius, radial_clearance_mm=.08,
                   removed_TPU_mm3=float((muscle-revised).volume()),
                   components=len(revised.decompose()),
                   insert_overlap_mm3=float(A.I.intersection_volume(revised, body)),
                   local_wire_envelope_overlap_mm3=float(A.I.intersection_volume(revised, relief^local)))
        pockets.append(row)
        if (not any('sources' in option for option in pockets)
                and row['components'] == 1 and row['insert_overlap_mm3'] <= .001
                and row['local_wire_envelope_overlap_mm3'] <= .001):
            try:
                _, tpu_source = export(revised, folder/'extensor_lattice_guide_pocket.obj', guide)
                _, guide_source = export(body, folder/'extensor_guide_insert.obj', guide)
            except (AssertionError, ValueError) as error:
                row['export_error'] = str(error)
            else:
                row['sources'] = dict(tpu=tpu_source, insert=guide_source)
    selected = next((row['sources'] for row in pockets if 'sources' in row), None)
    mounts = screen_mount(bridge, selected, guide, incoming, relief, folder) if selected else []
    mount = next((row['source'] for row in mounts if 'source' in row), None)
    entry = None
    cap_after_tpu = None
    detachable_guide = None
    if mount is not None:
        tpu = load_component(selected['tpu'])
        cap = load_component(mount)
        prior = bridge['paired_lattice_entry']['extensor_second']
        delta = np.array(prior['offset_mm'])
        failures = []
        for fraction in np.linspace(1., 0., 73):
            volume = float(A.I.intersection_volume(tpu.translate(delta*fraction), cap))
            if volume > .001:
                failures.append(dict(distance_from_seat_mm=float(np.linalg.norm(delta)*fraction),
                                     overlap_mm3=volume))
        entry = dict(offset_mm=delta.tolist(), samples=73, step_mm=.25,
                     prior_assembly_path_pass=prior['pass_sampled'],
                     pass_against_new_cap=not failures, failures=failures)
        cap_after_tpu = screen_cap_after_tpu(cap, tpu, bridge, incoming)
        reaction = next(row for row in bridge['proximal_reaction_seats'] if row['muscle'] == 'extensor')
        original_cap = load_component(next(flat['source'] for trial in reaction['trials']
                                            for flat in trial.get('flat_reaction_cap', []) if 'source' in flat))
        removable = (cap-original_cap).as_original().simplify(1e-6)
        detachable_guide = dict(components=len(removable.decompose()))
        if detachable_guide['components'] == 1:
            try:
                removable, detachable_guide['source'] = export(
                    removable, folder/'extensor_detachable_guide_trial.obj', guide)
            except (AssertionError, ValueError) as error:
                detachable_guide['export_error'] = str(error)
            else:
                detachable_guide['entry'] = screen_cap_after_tpu(
                    removable, tpu, bridge, incoming, {'extensor_fixed_cap': original_cap})
                detachable_guide['side_opening_trials'] = screen_side_opening(
                    removable, tpu, original_cap, bridge, detachable_guide['entry']['trials'], folder, guide)
                side = next((row for row in detachable_guide['side_opening_trials'] if 'source' in row), None)
                if side is not None:
                    detachable_guide['assembly_readback'] = verify_side_opening_assembly(
                        base, bridge, channels, removable, side, paths)
    save(base/'extensor_guide_body_audit.json', dict(
        bridge_report_sha256=sha(bridge_file), channels_sha256=sha(channel_file),
        guide_source_mm=guide.tolist(), samples=len(paths), angular_step_deg=1.,
        guide_lumen_sweep_radius_mm=.5, trials=trials, pocket_trials=pockets,
        mount_trials=mounts, existing_entry_readback=entry,
        cap_after_tpu=cap_after_tpu,
        detachable_guide=detachable_guide,
        scope='Neutral local guide mounting screen against extensor TPU, elbow adapters and local wire fan. Other solids and assembly path pending.'))
    print('mounts', len(mounts), 'selected', next((row for row in mounts if 'source' in row), None))


if __name__ == '__main__':
    main()
