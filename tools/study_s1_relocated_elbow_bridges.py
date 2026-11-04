"""Screen bridge stock between relocated original adapters and new bone seats."""
import numpy as np
from itertools import product

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha
from study_s1_elbow_proximal_stations import proximal_groups


def check_fasteners(seats, axes, solids, selected, compact=False, channel_entry=None, plate_hardware=None,
                    installed_flexor=None):
    mf = A.I._manifold()
    results = []
    for seat, axis in zip(seats, axes, strict=True):
        front, back = (np.array(seat[f'{side}_face_mm']) for side in ('front', 'back'))
        n = np.array(axis['axis'])
        narrow = compact and seat['muscle'] == 'extensor'
        washer_radius = 2.75 if narrow else 3.5
        washer_thickness = 1.5 if narrow and seat['anchor'] == 1 else .5
        head = front + n*3.5
        span = float((head-back) @ n)
        length = next(size for size in (25, 30, 35, 40, 45, 50)
                      if size >= span+washer_thickness+3.4)
        bore = rod(back-n*10, head+n*5, 1.7)
        side = np.eye(3)[np.argmin(abs(n))]
        u = A.D.unit(np.cross(n, side))
        nut_frame = np.c_[u, np.cross(-n, u), -n, back-n*washer_thickness]
        parts = {
            'bolt': rod(head, head-n*length, 1.5)+rod(head, head+n*3, 2.75),
            'front_washer': rod(front+n*3, head, washer_radius)-bore,
            'back_washer': rod(back, back-n*washer_thickness, washer_radius)-bore,
            'nut': mf.Manifold.cylinder(2.4, 5.5/np.sqrt(3), 5.5/np.sqrt(3), 6).transform(nut_frame)
                   - rod(back+n, back-n*6, 1.6),
        }
        tools = {
            'head_driver': rod(head+n*3.01, head+n*28, 2.),
            'rear_socket': rod(back-n*(washer_thickness+.01), back-n*28, 3.75 if narrow else 4.2)
                           -rod(back, back-n*29, 3.2),
            'bolt_insertion': rod(head, head+n*length, 2.75),
        }
        label = f"{seat['muscle']}_{seat['anchor']}"
        results.append(dict(label=label, nominal_length_mm=length, grip_to_back_face_mm=span,
                            washer_radius_mm=washer_radius, rear_washer_thickness_mm=washer_thickness,
                            adapter_selected=any(label.startswith(muscle) and key in selected
                                                 for muscle, key in [('flexor', 'elbow:FJ1486_bone_adapter'),
                                                                     ('extensor', 'elbow:FJ1477_bone_adapter')]),
                            parts=parts, tools=tools))
    hardware = {f"{row['label']}:{name}": part for row in results for name, part in row['parts'].items()}
    if installed_flexor is not None:
        installed_flexor.update((key, part) for key, part in hardware.items() if key.startswith('flexor_'))
    obstacles = {**solids, **hardware}
    if channel_entry is not None:
        fixed = {key: part for key, part in obstacles.items()
                 if key != 'channel:extensor' and not key.startswith(('extensor_0:', 'extensor_1:'))}
        fixed.update(plate_hardware or {})
        channel_entry.extend(check_channel_entry(solids['channel:extensor'], fixed, mf))
    for row in results:
        label = row['label']
        for group in ('parts', 'tools'):
            hits = {}
            for name, body in row.pop(group).items():
                hits[name] = {key: volume for key, solid in obstacles.items()
                              if not key.startswith(label+':')
                              and (volume := float(A.I.intersection_volume(body, solid))) > .001}
            row[f'{group}_hits'] = hits
    return results


def check_channel_entry(core, obstacles, mf):
    """Use a conservative convex sweep of the straight channel core."""
    rows = []
    radius = 18.
    for elevation in (-60, -30, 0, 30, 60):
        pitch = np.radians(elevation)
        for azimuth in range(0, 360, 15):
            yaw = np.radians(azimuth)
            direction = np.array([np.cos(pitch)*np.cos(yaw), np.cos(pitch)*np.sin(yaw), np.sin(pitch)])
            sweep = mf.Manifold.batch_hull([core, core.translate(radius*direction)])
            hits = {key: vol for key, obstacle in obstacles.items()
                    if (vol := float(A.I.intersection_volume(sweep, obstacle))) > .001}
            rows.append(dict(offset_mm=(radius*direction).tolist(), hits=hits))
    return rows


def build_extensor_lattice(channel, obstacles, folder, entry_paths):
    mf = A.I._manifold()
    start, end = np.array(channel['path_source_mm'])[[0, -1]]
    axis = A.D.unit(end-start)
    u = A.D.unit(np.cross(axis, np.eye(3)[np.argmin(abs(axis))]))
    v = np.cross(axis, u)
    core = obstacles['channel:extensor']
    lumen = A.D.outer_envelope(np.array([start-axis*4, end+axis*4]), channel['inner_radius_mm'])
    rows = []
    for peak_radius in (5.5, 4.5, 3.5):
        rings, spokes = [], []
        count = 8
        for index, t in enumerate(np.linspace(.025, .975, 17)):
            centre = start+t*(end-start)
            theta = 2*np.pi*(np.arange(count)+.5*(index%2))/count
            radius = 1.4+(peak_radius-1.4)*np.sin(np.pi*t)**.75
            directions = np.outer(np.cos(theta), u)+.8*np.outer(np.sin(theta), v)
            ring = centre+radius*directions
            rings.append(ring)
            if index % 2 == 0:
                spokes.extend(rod(centre+direction*1.2, point, .55)
                              for direction, point in zip(directions, ring, strict=True))
        struts = []
        for index, (lower, upper) in enumerate(zip(rings, rings[1:])):
            for j in range(count):
                for k in (j, (j+(1 if index%2 else -1))%count):
                    struts.append(rod(lower[j], upper[k], .55))
        nodes = [mf.Manifold.sphere(.62, 12).translate(point) for ring in rings for point in ring]
        lattice = (mf.Manifold.batch_boolean([core, *struts, *spokes, *nodes], mf.OpType.Add)
                   -lumen).as_original().simplify(1e-6)
        raw = lattice
        adapter = obstacles['elbow:FJ1486_bone_adapter']
        crop = rod(start-axis*3, end+axis*3, peak_radius+2.)
        clearance = (adapter^crop).minkowski_sum(mf.Manifold.cube([.2]*3, True))
        lattice = (raw-clearance).as_original().simplify(1e-6)
        hits = {key: vol for key, part in obstacles.items() if key != 'channel:extensor'
                and (vol := float(A.I.intersection_volume(lattice, part))) > .001}
        rows.append(dict(peak_radius_mm=peak_radius, hits=hits, components=len(lattice.decompose()),
                         relieved_mm3=float((raw-lattice).volume()),
                         channel_wall_loss_mm3=float((core-lattice).volume())))
        if not hits and rows[-1]['components'] == 1 and rows[-1]['channel_wall_loss_mm3'] <= .001:
            try:
                actual, source = export(lattice, folder/'extensor_lattice_midspan.obj', start)
            except (AssertionError, ValueError) as error:
                rows[-1]['export_error'] = str(error)
                continue
            rows[-1]['source'] = source
            rows[-1]['wire_overlap_mm3'] = float(A.I.intersection_volume(
                actual, A.D.outer_envelope(np.array([start, end]), .3)))
            terminals = []
            for distal_extension in (7., 5., 3., 1., 0.):
                terminal_start, terminal_end = start-axis*7., end+axis*distal_extension
                passage = A.D.outer_envelope(
                    np.array([terminal_start-axis*2, terminal_end+axis*2]), channel['inner_radius_mm'])
                sleeves = [rod(terminal_start, start+axis, 1.6),
                           rod(end-axis, terminal_end, 1.6)]
                body = (mf.Manifold.batch_boolean([actual, *sleeves], mf.OpType.Add)
                        -passage).as_original().simplify(1e-6)
                hits = {key: vol for key, part in obstacles.items() if key != 'channel:extensor'
                        and (vol := float(A.I.intersection_volume(body, part))) > .001}
                terminal = dict(proximal_extension_mm=7., distal_extension_mm=distal_extension,
                                remaining_nominal_endpoint_gap_mm=[1., 8.-distal_extension],
                                hits=hits, components=len(body.decompose()))
                terminals.append(terminal)
                if hits or terminal['components'] != 1:
                    continue
                body, terminal['source'] = export(
                    body, folder/f'extensor_lattice_distal_extension_{int(distal_extension)}mm.obj', start)
                terminal['entry_samples'] = []
                for path in entry_paths:
                    if path['hits']:
                        continue
                    delta = np.array(path['offset_mm'])
                    sweep = mf.Manifold.batch_hull([body, body.translate(delta)])
                    sweep_hits = {key: vol for key, part in obstacles.items() if key != 'channel:extensor'
                                  and (vol := float(A.I.intersection_volume(sweep, part))) > .001}
                    first_hit = None
                    for distance in np.linspace(0., np.linalg.norm(delta), 37):
                        pose = body.translate(A.D.unit(delta)*distance)
                        hits = {key: vol for key in sweep_hits
                                if (vol := float(A.I.intersection_volume(pose, obstacles[key]))) > .001}
                        if hits:
                            first_hit = dict(distance_from_seat_mm=float(distance), hits=hits)
                            break
                    terminal['entry_samples'].append(dict(
                        offset_mm=delta.tolist(), convex_sweep_hits=sweep_hits,
                        sample_step_mm=.5, first_sampled_hit=first_hit))
                protected = mf.Manifold.batch_boolean([core, *sleeves], mf.OpType.Add)-passage
                relief_trials = []
                for path in sorted(terminal['entry_samples'],
                                   key=lambda item: sum(item['convex_sweep_hits'].values())):
                    delta = np.array(path['offset_mm'])
                    sweep = mf.Manifold.batch_hull([body, body.translate(delta)])
                    travel = rod(np.zeros(3), -delta, .05)
                    cutters = [(obstacles[key]^sweep).minkowski_sum(travel)
                               for key in path['convex_sweep_hits']]
                    relief = mf.Manifold.batch_boolean(cutters, mf.OpType.Add)
                    candidate = (body-relief).as_original().simplify(1e-6)
                    loss = float((protected-candidate).volume())
                    trial = dict(offset_mm=delta.tolist(), protected_loss_mm3=loss,
                                 removed_mm3=float((body-candidate).volume()),
                                 components=len(candidate.decompose()))
                    relief_trials.append(trial)
                    if loss > .001 or trial['components'] != 1:
                        continue
                    try:
                        _, trial['source'] = export(candidate, folder/'extensor_lattice_entry_relief.obj', start)
                    except (AssertionError, ValueError) as error:
                        trial['export_error'] = str(error)
                        continue
                    break
                terminal['entry_relief_trials'] = relief_trials
                if any('source' in trial for trial in relief_trials):
                    break
            rows[-1]['terminal_trials'] = terminals
            break
    return rows


def verify_lattice_entry(trials, obstacles):
    selected = next((item for row in trials for terminal in row.get('terminal_trials', [])
                     for item in terminal.get('entry_relief_trials', []) if 'source' in item), None)
    if selected is None:
        return dict(pass_sampled=False, reason='No exported entry-relief candidate')
    body = load_component(selected['source'])
    delta = np.array(selected['offset_mm'])
    mf = A.I._manifold()
    hull = mf.Manifold.batch_hull([body, body.translate(delta)])
    nearby = {key: part for key, part in obstacles.items() if key != 'channel:extensor'
              and A.I.intersection_volume(hull, part) > .001}
    failures = []
    for fraction in np.linspace(1., 0., 73):
        pose = body.translate(delta*fraction)
        hits = {key: vol for key, part in nearby.items()
                if (vol := float(A.I.intersection_volume(pose, part))) > .001}
        if hits:
            failures.append(dict(distance_from_seat_mm=float(np.linalg.norm(delta)*fraction), hits=hits))
    return dict(pass_sampled=not failures, checked_obstacles=len(obstacles)-1,
                obstacles_inside_conservative_hull=sorted(nearby), samples=73,
                nominal_step_mm=.25, offset_mm=delta.tolist(), failures=failures)


def sample_flexor_entry(body, trials, obstacles):
    spacer = 'elbow:elbow_spacer_left'
    results = []
    for trial in trials:
        delta = np.array(trial['offset_mm'])
        nearby = {key: obstacles[key] for key in trial['hits'] if key != 'tpu:extensor_relocated'}
        spacer_hit = None
        blocked = None
        for step in range(37):
            pose = body.translate(delta*(step/36))
            if spacer in nearby and spacer_hit is None:
                overlap = float(A.I.intersection_volume(pose, nearby[spacer]))
                if overlap > .001:
                    spacer_hit = dict(distance_from_seat_mm=step*.5, overlap_mm3=overlap)
            hits = {key: overlap for key, part in nearby.items() if key != spacer
                    and (overlap := float(A.I.intersection_volume(pose, part))) > .001}
            if hits:
                blocked = dict(distance_from_seat_mm=step*.5, hits=hits)
                break
        results.append(dict(offset_mm=trial['offset_mm'], spacer_first_hit=spacer_hit,
                            other_first_hit=blocked))
    return dict(sample_step_mm=.5, travel_mm=18., candidates=results,
                passing_before_spacer=sum(row['other_first_hit'] is None for row in results),
                passing_after_spacer=sum(row['other_first_hit'] is None and row['spacer_first_hit'] is None
                                         for row in results))


def audit_terminal_seats(channel_data, obstacles, extensor_distal_extension):
    rows = []
    for channel in channel_data['channels']:
        muscle = channel['label']
        ends = np.array(channel['path_source_mm'])[[0, -1]]
        endpoints = np.array(channel_data['choice'][f'{muscle}_endpoints_source_mm'])
        for side, core_end, anchor in zip(('proximal', 'distal'), ends, endpoints, strict=True):
            axis = A.D.unit(anchor-core_end)
            extension = extensor_distal_extension if (muscle, side) == ('extensor', 'distal') else 7.
            tip = core_end+axis*extension
            # A collar around the free span needs space independently of the
            # running cable. Probe beyond the anchor for rigid seating stock.
            sleeve = rod(tip, anchor+axis*.25, 1.6)-rod(tip-axis, anchor+axis, .85)
            hits, seat_surfaces = {}, {}
            for key, solid in obstacles.items():
                if not key.startswith('elbow:'):
                    continue
                volume = float(A.I.intersection_volume(sleeve, solid))
                if volume > .001:
                    hits[key] = volume
                for hit in solid.ray_cast(tip, anchor+axis*8):
                    distance = float((np.array(hit.position)-tip)@axis)
                    seat_surfaces.setdefault(key, []).append(distance)
            rows.append(dict(muscle=muscle, side=side, tip_source_mm=tip.tolist(),
                             anchor_source_mm=anchor.tolist(), outward_axis=axis.tolist(),
                             nominal_gap_mm=float(np.linalg.norm(anchor-tip)), collar_hits=hits,
                             centreline_surfaces_from_tip_mm=seat_surfaces))
    return rows


def screen_split_reaction_seat(seat, tip, axis, obstacles, folder, muscle, cap_key):
    mf = A.I._manifold()
    u = A.D.unit(np.cross(axis, np.eye(3)[np.argmin(abs(axis))]))
    v = np.cross(axis, u)
    trials = []
    for angle in range(0, 180, 15):
        theta = np.radians(angle)
        lateral = u*np.cos(theta)+v*np.sin(theta)
        frame = np.c_[lateral, np.cross(axis, lateral), axis, tip]
        halves = [seat ^ mf.Manifold.cube([12., 24., 24.])
                  .translate([offset, -12., -12.]).transform(frame) for offset in (-12., 0.)]
        sweeps = []
        paths = []
        for sign, half in zip((-1, 1), halves, strict=True):
            direction = lateral*sign
            hull = mf.Manifold.batch_hull([half, half.translate(direction*8.)])
            sweeps.append(hull)
            nearby = {key: solid for key, solid in obstacles.items()
                      if A.I.intersection_volume(hull, solid) > .001}
            first = None
            for distance in np.linspace(0., 8., 33):
                moved = half.translate(direction*distance)
                hits = {key: vol for key, solid in nearby.items()
                        if (vol := float(A.I.intersection_volume(moved, solid))) > .001}
                if hits:
                    first = dict(distance_mm=float(distance), hits=hits)
                    break
            paths.append(dict(direction=direction.tolist(), first_hit=first,
                              components=len(half.decompose())))
        row = dict(split_deg=angle, travel_mm=8., sample_step_mm=.25, paths=paths)
        trials.append(row)
        if all(path['components'] == 1 and path['first_hit'] is None for path in paths):
            row['sources'] = [export(half.as_original().simplify(1e-6),
                                     folder/f'{muscle}_reaction_half_{i}.obj', tip)[1]
                              for i, half in enumerate(halves)]
            break
        if any(path['components'] != 1 for path in paths):
            continue
        other_hits = sorted({key for sweep in sweeps for key, part in obstacles.items()
                             if key != cap_key and A.I.intersection_volume(sweep, part) > .001})
        row['slot_other_envelope_hits'] = other_hits
        if other_hits:
            continue
        cap = obstacles[cap_key]
        entry_stock = mf.Manifold.batch_boolean(sweeps, mf.OpType.Add)
        removed = cap ^ entry_stock
        revised = (cap-entry_stock).as_original().simplify(1e-6)
        prefix = cap_key.removesuffix('_rigid_cap')
        guards = {}
        for key, part in obstacles.items():
            if not key.startswith(prefix) or not any(word in key for word in ('washer', 'bolt', 'nut')):
                continue
            lo, hi = np.array(part.bounding_box()).reshape(2, 3)
            guard = mf.Manifold.cube(hi-lo+2.).translate(lo-1.)
            volume = float(A.I.intersection_volume(removed, guard))
            if volume > .001:
                guards[key] = volume
        slot = dict(removed_mm3=float(removed.volume()), cap_components=len(revised.decompose()),
                    hardware_guard_hits=guards)
        row['cap_slot'] = slot
        if slot['cap_components'] != 1 or guards:
            continue
        try:
            _, slot['source'] = export(revised, folder/f'{muscle}_side_entry_cap.obj', tip)
            slot['halves'] = [export(half.as_original().simplify(1e-6),
                                     folder/f'{muscle}_side_entry_half_{i}.obj', tip)[1]
                              for i, half in enumerate(halves)]
        except (ValueError, AssertionError) as error:
            slot['export_error'] = str(error)
            continue
        break
    return trials


def screen_open_reaction_shim(seat, tip, axis, obstacles, folder, muscle, cap_key):
    mf = A.I._manifold()
    u = A.D.unit(np.cross(axis, np.eye(3)[np.argmin(abs(axis))]))
    v = np.cross(axis, u)
    trials = []
    for angle in range(0, 360, 15):
        theta = np.radians(angle)
        lateral = u*np.cos(theta)+v*np.sin(theta)
        frame = np.c_[lateral, np.cross(axis, lateral), axis, tip]
        keep = mf.Manifold.cube([12., 12., 12.]).translate([-6., -6., .06]).transform(frame)
        slot = mf.Manifold.cube([6., 1.2, 14.]).translate([-6., -.6, -1.]).transform(frame)
        shim = ((seat^keep)-slot).as_original().simplify(1e-6)
        row = dict(angle_deg=angle, components=len(shim.decompose()))
        trials.append(row)
        if row['components'] != 1:
            continue
        hull = mf.Manifold.batch_hull([shim, shim.translate(lateral*8)])
        nearby = {key: part for key, part in obstacles.items() if A.I.intersection_volume(hull, part) > .001}
        cap = obstacles[cap_key]
        for distance in np.linspace(0., 8., 33):
            pose = shim.translate(lateral*distance)
            hits = {key: vol for key, part in nearby.items() if key != cap_key
                    if (vol := float(A.I.intersection_volume(pose, part))) > .001}
            if hits:
                row['first_hit'] = dict(distance_mm=float(distance), hits=hits)
                break
        else:
            removed = cap ^ hull
            revised = (cap-hull).as_original().simplify(1e-6)
            prefix = cap_key.removesuffix('_rigid_cap')
            guards = {}
            for key, part in obstacles.items():
                if not key.startswith(prefix) or not any(word in key for word in ('washer', 'bolt', 'nut')):
                    continue
                lo, hi = np.array(part.bounding_box()).reshape(2, 3)
                guard = mf.Manifold.cube(hi-lo+2.).translate(lo-1.)
                volume = float(A.I.intersection_volume(removed, guard))
                if volume > .001:
                    guards[key] = volume
            support = float(A.I.intersection_volume(shim.translate(axis*.1), revised))/.1
            row.update(cap_removed_mm3=float(removed.volume()), cap_components=len(revised.decompose()),
                       hardware_guard_hits=guards, axial_bearing_proxy_mm2=support)
            if guards or row['cap_components'] != 1 or support < 1.:
                continue
            try:
                _, row['source'] = export(shim, folder/f'{muscle}_open_reaction_shim.obj', tip)
                _, row['cap_source'] = export(revised, folder/f'{muscle}_open_shim_entry_cap.obj', tip)
            except (AssertionError, ValueError) as error:
                row.pop('source', None)
                row['export_error'] = str(error)
                continue
            row['nominal_axial_gap_mm'] = .06
            break
    return dict(travel_mm=8., sample_step_mm=.25, slot_width_mm=1.2, trials=trials)


def build_proximal_reaction_seats(terminals, obstacles, folder, entry_paths):
    results = []
    for row in terminals:
        if row['side'] != 'proximal':
            continue
        muscle = row['muscle']
        head = {'flexor': 'FJ1486', 'extensor': 'FJ1477'}[muscle]
        cap_key = f'elbow:{head}_proximal_rigid_cap'
        cap = obstacles[cap_key]
        tip = np.array(row['tip_source_mm'])
        axis = np.array(row['outward_axis'])
        trials = []
        for radius in (2.4, 2., 1.8):
            blank = rod(tip-axis*1.5, tip+axis*8., radius)
            socket = rod(tip-axis*2, tip+axis*.05, 1.65)
            lumen = rod(tip-axis*2, tip+axis*9., 1.)
            raw = blank-socket-lumen
            shared = float(A.I.intersection_volume(raw, cap))
            seat = (raw-cap).as_original().simplify(1e-6)
            components = seat.decompose()
            mouth = rod(tip-axis*1.5, tip-axis*.5, radius)-socket
            retained = [part for part in components if A.I.intersection_volume(mouth, part) > .001]
            trial = dict(radius_mm=radius, cap_intersection_before_fit_mm3=shared,
                         components_before_mouth_selection=len(components),
                         socket_side_components=len(retained))
            trials.append(trial)
            if len(retained) != 1 or shared <= .001:
                continue
            seat = retained[0]
            trial['stock_beyond_cap_removed_mm3'] = float(sum(part.volume() for part in components)-seat.volume())
            hits = {key: vol for key, part in obstacles.items()
                    if (vol := float(A.I.intersection_volume(seat, part))) > .001}
            trial['hits'] = hits
            if hits:
                continue
            try:
                _, trial['source'] = export(seat, folder/f'{muscle}_proximal_reaction_seat.obj', tip)
            except (ValueError, AssertionError) as error:
                trial['export_error'] = str(error)
                continue
            trial['split_seat_trials'] = screen_split_reaction_seat(seat, tip, axis, obstacles, folder, muscle, cap_key)
            trial['open_shim'] = screen_open_reaction_shim(seat, tip, axis, obstacles, folder, muscle, cap_key)
            mf = A.I._manifold()
            body = obstacles.get(f'tpu:{muscle}_relocated')
            flat_trials = []
            for gap in (.06, .25, .5, .75, 1.):
                face = tip+axis*gap
                stock = rod(face, tip+axis*8., 2.4)-lumen
                # Keep the original cap material, including its original bore.
                flat = (cap+stock).as_original().simplify(1e-6)
                added = flat-cap
                hits = {key: vol for key, part in obstacles.items() if key != cap_key
                        and (vol := float(A.I.intersection_volume(added, part))) > .001}
                item = dict(axial_gap_mm=gap, hits=hits, components=len(flat.decompose()))
                flat_trials.append(item)
                if body is None or hits or item['components'] != 1:
                    continue
                item['paths'] = []
                for offset in entry_paths.get(muscle, []):
                    delta = np.array(offset)
                    first_hit = next((float(np.linalg.norm(delta)*t) for t in np.linspace(0., 1., 73)
                                      if A.I.intersection_volume(body.translate(delta*t), added) > .001), None)
                    item['paths'].append(dict(offset_mm=offset, first_hit_mm=first_hit))
                if any(path['first_hit_mm'] is None for path in item['paths']):
                    try:
                        _, item['source'] = export(flat, folder/f'{muscle}_flat_reaction_cap.obj', tip)
                    except (AssertionError, ValueError) as error:
                        item['export_error'] = str(error)
                        continue
                    item['bearing_area_proxy_mm2'] = float(
                        A.I.intersection_volume(body.translate(axis*(gap+.1)), added))/.1
                    break
            trial['flat_reaction_cap'] = flat_trials
            # Grow only toward the cap so the socket mouth stays in place.
            fused = (cap+seat+seat.translate(axis*.2)).as_original().simplify(1e-6)
            added = fused-cap
            integration_hits = {key: vol for key, part in obstacles.items() if key != cap_key
                                and (vol := float(A.I.intersection_volume(added, part))) > .001}
            integration = dict(components=len(fused.decompose()), hits=integration_hits,
                               added_mm3=float(added.volume()))
            trial['integral_cap'] = integration
            if integration['components'] != 1 or integration_hits:
                break
            try:
                actual, integration['source'] = export(fused, folder/f'{muscle}_cap_with_reaction_seat.obj', tip)
            except (ValueError, AssertionError) as error:
                integration['export_error'] = str(error)
                break
            body = obstacles.get(f'tpu:{muscle}_relocated')
            if body is not None:
                integration['prior_entry_paths'] = []
                for offset in entry_paths.get(muscle, []):
                    delta = np.array(offset)
                    first = next((dict(distance_mm=float(np.linalg.norm(delta)*t), overlap_mm3=vol)
                                  for t in np.linspace(0., 1., 73)
                                  if (vol := float(A.I.intersection_volume(body.translate(delta*t), actual))) > .001),
                                 None)
                    integration['prior_entry_paths'].append(dict(offset_mm=offset, first_hit=first))
            break
        results.append(dict(muscle=muscle, trials=trials))
    return results


def audit_combined_terminal_wires(reactions, channels, obstacles, folder):
    combined = dict(obstacles)
    original_caps = {}
    for reaction in reactions:
        muscle = reaction['muscle']
        head = {'flexor': 'FJ1486', 'extensor': 'FJ1477'}[muscle]
        key = f'elbow:{head}_proximal_rigid_cap'
        selected = next((flat['source'] for trial in reaction['trials']
                         for flat in trial.get('flat_reaction_cap', []) if 'source' in flat), None)
        if selected is None:
            return dict(complete=False, missing=muscle)
        original_caps[key] = combined[key]
        combined[key] = load_component(selected)
    rows = []
    for channel in channels['channels']:
        muscle = channel['label']
        head = {'flexor': 'FJ1486', 'extensor': 'FJ1477'}[muscle]
        endpoints = np.array(channels['choice'][f'{muscle}_endpoints_source_mm'])
        path = np.vstack([endpoints[0], channel['path_source_mm'], endpoints[1]])
        wire = A.D.outer_envelope(path, .3)
        wire_hits = {}
        for key, solid in combined.items():
            overlap = float(A.I.intersection_volume(wire, solid))
            if overlap > .001:
                wire_hits[key] = overlap
        added_hits = {key: float(A.I.intersection_volume(wire, combined[key]-part))
                      for key, part in original_caps.items()}
        # Report the load-transfer region separately from the free cable.
        distal = endpoints[1]
        outgoing = A.D.unit(distal-path[-2])
        crimp_key = f'elbow:{head}_active_wire_crimp_barrel'
        crimp = combined[crimp_key]
        crimp_probe = crimp.ray_cast(distal-outgoing*2., distal+outgoing*8.)
        crossings = [float((np.array(hit.position)-distal) @ outgoing) for hit in crimp_probe]
        far = max(crossings, default=0.)+.5
        tip = distal+outgoing*far
        feed = rod(distal, tip, .3)
        bore = rod(distal-outgoing, tip+outgoing, .5)
        revised = (crimp-bore).as_original().simplify(1e-6)
        other_hits = {key: vol for key, part in combined.items() if key != crimp_key
                      and (vol := float(A.I.intersection_volume(feed, part))) > .001}
        trial = dict(feed_length_mm=far, bore_diameter_mm=1., other_feed_hits=other_hits,
                     components=len(revised.decompose()), removed_mm3=float((crimp-revised).volume()),
                     residual_feed_overlap_mm3=float(A.I.intersection_volume(feed, revised)),
                     nominal_centreline_stock_length_mm=(max(crossings)-min(crossings)) if len(crossings) >= 2 else None)
        if crossings and not other_hits and trial['components'] == 1 and trial['residual_feed_overlap_mm3'] <= .001:
            try:
                _, trial['source'] = export(revised, folder/f'{muscle}_distal_feed_crimp.obj', distal)
            except (AssertionError, ValueError) as error:
                trial['export_error'] = str(error)
        if muscle == 'extensor' and len(crossings) >= 2:
            trial['aligned_sleeves'] = []
            for length in (4., 3., 2.5):
                for radius in (1.5, 1.25, 1.):
                    start = distal+outgoing*(min(crossings)-.25)
                    end = start+outgoing*length
                    outer = rod(start, end, radius)
                    passage = rod(distal-outgoing, end+outgoing, .5)
                    candidate = ((crimp+outer)-passage).as_original().simplify(1e-6)
                    added = candidate-crimp
                    extended_feed = rod(distal, end+outgoing*.5, .3)
                    hits = {key: vol for key, part in combined.items() if key != crimp_key
                            and (vol := float(A.I.intersection_volume(added+extended_feed, part))) > .001}
                    item = dict(sleeve_length_mm=length, radial_wall_mm=radius-.5,
                                components=len(candidate.decompose()), hits=hits,
                                feed_overlap_mm3=float(A.I.intersection_volume(extended_feed, candidate)))
                    trial['aligned_sleeves'].append(item)
                    if hits or item['components'] != 1 or item['feed_overlap_mm3'] > .001:
                        continue
                    try:
                        _, item['source'] = export(candidate, folder/'extensor_distal_aligned_crimp.obj', distal)
                    except (AssertionError, ValueError) as error:
                        item['export_error'] = str(error)
                        continue
                    break
                if 'source' in item:
                    break
        rows.append(dict(muscle=muscle, free_wire_hits=wire_hits, added_cap_wire_overlap_mm3=added_hits,
                         crimp_centreline_crossings_from_anchor_mm=crossings, distal_feed_trial=trial))
    return dict(complete=True, checked_solids=len(combined), rows=rows,
                scope='Nominal cable between anchors. Crimp bore and feed path beyond the anchors are unverified.')


def screen_extensor_detour(path, centre, axis, fixed):
    fractions = np.linspace(0., 1., 33)
    chord = A.D.unit(path[-1]-path[0])
    u = A.D.unit(np.cross(chord, np.eye(3)[np.argmin(abs(chord))]))
    v = np.cross(chord, u)
    neutral = path[0]+np.outer(fractions, path[-1]-path[0])
    rows = []
    for displacement in (2., 4., 6., 8.):
        for azimuth in range(0, 360, 45):
            theta = np.radians(azimuth)
            offset = displacement*(np.cos(theta)*u+np.sin(theta)*v)
            def check(angle):
                affine = A.I.rotation_about(centre, axis, angle)
                distal = affine[:, :3]@path[-1]+affine[:, 3]
                posed = (neutral+np.outer(fractions, distal-path[-1])
                         +np.outer(np.sin(np.pi*fractions)**2*(angle/120.), offset))
                wire = A.D.outer_envelope(posed, .3)
                hits = {key: vol for key, part in fixed.items()
                        if (vol := float(A.I.intersection_volume(wire, part))) > .001}
                return dict(angle_deg=angle, hits=hits,
                            length_mm=float(np.linalg.norm(np.diff(posed, axis=0), axis=1).sum()))
            row = dict(max_offset_mm=displacement, azimuth_deg=azimuth, offset_source_mm=offset.tolist(),
                       at_120_deg=check(120))
            rows.append(row)
            if row['at_120_deg']['hits']:
                continue
            row['first_hit'] = next((result for angle in range(121)
                                     if (result := check(angle))['hits']), None)
            if row['first_hit'] is None:
                return dict(selected=row, trials=rows, scope='Prescribed smooth path only. Guide hardware and TPU excluded.')
    return dict(selected=None, trials=rows, scope='Prescribed smooth path only. Guide hardware and TPU excluded.')


def screen_fixed_extensor_guide(path, centre, axis, fixed, offset):
    proximal, distal = path[[0, -1]]
    direction = A.D.unit(offset)
    trials = []
    for distance in (8., 16., 24., 32., 48.):
        for shift in (1., 2., 4., 6.):
            guide = proximal+A.D.unit(distal-proximal)*distance+direction*shift
            def check(angle):
                transform = A.I.rotation_about(centre, axis, angle)
                far = transform[:, :3]@distal+transform[:, 3]
                route = np.stack([proximal, guide, far])
                wire = A.D.outer_envelope(route, .3)
                hits = {key: vol for key, part in fixed.items()
                        if (vol := float(A.I.intersection_volume(wire, part))) > .001}
                legs = np.diff(route, axis=0)
                turn = float(np.degrees(np.arccos(np.clip(A.D.unit(legs[0])@A.D.unit(legs[1]), -1., 1.))))
                return dict(angle_deg=angle, hits=hits, turn_deg=turn,
                            length_mm=float(np.linalg.norm(legs, axis=1).sum()))
            row = dict(guide_source_mm=guide.tolist(), distance_mm=distance, shift_mm=shift,
                       end_poses=[check(angle) for angle in (0, 120)])
            trials.append(row)
            if any(pose['hits'] for pose in row['end_poses']):
                continue
            samples = [check(angle) for angle in range(121)]
            row['first_hit'] = next((sample for sample in samples if sample['hits']), None)
            row['maximum_turn_deg'] = max(sample['turn_deg'] for sample in samples)
            if row['first_hit'] is None:
                row['length_change_mm'] = samples[-1]['length_mm']-samples[0]['length_mm']
                return dict(selected=row, trials=trials)
    return dict(selected=None, trials=trials)


def sweep_extensor_crimp(wire_report, obstacles, channels, reactions, folder):
    source = next((item['source'] for row in wire_report.get('rows', []) if row['muscle'] == 'extensor'
                   for item in row['distal_feed_trial'].get('aligned_sleeves', []) if 'source' in item), None)
    if source is None:
        return dict(checked=False)
    elbow_manifest = (ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces'
                      /'mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
                      /'continuous_channel_candidate/assembly.json')
    owners = {f"elbow:{row['key']}": row['owner'] for row in read(elbow_manifest)['instances']}
    moving = {key for key, owner in owners.items() if owner == 'elbow' or owner.endswith('_distal_pivot')}
    moving.add('bone:forearm')
    key = 'elbow:FJ1477_active_wire_crimp_barrel'
    moving.add(key)
    moving.add('elbow:FJ1486_active_wire_crimp_barrel')
    candidate = load_component(source)
    added = candidate-obstacles[key]
    centre, axis = A.frames.elbow_axis()
    fixed = {name: part for name, part in obstacles.items()
             if name not in moving and not name.startswith(('tpu:', 'channel:'))}
    for reaction in reactions:
        head = {'flexor': 'FJ1486', 'extensor': 'FJ1477'}[reaction['muscle']]
        source = next((flat['source'] for trial in reaction['trials']
                       for flat in trial.get('flat_reaction_cap', []) if 'source' in flat), None)
        if source is not None:
            fixed[f'elbow:{head}_proximal_rigid_cap'] = load_component(source)
    rows = []
    wire_rows = {}
    adapter_relief = None
    detour = None
    fixed_guide = None
    for channel in channels['channels']:
        muscle = channel['label']
        endpoints = np.array(channels['choice'][f'{muscle}_endpoints_source_mm'])
        path = np.vstack([endpoints[0], channel['path_source_mm'], endpoints[1]])
        distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
        fraction = distance/distance[-1]
        def sample_wire(angle):
            affine = A.I.rotation_about(centre, axis, angle)
            distal = affine[:, :3]@path[-1]+affine[:, 3]
            posed = path+np.outer(fraction, distal-path[-1])
            length = float(np.linalg.norm(np.diff(posed, axis=0), axis=1).sum())
            sample = dict(angle_deg=angle, path_length_mm=length,
                          change_from_neutral_mm=length-float(distance[-1]))
            cable = A.D.outer_envelope(posed, .3)
            sample['rigid_hits'] = {name: vol for name, part in fixed.items()
                                    if (vol := float(A.I.intersection_volume(cable, part))) > .001}
            return sample
        samples = [sample_wire(angle) for angle in range(121)]
        first = next((sample for sample in samples if sample['rigid_hits']), None)
        fine = []
        if first is not None and first['angle_deg'] > 0:
            end = first['angle_deg']*10
            fine = [sample_wire(step/10) for step in range(end-10, end+1)]
        wire_rows[muscle] = dict(
            first_sampled_fixed_hit=next((sample for sample in fine if sample['rigid_hits']), first),
            boundary_samples=fine,
            poses=samples)
        if muscle == 'extensor':
            detour = screen_extensor_detour(path, centre, axis, fixed)
            if detour['selected'] is not None:
                fixed_guide = screen_fixed_extensor_guide(
                    path, centre, axis, fixed, detour['selected']['offset_source_mm'])
            mf = A.I._manifold()
            adapter_key = 'elbow:FJ1486_bone_adapter'
            adapter = fixed[adapter_key]
            cutters = []
            for step in range(1160, 1201):
                transform = A.I.rotation_about(centre, axis, step/10)
                distal = transform[:, :3]@path[-1]+transform[:, 3]
                posed = path+np.outer(fraction, distal-path[-1])
                cutters.append(A.D.outer_envelope(posed, .55))
            clearance = mf.Manifold.batch_boolean(cutters, mf.OpType.Add)
            revised = (adapter-clearance).as_original().simplify(1e-6)
            removed = adapter-revised
            guards = {}
            for name, part in fixed.items():
                if not (name.startswith('flexor_') or name.startswith('elbow:FJ1486_proximal_')):
                    continue
                lo, hi = np.array(part.bounding_box()).reshape(2, 3)
                guard = mf.Manifold.cube(hi-lo+1.).translate(lo-.5)
                volume = float(A.I.intersection_volume(removed, guard))
                if volume > .001:
                    guards[name] = volume
            adapter_relief = dict(sample_step_deg=.1, angle_range_deg=[116., 120.],
                                  wire_radius_mm=.3, cutting_radius_mm=.55,
                                  removed_mm3=float(removed.volume()), components=len(revised.decompose()),
                                  hardware_guard_hits=guards,
                                  other_obstacles=sorted({key for sample in samples for key in sample['rigid_hits']
                                                          if key != adapter_key}))
            if adapter_relief['components'] == 1 and not guards:
                try:
                    _, adapter_relief['source'] = export(revised, folder/'flexor_adapter_wire_relief.obj', centre)
                except (AssertionError, ValueError) as error:
                    adapter_relief['export_error'] = str(error)
    for angle in range(121):
        posed = added.transform(A.I.rotation_about(centre, axis, angle))
        hits = {name: vol for name, part in fixed.items()
                if (vol := float(A.I.intersection_volume(posed, part))) > .001}
        rows.append(dict(angle_deg=angle, hits=hits))
    return dict(checked=True, assembly_owner_sha256=sha(elbow_manifest), checked_fixed_solids=len(fixed),
                step_deg=1., first_hit=next((row for row in rows if row['hits']), None), poses=rows,
                wire_proxy=dict(wire_rows=wire_rows, length_step_deg=1., collision_step_deg=1., boundary_step_deg=.1,
                                model='Linear blend of distal anchor displacement along neutral cable. Fixed rigid obstacles only.'),
                flexor_adapter_relief=adapter_relief,
                extensor_detour_screen=detour,
                extensor_fixed_guide_screen=fixed_guide,
                scope='Added stock only, nominal forearm ownership. Shoulder held at zero. TPU and actual swivel motion excluded.')


def build_flexor_lattice(channel, obstacles, folder, endpoints):
    mf = A.I._manifold()
    path = np.array(channel['path_source_mm'])
    stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
    t = np.linspace(.025, .975, 19)
    centres = np.column_stack([np.interp(t*stations[-1], stations, coord) for coord in path.T])
    tangents = np.gradient(centres, axis=0)
    chord = A.D.unit(path[-1]-path[0])
    reference = np.eye(3)[np.argmin(abs(chord))]
    core = obstacles['channel:flexor']
    lumen_path = np.vstack([path[0]-4*A.D.unit(path[1]-path[0]), path,
                            path[-1]+4*A.D.unit(path[-1]-path[-2])])
    lumen = A.D.outer_envelope(lumen_path, channel['inner_radius_mm'])
    trials = []
    for radius in (5.5, 4.5, 3.5):
        rings, members = [], []
        for i, (point, tangent) in enumerate(zip(centres, tangents, strict=True)):
            tangent = A.D.unit(tangent)
            u = A.D.unit(np.cross(tangent, reference))
            v = np.cross(tangent, u)
            theta = (np.arange(8)+.5*(i%2))*np.pi/4
            radial = np.outer(np.cos(theta), u)+.8*np.outer(np.sin(theta), v)
            ring = point+(1.4+(radius-1.4)*np.sin(np.pi*t[i])**.75)*radial
            rings.append(ring)
            members.extend(mf.Manifold.sphere(.62, 12).translate(p) for p in ring)
            if i%2 == 0:
                members.extend(rod(point+1.2*d, p, .55) for d, p in zip(radial, ring, strict=True))
        for i, (lower, upper) in enumerate(zip(rings, rings[1:])):
            for j in range(8):
                members.extend(rod(lower[j], upper[k], .55) for k in (j, (j+(1 if i%2 else -1))%8))
        raw = mf.Manifold.batch_boolean([core, *members], mf.OpType.Add)-lumen
        nearby = {key: part for key, part in obstacles.items() if key != 'channel:flexor'
                  and A.I.intersection_volume(raw, part) > .001}
        crop = A.D.outer_envelope(path, radius+2.)
        margin = mf.Manifold.cube([.2]*3, True)
        relief = mf.Manifold.batch_boolean([(part^crop).minkowski_sum(margin)
                                            for part in nearby.values()], mf.OpType.Add)
        body = (raw-relief).as_original().simplify(1e-6)
        trial = dict(radius_mm=radius, relieved_against=sorted(nearby),
                     removed_mm3=float((raw-body).volume()), components=len(body.decompose()),
                     wall_loss_mm3=float((core-body).volume()))
        trials.append(trial)
        if trial['components'] != 1 or trial['wall_loss_mm3'] > .001:
            continue
        try:
            actual, trial['source'] = export(body, folder/'flexor_lattice_midspan.obj', path[0])
        except (AssertionError, ValueError) as error:
            trial['export_error'] = str(error)
            continue
        trial['readback_hits'] = {key: vol for key, part in obstacles.items() if key != 'channel:flexor'
                                 and (vol := float(A.I.intersection_volume(actual, part))) > .001}
        trial['wire_overlap_mm3'] = float(A.I.intersection_volume(actual, A.D.outer_envelope(path, .3)))
        if not trial['readback_hits'] and trial['wire_overlap_mm3'] <= .001:
            trial['terminal_trials'] = []
            ends = path[[0, -1]]
            leads = np.array(endpoints)-ends
            lengths = np.linalg.norm(leads, axis=1)
            directions = leads / lengths[:, None]
            for extension in (7., 5., 3.):
                tips = ends+directions*extension
                full_path = np.vstack([tips[0]+directions[0]*2, path,
                                       tips[1]+directions[1]*2])
                passage = A.D.outer_envelope(full_path, channel['inner_radius_mm'])
                sleeves = [rod(tip, end-direction, channel['outer_radius_mm'])
                           for tip, end, direction in zip(tips, ends, directions, strict=True)]
                extended = (mf.Manifold.batch_boolean([actual, *sleeves], mf.OpType.Add)
                            -passage).as_original().simplify(1e-6)
                terminal = dict(extension_mm=extension, remaining_endpoint_gap_mm=(lengths-extension).tolist(),
                                components=len(extended.decompose()))
                trial['terminal_trials'].append(terminal)
                if terminal['components'] != 1:
                    continue
                terminal['hits'] = {key: vol for key, part in obstacles.items() if key != 'channel:flexor'
                                    and (vol := float(A.I.intersection_volume(extended, part))) > .001}
                terminal['wire_overlap_mm3'] = float(A.I.intersection_volume(
                    extended, A.D.outer_envelope(full_path, .3)))
                if terminal['hits'] or terminal['wire_overlap_mm3'] > .001:
                    continue
                try:
                    _, terminal['source'] = export(extended, folder/'flexor_lattice_extended.obj', path[0])
                except (AssertionError, ValueError) as error:
                    terminal['export_error'] = str(error)
                    continue
                break
            break
    return trials


def check_extensor_thread_plate(seats, axes, solids, folder, installed):
    """Screen a shared nut plate and the two independent bearing spacers."""
    mf = A.I._manifold()
    pairs = [(seat, axis) for seat, axis in zip(seats, axes, strict=True)
             if seat['muscle'] == 'extensor']
    n = np.array(pairs[0][1]['axis'])
    backs = np.array([seat['back_face_mm'] for seat, _ in pairs])
    if any(np.linalg.norm(n-np.array(axis['axis'])) > 1e-4 for _, axis in pairs):
        raise ValueError('Shared plate requires parallel screw axes')
    bone = solids['bone:humerus']
    candidates = []
    selected = None
    for offset in (.5, 1., 1.5, 2., 3.):
        near = backs - np.outer(backs@n - min(backs@n) + offset, n)
        blanks = [rod(p, p-n*3., 3.) for p in near]
        plate = mf.Manifold.batch_hull(blanks)
        holes = [rod(p+n*4., p-n*5., 1.5) for p in near]
        plate -= mf.Manifold.batch_boolean(holes, mf.OpType.Add)
        spacers = [(rod(back+n*.2, point, 2.75)-rod(back+n, point-n, 1.7))-bone
                   for back, point in zip(backs, near, strict=True)]
        parts = {'plate': plate, **{f'spacer_{i}': part for i, part in enumerate(spacers)}}
        hits = {name: {key: vol for key, solid in solids.items()
                       if (vol := float(A.I.intersection_volume(part, solid))) > .001}
                for name, part in parts.items()}
        # Check translation of the plate as a solid before the screws go in.
        tangent = A.D.unit(near[1]-near[0])
        entry = {}
        for side in (-1, 1):
            sweep = plate.minkowski_sum(rod(np.zeros(3), tangent*side*18., .025))
            entry[str(side)] = {key: vol for key, solid in solids.items()
                                if (vol := float(A.I.intersection_volume(sweep, solid))) > .001}
        counts = {key: len(part.decompose()) for key, part in parts.items()}
        candidates.append(dict(offset_mm=offset, plate_face_centres_mm=near.tolist(),
                               component_counts=counts, part_hits=hits, entry_hits=entry))
        if selected is None and all(count == 1 for count in counts.values()) and not any(hits.values()) and not entry['1']:
            sources = {}
            for name, part in parts.items():
                actual, source = export(part.as_original().simplify(1e-6),
                                        folder / f'extensor_thread_{name}.obj', near[0])
                assert all(A.I.intersection_volume(actual, obstacle) <= .001 for obstacle in solids.values())
                sources[name] = source
                installed[f'extensor_thread:{name}'] = actual
            selected = dict(offset_mm=offset, sources=sources, entry_direction=tangent.tolist())
            placed = dict(solids)
            sequence = []
            for index, spacer in enumerate(spacers):
                # A swept solid cannot use radial inflation at its final
                # conformal bearing surface, where contact is intentional.
                sweep = mf.Manifold.batch_boolean([
                    spacer.translate(-n*distance) for distance in np.linspace(15., 0., 61)
                ], mf.OpType.Add)
                collisions = {key: vol for key, obstacle in placed.items()
                              if (vol := float(A.I.intersection_volume(sweep, obstacle))) > .001}
                sequence.append(dict(spacer=index, travel_mm=15., sample_step_mm=.25,
                                     hits=collisions))
                placed[f'spacer_{index}'] = spacer
            selected['spacer_insertion'] = sequence
            selected['plate_spacer_static_overlap_mm3'] = [
                float(A.I.intersection_volume(plate, spacer)) for spacer in spacers]
            bolts = []
            for index, ((seat, _), point) in enumerate(zip(pairs, near, strict=True)):
                front = np.array(seat['front_face_mm'])
                head = front+n*3.5
                grip = float((head-point) @ n)
                length = next(size for size in (25, 30, 35, 40, 45, 50) if size >= grip+4.)
                # Nominal screw and plate geometry. Thread flank geometry
                # is omitted, so verify shank fit against the surrounding solids.
                bolt = rod(head, head-n*length, 1.5)+rod(head, head+n*3, 2.75)
                washer = rod(front+n*3, head, 2.75)-rod(front, head+n, 1.7)
                installed[f'extensor_thread:bolt_{index}'] = bolt
                installed[f'extensor_thread:washer_{index}'] = washer
                hits = {key: vol for key, obstacle in {**placed, 'plate': plate}.items()
                        if (vol := float(A.I.intersection_volume(bolt, obstacle))) > .001}
                bolts.append(dict(anchor=index, nominal_length_mm=length, grip_mm=grip,
                                  plate_engagement_mm=min(3., length-grip),
                                  protrusion_mm=length-grip-3., static_hits=hits))
            selected['bolts'] = bolts
    return dict(candidates=candidates, selected=selected)


def main():
    base = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = base / 'assembly_tpu_readback_candidate.json'
    seat_file = base / 'elbow_relocated_bone_seats/bone_seat_check.json'
    axes_file = base / 'elbow_relocated_seat_audit.json'
    channels = read(base / 'elbow_relocated_channel_cores/channel_core_check.json')
    stations = read(base / 'elbow_proximal_rotated_station_screen.json')
    seats = read(seat_file)
    assert seats['assembly_sha256'] == sha(manifest)
    assert seats['seat_report_sha256'] == sha(axes_file)
    axes = read(axes_file)['anchors']
    frame = read(A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json')['frame_source']
    axial = np.array(frame['basis_columns'])[:, 2]
    solids = {key: load_component(source) for key, source in read(manifest)['sources'].items()
              if key not in ('tpu:elbow_flexor', 'tpu:elbow_extensor', 'bone:humerus')}
    solids['bone:humerus'] = load_component(seats['source'])
    for row in channels['channels']:
        solids[f"channel:{row['label']}"] = load_component(row['source'])
    targets = {}
    corrected_groups = proximal_groups(stations['moving_groups'])
    for muscle, head in [('flexor', 'FJ1486'), ('extensor', 'FJ1477')]:
        pair = channels['choice']['pair']
        pose = A.I.rotation_about(frame['origin_mm'], axial, pair[f'{muscle}_angle_deg'])
        pose[:, 3] += axial * pair[f'{muscle}_shift_mm']
        for key in corrected_groups[head]:
            if '_bone_' in key and not key.endswith('_bone_adapter'):
                del solids[key]
            else:
                solids[key] = solids[key].transform(pose)
        targets[muscle] = f'elbow:{head}_bone_adapter'
    mf = A.I._manifold()
    results = []
    options = {}
    for seat, axis in zip(seats['seats'], axes, strict=True):
        key = targets[seat['muscle']]
        old = solids[key]
        front = np.array(seat['front_face_mm'])
        normal = np.array(axis['axis'])
        start = np.array(axis['front_washer_source_mm']) - normal*1.8
        end = front + normal*1.5
        bore = rod(front-normal, front+normal*9, 1.7)
        seat_ear = rod(front, front+normal*3, 4.5)
        relief = mf.Manifold()
        if seat['muscle'] == 'extensor' and seat['anchor'] == 0:
            crop = mf.Manifold.cube([18., 18., 18.], True).translate(front)
            margin = mf.Manifold.cube([.1, .1, .1], True)
            relief = mf.Manifold.batch_boolean([
                (solids[name]^crop).minkowski_sum(margin)
                for name in ('bone:humerus', 'distal:FJ1468/backing_plate')
            ], mf.OpType.Add)
        washer_land = rod(front+normal*2.5, front+normal*3, 3.5)-bore
        tangent = A.D.unit(np.cross(normal, axial))
        rows = []
        for side in (-1, 1):
            entry = start + side*tangent*3.4
            for bow, arrival in [(b, a) for b in (0., 4., 8., 12.)
                                 for a in (0., -35., 35., -70., 70., -105., 105.)]:
                angle = np.radians(arrival)
                radial = side*tangent*np.cos(angle) + axial*np.sin(angle)
                exit_point = end + radial*3
                approach = end + radial*5.5 + normal*4
                mid = (entry+approach)/2 + normal*bow
                points = (entry, mid, approach, exit_point)
                pieces = [seat_ear, *[rod(a, b, 2.5) for a, b in zip(points, points[1:])],
                          *[mf.Manifold.sphere(2.5, 24).translate(p) for p in points[1:-1]]]
                raw = mf.Manifold.batch_boolean(pieces, mf.OpType.Add)-bore
                bridge = raw-relief
                merged = old+bridge
                connected = len(merged.decompose()) == 1
                hits = []
                for other, part in solids.items():
                    if other == key:
                        continue
                    volume = float(A.I.intersection_volume(bridge, part))
                    if volume > .001:
                        hits.append(dict(key=other, overlap_mm3=volume))
                candidate = dict(side=side, bow_mm=bow, arrival_deg=arrival, connected=connected, hits=hits,
                                 relief_removed_mm3=float((raw-bridge).volume()),
                                 washer_land_loss_mm3=float((washer_land-bridge).volume()),
                                 points_mm=[point.tolist() for point in points])
                rows.append(candidate)
                if connected and not hits and candidate['washer_land_loss_mm3'] <= .001:
                    options.setdefault((seat['muscle'], seat['anchor']), []).append((candidate, bridge))
        ear_hits = {name: volume for name, part in solids.items() if name != key
                    and (volume := float(A.I.intersection_volume((seat_ear-bore)-relief, part))) > .001}
        results.append(dict(muscle=seat['muscle'], anchor=seat['anchor'], candidates=rows,
                            seat_ear_hits=ear_hits,
                            passing=len(options.get((seat['muscle'], seat['anchor']), []))))
    selected = {}
    folder = base / 'elbow_relocated_bone_seats'
    for muscle, key in targets.items():
        pockets = []
        for seat, axis in zip(seats['seats'], axes, strict=True):
            if seat['muscle'] == muscle:
                front = np.array(seat['front_face_mm'])
                normal = np.array(axis['axis'])
                pockets.append(rod(front+normal*3, front+normal*40, 3.6))
        access = mf.Manifold.batch_boolean(pockets, mf.OpType.Add)
        for (first, left), (second, right) in product(options.get((muscle, 0), []),
                                                      options.get((muscle, 1), [])):
            uncut = solids[key]+left+right
            merged = (uncut-access).as_original().simplify(1e-6)
            if len(merged.decompose()) != 1:
                continue
            if any(A.I.intersection_volume(merged, part) > .001 for name, part in solids.items()
                   if name != key and name not in selected):
                continue
            try:
                actual, source = export(merged, folder / f'{muscle}_bridge_adapter.obj', np.zeros(3))
            except (AssertionError, ValueError):
                continue
            if any(A.I.intersection_volume(actual, part) > .001 for name, part in solids.items()
                   if name != key):
                continue
            solids[key] = actual
            selected[key] = dict(source=source, bridges=[first, second],
                                 access_removed_mm3=float((uncut-actual).volume()))
            break
    channel_entry = []
    plate_hardware = {}
    flexor_hardware = {}
    plate_report = check_extensor_thread_plate(seats['seats'], axes, solids, folder, plate_hardware)
    fasteners = check_fasteners(seats['seats'], axes, solids, selected, channel_entry=channel_entry,
                                plate_hardware=plate_hardware, installed_flexor=flexor_hardware)
    lattice = build_extensor_lattice(next(row for row in channels['channels'] if row['label'] == 'extensor'),
                                    {**solids, **plate_hardware}, folder, channel_entry)
    obstacles = {**solids, **plate_hardware, **flexor_hardware}
    for row in lattice:
        for terminal in row.get('terminal_trials', []):
            for trial in terminal.get('entry_relief_trials', []):
                if 'source' in trial:
                    obstacles['tpu:extensor_relocated'] = load_component(trial['source'])
    flexor = build_flexor_lattice(next(row for row in channels['channels'] if row['label'] == 'flexor'),
                                 obstacles, folder, channels['choice']['flexor_endpoints_source_mm'])
    flexor_source = next((trial['source'] for row in flexor for trial in row.get('terminal_trials', [])
                          if 'source' in trial), None)
    extensor_terminal = next((terminal for row in lattice for terminal in row.get('terminal_trials', [])
                              if any('source' in trial for trial in terminal.get('entry_relief_trials', []))), None)
    terminal_audit = audit_terminal_seats(
        channels, obstacles, extensor_terminal['distal_extension_mm'] if extensor_terminal else 7.)
    terminal_obstacles = dict(obstacles)
    if flexor_source is not None:
        terminal_obstacles['tpu:flexor_relocated'] = load_component(flexor_source)
    pair_entry = None
    entry_paths = {}
    if flexor_source is not None and 'tpu:extensor_relocated' in obstacles:
        flexor_body = load_component(flexor_source)
        fixed = {key: solid for key, solid in obstacles.items()
                 if key not in ('channel:flexor', 'channel:extensor')}
        extensor_second = verify_lattice_entry(
            lattice, {**solids, **plate_hardware, **flexor_hardware, 'tpu:flexor_relocated': flexor_body})
        # The conservative convex envelope can reject a feasible nonconvex
        # path. Keep such rejections distinct from sampled solid collisions.
        flexor_second = check_channel_entry(flexor_body, fixed, mf)
        flexor_first = [
            row for row in flexor_second
            if not (set(row['hits']) - {'tpu:extensor_relocated'})
        ]
        pair_entry = dict(extensor_second=extensor_second, flexor_trials=flexor_second,
                          flexor_first_clear_envelopes=len(flexor_first),
                          flexor_second_clear_envelopes=sum(not row['hits'] for row in flexor_second),
                          flexor_first_sampled=sample_flexor_entry(flexor_body, flexor_second, fixed))
        entry_paths['extensor'] = [extensor_second['offset_mm']] if extensor_second['pass_sampled'] else []
        entry_paths['flexor'] = [row['offset_mm'] for row in pair_entry['flexor_first_sampled']['candidates']
                                 if row['other_first_hit'] is None and row['spacer_first_hit'] is None]
    reactions = build_proximal_reaction_seats(terminal_audit, terminal_obstacles, folder, entry_paths)
    wire_report = audit_combined_terminal_wires(reactions, channels, terminal_obstacles, folder)
    save(base / 'elbow_relocated_bridge_screen.json', dict(
        assembly_sha256=sha(manifest), seats_sha256=sha(seat_file), results=results,
        restored_distal_keys=sorted(set(sum(stations['moving_groups'].values(), []))
                                    - set(sum(corrected_groups.values(), []))),
        selected_adapters=selected,
        fasteners=fasteners,
        extensor_core_entry_after_flexor_hardware=dict(
            travel_mm=18., convex_outer_sweep=True, candidates=channel_entry,
            passing=sum(not row['hits'] for row in channel_entry),
            installed_extensor_hardware=sorted(plate_hardware),
            scope='Core only, includes both revised fastener sets when selected. Regenerated TPU lattice excluded.'),
        compact_fasteners=check_fasteners(seats['seats'], axes, solids, selected, compact=True),
        extensor_thread_plate=plate_report,
        extensor_lattice_midspan=dict(trials=lattice, strut_radius_mm=.55,
                                    scope='Midspan lattice and continuous channel. Terminal interfaces, flexor hardware and assembly path pending.'),
        extensor_lattice_entry_readback=verify_lattice_entry(
            lattice, {**solids, **plate_hardware, **flexor_hardware}),
        flexor_lattice_midspan=flexor,
        paired_lattice_entry=pair_entry,
        terminal_seat_audit=terminal_audit,
        proximal_reaction_seats=reactions,
        combined_terminal_wire_audit=wire_report,
        extensor_crimp_elbow_sweep=sweep_extensor_crimp(wire_report, terminal_obstacles, channels, reactions, folder),
        scope='Neutral bridge checks with nominal M3 hardware and straight tool envelopes. No thread, installation sequence, wire or motion verification.'))
    print([(row['muscle'], row['anchor'], row['passing']) for row in results], list(selected))


if __name__ == '__main__':
    main()
