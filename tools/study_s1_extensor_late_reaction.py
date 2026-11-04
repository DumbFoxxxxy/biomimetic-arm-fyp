"""Test post-installation of the extensor reaction face on the corrected elbow."""
import numpy as np
from scipy.optimize import brentq

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha
from audit_s1_extensor_guide_body import load_relocated_assembly, screen_cap_after_tpu


def verify_sequence(solids, sequence):
    mf = A.I._manifold()
    rows = []
    for key, part, offset in sequence:
        delta = np.array(offset)
        hull = mf.Manifold.batch_hull([part, part.translate(delta)])
        nearby = {name: other for name, other in solids.items()
                  if A.I.intersection_volume(hull, other) > .001}
        count = int(round(np.linalg.norm(delta)/.25))+1
        failures = []
        for t in np.linspace(1., 0., count):
            pose = part.translate(delta*t)
            hits = {name: vol for name, other in nearby.items()
                    if (vol := float(A.I.intersection_volume(pose, other))) > .001}
            if hits:
                failures.append(dict(distance_mm=float(np.linalg.norm(delta)*t), hits=hits))
        rows.append(dict(part=key, obstacles=len(solids), samples=count, failures=failures))
        solids[key] = part
    return rows


def check_revised_elbow_motion(solids, channels, guide_point, assembly_sources):
    owner_file = (ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces'
                  /'mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
                  /'continuous_channel_candidate/assembly.json')
    moving_keys = {'bone:forearm'}
    moving_keys.update(f"elbow:{row['key']}" for row in read(owner_file)['instances']
                       if row['owner'] == 'elbow' or row['owner'].endswith('_distal_pivot'))
    moving = {key: part for key, part in solids.items() if key in moving_keys}
    local_keys = ('elbow:FJ1477_proximal_rigid_cap', 'extensor:dual_key_guide',
                  'extensor:late_reaction_seat', 'extensor:stop_washer', 'extensor:retention_pin')
    local = {key: solids[key] for key in local_keys}
    centre, axis = A.frames.elbow_axis()
    bare = read(A.DESIGN/'s1_y_muscle_pair/bare_motion_review/bare_joint_motion_checks.json')
    owner_reference = read(A.DESIGN/'s1_integrated_humerus_anchors/assembly_check.json')['intended_body_owners']
    stationary_owners = {'fixed_context', 'clavicle_ideal_fastened', 'scapula_ideal_fastened',
                         'nominal_scapula_clamp', 'nominal_clavicle_bench_clamp', 'nominal_bench_fastened',
                         'assumed_external_fixed_bench_region', 'front_liner_nominal_retainer'}
    shoulder = {key: solid for key, solid in solids.items() if owner_reference.get(key) in stationary_owners}
    routes = {}
    swivels = {}
    for channel in channels['channels']:
        label = channel['label']
        start, end = np.array(channels['choice'][f'{label}_endpoints_source_mm'])
        if label == 'extensor':
            routes[label] = np.stack([start, guide_point, end]), np.array([0., 0., 1.])
        else:
            points = np.vstack([start, channel['path_source_mm'], end])
            lengths = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
            routes[label] = points, lengths/lengths[-1]
        head = {'flexor': 'FJ1486', 'extensor': 'FJ1477'}[label]
        key = f'elbow:{head}_cap_rotor'
        source = assembly_sources[key]
        definition = read((ROOT/source['file']).parent/'screen.json')
        pivot = np.array(definition['pivot_zero_display_mm'])+np.array(source['translation_source_mm'])
        hinge = np.array(definition['trunnion_axis_zero_display'])
        pack = [f"elbow:{item['key']}" for item in read(owner_file)['instances']
                if item['owner'] == f'{head}_distal_pivot']
        swivels[label] = (pivot, hinge, pack)
    coupled = {label: [] for label in ('flexion', 'abduction')}
    fixed_elbow = {key: part for key, part in solids.items() if key not in moving
                   and (key.startswith(('elbow:', 'extensor:', 'extensor_plate:', 'flexor_'))
                        or key == 'bone:humerus')}
    internal_wire_hits = []
    swivel_trials = {label: [] for label in routes}
    simultaneous_hits = []
    equilibrium_hits = []
    equilibrium_samples = 0
    rows = []
    for angle in range(121):
        transform = A.I.rotation_about(centre, axis, angle)
        moving_posed = {key: part.transform(transform) for key, part in moving.items()}
        hits = []
        for key, posed in moving_posed.items():
            hits.extend(dict(moving=key, fixed=target, overlap_mm3=volume)
                        for target, obstacle in local.items()
                        if (volume := float(A.I.intersection_volume(posed, obstacle))) > .001)
        if hits:
            rows.append(dict(angle_deg=angle, hits=hits))
        posed_wires = {}
        for muscle, (points, fraction) in routes.items():
            far = transform[:, :3]@points[-1]+transform[:, 3]
            posed_wires[f'wire:{muscle}'] = A.D.outer_envelope(points+np.outer(fraction, far-points[-1]), .3)
        for name, cable in posed_wires.items():
            collisions = {key: volume for key, part in {**fixed_elbow, **moving_posed}.items()
                          if (volume := float(A.I.intersection_volume(cable, part))) > .001}
            if collisions:
                internal_wire_hits.append(dict(elbow_deg=angle, cable=name, hits=collisions))
        tracked_packs, tracked_wires = {}, {}
        equilibrium_packs, equilibrium_wires = {}, {}
        for muscle, (pivot, hinge, pack) in swivels.items():
            points, fraction = routes[muscle]
            far = transform[:, :3]@points[-1]+transform[:, 3]
            preceding = points[-2]+fraction[-2]*(far-points[-1])
            reference = points[-2]-points[-1]
            direction = transform[:, :3].T@(preceding-far)
            reference -= hinge*(reference@hinge)
            direction -= hinge*(direction@hinge)
            swing = float(np.degrees(np.arctan2(hinge@np.cross(reference, direction), reference@direction)))
            swivel = A.I.rotation_about(pivot, hinge, swing)
            end_local = swivel[:, :3]@points[-1]+swivel[:, 3]
            end = transform[:, :3]@end_local+transform[:, 3]
            route = points+np.outer(fraction, end-points[-1])
            cable = A.D.outer_envelope(route, .3)
            length = float(np.linalg.norm(np.diff(route, axis=0), axis=1).sum())
            neutral_length = float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())
            pull = transform[:, :3].T@A.D.unit(route[-2]-route[-1])
            swivel_arm = float(np.cross(end_local-pivot, pull)@hinge)
            def torque_at(value):
                candidate = A.I.rotation_about(pivot, hinge, value)
                anchor = candidate[:, :3]@points[-1]+candidate[:, 3]
                far_point = transform[:, :3]@anchor+transform[:, 3]
                previous = points[-2]+fraction[-2]*(far_point-points[-1])
                load = transform[:, :3].T@A.D.unit(previous-far_point)
                return float(np.cross(anchor-pivot, load)@hinge)
            lo, hi = swing-20., swing+20.
            equilibrium = None
            if torque_at(lo) > 0. > torque_at(hi):
                stable = float(brentq(torque_at, lo, hi, xtol=1e-7))
                turn = A.I.rotation_about(pivot, hinge, stable)
                anchor = turn[:, :3]@points[-1]+turn[:, 3]
                far_point = transform[:, :3]@anchor+transform[:, 3]
                stable_route = points+np.outer(fraction, far_point-points[-1])
                stable_wire = A.D.outer_envelope(stable_route, .3)
                rotor_key = next(key for key in pack if key.endswith('_cap_rotor'))
                rotor = solids[rotor_key].transform(turn).transform(transform)
                equilibrium = dict(swivel_deg=stable, deviation_deg=stable-swing,
                                   rotor_wire_overlap_mm3=float(A.I.intersection_volume(stable_wire, rotor)))
                if angle % 5 == 0:
                    sensitivity = []
                    for error in (-5., -3., -1., 1., 3., 5.):
                        perturbed = A.I.rotation_about(pivot, hinge, stable+error)
                        anchor = perturbed[:, :3]@points[-1]+perturbed[:, 3]
                        far_point = transform[:, :3]@anchor+transform[:, 3]
                        varied_wire = A.D.outer_envelope(
                            points+np.outer(fraction, far_point-points[-1]), .3)
                        varied_rotor = solids[rotor_key].transform(perturbed).transform(transform)
                        sensitivity.append(dict(offset_deg=error,
                                                balancing_moment_per_newton_Nmm=-torque_at(stable+error),
                                                wire_overlap_mm3=float(A.I.intersection_volume(
                                                    varied_wire, varied_rotor))))
                    equilibrium['sensitivity'] = sensitivity
                equilibrium_wires[f'wire:{muscle}'] = stable_wire
                equilibrium_packs[muscle] = {key: solids[key].transform(turn).transform(transform)
                                            for key in pack if key in solids}
            posed_pack = {key: solids[key].transform(swivel).transform(transform) for key in pack if key in solids}
            tracked_packs[muscle] = posed_pack
            tracked_wires[muscle] = cable
            obstacles = {**fixed_elbow, **moving_posed, **posed_pack}
            wire_hits = {key: vol for key, part in obstacles.items()
                         if (vol := float(A.I.intersection_volume(cable, part))) > .001}
            first_rigid = None
            if not any(row.get('rigid_first_hit') for row in swivel_trials[muscle]):
                first_rigid = next((dict(part=key, obstacle=target, overlap_mm3=vol)
                                    for key, part in posed_pack.items() for target, obstacle in obstacles.items()
                                    if target not in posed_pack
                                    and (vol := float(A.I.intersection_volume(part, obstacle))) > .001), None)
            swivel_trials[muscle].append(dict(elbow_deg=angle, prescribed_swivel_deg=swing,
                                              wire_hits=wire_hits, rigid_first_hit=first_rigid,
                                              cable_length_mm=length, cable_stroke_mm=length-neutral_length,
                                              tension_moment_per_newton_Nmm=swivel_arm,
                                              local_cable_only_equilibrium=equilibrium))
        flexor = {**tracked_packs['flexor'], 'wire:flexor': tracked_wires['flexor']}
        extensor = {**tracked_packs['extensor'], 'wire:extensor': tracked_wires['extensor']}
        cross_swivel = [dict(flexor=key, extensor=target, overlap_mm3=vol)
                        for key, part in flexor.items() for target, obstacle in extensor.items()
                        if (vol := float(A.I.intersection_volume(part, obstacle))) > .001]
        if cross_swivel:
            simultaneous_hits.append(dict(elbow_deg=angle, pairs=cross_swivel))
        if len(equilibrium_packs) == 2:
            equilibrium_samples += 1
            flex_pack, ext_pack = equilibrium_packs['flexor'], equilibrium_packs['extensor']
            all_packs = {**flex_pack, **ext_pack}
            surroundings = {key: part for key, part in {**fixed_elbow, **moving_posed, **shoulder}.items()
                            if key not in all_packs}
            pair_groups = [(all_packs, surroundings), (flex_pack, ext_pack),
                           (equilibrium_wires, {**surroundings, **all_packs}),
                           ({'wire:flexor': equilibrium_wires['wire:flexor']},
                            {'wire:extensor': equilibrium_wires['wire:extensor']})]
            contacts = [dict(part=key, obstacle=target, overlap_mm3=vol)
                        for left, right in pair_groups for key, part in left.items() for target, other in right.items()
                        if (vol := float(A.I.intersection_volume(part, other))) > .001]
            if contacts:
                equilibrium_hits.append(dict(elbow_deg=angle, pairs=contacts))
        for label, results in coupled.items():
            gh = A.I.rotation_about(bare['joint_centre_source_mm'], bare['axes_source'][label], angle/4)
            candidates = {key: part.transform(gh) for key, part in moving_posed.items()}
            candidates.update({key: part.transform(gh) for key, part in {**local, **posed_wires}.items()})
            cross = [dict(moving=key, fixed=target, overlap_mm3=volume)
                     for key, part in candidates.items() for target, obstacle in shoulder.items()
                     if (volume := float(A.I.intersection_volume(part, obstacle))) > .001]
            if cross:
                results.append(dict(gh_deg=angle/4, elbow_deg=angle, hits=cross))
    return dict(step_deg=1., samples=121, moving_keys=sorted(moving), local_keys=list(local),
                owner_reference_sha256=sha(owner_file), hits=rows,
                internal_wires=dict(wire_radius_mm=.3, fixed_keys=sorted(fixed_elbow),
                                    moving_keys=sorted(moving), hits=internal_wire_hits,
                                    scope='Prescribed cable paths between nominal anchors. No TPU or real swivel response.'),
                tracked_swivels=dict(rows=swivel_trials,
                                     simultaneous_equilibrium=dict(samples=equilibrium_samples, step_deg=1.,
                                                                   hits=equilibrium_hits,
                                                                   scope='Shoulder neutral. Rigid parts and taut cables; friction and TPU loads omitted.'),
                                     simultaneous_cross_checks=dict(samples=121, step_deg=1.,
                                                                    hits=simultaneous_hits),
                                     scope='Each swivel tracks projected incoming direction while the other is neutral. Rigid checks stop at first hit per pack. Equilibrium unverified.'),
                coupled=dict(gh_step_deg=.25, elbow_step_deg=1., shoulder_keys=sorted(shoulder),
                             paths=coupled, samples_per_path=121,
                             scope='Rigid forearm and revised proximal fittings, plus prescribed elbow cables versus stationary shoulder solids.'),
                scope='Nominal rigid elbow rotation, distal swivels held neutral, shoulder fixed. No TPU deformation.')


def screen_seat_pin(seat, direction, cap, solids, socket, bearing_axis, folder):
    """Probe a transverse pin supported by the cap outside the seat footprint."""
    mf = A.I._manifold()
    vertices = np.asarray(seat.to_mesh64().vert_properties)[:, :3]
    front = vertices[np.argmax(vertices@direction)]
    transverse = A.D.unit(np.cross(direction, np.eye(3)[np.argmin(abs(direction))]))
    other = np.cross(direction, transverse)
    rows = []
    selected = None
    for gap in (.75, 1.25, 2.):
        point = front+direction*gap
        for angle in range(0, 180, 30):
            theta = np.radians(angle)
            axis = transverse*np.cos(theta)+other*np.sin(theta)
            pin = rod(point-axis*4, point+axis*4, .6)
            bore = rod(point-axis*5, point+axis*5, .68)
            revised = (cap-bore).as_original().simplify(1e-6)
            support = float(A.I.intersection_volume(pin, cap))
            if support < 1. or len(revised.decompose()) != 1:
                continue
            fixed = {key: part for key, part in solids.items() if key != 'elbow:FJ1477_proximal_rigid_cap'}
            fixed['revised_cap'] = revised
            fixed['tool_envelope'] = socket
            for sign in (-1, 1):
                delta = axis*(sign*8.)
                envelope = mf.Manifold.batch_hull([pin, pin.translate(delta)])
                hits = {key: volume for key, part in fixed.items()
                        if (volume := float(A.I.intersection_volume(envelope, part))) > .001}
                stop = next((float(distance) for distance in np.arange(.05, 3.01, .05)
                             if A.I.intersection_volume(seat.translate(direction*distance), pin) > .001), None)
                row = dict(point_mm=point.tolist(), axis=axis.tolist(), sign=sign,
                                 bore_support_proxy_mm3=support, entry_envelope_hits=hits,
                                 seat_exit_stop_mm=stop)
                rows.append(row)
                if hits or stop is None or selected is not None:
                    continue
                bearing = float(A.I.intersection_volume(seat.translate(bearing_axis*.1), revised))/.1
                if bearing < 1.:
                    continue
                insertion = axis*sign
                cap_hits = cap.ray_cast(point+insertion*14., point-insertion*14.)
                crossings = [float((np.array(hit.position)-point)@insertion) for hit in cap_hits]
                axial_fit = dict(cap_crossings_mm=crossings, pin_end_coordinates_mm=[-4., 4.])
                if len(crossings) >= 2:
                    entry_face = crossings[0]
                    head_position = point+insertion*entry_face
                    head = rod(head_position, head_position+insertion, 1.1)
                    stem = rod(point-insertion*4, head_position+insertion*.15, .6)
                    headed = (stem+head).as_original().simplify(1e-6)
                    sweep = mf.Manifold.batch_hull([headed, headed.translate(insertion*8)])
                    axial_fit.update(head_face_mm=entry_face,
                                     head_entry_hits={key: vol for key, part in fixed.items()
                                                      if (vol := float(A.I.intersection_volume(sweep, part))) > .001},
                                     near_cap_stock_mm=entry_face-crossings[1])
                    if not axial_fit['head_entry_hits'] and len(headed.decompose()) == 1:
                        try:
                            _, axial_fit['headed_pin_source'] = export(
                                headed, folder/'extensor_late_headed_pin.obj', point)
                        except (ValueError, AssertionError) as error:
                            axial_fit['export_error'] = str(error)
                    if entry_face < 0 and crossings[1] < -4.2:
                        blind_bore = rod(point-insertion*4.2, point+insertion*5., .68)
                        blind_cap = (cap-blind_bore).as_original().simplify(1e-6)
                        tool = rod(point+insertion*4.01, point+insertion*20., .5)
                        tool_hits = {key: vol for key, part in fixed.items() if key != 'revised_cap'
                                     and (vol := float(A.I.intersection_volume(tool, part))) > .001}
                        tool_hits['blind_cap_mm3'] = float(A.I.intersection_volume(tool, blind_cap))
                        engagement = min(0., entry_face)+4.
                        thread = dict(nominal_engagement_mm=engagement,
                                      axial_bottom_stock_mm=-4.2-crossings[1],
                                      axial_bottom_clearance_mm=.2,
                                      cap_components=len(blind_cap.decompose()), driver_hits=tool_hits,
                                      nominal_pin_cap_overlap_mm3=float(A.I.intersection_volume(pin, blind_cap)))
                        axial_fit['blind_thread_trial'] = thread
                        if (engagement > 0. and thread['cap_components'] == 1
                                and max(tool_hits.values()) <= .001
                                and thread['nominal_pin_cap_overlap_mm3'] <= .001):
                            _, thread['cap_source'] = export(
                                blind_cap, folder/'extensor_late_blind_pin_cap.obj', point)
                try:
                    _, cap_source = export(revised, folder/'extensor_late_pin_cap.obj', point)
                    _, pin_source = export(pin, folder/'extensor_late_retention_pin.obj', point)
                except (ValueError, AssertionError) as error:
                    row['export_error'] = str(error)
                    continue
                selected = dict(**row, bearing_proxy_mm2=bearing, axial_fit=axial_fit,
                                sources=dict(cap=cap_source, pin=pin_source))
    return dict(trials=rows, selected=selected)


def main():
    base = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    bridge_file = base/'elbow_relocated_bridge_screen.json'
    bridge = read(bridge_file)
    assembly = base/'assembly_tpu_readback_candidate.json'
    assert bridge['assembly_sha256'] == sha(assembly)
    assert len(bridge['restored_distal_keys']) == 20
    channels = read(base/'elbow_relocated_channel_cores/channel_core_check.json')
    pair = channels['choice']['pair']
    frame = read(A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json')['frame_source']
    shaft = np.array(frame['basis_columns'])[:, 2]
    pose = A.I.rotation_about(frame['origin_mm'], shaft, pair['extensor_angle_deg'])
    pose[:, 3] += shaft*pair['extensor_shift_mm']
    cap = load_component(read(assembly)['sources']['elbow:FJ1477_proximal_rigid_cap']).transform(pose)
    lattice = next(item['source'] for row in bridge['extensor_lattice_midspan']['trials']
                   for terminal in row.get('terminal_trials', [])
                   for item in terminal.get('entry_relief_trials', []) if 'source' in item)
    body = load_component(lattice)
    tip_row = next(row for row in bridge['terminal_seat_audit']
                   if (row['muscle'], row['side']) == ('extensor', 'proximal'))
    tip, axis = (np.array(tip_row[key]) for key in ('tip_source_mm', 'outward_axis'))
    rows = []
    for gap in (.06, .25, .5):
        raw = rod(tip+axis*gap, tip+axis*7.8, 2.4)-rod(tip-axis, tip+axis*9., 1.)
        pieces = (raw-cap).decompose()
        mouth = rod(tip+axis*gap, tip+axis*(gap+.4), 2.4)
        retained = [part for part in pieces if A.I.intersection_volume(part, mouth) > .001]
        row = dict(gap_mm=gap, pieces=len(pieces), bearing_side_pieces=len(retained))
        rows.append(row)
        if len(retained) != 1:
            continue
        seat = retained[0]
        row['entry'] = screen_cap_after_tpu(seat, body, bridge, axis, {'extensor_cap': cap})
        if not row['entry']['passing']:
            continue
        _, row['source'] = export(seat.as_original().simplify(1e-6),
                                  base/'elbow_relocated_bone_seats/extensor_late_reaction_face.obj', tip)
        break
    mf = A.I._manifold()
    bearing_trials = []
    viable = []
    for back in (3.5, 4., 4.5, 5.):
        seat = rod(tip+axis*.25, tip+axis*back, 2.4)-rod(tip-axis, tip+axis*9., 1.)
        entry = screen_cap_after_tpu(seat, body, bridge, axis)
        for path in entry['trials']:
            if path['first_hit']:
                continue
            delta = np.array(path['offset_mm'])
            if abs(float(A.D.unit(delta)@axis)) > .01:
                continue
            # Solid outer sweep keeps the entire lateral slot clear. Its
            # axial end plane remains available as a bearing shoulder.
            stock = rod(tip+axis*.25, tip+axis*back, 2.4)
            sweep = mf.Manifold.batch_hull([stock, stock.translate(delta)])
            revised = (cap-sweep).as_original().simplify(1e-6)
            removed = cap-revised
            sources = read(assembly)['sources']
            hardware_hits = {}
            for key in ('left', 'right'):
                bolt = load_component(sources[f'elbow:FJ1477_proximal_bolt_{key}']).transform(pose)
                lo, hi = np.array(bolt.bounding_box()).reshape(2, 3)
                guard = mf.Manifold.cube(hi-lo+1).translate(lo-.5)
                overlap = float(A.I.intersection_volume(removed, guard))
                if overlap > .001:
                    hardware_hits[key] = overlap
            bearing = float(A.I.intersection_volume(seat.translate(axis*.1), revised))/.1
            row = dict(depth_mm=back, offset_mm=delta.tolist(),
                                       removed_cap_mm3=float(removed.volume()), bearing_proxy_mm2=bearing,
                                       cap_components=len(revised.decompose()), bolt_guard_hits=hardware_hits)
            bearing_trials.append(row)
            if row['cap_components'] == 1 and not hardware_hits and bearing >= 1.:
                viable.append((row, seat, revised))
    for row, seat, revised in sorted(viable, key=lambda item: item[0]['removed_cap_mm3']):
        folder = base/'elbow_relocated_bone_seats'
        try:
            seat, seat_source = export(seat, folder/'extensor_late_planar_seat.obj', tip)
            revised, cap_source = export(revised, folder/'extensor_late_planar_cap.obj', tip)
        except (AssertionError, ValueError):
            continue
        delta = np.array(row['offset_mm'])
        overlap = max(float(A.I.intersection_volume(seat.translate(delta*t), revised))
                      for t in np.linspace(0., 1., 49))
        if overlap <= .001:
            row['sources'] = dict(seat=seat_source, cap=cap_source)
            row['readback_path_peak_cap_overlap_mm3'] = overlap
            guide_file = base/'extensor_guide_body_audit.json'
            legacy = read(guide_file)['detachable_guide']
            prior = next(candidate for candidate in legacy['assembly_readback']['guide_localization']
                         ['dual_key_trials'] if 'sources' in candidate)
            guide = load_component(prior['sources']['guide'])
            guide_offset = np.array(next(item['offset_mm'] for item in legacy['side_opening_trials']
                                         if 'source' in item))
            def peak(moving, offset, fixed):
                return max(float(A.I.intersection_volume(moving.translate(offset*t), fixed))
                           for t in np.linspace(0., 1., 49))
            row['legacy_guide_compatibility'] = dict(
                report_sha256=sha(guide_file), static_overlap_mm3=float(A.I.intersection_volume(guide, seat)),
                seat_after_guide_peak_mm3=peak(seat, delta, guide),
                guide_after_seat_peak_mm3=peak(guide, guide_offset, seat),
                guide_path_vs_revised_TPU_peak_mm3=peak(guide, guide_offset, body),
                guide_path_vs_slotted_cap_before_keyways_peak_mm3=peak(guide, guide_offset, revised),
                step_mm=.25, scope='Prior guide geometry used as a fit probe. This cap has not been cut for the old keys.')
            # In seat coordinates, the installed guide sweeps opposite to the
            # insertion direction. Limit the cut to this occupied region.
            local = mf.Manifold.batch_hull([seat, seat.translate(delta)])
            guide_stock = guide^local
            clearance = guide_stock.minkowski_sum(rod(np.zeros(3), -delta, .06))
            relieved = (seat-clearance).as_original().simplify(1e-6)
            compat = row['legacy_guide_compatibility']
            compat['seat_relief'] = dict(
                removed_mm3=float((seat-relieved).volume()),
                components=len(relieved.decompose()),
                bearing_proxy_mm2=float(A.I.intersection_volume(relieved.translate(axis*.1), revised))/.1,
                seat_after_guide_peak_mm3=peak(relieved, delta, guide))
            relief_result = compat['seat_relief']
            if (relief_result['components'] == 1 and relief_result['bearing_proxy_mm2'] >= 1.
                    and relief_result['seat_after_guide_peak_mm3'] <= .001):
                _, relief_result['source'] = export(relieved, folder/'extensor_late_seat_guide_relief.obj', tip)
            # The old keyed cap was cut from an integral flat cap. Recover
            # only the groove removal that belongs to the original rigid cap.
            old_key_cap = load_component(prior['sources']['cap'])
            keyed_cap = (revised-(cap-old_key_cap)).as_original().simplify(1e-6)
            witness_stock = mf.Manifold.batch_boolean(
                [keyed_cap^guide.translate(guide_offset*t) for t in np.linspace(0., 1., 97)],
                mf.OpType.Add)
            witness_volume = float(witness_stock.volume())
            if witness_volume > .001:
                keyed_cap = (keyed_cap-witness_stock.minkowski_sum(
                    mf.Manifold.cube([.03]*3, True))).as_original().simplify(1e-6)
            guide_hull = mf.Manifold.batch_hull([guide, guide.translate(guide_offset)])
            opening = guide_hull.minkowski_sum(mf.Manifold.cube([.12]*3, True))
            open_tpu = (body-opening).as_original().simplify(1e-6)
            joint = dict(cap_components=len(keyed_cap.decompose()), tpu_components=len(open_tpu.decompose()),
                         added_cap_cut_mm3=float((revised-keyed_cap).volume()),
                         supplemental_entry_witness_mm3=witness_volume,
                         removed_tpu_mm3=float((body-open_tpu).volume()),
                         guide_vs_cap_peak_mm3=peak(guide, guide_offset, keyed_cap),
                         guide_vs_tpu_peak_mm3=peak(guide, guide_offset, open_tpu),
                         seat_vs_cap_peak_mm3=peak(relieved, delta, keyed_cap),
                         seat_vs_tpu_peak_mm3=peak(relieved, delta, open_tpu),
                         seat_bearing_proxy_mm2=float(A.I.intersection_volume(
                             relieved.translate(axis*.1), keyed_cap))/.1)
            row['combined_key_and_seat'] = joint
            if (joint['cap_components'] == joint['tpu_components'] == 1
                    and max(value for key, value in joint.items() if key.endswith('peak_mm3')) <= .001
                    and joint['seat_bearing_proxy_mm2'] >= 1.):
                joint['sources'] = {
                    'cap': export(keyed_cap, folder/'extensor_late_keyed_cap.obj', tip)[1],
                    'tpu': export(open_tpu, folder/'extensor_late_keyed_lattice.obj', tip)[1],
                    'guide': prior['sources']['guide'], 'seat': relief_result['source']}
                sources = joint['sources']
                solids = load_relocated_assembly(base, bridge, channels, sources['tpu'], sources['cap'])
                revised_tpu = solids.pop('tpu:extensor_relocated')
                delayed = ('elbow:FJ1477_proximal_nut_washer_right', 'elbow:FJ1477_proximal_nut_right')
                washer, nut = (solids.pop(key) for key in delayed)
                washer_centre = np.array(washer.bounding_box()).reshape(2, 3).mean(axis=0)
                nut_centre = np.array(nut.bounding_box()).reshape(2, 3).mean(axis=0)
                fastener_offset = A.D.unit(nut_centre-washer_centre)*12.
                stop = next(item['source'] for item in legacy['assembly_readback']['washer_tongue_trials']
                            if 'source' in item)
                joint['assembly_sequence'] = verify_sequence(solids, [
                    ('tpu:extensor_relocated', revised_tpu, bridge['paired_lattice_entry']['extensor_second']['offset_mm']),
                    ('extensor:dual_key_guide', guide, guide_offset),
                    ('extensor:late_reaction_seat', load_component(sources['seat']), delta),
                    ('extensor:stop_washer', load_component(stop), fastener_offset),
                    (delayed[1], nut, fastener_offset),
                ])
                fixed = {key: solid for key, solid in solids.items()
                         if key != 'extensor:late_reaction_seat'}
                seat_actual = solids['extensor:late_reaction_seat']
                seat_travel = A.D.unit(delta)
                def first_exit_contact():
                    for distance in np.linspace(0., 12., 241):
                        moved = seat_actual.translate(seat_travel*distance)
                        hits = {key: volume for key, obstacle in fixed.items()
                                if (volume := float(A.I.intersection_volume(moved, obstacle))) > .001}
                        if hits:
                            return dict(distance_mm=float(distance), hits=hits)
                    return None
                normal = A.D.unit(fastener_offset)
                socket = rod(nut_centre+normal*1.5, nut_centre+normal*22., 4.2)
                socket_hits = {key: vol for key, solid in solids.items()
                               if key not in (delayed[1], 'elbow:FJ1477_proximal_bolt_right')
                               and (vol := float(A.I.intersection_volume(socket, solid))) > .001}
                joint['retention_and_tool'] = dict(
                    seat_exit_step_mm=.05, seat_exit_first_contact=first_exit_contact(),
                    socket_proxy_hits=socket_hits,
                    scope='Straight reverse path and nominal socket envelope. Other escape paths and torque unverified.')
                joint['retention_and_tool']['transverse_pin_trials'] = screen_seat_pin(
                    seat_actual, seat_travel, solids['elbow:FJ1477_proximal_rigid_cap'], solids, socket, axis, folder)
                pin = joint['retention_and_tool']['transverse_pin_trials']['selected']
                if pin is not None:
                    blind = pin['axial_fit'].get('blind_thread_trial', {})
                    cap_source = blind.get('cap_source', pin['sources']['cap'])
                    motion_solids = {**solids,
                                     'elbow:FJ1477_proximal_rigid_cap': load_component(cap_source),
                                     'extensor:retention_pin': load_component(pin['sources']['pin'])}
                    joint['revised_elbow_motion'] = check_revised_elbow_motion(
                        motion_solids, channels, np.array(read(guide_file)['guide_source_mm']), read(assembly)['sources'])
            print('selected', row['depth_mm'], row['bearing_proxy_mm2'], row['removed_cap_mm3'])
            break
    save(base/'extensor_late_reaction_audit.json', dict(
        bridge_sha256=sha(bridge_file), trials=rows, planar_bearing_trials=bearing_trials,
        scope='Local fit and rigid insertion with corrected TPU, cap and two elbow adapters. Fixation and remaining parts excluded.'))
    print([(row['gap_mm'], row.get('entry', {}).get('passing')) for row in rows])


if __name__ == '__main__':
    main()
