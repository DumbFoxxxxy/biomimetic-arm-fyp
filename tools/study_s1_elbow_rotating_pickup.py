"""Screen a keyed rotating elbow pin as a synchronous angle pickup."""
import argparse

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--d-flat', action='store_true')
    args = parser.parse_args()
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    package_file = folder/'anatomical_elbow_cam_package_screen.json'
    cam_file = folder/'anatomical_elbow_cam_pitch_geometry.json'
    assembly_file = folder/'assembly_anatomical_elbow_angled_candidate.json'
    motion_file = folder/'anatomical_elbow_rigid_motion.json'
    motion = read(motion_file)
    assert motion['assembly_sha256'] == sha(assembly_file)
    package = read(package_file)
    assert package['cam_pitch_report_sha256'] == sha(cam_file)
    cam = read(cam_file)
    assert cam['rigid_motion_sha256'] == sha(motion_file)
    belt = next(row for row in package['synchronization_belt']['trials'] if 'passed_with_min_gap_mm' in row)
    centre, axis = np.asarray(motion['centre_source_mm']), np.asarray(motion['axis_source'])
    sources = read(assembly_file)['sources']
    keys = {'pin': 'elbow:elbow_shoulder_pin', 'bush': 'elbow:elbow_bearing_bush', 'hub': 'bone:forearm'}
    originals = {name: load_component(sources[key]) for name, key in keys.items()}
    # The rotating extension is joined to the existing head. The new cross
    # key makes the shank, bushing and forearm move together inside the fork.
    end_station = package['synchronization_belt']['hinge_pin_end_station_mm']
    extension = rod(centre+axis*(end_station+.2), centre+axis*(belt['belt_axial_station_mm']-2.5), 2.5)
    extended_pin = (originals['pin']+extension).as_original().simplify(1e-6)
    fixed = {key: load_component(source) for key, source in sources.items()
             if key not in set(motion['moving_keys']) | set(motion['excluded_unposed_elbow_tpu'])
             and key != keys['pin']}
    neutral_obstacles = {**fixed, **{key: load_component(sources[key]) for key in motion['moving_keys']
                                    if key not in keys.values()}}
    moving = {key: neutral_obstacles[key] for key in motion['moving_keys'] if key not in keys.values()}
    swivels = {}
    for label, head in (('flexor', 'FJ1486'), ('extensor', 'FJ1477')):
        key = f'elbow:{head}_cap_rotor'
        source = sources[key]
        screen_file = (ROOT/source['file']).parent/'screen.json'
        assert sha(screen_file) == cam['swivel_screen_sha256'][label]
        screen = read(screen_file)
        pack = {name for name in moving if name.startswith(f'elbow:{head}_distal_')
                or name in {key, f'elbow:{head}_active_wire_crimp_barrel'}}
        swivels[label] = (np.array(screen['pivot_zero_display_mm'])+source['translation_source_mm'],
                          screen['trunnion_axis_zero_display'], pack)
    u = A.D.unit(np.cross(axis, [0., 0., 1.]))
    v = np.cross(axis, u)
    mf = A.I._manifold()
    trials = []
    winner = None
    for station in (0., -3., 3.):
        location = centre+axis*station
        for azimuth in range(0, 360, 30):
            theta = np.radians(azimuth)
            direction = np.cos(theta)*u+np.sin(theta)*v
            bore_start, dowel_start = (3.20, 3.23) if args.d_flat else (-1.8, -1.65)
            bore = rod(location+direction*bore_start, location+direction*12., .95)
            dowel = rod(location+direction*dowel_start, location+direction*6.7, .9)
            engage = {name: A.I.intersection_volume(dowel, solid) for name, solid in originals.items()}
            row = dict(axial_station_mm=station, azimuth_deg=azimuth, direction_source=direction.tolist(),
                       engagement_stock_mm3=engage)
            trials.append(row)
            if min(engage['bush'], engage['hub']) < 1. or (not args.d_flat and engage['pin'] < 1.):
                continue
            parts = {name: ((extended_pin if name == 'pin' else solid)-bore).as_original().simplify(1e-6)
                     for name, solid in originals.items()}
            if args.d_flat:
                normal = -direction
                frame = np.stack([normal, np.cross(axis, normal), axis], axis=1)

                def flat_cutter(distance):
                    return mf.Manifold.cube([8., 16., 70.]).translate(
                        [distance, -8., -8.2]).transform(np.c_[frame, centre])

                parts['pin'] = (extended_pin-flat_cutter(2.55)).as_original().simplify(1e-6)
                liner_fill = rod(centre-axis*7.8, centre+axis*7.8, 3.25)
                liner_void = rod(centre-axis*9., centre+axis*9., 3.15)-flat_cutter(2.60)
                parts['bush'] = ((originals['bush']+liner_fill)-liner_void-bore).as_original().simplify(1e-6)
                if A.I.intersection_volume(parts['pin'], parts['bush']) > .001:
                    continue
            row['components'] = {name: sum(part.volume() > 1e-8 for part in solid.decompose())
                                 for name, solid in parts.items()}
            if set(row['components'].values()) != {1}:
                continue
            row['hub_backing_at_dowel_tip'] = A.I.intersection_volume(
                rod(location+direction*5.9, location+direction*6.7, .9), originals['hub'])
            if row['hub_backing_at_dowel_tip'] < 1.:
                continue
            insertion = rod(location+direction*dowel_start, location+direction*26., .9)
            row['neutral_insertion_hits_mm3'] = {
                key: volume for key, solid in {**neutral_obstacles, **parts}.items()
                if (volume := A.I.intersection_volume(insertion, solid)) > .001}
            if row['neutral_insertion_hits_mm3']:
                continue
            winner = parts, dowel, row
            break
        if winner:
            break
    result = dict(package_sha256=sha(package_file), original_assembly_sha256=sha(assembly_file), trials=trials,
                  scope='Nominal keyed-angle pickup. Cross-hole and key affect the loaded hinge. '
                        'Strength, wear, bearing material, and installation sequence remain unverified.')
    if winner:
        parts, dowel, row = winner
        readback = {}
        exported = {}
        prefix = 'anatomical_elbow_d_flat_pickup' if args.d_flat else 'anatomical_elbow_pickup'
        for name, solid in {**parts, 'dowel': dowel}.items():
            readback[name], exported[name] = export(solid, folder/f'{prefix}_{name}.obj', centre)
        hits = []
        near_pairs = {}
        checked = 0
        for pose in cam['poses']:
            angle = pose['elbow_deg']
            turn = A.I.rotation_about(centre, axis, angle)
            posed_moving = {key: part.transform(turn) for key, part in moving.items()}
            for label, (pivot, hinge, pack) in swivels.items():
                swivel = A.I.rotation_about(pivot, hinge, pose['muscles'][label]['swivel_deg'])
                posed_moving.update({key: moving[key].transform(swivel).transform(turn) for key in pack})
            obstacles = {**fixed, **posed_moving}
            obstacle_keys = list(obstacles)
            bounds = np.array([solid.bounding_box() for solid in obstacles.values()])
            for name, part in readback.items():
                posed = part.transform(turn)
                lo, hi = np.array(posed.bounding_box()).reshape(2, 3)
                separation = np.maximum(lo-bounds[:, 3:], bounds[:, :3]-hi)
                for index in np.flatnonzero(np.all(separation < .5, axis=1)):
                    key = obstacle_keys[index]
                    other = obstacles[key]
                    checked += 1
                    volume = A.I.intersection_volume(posed, other)
                    if volume > .001:
                        hits.append(dict(elbow_deg=angle, pickup_part=name, obstacle=key,
                                         obstacle_moves=key in moving, overlap_mm3=volume))
                    gap = A.I.gap(posed, other, .5)
                    pair = name, key
                    if gap < near_pairs.get(pair, (0.5, 0))[0]:
                        near_pairs[pair] = (gap, angle)
        internal = {name: A.I.intersection_volume(readback['dowel'], solid)
                    for name, solid in readback.items() if name != 'dowel'}
        insertion_direction = np.array(row['direction_source'])
        start = centre+axis*row['axial_station_mm']
        insertion = rod(start+insertion_direction*dowel_start, start+insertion_direction*26., .9)
        insertion_hits = {key: volume for key, solid in {**neutral_obstacles, **readback}.items()
                          if key != 'dowel' and (volume := A.I.intersection_volume(insertion, solid)) > .001}
        internal.update({f'{left}:{right}': A.I.intersection_volume(readback[left], readback[right])
                         for left, right in (('pin', 'bush'), ('pin', 'hub'), ('bush', 'hub'))})
        slab_width = .05
        slab = rod(start-axis*slab_width/2, start+axis*slab_width/2, 3.1)
        original_area = A.I.intersection_volume(originals['pin'], slab)/slab_width
        bored_area = A.I.intersection_volume(readback['pin'], slab)/slab_width
        result.update(selected_trial=row, sources=exported,
                      fixed_sweep_hits=[hit for hit in hits if not hit['obstacle_moves']],
                      moving_sweep_hits=[hit for hit in hits if hit['obstacle_moves']],
                      cam_pitch_report_sha256=sha(cam_file), nearby_pair_checks=checked,
                      close_pairs=[dict(pickup_part=name, obstacle=key, gap_mm=gap, elbow_deg=angle)
                                   for (name, key), (gap, angle) in sorted(near_pairs.items(), key=lambda item: item[1])],
                      nominal_pin_section=dict(slice_width_mm=slab_width, original_area_mm2=original_area,
                                               bored_area_mm2=bored_area,
                                               area_reduction_percent=100.*(1.-bored_area/original_area)),
                      dowel_to_bored_parts_mm3=internal, sample_step_deg=1,
                      neutral_dowel_insertion_hits_mm3=insertion_hits,
                      geometric_screen_passed=not hits and not insertion_hits and max(internal.values()) <= .001)
        if args.d_flat:
            def overlap_at(angle):
                return A.I.intersection_volume(
                    readback['pin'].transform(A.I.rotation_about(centre, axis, angle)), readback['bush'])

            angular_slack = {}
            for sign in (1, -1):
                lo, hi = 0., 30.
                assert overlap_at(sign*hi) > .001
                for _ in range(12):
                    mid = (lo+hi)/2
                    if overlap_at(sign*mid) <= .001:
                        lo = mid
                    else:
                        hi = mid
                angular_slack[str(sign)] = hi
            result['d_flat'] = dict(flat_distance_from_axis_mm=2.55,
                                    nominal_flat_clearance_mm=.05, first_interference_deg=angular_slack,
                                    short_dowel_min_axis_distance_mm=dowel_start)
            insertion_obstacles = {**neutral_obstacles, 'bush': readback['bush'], 'hub': readback['hub']}
            axial_hits = []
            for shift in np.arange(-110., .01, .5):
                displaced = readback['pin'].translate(axis*shift)
                for key, obstacle in insertion_obstacles.items():
                    volume = A.I.intersection_volume(displaced, obstacle)
                    if volume > .001:
                        axial_hits.append(dict(shift_mm=float(shift), obstacle=key, overlap_mm3=volume))
            result['d_flat'].update(pin_axial_insertion_step_mm=.5, pin_axial_insertion_hits=axial_hits)
    output = 'anatomical_elbow_d_flat_pickup.json' if args.d_flat else 'anatomical_elbow_rotating_pickup.json'
    save(folder/output, result)
    print('selected', result.get('selected_trial'), 'passes', result.get('geometric_screen_passed'),
          'first hits', result.get('fixed_sweep_hits', [])[:3])


if __name__ == '__main__':
    main()
