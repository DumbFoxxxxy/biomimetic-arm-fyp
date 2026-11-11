"""Find nominal axial package space for the elbow cam and output drums."""
import numpy as np
from numpy.polynomial import Polynomial

from build_s1_teres_origin_mount import A, ROOT, load_component, mesh_geom, read, rod, save, sha


def check_envelope(envelope, poses, cap=.5, mating_face=None):
    lo, hi = np.array(envelope.bounding_box()).reshape(2, 3)
    nearest = dict(gap_mm=cap)
    for angle, keys, bounds, solids in poses:
        separation = np.maximum(lo-bounds[:, 3:], bounds[:, :3]-hi)
        for index in np.flatnonzero(np.all(separation < cap, axis=1)):
            part = solids[index]
            if np.all(separation[index] < 0):
                volume = A.I.intersection_volume(envelope, part)
                if volume > .001:
                    return dict(first_overlap=dict(elbow_deg=angle, part=keys[index], overlap_mm3=volume),
                                closest=None)
            if keys[index] == mating_face:
                continue
            gap = A.I.gap(envelope, part, cap)
            if gap < nearest['gap_mm']:
                nearest = dict(gap_mm=gap, elbow_deg=angle, part=keys[index])
    return dict(first_overlap=None, closest=nearest)


def scan_envelope(centre, axis, outer_radius, width, offsets, poses):
    return [dict(axial_offset_mm=float(offset), **check_envelope(
        rod(centre+axis*(offset-width/2), centre+axis*(offset+width/2), outer_radius), poses))
            for offset in offsets]


def pack_pose(angle, solids):
    keys = list(solids)
    parts = list(solids.values())
    return angle, keys, np.array([part.bounding_box() for part in parts]), parts


def screen_parallel_axis(centre, axis, radius, width, poses):
    first = A.D.unit(np.cross(axis, [0., 0., 1.]))
    second = np.cross(axis, first)
    tested = []
    for radial_offset in (12., 16., 20., 24., 28., 32., 36., 40.):
        for azimuth in range(0, 360, 30):
            angle = np.radians(azimuth)
            radial = radial_offset*(np.cos(angle)*first+np.sin(angle)*second)
            trials = scan_envelope(centre+radial, axis, radius, width,
                                   [0., 8., -8., 16., -16., 24., -24.], poses)
            tested.extend(dict(radial_offset_mm=radial_offset, azimuth_deg=azimuth,
                               centre_source_mm=(centre+radial+axis*row['axial_offset_mm']).tolist(),
                               **row) for row in trials)
        viable = [row for row in tested if row['closest'] and row['closest']['gap_mm'] >= .2]
        if viable:
            break
    else:
        viable = []
    viable.sort(key=lambda row: (abs(row['axial_offset_mm']), -row['closest']['gap_mm'], row['azimuth_deg']))
    return dict(trials=tested, shortlist=viable,
                scope='Parallel-axis envelope only. Existing hinge-to-cam synchronization '
                      'and mounts remain to be designed.')


def screen_sync_belt(centre, axis, cam_width, shortlist, pin, poses):
    mf = A.I._manifold()
    vertices, _ = mesh_geom.manifold_mesh_to_arrays(pin.to_mesh64())
    pin_end = float(((vertices-centre)@axis).min())
    pin_face = centre+axis*pin_end
    rows = []
    for item in shortlist:
        cam = np.asarray(item['centre_source_mm'])
        start_cam = cam-axis*cam_width/2
        for station in range(-32, -57, -2):
            drive = centre+axis*station
            driven = cam+axis*(station-item['axial_offset_mm'])
            hubs = [rod(point-axis*2.5, point+axis*2.5, 5.5) for point in (drive, driven)]
            belt = mf.Manifold.batch_hull(hubs)
            # Envelope includes pulleys and both straight spans, including
            # otherwise empty space between spans for a conservative first pass.
            checks = {'belt_and_pulleys': check_envelope(belt, poses)}
            if checks['belt_and_pulleys']['first_overlap'] is None:
                checks['hinge_extension'] = check_envelope(
                    rod(pin_face, drive, 2.5), poses, mating_face='elbow:elbow_shoulder_pin')
                checks['cam_extension'] = check_envelope(rod(start_cam, driven, 2.5), poses)
            rows.append(dict(cam_azimuth_deg=item['azimuth_deg'], belt_axial_station_mm=station,
                             tested_components=checks))
            if len(checks) == 3 and all(check['closest'] is not None and check['closest']['gap_mm'] >= .2
                                        for check in checks.values()):
                rows[-1]['passed_with_min_gap_mm'] = min(check['closest']['gap_mm'] for check in checks.values())
                break
    return dict(hinge_pin_end_station_mm=pin_end, trials=rows,
                nominal_mating_face='elbow:elbow_shoulder_pin to hinge extension',
                scope='1:1 synchronization packaging, with nominal 11 mm pulley outer diameter, '
                      '5 mm belt plane width and 5 mm shaft envelopes. The shaft pickup, hub bores, '
                      'bearings and bracket are unresolved.')


def screen_bias_drum(cam, axis, width, shortlist, belt, poses):
    """Package a constant-radius spring drum and bound the cable-only bias load."""
    pitch_radius, outer_radius, drum_width = 4.5, 6., 4.
    angles = np.radians([row['elbow_deg'] for row in cam['poses']])
    reflected = {}
    for label, profile in cam['cam_pitch_profiles'].items():
        ratio = -Polynomial(profile['fit_coefficients_rad']).deriv()(angles)
        torque = np.array([row['torque_per_N_Nmm'][label] for row in cam['poses']])
        reflected[label] = ratio*torque
    # For independent nonnegative cable tensions, negative bias only needs
    # to exceed each cable's positive contribution at its upper tension bound.
    adverse = {label: float(max(0., values.max())) for label, values in reflected.items()}
    rows = []
    for candidate in shortlist:
        match = next((row for row in belt['trials'] if row['cam_azimuth_deg'] == candidate['azimuth_deg']
                      and 'passed_with_min_gap_mm' in row), None)
        if match is None:
            continue
        cam_station = candidate['axial_offset_mm']
        belt_station = match['belt_axial_station_mm']
        cam_outer_face = cam_station-width/2
        belt_inner_face = belt_station+2.5
        start = int(np.ceil(belt_inner_face+drum_width/2+.5))
        stop = int(np.floor(cam_outer_face-drum_width/2-.5))
        cam_centre = np.asarray(candidate['centre_source_mm'])
        trials = scan_envelope(cam_centre, axis, outer_radius, drum_width,
                               [station-cam_station for station in range(start, stop+1)], poses)
        rows.extend(dict(cam_azimuth_deg=candidate['azimuth_deg'],
                         drum_station_mm=row['axial_offset_mm']+cam_station, **row) for row in trials)
    feasible = [row for row in rows if row['closest'] and row['closest']['gap_mm'] >= .2]
    return dict(pitch_radius_mm=pitch_radius, envelope_radius_mm=outer_radius, width_mm=drum_width,
                spring_cable_travel_mm=float(pitch_radius*(angles[-1]-angles[0])),
                reflected_torque_per_N_Nmm={label: [float(values.min()), float(values.max())]
                                            for label, values in reflected.items()},
                negative_bias_adverse_torque_per_N_Nmm=adverse,
                spring_force_per_N_of_max_muscle_tension={label: value/pitch_radius
                                                         for label, value in adverse.items()},
                trials=rows, feasible=feasible,
                scope='Ideal 1:1 belt and negative flank bias. Cable-only reflected torque lower bound. '
                      'TPU reaction, friction, spring housing, actual spring and supports remain unresolved.')


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    cam_file = folder/'anatomical_elbow_cam_pitch_geometry.json'
    motion_file = folder/'anatomical_elbow_rigid_motion.json'
    assembly_file = folder/'assembly_anatomical_elbow_angled_candidate.json'
    combined_file = folder/'anatomical_elbow_combined_bypass.json'
    cam, motion, assembly, combined = map(read, (cam_file, motion_file, assembly_file, combined_file))
    assert motion['assembly_sha256'] == sha(assembly_file)
    assert cam['rigid_motion_sha256'] == sha(motion_file)
    assert cam['combined_bypass_sha256'] == sha(combined_file)
    sources = {**assembly['sources'], 'tpu:deltoid': combined['deltoid_source']}
    solids = {key: load_component(source) for key, source in sources.items()
              if key not in motion['excluded_unposed_elbow_tpu']}
    moving = set(motion['moving_keys'])
    fixed = {key: solid for key, solid in solids.items() if key not in moving}
    rotors = {}
    for label, head in (('flexor', 'FJ1486'), ('extensor', 'FJ1477')):
        key = f'elbow:{head}_cap_rotor'
        screen = read((ROOT/sources[key]['file']).parent/'screen.json')
        rotors[label] = dict(
            centre=np.asarray(screen['pivot_zero_display_mm'])+sources[key]['translation_source_mm'],
            axis=np.asarray(screen['trunnion_axis_zero_display']),
            keys={name for name in moving if name.startswith(f'elbow:{head}_distal_')
                  or name in {key, f'elbow:{head}_active_wire_crimp_barrel'}})
    centre, axis = np.array(motion['centre_source_mm']), np.array(motion['axis_source'])
    cam_poses = []
    drum_poses = []
    for row in cam['poses']:
        angle = row['elbow_deg']
        bend = A.I.rotation_about(centre, axis, angle)
        posed = {**fixed, **{key: solids[key].transform(bend) for key in moving}}
        for label, rotor in rotors.items():
            turn = A.I.rotation_about(rotor['centre'], rotor['axis'], row['muscles'][label]['swivel_deg'])
            posed.update({key: solids[key].transform(turn).transform(bend) for key in rotor['keys']})
        cam_poses.append(pack_pose(angle, posed))
        # Express the posed assembly in the moving forearm frame so each
        # axisymmetric output envelope can stay fixed for the sweep.
        reverse = A.I.rotation_about(centre, axis, -angle)
        local = {key: part.transform(reverse) for key, part in posed.items()}
        drum_poses.append(pack_pose(angle, local))
    profiles = cam['cam_pitch_profiles']
    cam_radius = max(row['pitch_curve_max_radial_extent_mm'] for row in profiles.values())+1.25
    cam_width = 8.
    cam_offsets = [sign*offset for offset in range(8, 33) for sign in (1, -1)]
    drum_offsets = [sign*offset for offset in range(4, 25) for sign in (1, -1)]
    cam_trials = scan_envelope(centre, axis, cam_radius, cam_width, cam_offsets, cam_poses)
    parallel = screen_parallel_axis(centre, axis, cam_radius, cam_width, cam_poses)
    belt = screen_sync_belt(centre, axis, cam_width, parallel['shortlist'],
                            solids['elbow:elbow_shoulder_pin'], cam_poses)
    bias_drum = screen_bias_drum(cam, axis, cam_width, parallel['shortlist'], belt, cam_poses)
    drums = {}
    for label, rotor in rotors.items():
        radius = profiles[label]['output_drum_pitch_radius_mm']+1.25
        width = 4.
        drums[label] = dict(axis_source=rotor['axis'].tolist(), pivot_source_mm=rotor['centre'].tolist(),
                            outer_radius_mm=radius, width_mm=width,
                            trials=scan_envelope(rotor['centre'], rotor['axis'], radius, width,
                                                 drum_offsets, drum_poses))
    package_file = folder/'anatomical_elbow_cam_package_screen.json'
    save(package_file, dict(
        cam_pitch_report_sha256=sha(cam_file), sample_step_deg=1, screened_instances=len(solids),
        cam_stack=dict(axis_source=axis.tolist(), pivot_source_mm=centre.tolist(),
                       outer_radius_mm=cam_radius, width_mm=cam_width, trials=cam_trials),
        output_drums=drums, excluded_unposed_elbow_tpu=motion['excluded_unposed_elbow_tpu'],
        parallel_axis_cam=parallel,
        synchronization_belt=belt,
        scope='Packaging envelopes against existing posed parts. Support, cable, '
              'flange details, and flexible elbow muscle envelopes remain unmodelled.'))
    bias_drum['package_sha256'] = sha(package_file)
    save(folder/'anatomical_elbow_bias_drum_screen.json', bias_drum)
    for label, rows in [('cam', cam_trials), *[(key, row['trials']) for key, row in drums.items()]]:
        for sign in (1, -1):
            print(label, sign, next((row for row in rows if row['axial_offset_mm']*sign > 0
                                     and row['closest'] and row['closest']['gap_mm'] >= .2), None))
    print('parallel-axis tested', len(parallel['trials']), 'shortlist', parallel['shortlist'][:4])
    print('belt screen', len(belt['trials']),
          [row for row in belt['trials'] if 'passed_with_min_gap_mm' in row])
    print('bias drum', len(bias_drum['feasible']), bias_drum['negative_bias_adverse_torque_per_N_Nmm'])


if __name__ == '__main__':
    main()
