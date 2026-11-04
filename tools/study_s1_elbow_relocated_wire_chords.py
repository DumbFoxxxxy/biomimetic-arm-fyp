"""Reject elbow interface placements with obstructed or reversed tension chords."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, rod, save, sha


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_tpu_readback_candidate.json'
    screen_file = folder / 'elbow_proximal_rotated_station_screen.json'
    repair_file = folder / 'elbow_tpu_readback_repair.json'
    model, screen, repair = map(read, (manifest, screen_file, repair_file))
    assert screen['assembly_sha256'] == sha(manifest)
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    assert screen['frame_sha256'] == sha(frame_file)
    frame = read(frame_file)['frame_source']
    axis = np.array(frame['basis_columns'])[:, 2]
    centre, elbow_axis = map(np.array, A.frames.elbow_axis())
    solids = {key: load_component(row) for key, row in model['sources'].items()}
    mobile = {key for group in screen['moving_groups'].values() for key in group}
    fixed = set(solids) - mobile - {'tpu:elbow_flexor', 'tpu:elbow_extensor'}
    rows = {}
    probes = {}
    bypasses = {}
    for head, label in [('FJ1486', 'flexor'), ('FJ1477', 'extensor')]:
        source = next(row['original_source'] for row in repair['rows'] if row['key'] == f'tpu:elbow_{label}')
        design_file = (ROOT/source['file']).parent/'design.json'
        points = np.array(read(design_file)['pad_centres_display_mm']) + np.array(source['translation_source_mm'])
        original, distal = points
        original_moment = float(np.cross(distal-centre, A.D.unit(original-distal)) @ elbow_axis)
        rows[label], probes[label] = {}, {}
        for pair in screen['collision_free_hardware_pairs']:
            angle, shift = pair[f'{label}_angle_deg'], pair[f'{label}_shift_mm']
            key = (angle, shift)
            if key in rows[label]:
                continue
            transform = A.I.rotation_about(frame['origin_mm'], axis, angle)
            proximal = transform[:, :3] @ original + transform[:, 3] + axis*shift
            direction = A.D.unit(distal-proximal)
            length = float(np.linalg.norm(distal-proximal))
            if length <= 16:
                raise ValueError('Relocated pads leave no interior chord')
            probe = rod(proximal+direction*8, distal-direction*8, 2.)
            moment = float(np.cross(distal-centre, -direction) @ elbow_axis)
            hits = []
            for target in sorted(fixed):
                volume = float(A.I.intersection_volume(probe, solids[target]))
                if volume > .001:
                    hits.append(dict(target=target, overlap_mm3=volume))
            rows[label][key] = dict(angle_deg=angle, shift_mm=shift, proximal_source_mm=proximal.tolist(),
                                    distal_source_mm=distal.tolist(), chord_mm=length,
                                    original_moment_arm_mm=original_moment, chord_moment_arm_mm=moment,
                                    moment_sign_preserved=original_moment*moment > 0, hits=hits)
            probes[label][key] = probe
    obstacle_order = [key for key in ('bone:humerus', 'bone:forearm', 'elbow:elbow_spacer_left',
                                      'tpu:deltoid') if key in fixed]
    obstacle_order += sorted(fixed - set(obstacle_order))
    transverse = np.array(frame['basis_columns'])[:, :2]
    for key, row in rows['flexor'].items():
        proximal = np.array(row['proximal_source_mm'])
        distal = np.array(row['distal_source_mm'])
        options = []
        bypasses[key] = options
        for radius in [8., 16., 24., 32., 40.]:
            for fraction in [.35, .55, .75]:
                for degrees in range(0, 360, 30):
                    angle = np.deg2rad(degrees)
                    waypoint = proximal+(distal-proximal)*fraction + transverse @ (
                        radius*np.array([np.cos(angle), np.sin(angle)]))
                    toward_distal = A.D.unit(distal-waypoint)
                    moment = float(np.cross(distal-centre, -toward_distal) @ elbow_axis)
                    if moment * row['original_moment_arm_mm'] <= 0:
                        continue
                    start = proximal+A.D.unit(waypoint-proximal)*8
                    finish = distal-toward_distal*8
                    probe = A.D.outer_envelope(np.stack([start, waypoint, finish]), 2.)
                    if any(A.I.intersection_volume(probe, solids[target]) > .001 for target in obstacle_order):
                        continue
                    options.append((dict(waypoint_source_mm=waypoint.tolist(), detour_radius_mm=radius,
                                         moment_arm_mm=moment,
                                         length_mm=float(np.linalg.norm(proximal-waypoint)
                                                         + np.linalg.norm(distal-waypoint))), probe))
                if options:
                    break
            if options:
                break
    pairs = []
    bypass_pairs = []
    posed_cache = {}
    for pair in screen['collision_free_hardware_pairs']:
        a = (pair['flexor_angle_deg'], pair['flexor_shift_mm'])
        b = (pair['extensor_angle_deg'], pair['extensor_shift_mm'])
        if not rows['extensor'][b]['hits'] and rows['extensor'][b]['moment_sign_preserved']:
            for detour, probe in bypasses[a]:
                if A.I.intersection_volume(probe, probes['extensor'][b]) <= .001:
                    hardware_hits = []
                    for head, placement in [('FJ1486', a), ('FJ1477', b)]:
                        cache_key = (head, placement)
                        if cache_key not in posed_cache:
                            angle, shift = placement
                            transform = A.I.rotation_about(frame['origin_mm'], axis, angle)
                            transform[:, 3] += axis*shift
                            posed_cache[cache_key] = {
                                key: solids[key].transform(transform) for key in screen['moving_groups'][head]
                            }
                        for key, body in posed_cache[cache_key].items():
                            for label, path in [('flexor', probe), ('extensor', probes['extensor'][b])]:
                                volume = float(A.I.intersection_volume(path, body))
                                if volume > .001:
                                    hardware_hits.append(dict(route=label, target=key, overlap_mm3=volume))
                    bypass_pairs.append(dict(**pair, flexor_waypoint=detour, relocated_hardware_hits=hardware_hits))
        if any(row['hits'] or not row['moment_sign_preserved'] for row in (rows['flexor'][a], rows['extensor'][b])):
            continue
        overlap = float(A.I.intersection_volume(probes['flexor'][a], probes['extensor'][b]))
        if overlap <= .001:
            pairs.append(dict(**pair, chord_pair_overlap_mm3=overlap))
    save(folder/'elbow_relocated_wire_chord_screen.json', dict(
        assembly_sha256=sha(manifest), station_screen_sha256=sha(screen_file),
        script_sha256=sha(Path(__file__)), radius_mm=2, excluded_terminal_lengths_mm=8,
        rows={label: list(values.values()) for label, values in rows.items()}, pairs=pairs,
        waypoint_pairs=bypass_pairs,
        scope='Neutral straight and single-waypoint screening. Actual guides, terminal lengths, lattice and motion pending.'))
    clear_bypass = [row for row in bypass_pairs if not row['relocated_hardware_hits']]
    print('unobstructed sign-preserving chords',
          {label: sum(not row['hits'] and row['moment_sign_preserved'] for row in values.values())
           for label, values in rows.items()}, 'pairs', len(pairs), 'waypoint pairs', len(bypass_pairs),
          'clear of relocated hardware', len(clear_bypass), 'first', clear_bypass[:1], flush=True)


if __name__ == '__main__':
    main()
