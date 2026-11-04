"""Screen axial relocation of the two proximal elbow assemblies."""
import argparse
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def proximal_groups(groups):
    """Keep distal swivel and forearm hardware in its source frame."""
    manifest = (ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces'
                /'mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
                /'continuous_channel_candidate/assembly.json')
    owners = {f"elbow:{row['key']}": row['owner'] for row in read(manifest)['instances']}
    return {head: [key for key in keys if owners[key] != 'elbow'
                   and not owners[key].endswith('_distal_pivot')] for head, keys in groups.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rotate', action='store_true')
    rotate = parser.parse_args().rotate
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_tpu_readback_candidate.json'
    model = read(manifest)
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    frame = read(frame_file)['frame_source']
    axis = np.array(frame['basis_columns'])[:, 2]
    solids = {key: load_component(row) for key, row in model['sources'].items()}
    groups = {}
    for head in ('FJ1486', 'FJ1477'):
        prefix = f'elbow:{head}_'
        groups[head] = sorted(key for key in solids if key.startswith(prefix) and key[len(prefix):].startswith(
            ('proximal_', 'bone_', 'cap_rotor', 'pin_', 'washer_', 'keeper')))
    groups = proximal_groups(groups)
    moving = set(groups['FJ1486'] + groups['FJ1477'])
    pending = {'bone:humerus', 'tpu:elbow_flexor', 'tpu:elbow_extensor'}
    fixed = sorted(set(solids) - moving - pending)
    bone = solids['bone:humerus']
    rows = {}
    passing = {}
    for head, keys in groups.items():
        rows[head] = []
        passing[head] = []
        stations = [(angle, shift) for angle in (range(0, 360, 30) if rotate else [0])
                    for shift in (range(0, 61, 6) if rotate else range(0, 61, 2))]
        for angle, shift in stations:
            transform = A.I.rotation_about(frame['origin_mm'], axis, angle)
            transform[:, 3] += axis * shift
            posed = {key: solids[key].transform(transform) for key in keys}
            hits = []
            for key, part in posed.items():
                for target in fixed:
                    volume = float(A.I.intersection_volume(part, solids[target]))
                    if volume > .001:
                        hits.append(dict(moving=key, target=target, volume_mm3=volume))
                        if rotate:
                            break
                if rotate and hits:
                    break
            adapter = posed[f'elbow:{head}_bone_adapter']
            row = dict(angle_deg=angle, shift_mm=shift, fixed_hits=hits)
            if not rotate or not hits:
                row.update(adapter_bone_overlap_mm3=float(A.I.intersection_volume(adapter, bone)),
                           adapter_bone_min_gap_capped_2mm=float(A.I.gap(adapter, bone, 2)))
            rows[head].append(row)
            if not hits:
                passing[head].append((row, posed))
    pairs = []
    for flex, a in passing['FJ1486']:
        for ext, b in passing['FJ1477']:
            peak = max(float(A.I.intersection_volume(left, right)) for left in a.values() for right in b.values())
            if peak <= .001:
                pairs.append(dict(flexor_shift_mm=flex['shift_mm'], extensor_shift_mm=ext['shift_mm'],
                                  flexor_angle_deg=flex['angle_deg'], extensor_angle_deg=ext['angle_deg'],
                                  cross_group_peak_mm3=peak,
                                  adapter_bone_overlap_mm3=[flex['adapter_bone_overlap_mm3'],
                                                            ext['adapter_bone_overlap_mm3']]))
    pairs.sort(key=lambda row: (row['flexor_shift_mm']+row['extensor_shift_mm'],
                                sum(row['adapter_bone_overlap_mm3'])))
    report = 'elbow_proximal_rotated_station_screen.json' if rotate else 'elbow_proximal_station_screen.json'
    save(folder/report, dict(
        assembly_sha256=sha(manifest), frame_sha256=sha(frame_file), script_sha256=sha(Path(__file__)),
        moving_groups=groups, pending_reconstruction_or_fit=sorted(pending),
        axis_source=axis.tolist(), rows=rows, collision_free_hardware_pairs=pairs,
        only_first_fixed_hit_recorded=rotate,
        scope='Rigid neutral placement screen. Shared-bone seats and holes require remodelling, both muscles require regeneration.'))
    print('clear station counts', {head: len(options) for head, options in passing.items()},
          'clear pairs', len(pairs), 'first', pairs[:3], flush=True)


if __name__ == '__main__':
    main()
