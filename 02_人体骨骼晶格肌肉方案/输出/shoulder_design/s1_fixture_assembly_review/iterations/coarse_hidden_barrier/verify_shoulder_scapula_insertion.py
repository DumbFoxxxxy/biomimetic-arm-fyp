"""Dense independent geometric samples of the selected scapula insertion step.

Explicitly requires preplaced front liners, no GH/humerus/clavicle bone or
muscles, and enumerates present/omitted objects. Temporary liner holding,
tools, strength and the complete assembly sequence remain unverified.
"""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import audit_shoulder_scapula_insertion as S


def gap_capped(body, obstacle, cap=2.):
    a, b = np.array(body.bounding_box()), np.array(obstacle.bounding_box())
    separation = np.maximum(np.maximum(a[:3]-b[3:], b[:3]-a[3:]), 0.)
    if np.linalg.norm(separation) >= cap:
        return cap
    return A.I.gap(body, obstacle, cap)


def run():
    audit_file = S.DEST/'scapula_insertion_audit.json'
    audit = json.loads(audit_file.read_text(encoding='utf8'))
    assembly = A.ROOT/audit['assembly_source']
    if audit['assembly_sha256'] != A.hashed(assembly):
        raise ValueError('Changed assembly since insertion audit')
    if audit['script_sha256'] != A.hashed(Path(S.__file__).resolve()):
        raise ValueError('Changed insertion audit script')
    stage = next(r for r in audit['stages'] if r['name'] == 'bone_against_preplaced_front_liners')
    candidates = [(i, r) for i, r in enumerate(stage['trials']) if r['sampled_path_clear']]
    if not candidates:
        raise ValueError('No sampled candidate in stated liner-preplacement stage')
    index, trial = min(candidates, key=lambda row:np.linalg.norm(
        np.array(row[1]['direction_source'])-np.array([0., -1, 0])))
    direction = np.array(trial['direction_source'])
    for key, row in audit['source_manifest'].items():
        if A.hashed(A.ROOT/row['file']) != row['sha256']:
            raise ValueError('Changed source: '+key)
    present = stage['moving_keys']+stage['fixed_keys']
    solids = {}
    for key in present:
        row = audit['source_manifest'][key]
        solids[key], _, _ = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
    body = solids['bone:scapula']
    fixed = {k:solids[k] for k in stage['fixed_keys']}
    maximum = audit['maximum_translation_mm']
    distances = np.linspace(0., maximum, int(round(maximum/.1))+1)
    phases = []
    for phase, positions in [('withdrawal', distances), ('inverse_insertion', distances[::-1])]:
        rows = []
        for i, distance in enumerate(positions):
            # Recompute the transform and all pair checks in each direction;
            # don't reuse cached forward results for the inverse path.
            posed = body.translate(direction*distance)
            volumes = {k:A.I.intersection_volume(posed, s) for k, s in fixed.items()}
            gaps = {k:gap_capped(posed, s) for k, s in fixed.items()}
            nearest = min(gaps, key=gaps.get)
            rows.append(dict(displacement_mm=round(float(distance), 6),
                             source_translation_mm=(direction*distance).tolist(),
                             pair_volumes_mm3=volumes,
                             maximum_pair_overlap_mm3=max(volumes.values()),
                             hits_over_threshold=[dict(part=k, volume_mm3=v) for k, v in volumes.items() if v > .001],
                             minimum_gap_capped_2mm=gaps[nearest], nearest_part=nearest))
            if i % 50 == 0:
                print(phase, i+1, '/', len(positions), 'pose', round(float(distance), 2), flush=True)
        phases.append(dict(name=phase, sample_count=len(rows), samples=rows,
                           all_samples_below_interference_threshold=all(not r['hits_over_threshold'] for r in rows)))
    selected = dict(scope='Selected partial assembly scene; one scapula slide, not shoulder actuation',
                    n_instances=len(present), sources={k:audit['source_manifest'][k] for k in present},
                    moving_keys=stage['moving_keys'], fixed_keys=stage['fixed_keys'],
                    direction_source=direction.tolist(), maximum_translation_mm=maximum,
                    not_yet_installed_keys=stage['not_yet_installed_keys'],
                    stage_pass=False, actual_drive_range_deg=None)
    (S.DEST/'selected_stage_scene.json').write_text(json.dumps(selected, indent=2)+'\n', encoding='utf8', newline='\n')
    result = dict(scope=__doc__, audit_sha256=A.hashed(audit_file), assembly_sha256=A.hashed(assembly),
                  source_assembly=audit['assembly_source'], script_sha256=A.hashed(Path(__file__).resolve()),
                  audit_script_sha256=A.hashed(Path(S.__file__).resolve()),
                  all_source_hashes_match=True, source_object_count=len(audit['source_manifest']),
                  stage_name=stage['name'], selected_coarse_trial_index=index,
                  selected_direction_source=direction.tolist(), maximum_translation_mm=maximum,
                  moving_keys=stage['moving_keys'], fixed_keys=stage['fixed_keys'],
                  not_yet_installed_keys=stage['not_yet_installed_keys'],
                  pose_step_mm=.1, interference_threshold_mm3=.001, gap_query_cap_mm=2., phases=phases,
                  selected_partial_path_samples_clear=all(p['all_samples_below_interference_threshold'] for p in phases),
                  required_but_unverified=['three front liners must remain at registered positions before bone insertion',
                                          'actual liner temporary retention/placement and hand/tool clearance',
                                          'subsequent rear liners, jaws, nuts/washers, GH and muscle installation'],
                  stage_pass=False, actual_drive_range_deg=None)
    (S.DEST/'selected_path_dense_check.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Dense partial path samples clear', result['selected_partial_path_samples_clear'],
          'maximum overlap', max(r['maximum_pair_overlap_mm3'] for p in phases for r in p['samples']), flush=True)
    return result


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
