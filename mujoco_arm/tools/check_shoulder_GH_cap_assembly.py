"""Finite staged GH cap/hardware withdrawal and inverse-placement probes.

Bare fixture and actual humerus only; muscles and other absent modules are
listed. A clear selected path is not a full assembly, tool or strength pass.
"""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from check_shoulder_front_liner_keepers import dense_slide
from relieve_shoulder_GH_retaining_caps import DEST


def run():
    source = DEST/'assembly_check.json'
    data = json.loads(source.read_text(encoding='utf8'))
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k,r in data['sources'].items()}
    joint_path = A.G1/'assembly.json'
    n = np.array(json.loads(joint_path.read_text(encoding='utf8'))['neck_axis'])
    bare = {k for k in solids if k.startswith(('fixture:', 'bench:', 'clamp:', 'liner_keeper:'))
            or k in ['bone:scapula', 'bone:humerus', 'clavicle:fixed_support'] or k.startswith('clavicle:M4_')}
    gh = {k for k in solids if k.startswith('GH_')}
    cases = []
    for i in range(2):
        # All hardware is absent when placing the cap. The other cap remains.
        jobs = [(f'GH_cap_{i}', n, 20., bare|{f'GH_cap_{1-i}'}),
                (f'GH_bolt_{i}', n, 15., bare|(gh-{f'GH_bolt_{i}', f'GH_nut_{i}'})),
                (f'GH_nut_{i}', -n, 15., bare|(gh-{f'GH_bolt_{i}', f'GH_nut_{i}'})),
                (f'GH_top_washer_{i}', n, 15., bare|(gh-{f'GH_bolt_{i}', f'GH_nut_{i}', f'GH_top_washer_{i}'})),
                (f'GH_bottom_washer_{i}', -n, 15., bare|(gh-{f'GH_bolt_{i}', f'GH_nut_{i}', f'GH_bottom_washer_{i}'}))]
        for key, direction, maximum, fixed_keys in jobs:
            result = dense_slide({key:solids[key]}, {k:solids[k] for k in sorted(fixed_keys)}, direction, maximum)
            firsts = []
            for phase in result['phases']:
                hit = next((s for s in phase['samples'] if s['hits_over_threshold']), None)
                firsts.append(dict(phase=phase['name'], first_hit=hit))
            row = dict(part=key, path=result, first_interferences=firsts,
                       not_yet_installed_keys=sorted(set(solids)-fixed_keys-{key}))
            cases.append(row)
            print(key, result['sampled_path_clear'], firsts if not result['sampled_path_clear'] else '', flush=True)
    report = dict(scope=__doc__, assembly_source=source.relative_to(A.ROOT).as_posix(),
                  assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
                  joint_reference_sha256=A.hashed(joint_path),
                  dense_slide_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/check_shoulder_front_liner_keepers.py'),
                  cases=cases, step_mm=.1, all_source_hashes_match=True,
                  all_selected_paths_clear=all(r['path']['sampled_path_clear'] for r in cases),
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['full assembly ordering, temporary holding and tool/head access',
                                'thread engagement, preload, friction and strength',
                                'muscle placement/deformation and motion-range acceptance'])
    (DEST/'staged_cap_hardware_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
