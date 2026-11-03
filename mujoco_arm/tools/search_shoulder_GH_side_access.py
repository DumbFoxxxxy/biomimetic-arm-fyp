"""Search finite straight GH subassembly entrances without changing geometry.

Bare fixture, scapula and humerus. Search failures are directional evidence,
not an impossibility proof. Selected clear paths require reverse dense checks.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from relieve_shoulder_GH_retaining_caps import DEST as SOURCE_DIR
from audit_shoulder_GH_assembly import scan
from check_shoulder_front_liner_keepers import dense_slide

DEST = A.DESIGN/'s1_GH_assembly_review/side_access_search'


def run():
    DEST.mkdir(parents=True, exist_ok=True)
    source = SOURCE_DIR/'assembly_check.json'
    data = json.loads(source.read_text(encoding='utf8'))
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k,r in data['sources'].items()}
    reference = A.G1/'assembly.json'
    joint = json.loads(reference.read_text(encoding='utf8'))
    n = np.array(joint['neck_axis']); x = A.D.unit(np.cross(n,[0.,-1.,0.]))
    F = np.column_stack([x,np.cross(n,x),n])
    directions = []
    for basis in [np.eye(3), F]:
        for v in itertools.product([-1,0,1], repeat=3):
            if not any(v):continue
            d = A.D.unit(basis@v)
            if not any(np.linalg.norm(d-other)<1e-6 for other in directions):directions.append(d)
    bare = {k for k in solids if k.startswith(('fixture:', 'bench:', 'clamp:', 'liner_keeper:'))
            or k in ['bone:scapula','bone:humerus','clavicle:fixed_support'] or k.startswith('clavicle:M4_')}
    gh = {k for k in solids if k.startswith('GH_')}
    cases = []
    for i in range(2):
        jobs = [(f'cap_and_top_washer_{i}', [f'GH_cap_{i}',f'GH_top_washer_{i}'], bare|{f'GH_cap_{1-i}'}),
                (f'nut_{i}', [f'GH_nut_{i}'], bare|(gh-{f'GH_bolt_{i}',f'GH_nut_{i}'})),
                (f'bottom_washer_{i}', [f'GH_bottom_washer_{i}'], bare|(gh-{f'GH_bolt_{i}',f'GH_nut_{i}',f'GH_bottom_washer_{i}'}))]
        for label, moving_keys, fixed_keys in jobs:
            moving = A.I._manifold().Manifold.batch_boolean([solids[k] for k in moving_keys], A.I._manifold().OpType.Add)
            fixed = {k:solids[k] for k in sorted(fixed_keys)}
            attempts=[]; selected=None
            for direction in directions:
                path=scan(moving,direction,fixed,maximum=25.,step=.1)
                first=path['first_interference']
                attempts.append(dict(direction_source=direction.tolist(), entire_sampled_path_clear=path['entire_sampled_path_clear'],
                                     first_hit_mm=None if first is None else first['distance_mm'],
                                     first_hits=[] if first is None else first['hits_over_threshold']))
                if path['entire_sampled_path_clear']:
                    selected=dense_slide({k:solids[k] for k in moving_keys},fixed,direction,25.)
                    break
            row=dict(name=label,moving_keys=moving_keys,fixed_keys=sorted(fixed_keys),
                     not_yet_installed_keys=sorted(set(solids)-fixed_keys-set(moving_keys)),
                     attempts=attempts,selected_path=selected)
            cases.append(row)
            print(label,'attempts',len(attempts),'selected',None if selected is None else selected['direction_source'],flush=True)
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),joint_reference_sha256=A.hashed(reference),
                step_mm=.1,maximum_translation_mm=25.,direction_count=len(directions),cases=cases,
                stage_pass=False,actual_drive_range_deg=None,
                not_verified=['tool access, temporary holding and full assembly order',
                              'coupled translation/rotation and unsearched directions',
                              'preload, thread engagement, print fit, strength and muscle drive'])
    (DEST/'search_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
