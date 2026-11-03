"""Finite diagnostic-ball translation probes after GH entry openings.

World grid, recorded joint basis grid and true neck axis are sampled.
No general escape, strength, dynamic stability or muscle-drive proof.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from open_shoulder_GH_assembly_entries import DEST
from audit_shoulder_GH_assembly import scan


def run():
    source=DEST/'assembly_check.json';data=json.loads(source.read_text(encoding='utf8'))
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in data['sources'].items()}
    joint=json.loads((A.G1/'assembly.json').read_text(encoding='utf8'));n=np.array(joint['neck_axis']);G=np.array(joint['centre_source_mm'])
    x=A.D.unit(np.cross(n,[0.,-1.,0.]));F=np.column_stack([x,np.cross(n,x),n])
    directions=[]
    for basis in [np.eye(3),F]:
        for vector in itertools.product([-1.,0.,1.],repeat=3):
            if not any(vector):continue
            direction=A.D.unit(basis@vector)
            if not any(np.linalg.norm(direction-other)<1e-6 for other in directions):directions.append(direction)
    fixed={k:S[k] for k in ['bone:scapula','GH_cap_0','GH_cap_1']}
    ball=A.I._manifold().Manifold.sphere(joint['ball_radius_mm'],40).translate(G)
    rows=[]
    for direction in directions:
        path=scan(ball,direction,fixed,25.,.1)
        first=path['first_interference']
        rows.append(dict(direction_source=direction.tolist(),first_hit_mm=None if first is None else first['distance_mm'],
                         diagnostic_path=path))
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),fixed_keys=list(fixed),ball_radius_mm=joint['ball_radius_mm'],
                step_mm=.1,maximum_requested_translation_mm=25.,directions=rows,
                all_sampled_directions_hit_before_25mm=all(r['first_hit_mm'] is not None for r in rows),
                stage_pass=False,actual_drive_range_deg=None)
    (DEST/'finite_capture_direction_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Directions',len(rows),'all hit',report['all_sampled_directions_hit_before_25mm'],
          'first hit range',[min(r['first_hit_mm'] for r in rows if r['first_hit_mm'] is not None),
                             max(r['first_hit_mm'] for r in rows if r['first_hit_mm'] is not None)],flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
