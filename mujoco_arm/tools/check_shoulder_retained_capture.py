"""Finite straight head translations against actual retained GH meshes.

This checks only 26 straight translations at the neutral neck orientation.
It does not prove arbitrary escape paths, retention strength or assembly.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A


def run(source):
    report=json.loads(source.read_text(encoding='utf8'))
    fixed={}
    for key in ['bone:scapula','GH_cap_0','GH_cap_1']:
        row=report['sources'][key]
        fixed[key]=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])[0]
    joint_file=A.G1/'assembly.json'
    radius=json.loads(joint_file.read_text(encoding='utf8'))['ball_radius_mm']
    centre=A.frames.joint_centres('source')['GH']
    ball=A.I._manifold().Manifold.sphere(radius,40).translate(centre)
    rows=[]
    for raw in itertools.product([-1,0,1],repeat=3):
        if not any(raw):continue
        axis=np.array(raw,dtype=float);axis/=np.linalg.norm(axis);first=None
        for step in range(1,31):
            moved=ball.translate(axis*step*.1)
            hits=[dict(part=k,volume_mm3=A.I.intersection_volume(moved,s)) for k,s in fixed.items()]
            hits=[h for h in hits if h['volume_mm3']>.001]
            if hits:
                first=dict(distance_mm=step*.1,hits=hits);break
        rows.append(dict(axis=axis.tolist(),first=first))
    result=dict(scope=__doc__,assembly_sha256=A.hashed(source),
        script_sha256=A.hashed(Path(__file__).resolve()),joint_definition_sha256=A.hashed(joint_file),
        radius_mm=radius,step_mm=.1,maximum_mm=3,intersection_threshold_mm3=.001,
        rows=rows,all_directions_encounter_material=all(x['first'] is not None for x in rows),
        stage_pass=False)
    (source.parent/'finite_retention_probes.json').write_text(json.dumps(result,indent=2)+'\n',
                                                            encoding='utf8',newline='\n')
    print('Finite capture:',sum(x['first'] is not None for x in rows),'of',len(rows),
          'max first',max((x['first']['distance_mm'] for x in rows if x['first']),default=None))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assembly',type=Path,default=A.DEST/'GH_only_candidate/assembly_check.json')
    run(parser.parse_args().assembly)
