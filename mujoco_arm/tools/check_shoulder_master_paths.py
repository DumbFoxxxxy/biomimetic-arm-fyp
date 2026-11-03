"""Audit the complete static steel routes, including external pull leads.

Does not infer deformable motion, channel friction or a physical drive range.
Sections are overlapping diagnostic probes; their volumes must not be summed.
"""
import hashlib
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import build_shoulder_tendon_segments as T
import interference as I


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--routes',type=Path,default=D.DEST/'routes_separated.json')
    args=ap.parse_args();source=args.routes
    spec=json.loads(source.read_text(encoding='utf8'))
    bone_files={k:T.BONE_DIR/v for k,v in D.SP.BONES.items()}
    bones={k:D.SP.load_solid(path) for k,path in bone_files.items()}
    rows=[]
    for route in spec['routes']:
        p=route['ports']['origin']
        sections={
            'complete':route['full_wire_path_source_mm'],
            'external_pull_lead':[p['hand_pull_lead_source_mm'],p['outer_face_source_mm']],
            'origin_connector':route['connectors']['origin']['path_source_mm'],
            'retained_core':route['core_path_source_mm'],
            'distal_connector':route['connectors']['distal']['path_source_mm'],
        }
        checks={}
        for name,points in sections.items():
            probes={}
            for label,radius in [('steel',.3),('nominal_lumen',.86)]:
                body=D.outer_envelope(np.asarray(points),radius)
                probes[label]={k:I.intersection_volume(body,bone) for k,bone in bones.items()}
            checks[name]=probes
        row=dict(head=route['head'],module=route['module'],sections=checks)
        rows.append(row)
        print(route['head'],'steel/bone',checks['complete']['steel'],flush=True)
    failed=[r['head'] for r in rows if any(v>.001 for v in r['sections']['complete']['steel'].values())]
    report=dict(scope='six complete static steel paths versus three actual bones; not complete assembly or motion',
        source=str(source.relative_to(D.frames.PROJECT_ROOT)),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        bone_sha256={k:hashlib.sha256(path.read_bytes()).hexdigest() for k,path in bone_files.items()},
        steel_radius_mm=.3,nominal_lumen_probe_radius_mm=.86,
        intersection_reporting_threshold_mm3=.001,rows=rows,
        steel_bone_clear_at_reported_threshold=not failed,failed_heads=failed,
        section_volumes_are_not_additive=True,stage_pass=False,
        missing=['external lead guide and conduit reaction fixation',
                 'complete installed materials and hardware',
                 'deformation, motion, friction and physical tests'])
    name='master_paths_bone_check.json' if source.name=='routes_separated.json' else source.stem+'_bone_check.json'
    (D.DEST/name).write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Static steel/bone failures:',failed,flush=True)


if __name__=='__main__':main()
