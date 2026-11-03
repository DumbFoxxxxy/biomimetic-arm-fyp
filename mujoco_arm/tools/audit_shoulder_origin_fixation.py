"""Measure actual proximal M3 axes before designing bone reaction mounts.

RayHit.distance is a normalized segment parameter, not millimetres. The
positions are independently projected to verify depth. No holes are drilled
and no unbuilt spacer/backplate is counted as an existing connection.
"""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_tendon_link as L
import audit_shoulder_retained_joint_fit as A
import shoulder_pads


def run():
    source=L.DEST/'assembly_check.json';assembly=json.loads(source.read_text(encoding='utf8'))
    meta=assembly['sources']['bone:scapula']
    bone,_,_=A.load_checked(A.ROOT/meta['file'],meta['translation_source_mm'],meta['sha256'])
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(route_file.read_text(encoding='utf8'))['routes'];rows=[]
    for head in ['FJ1500','FJ1504']:
        port=next(r for r in routes if r['head']==head)['ports']['origin']
        inner=np.array(port['inner_face_source_mm']);outer=np.array(port['outer_face_source_mm'])
        n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
        holes=[]
        for sign in [-1,1]:
            start=outer+F[:,0]*(7*sign)+n*3
            hits=bone.ray_cast(start,start-n*60.)
            crossings=[]
            for hit in hits:
                depth=float(np.dot(start-np.array(hit.position),n))
                if abs(depth-hit.distance*60.)>1e-6:raise ValueError('Ray depth convention mismatch')
                crossings.append(dict(source_position_mm=list(hit.position),depth_from_ray_start_mm=depth,
                    position_along_cap_normal_mm=float(np.dot(np.array(hit.position)-inner,n)),
                    surface_normal=list(hit.normal)))
            if len(crossings)!=2:raise ValueError('Expected two blade crossings; review separately: '+head)
            holes.append(dict(sign=sign,axis_point_source_mm=(inner+F[:,0]*7*sign).tolist(),
                normal_source=n.tolist(),bone_crossings=crossings,
                cap_inner_to_bone_front_mm=-crossings[0]['position_along_cap_normal_mm'],
                bone_stock_mm=crossings[1]['depth_from_ray_start_mm']-crossings[0]['depth_from_ray_start_mm']))
        rows.append(dict(head=head,origin_cap_key='cap:'+head+'/origin',holes=holes,
            existing_fixation_verified=False,
            required='actual front support and blade load-spreading/backing connection; nominal M3 holes alone do not close the reaction path'))
    out=L.BASE/'origin_fixation_review';out.mkdir(parents=True,exist_ok=True)
    report=dict(scope=__doc__,source_sha256=A.hashed(source),bone_source=meta,
        routes_sha256=A.hashed(route_file),script_sha256=A.hashed(Path(__file__).resolve()),
        bolt_hole_radius_mm=1.7,bolt_pitch_mm=14.,ray_length_mm=60.,rows=rows,
        not_verified=['support geometry, fastener access and complete reaction path',
            'blade bearing strength, preload, printed materials and pull tests'],stage_pass=False)
    (out/'fixation_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    for row in rows:print(row['head'],[(round(h['cap_inner_to_bone_front_mm'],3),round(h['bone_stock_mm'],3)) for h in row['holes']],flush=True)
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
