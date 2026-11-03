"""Materialize selected compact terminal routes and actual local parts."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import compact_shoulder_terminal as CT
import build_shoulder_tendon_segments as T
import mesh_geom
import interference as I
from build_shoulder_unified_muscle import local_export


def main():
    source=D.DEST/'routes_upper_origin_edge14p0_candidate.json';spec=json.loads(source.read_text(encoding='utf8'))
    search_file=D.DEST/'compact_distal_revision/search.json';search=json.loads(search_file.read_text(encoding='utf8'))
    if hashlib.sha256(source.read_bytes()).hexdigest()!=search['master_sha256']:raise ValueError('Changed search master')
    updated=copy.deepcopy(spec);bones={k:D.SP.load_solid(T.BONE_DIR/v) for k,v in D.SP.BONES.items()};tris=[]
    for name in D.SP.BONES.values():
        v,f=mesh_geom.load_obj(str(T.BONE_DIR/name));v,f=mesh_geom.weld(v,f);tris.append(v[f])
    folder=D.DEST/'compact_distal_revision';parts={};checks={}
    for head,study in search['results'].items():
        selected=study['selected'];row=next(r for r in updated['routes'] if r['head']==head)
        geometry=selected['geometry'];n=D.unit(geometry['normal_source']);reduction=selected['reduction_mm'];delta=-n*reduction
        inner=np.array(geometry['inner_face_source_mm']);outer=np.array(geometry['outer_face_source_mm'])
        core=np.array(row['core_path_source_mm'])
        # Approach from the pad side on a real 1mm axial inlet, instead of
        # warping the old high cap route back through the stop counterbore.
        connector=np.vstack([inner,T.relax_off_bone([inner-n,core[-1]],tris,1.8,step_mm=.5,iterations=80,lam=.1)])
        row['ports']['distal'].update(inner_face_source_mm=inner.tolist(),outer_face_source_mm=outer.tolist())
        row['connectors']['distal']=dict(path_source_mm=connector.tolist(),bone_clear_at_nominal_geometry=False)
        row['distal_stack_reduction_mm']=reduction;row['distal_terminal_override']=geometry
        full=np.vstack([row['ports']['origin']['hand_pull_lead_source_mm'],row['ports']['origin']['outer_face_source_mm'],
            row['connectors']['origin']['path_source_mm'],core[1:],connector[-2::-1],outer])
        row['full_wire_path_source_mm']=full.tolist();row['status']='compact distal and nominal stop; mounting and deformation unverified'
        cap,stop,_=CT.build(inner,n,geometry['tangent_source']);cap=cap-D.outer_envelope(full,.86)
        actual={}
        for kind,solid in [('cap',cap),('stop_blank',stop)]:
            actual[kind],parts[head+'_'+kind]=local_export(solid,folder/(head+'_'+kind+'.stl'),inner)
        wire=D.outer_envelope(full,.3)
        checks[head]=dict(wire_bone_intersections_mm3={k:I.intersection_volume(wire,s) for k,s in bones.items()},
            wire_part_intersections_mm3={k:I.intersection_volume(wire,s) for k,s in actual.items()},
            stop_capture_inward_mm3={str(d):I.intersection_volume(actual['stop_blank'].translate(-n*d),actual['cap']) for d in [0,.04,.06]},
            connector_length_mm=float(np.linalg.norm(np.diff(connector,axis=0),axis=1).sum()))
    updated['scope']='compact distal plate and stop candidates; full fixed-GH/TPU checks pending'
    route_file=D.DEST/'routes_compact_distal_candidate.json';route_file.write_text(json.dumps(updated,indent=2),encoding='utf8')
    report=dict(scope='local compact caps and uncrimped stops; bone mounting holes not yet implemented',
        master_sha256=hashlib.sha256(route_file.read_bytes()).hexdigest(),parts=parts,checks=checks,
        search_sha256=hashlib.sha256(search_file.read_bytes()).hexdigest(),stage_pass=False,
        not_verified=['bone mounting hole, bolt head and fixation/anti-rotation','crimp grip and strength',
            'complete static, moving and flexible assembly'])
    (folder/'parts_check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Local compact meshes',{k:r['mesh_checks_pass'] for k,r in parts.items()},'checks',checks,flush=True)


if __name__=='__main__':main()
