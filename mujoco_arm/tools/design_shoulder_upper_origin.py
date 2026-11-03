"""Screen real upper acromial surface candidates for the FJ1467 interface.

No bone is moved or cut. This changes the engineering attachment patch and
rebuilds its cap; full fastening, flexible motion and drive remain unverified.
"""
import copy
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import check_shoulder_assembly as C
import build_shoulder_ends as BE
import build_shoulder_tendon_segments as T
import interference as I
import mesh_geom
from build_shoulder_unified_muscle import local_export


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--edge-offset',type=float,default=13.5)
    args=ap.parse_args()
    if not 12<=args.edge_offset<=16:raise ValueError('Inspect before extending the local transition further')
    assembly_file=D.DEST/'whole_shoulder_fixture_revision/actual_assembly_check.json'
    assembly=json.loads(assembly_file.read_text(encoding='utf8'));rigid={}
    for key,row in assembly['sources'].items():
        if key.startswith('tpu:'):continue
        file=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(file.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed rigid input')
        s=D.SP.load_solid(file) if file.suffix=='.obj' else C.load_stl(file)
        if row['translation_source_mm']:s=s.translate(row['translation_source_mm'])
        rigid[key]=s
    rigid.pop('cap:FJ1467/origin')
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    moving_keys={'bone:humerus'}|{'cap:'+s['head']+'/'+s['end'] for s in seats if s['nearest_bone']=='humerus'}
    moving_keys|={k for k in rigid if k.startswith('stop:')}
    moving={k:v for k,v in rigid.items() if k in moving_keys};fixed={k:v for k,v in rigid.items() if k not in moving_keys}
    source=D.DEST/'routes_outlet_candidate.json';spec=json.loads(source.read_text(encoding='utf8'))
    old=next(r for r in spec['routes'] if r['head']=='FJ1467');old_c=np.array(old['ports']['origin']['bone_seat_source_mm'])
    V,F=mesh_geom.load_obj(str(T.BONE_DIR/D.SP.BONES['scapula']));V,F=mesh_geom.weld(V,F)
    triangles=V[F];norm=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);normals=np.zeros_like(V)
    for k in range(3):np.add.at(normals,F[:,k],norm)
    normals/=np.maximum(np.linalg.norm(normals,axis=1)[:,None],1e-20)
    distances=np.linalg.norm(V-old_c,axis=1)
    ids=np.flatnonzero((distances<10)&(normals[:,2]>.9)&(V[:,0]<old_c[0]))
    ids=ids[np.argsort(distances[ids])]
    bare=json.loads((D.DEST/'origin_stack_revision/search.json').read_text(encoding='utf8'))['bare_paths']
    gh=D.frames.joint_centres('source')['GH'];axes=D.frames.gh_probe_axes();trials=[];best=None
    for index in ids[:25]:
        c=V[index];n=normals[index]
        for pad_t in [2.5,4.]:
            inner=c+n*(pad_t+.2);cap=BE.disc(inner,n,1.5,radius=9.5,holes=True)
            hits={k:I.intersection_volume(cap,s) for k,s in fixed.items()}
            trial=dict(source_vertex_index=int(index),centre_source_mm=c.tolist(),normal_source=n.tolist(),
                displacement_mm=float(distances[index]),pad_thickness_mm=pad_t,cap_radius_mm=9.5,cap_thickness_mm=1.5,
                fixed_intersections_mm3=hits,cap_only_paths=[],score=0.)
            trials.append(trial)
            if any(v>.001 for v in hits.values()):continue
            scores=[]
            for path in bare:
                last=0;first=None
                for angle in path['angles_deg']:
                    tr=I.rotation_about(gh,axes[path['axis']],angle)
                    if any(I.intersection_volume(cap,s.transform(tr))>.001 for s in moving.values()):first=angle;break
                    last=angle
                trial['cap_only_paths'].append(dict(axis=path['axis'],sign=path['sign'],last_sample_clear_deg=last,first_hit_deg=first))
                scores.append(abs(last)/max(1,abs(path['last_sample_clear_deg'])))
            trial['score']=min(scores)
            if best is None or (trial['score'],-trial['displacement_mm'])>(best['trial']['score'],-best['trial']['displacement_mm']):
                best=dict(trial=trial,c=c,n=n,inner=inner,cap=cap)
                print('Upper origin best',trial['displacement_mm'],pad_t,trial['score'],[(p['axis'],p['sign'],p['last_sample_clear_deg']) for p in trial['cap_only_paths']],flush=True)
            if trial['score']==1:break
        if best and best['trial']['score']==1:break
    tag='' if args.edge_offset==13.5 else '_edge'+str(args.edge_offset).replace('.','p')
    folder=D.DEST/('upper_origin_revision'+tag);folder.mkdir(exist_ok=True)
    report=dict(scope='upper acromial cap-only screening, not full rigid motion or attachment acceptance',
        source_routes_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_assembly_sha256=hashlib.sha256(assembly_file.read_bytes()).hexdigest(),trials=trials,selected=None,stage_pass=False)
    if best:
        c=best['c'];n=best['n'];inner=best['inner'];pad_t=best['trial']['pad_thickness_mm']
        core=np.array(old['core_path_source_mm']);side=D.unit(np.array([-1.,0,0])-n*np.dot(n,[-1.,0,0]))
        bend=inner-n*.8;edge=c+side*args.edge_offset+n*(pad_t-.6)
        tris=[]
        for file in D.SP.BONES.values():
            v,f=mesh_geom.load_obj(str(T.BONE_DIR/file));v,f=mesh_geom.weld(v,f);tris.append(v[f])
        connector=T.relax_off_bone([inner,bend,edge,core[0]],tris,1.8,step_mm=.5,iterations=80,lam=.1)
        outer=inner+n*1.5;lead=outer+n*10
        updated=copy.deepcopy(spec);row=next(r for r in updated['routes'] if r['head']=='FJ1467')
        row.pop('rigid_origin_outlet',None)
        if tag:row['origin_connector_edge_offset_mm']=args.edge_offset
        row['origin_seat_override']=dict(centre_source_mm=c.tolist(),normal_source=n.tolist(),
            pad_radius_mm=11.,pad_thickness_mm=pad_t,lift_mm=0.,pad_boss_mm=0.,
            cap_radius_mm=9.5,cap_thickness_mm=1.5)
        row['ports']['origin']=dict(bone_seat_source_mm=c.tolist(),offset_from_seat_mm=[0,0,0],normal_source=n.tolist(),
            inner_face_source_mm=inner.tolist(),outer_face_source_mm=outer.tolist(),hand_pull_lead_source_mm=lead.tolist())
        row['connectors']['origin']=dict(path_source_mm=np.array(connector).tolist(),bone_clear_at_nominal_geometry=False)
        full=np.vstack([lead,outer,connector,core[1:],np.array(old['connectors']['distal']['path_source_mm'])[-2::-1],old['ports']['distal']['outer_face_source_mm']])
        row['full_wire_path_source_mm']=full.tolist();row['status']='upper seat candidate requires actual TPU and full rigid checks'
        updated['scope']='upper acromial FJ1467 engineering seat candidate; no bone geometry change'
        deltoid=[r for r in updated['routes'] if r['module']=='deltoid']
        updated['deltoid_wire_pair_minimum_mm']=[dict(a=a['head'],b=b['head'],distance_mm=D.path_distance(
            np.array(a['full_wire_path_source_mm']),np.array(b['full_wire_path_source_mm']))) for a,b in itertools.combinations(deltoid,2)]
        new_source=D.DEST/('routes_upper_origin'+tag+'_candidate.json');new_source.write_text(json.dumps(updated,indent=2),encoding='utf8')
        cap=best['cap']-D.outer_envelope(full,.86)
        actual,part=local_export(cap,folder/'FJ1467_origin_cap.stl',c)
        wire=D.outer_envelope(full,.3)
        report['selected']=dict(**best['trial'],part=part,new_routes_sha256=hashlib.sha256(new_source.read_bytes()).hexdigest(),
            wire_bone_intersections_mm3={k:I.intersection_volume(wire,rigid['bone:'+k]) for k in D.SP.BONES},
            cap_readback_fixed_intersections_mm3={k:I.intersection_volume(actual,s) for k,s in fixed.items()},
            cap_readback_steel_intersection_mm3=I.intersection_volume(wire,actual),
            connector_length_mm=float(np.linalg.norm(np.diff(connector,axis=0),axis=1).sum()))
    (folder/'check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Upper origin candidate generated',bool(best),flush=True)


if __name__=='__main__':main()
