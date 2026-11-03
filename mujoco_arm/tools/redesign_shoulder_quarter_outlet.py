"""Motion-screened FJ1467 quarter-turn outlet candidate, keeping the bone seat.

Unlike the previous fillet, the bend starts 2mm below the cap face and can
exit along the face. Curved tube material is modelled and checked, not an
invisible fixed wire waypoint. Screening is not full shoulder acceptance.
"""
import copy
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
import interference as I
import shoulder_pads
from build_shoulder_unified_muscle import local_export


def main():
    assembly_file=D.DEST/'whole_shoulder_fixture_revision/actual_assembly_check.json'
    assembly=json.loads(assembly_file.read_text(encoding='utf8'));rigid={}
    for key,row in assembly['sources'].items():
        if key.startswith('tpu:'):continue
        path=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed rigid input')
        s=D.SP.load_solid(path) if path.suffix=='.obj' else C.load_stl(path)
        if row['translation_source_mm']:s=s.translate(row['translation_source_mm'])
        rigid[key]=s
    rigid.pop('cap:FJ1467/origin')
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    moving_keys={'bone:humerus'}|{'cap:'+s['head']+'/'+s['end'] for s in seats if s['nearest_bone']=='humerus'}
    moving_keys|={k for k in rigid if k.startswith('stop:')}
    moving={k:v for k,v in rigid.items() if k in moving_keys};fixed={k:v for k,v in rigid.items() if k not in moving_keys}
    bare=json.loads((D.DEST/'origin_stack_revision/search.json').read_text(encoding='utf8'))['bare_paths']
    source=D.DEST/'routes_outlet_candidate.json';spec=json.loads(source.read_text(encoding='utf8'))
    original=next(r for r in spec['routes'] if r['head']=='FJ1467')
    initial=json.loads((D.DEST/'routes_separated.json').read_text(encoding='utf8'))
    old=next(r for r in initial['routes'] if r['head']=='FJ1467')
    n=D.unit(old['ports']['origin']['normal_source']);basis=shoulder_pads.pad_basis(n)
    base0=np.array(old['ports']['origin']['inner_face_source_mm'])
    plate0=C.load_stl(D.frames.OUTPUT_DIR/'shoulder_parts/FJ1467_origin_cap.stl')
    gh=D.frames.joint_centres('source')['GH'];axes=D.frames.gh_probe_axes()
    trials=[];best=None
    phases=[90,75,105,60,120,45,135,30,150,15,165,0,180,195,210,225,240,255,270,285,300,315,330,345]
    for reduction in [1.5,1.,.5,0.]:
        delta=-n*reduction;base=base0+delta
        for tilt in [90,105,120]:
            for phase in phases:
                u=basis[:,0]*np.cos(np.radians(phase))+basis[:,1]*np.sin(np.radians(phase))
                angles=np.r_[np.linspace(0,np.radians(tilt),tilt+1),np.arcsin(2/3)]
                angles=np.unique(angles)
                p0=base-n*2
                arc=p0+3*np.sin(angles)[:,None]*n+3*(1-np.cos(angles))[:,None]*u
                direction=n*np.cos(np.radians(tilt))+u*np.sin(np.radians(tilt))
                end=arc[-1]+direction*5;guide_path=np.vstack([arc,end])
                full_lumen=np.vstack([p0-n*3,guide_path,end+direction*3])
                lumen=D.outer_envelope(full_lumen,.86)
                cap=(plate0.translate(delta)+D.outer_envelope(guide_path,1.8))-lumen
                fixed_hits={k:I.intersection_volume(cap,s) for k,s in fixed.items()}
                trial=dict(reduction_mm=reduction,tilt_deg=tilt,phase_deg=phase,fixed_hits_mm3=fixed_hits,paths=[],score=0.)
                trials.append(trial)
                if any(v>.001 for v in fixed_hits.values()):continue
                if sum(p.volume()>1e-9 for p in cap.decompose())!=1:continue
                scores=[]
                for path in bare:
                    last=0;first=None
                    for angle in path['angles_deg']:
                        tr=I.rotation_about(gh,axes[path['axis']],angle)
                        if any(I.intersection_volume(cap,s.transform(tr))>.001 for s in moving.values()):first=angle;break
                        last=angle
                    trial['paths'].append(dict(axis=path['axis'],sign=path['sign'],last_sample_clear_deg=last,first_hit_deg=first))
                    scores.append(abs(last)/max(1,abs(path['last_sample_clear_deg'])))
                trial['score']=min(scores)
                if best is None or (trial['score'],-reduction)>(best['trial']['score'],-best['trial']['reduction_mm']):
                    best=dict(trial=trial,arc=arc,end=end,direction=direction,cap=cap,guide_path=guide_path,base=base,delta=delta)
                    print('Quarter outlet best',reduction,tilt,phase,'score',trial['score'],[(p['axis'],p['sign'],p['last_sample_clear_deg']) for p in trial['paths']],flush=True)
                if trial['score']==1:break
            if best and best['trial']['score']==1:break
        if best and best['trial']['score']==1:break
    folder=D.DEST/'quarter_outlet_revision';folder.mkdir(exist_ok=True)
    report=dict(scope='cap-only motion screening; full assembly and TPU not accepted',
        source_routes_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_assembly_sha256=hashlib.sha256(assembly_file.read_bytes()).hexdigest(),trials=trials,selected=None,stage_pass=False)
    if best:
        cap,metadata=local_export(best['cap'],folder/'FJ1467_origin_cap.stl',best['base'])
        trial=best['trial'];arc=best['arc'];delta=best['delta'];base=best['base']
        j=int(np.argmin(abs((arc-base)@n)))
        oldpath=np.array(old['connectors']['origin']['path_source_mm'])
        distance=np.r_[0,np.cumsum(np.linalg.norm(np.diff(oldpath,axis=0),axis=1))]
        # Blend the small stack displacement back to the unchanged active core.
        shifted=oldpath+delta*(1-distance/distance[-1])[:,None]
        connector=np.vstack([arc[:j+1][::-1],arc[0]-n,shifted[2:]])
        new_spec=copy.deepcopy(spec);row=next(r for r in new_spec['routes'] if r['head']=='FJ1467')
        row['ports']['origin'].update(inner_face_source_mm=arc[j].tolist(),offset_from_seat_mm=(arc[j]-base).tolist(),
            outer_face_source_mm=arc[-1].tolist(),hand_pull_lead_source_mm=(best['end']+best['direction']*3).tolist())
        row['ports']['origin']['outlet_face_note']='tangential/side outlet, not a straight normal-face port'
        row['origin_stack_reduction_mm']=trial['reduction_mm']
        row['connectors']['origin']['path_source_mm']=connector.tolist()
        # Old geometric connector diagnostics no longer describe this connector.
        row['connectors']['origin'].pop('outer_envelope_bone_intersection_mm3',None)
        row['connectors']['origin']['bone_clear_at_nominal_geometry']=False
        row['rigid_origin_outlet']=dict(path_source_mm=best['guide_path'].tolist(),outer_radius_mm=1.8,lumen_radius_mm=.86,
            bend_radius_mm=3.,tilt_deg=trial['tilt_deg'],inward_start_depth_mm=2.,material_verified=False)
        row['full_wire_path_source_mm']=np.vstack([best['end']+best['direction']*3,best['end'],arc[j+1:][::-1],connector,
            np.array(old['core_path_source_mm'])[1:],np.array(old['connectors']['distal']['path_source_mm'])[-2::-1],old['ports']['distal']['outer_face_source_mm']]).tolist()
        row['status']='candidate_requires_actual_TPU_regeneration_and_full_rigid_motion_check'
        new_spec['scope']='quarter-turn FJ1467 candidate, not stage acceptance'
        deltoid=[r for r in new_spec['routes'] if r['module']=='deltoid']
        new_spec['deltoid_wire_pair_minimum_mm']=[dict(a=a['head'],b=b['head'],distance_mm=D.path_distance(
            np.array(a['full_wire_path_source_mm']),np.array(b['full_wire_path_source_mm']))) for a,b in itertools.combinations(deltoid,2)]
        route_file=D.DEST/'routes_quarter_outlet_candidate.json'
        route_file.write_text(json.dumps(new_spec,indent=2),encoding='utf8')
        report['selected']=dict(**trial,part=metadata,routes_sha256=hashlib.sha256(route_file.read_bytes()).hexdigest(),
            cap_readback_fixed_intersections_mm3={k:I.intersection_volume(cap,s) for k,s in fixed.items()})
    (folder/'check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Quarter candidate generated',bool(best),flush=True)


if __name__=='__main__':main()
