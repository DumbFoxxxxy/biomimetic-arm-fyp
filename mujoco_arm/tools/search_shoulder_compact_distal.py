"""Moving-geometry search for compact FJ1500/FJ1504 terminal plates."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import check_shoulder_assembly as C
import interference as I
import compact_shoulder_terminal as CT


def main():
    source=D.DEST/'whole_shoulder_upper_edge14/actual_assembly_check.json'
    assembly=json.loads(source.read_text(encoding='utf8'));rigid={}
    for key,row in assembly['sources'].items():
        if key.startswith('tpu:'):continue
        file=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(file.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed input '+key)
        s=D.SP.load_solid(file) if file.suffix=='.obj' else C.load_stl(file)
        if row['translation_source_mm']:s=s.translate(row['translation_source_mm'])
        rigid[key]=s
    route_file=D.DEST/'routes_upper_origin_edge14p0_candidate.json'
    routes={r['head']:r for r in json.loads(route_file.read_text(encoding='utf8'))['routes']}
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    moving_keys={'bone:humerus'}|{'cap:'+s['head']+'/'+s['end'] for s in seats if s['nearest_bone']=='humerus'}
    moving_keys|={k for k in rigid if k.startswith('stop:')}
    fixed={k:v for k,v in rigid.items() if k not in moving_keys}
    stack_report=json.loads((D.DEST/'distal_stack_revision/search.json').read_text(encoding='utf8'))
    if stack_report['assembly_sha256']!=hashlib.sha256(source.read_bytes()).hexdigest():raise ValueError('Changed bare path baseline')
    bare=stack_report['bare_paths'];gh=D.frames.joint_centres('source')['GH'];axes=D.frames.gh_probe_axes()
    results={}
    for head in ['FJ1500','FJ1504']:
        key='cap:'+head+'/distal';row=routes[head];port=row['ports']['distal'];n=D.unit(port['normal_source'])
        along=D.unit(axes['longitudinal']-n*np.dot(axes['longitudinal'],n));side=np.cross(n,along)
        maximum=stack_report['results'][head]['max_reduction_mm'];trials=[];best=None
        for reduction in np.r_[np.arange(0,maximum,.5),maximum]:
            inner=np.array(port['inner_face_source_mm'])-n*reduction
            for phase in [0,45,-45,90,-90,135,-135,180]:
                tangent=along*np.cos(np.radians(phase))+side*np.sin(np.radians(phase))
                cap,stop,geometry=CT.build(inner,n,tangent);candidate=cap+stop
                same={k:I.intersection_volume(candidate,s) for k,s in rigid.items() if k in moving_keys and k!=key}
                trial=dict(reduction_mm=float(reduction),phase_deg=phase,geometry=geometry,
                    same_group_intersections_mm3=same,paths=[],score=0.)
                trials.append(trial)
                entry=D.outer_envelope(np.array([inner,inner-n]),1.61)
                entry_hit=I.intersection_volume(entry,rigid['bone:humerus'])
                trial['entry_0_75wall_guard_bone_intersection_mm3']=entry_hit
                if entry_hit>.001:continue
                if any(v>.001 for v in same.values()):continue
                scores=[]
                for path in bare:
                    last=0;first=None;hits=[]
                    for angle in path['angles_deg']:
                        posed=candidate.transform(I.rotation_about(gh,axes[path['axis']],angle))
                        hits=[dict(fixed=k,volume_mm3=I.intersection_volume(posed,s)) for k,s in fixed.items()]
                        hits=[h for h in hits if h['volume_mm3']>.001]
                        if hits:first=angle;break
                        last=angle
                    trial['paths'].append(dict(axis=path['axis'],sign=path['sign'],last_sample_clear_deg=last,first_hit_deg=first,hits=hits))
                    scores.append(abs(last)/max(1,abs(path['last_sample_clear_deg'])))
                trial['score']=min(scores)
                if best is None or (trial['score'],-trial['reduction_mm'])>(best['score'],-best['reduction_mm']):
                    best=trial
                    print(head,'compact best',reduction,phase,'score',trial['score'],'limits',[p['last_sample_clear_deg'] for p in trial['paths']],flush=True)
                if best and best['score']==1:break
            if best and best['score']==1:break
        results[head]=dict(trials=trials,selected=best)
    folder=D.DEST/'compact_distal_revision';folder.mkdir(exist_ok=True)
    report=dict(scope='compact plates and nominal uncrimped stops screened through bare paths; not complete assembly',
        assembly_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),master_sha256=hashlib.sha256(route_file.read_bytes()).hexdigest(),
        results=results,angle_step_deg=1.,stage_pass=False,
        not_verified=['actual mounting bolt/bone interface and anti-rotation','TPU deformation, grip, friction and strength',
            'regenerated modules and combined rigid assembly'])
    (folder/'search.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Compact selected',{h:(r['selected']['reduction_mm'],r['selected']['phase_deg'],r['selected']['score']) for h,r in results.items()},flush=True)


if __name__=='__main__':main()
