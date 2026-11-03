"""Screen lower FJ1500/FJ1504 distal cap stacks without moving their bone seats.

Both same-rigid-group occupancy and cap versus fixed parts through bare-bone
paths are checked. This screen is followed by actual TPU and full assembly checks.
"""
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


def main():
    source=D.DEST/'whole_shoulder_upper_edge14/actual_assembly_check.json'
    assembly=json.loads(source.read_text(encoding='utf8'));rigid={}
    for key,row in assembly['sources'].items():
        if key.startswith('tpu:'):continue
        file=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(file.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed rigid input '+key)
        s=D.SP.load_solid(file) if file.suffix=='.obj' else C.load_stl(file)
        if row['translation_source_mm']:s=s.translate(row['translation_source_mm'])
        rigid[key]=s
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    moving_keys={'bone:humerus'}|{'cap:'+s['head']+'/'+s['end'] for s in seats if s['nearest_bone']=='humerus'}
    moving_keys|={k for k in rigid if k.startswith('stop:')}
    fixed={k:v for k,v in rigid.items() if k not in moving_keys}
    gh=D.frames.joint_centres('source')['GH'];axes=D.frames.gh_probe_axes()
    bare=[]
    for axis in ['flexion','abduction','longitudinal']:
        for sign in [1,-1]:
            angles=[]
            for magnitude in range(121):
                angle=magnitude*sign;posed=rigid['bone:humerus'].transform(I.rotation_about(gh,axes[axis],angle))
                if any(I.intersection_volume(posed,rigid['bone:'+b])>.001 for b in ['scapula','clavicle']):break
                angles.append(angle)
            bare.append(dict(axis=axis,sign=sign,angles_deg=angles,last_sample_clear_deg=angles[-1]))
    results={}
    for head in ['FJ1500','FJ1504']:
        key='cap:'+head+'/distal';original=rigid[key]
        seat=next(s for s in seats if s['head']==head and s['end']=='distal')
        n=D.unit(seat['normal_source']);maximum=3.+float(seat['lift_mm'])
        trials=[]
        for reduction in np.r_[np.arange(0,maximum,.5),maximum]:
            cap=original.translate(-n*reduction)
            same={k:I.intersection_volume(cap,v) for k,v in rigid.items() if k in moving_keys and k!=key}
            trial=dict(reduction_mm=float(reduction),same_group_intersections_mm3=same,paths=[],score=0.)
            trials.append(trial)
            if any(v>.001 for v in same.values()):
                print(head,'reduce',reduction,'same-group hit',{k:v for k,v in same.items() if v>.001},flush=True);continue
            scores=[]
            for path in bare:
                last=0;first=None;hits=[]
                for angle in path['angles_deg']:
                    posed=cap.transform(I.rotation_about(gh,axes[path['axis']],angle))
                    hits=[dict(fixed=k,volume_mm3=I.intersection_volume(posed,v)) for k,v in fixed.items()]
                    hits=[h for h in hits if h['volume_mm3']>.001]
                    if hits:first=angle;break
                    last=angle
                trial['paths'].append(dict(axis=path['axis'],sign=path['sign'],last_sample_clear_deg=last,first_hit_deg=first,hits=hits))
                scores.append(abs(last)/max(1,abs(path['last_sample_clear_deg'])))
            trial['score']=min(scores)
            print(head,'reduce',reduction,'score',trial['score'],'limits',[p['last_sample_clear_deg'] for p in trial['paths']],flush=True)
        selected=max(trials,key=lambda t:(t['score'],-t['reduction_mm']))
        results[head]=dict(max_reduction_mm=maximum,trials=trials,selected=selected)
    folder=D.DEST/'distal_stack_revision';folder.mkdir(exist_ok=True)
    report=dict(scope='two moving distal cap stack screens, not full assembly or drive acceptance',
        assembly_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),bare_paths=bare,results=results,
        stack_step_mm=.5,angle_step_deg=1.,intersection_threshold_mm3=.001,stage_pass=False)
    (folder/'search.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Selected',{h:r['selected']['reduction_mm'] for h,r in results.items()},flush=True)


if __name__=='__main__':main()
