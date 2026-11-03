"""Compare shorter FJ1467 cap stacks against contiguous bare-bone paths.

This is cap-only design screening. Full rigid assembly, TPU and steel must
be regenerated and checked before a selected candidate can be adopted.
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
    source=D.DEST/'whole_shoulder_fixture_revision/actual_assembly_check.json'
    assembly=json.loads(source.read_text(encoding='utf8'));rigid={}
    for key,row in assembly['sources'].items():
        if key.startswith('tpu:'):continue
        path=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed input '+key)
        s=D.SP.load_solid(path) if path.suffix=='.obj' else C.load_stl(path)
        if row['translation_source_mm']:s=s.translate(row['translation_source_mm'])
        rigid[key]=s
    cap_key='cap:FJ1467/origin';original=rigid.pop(cap_key)
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    moving_keys={'bone:humerus'}|{'cap:'+s['head']+'/'+s['end'] for s in seats if s['nearest_bone']=='humerus'}
    moving_keys|={k for k in rigid if k.startswith('stop:')}
    moving={k:v for k,v in rigid.items() if k in moving_keys}
    fixed={k:v for k,v in rigid.items() if k not in moving_keys}
    gh=D.frames.joint_centres('source')['GH'];axes=D.frames.gh_probe_axes()
    bare=[]
    for axis in ['flexion','abduction','longitudinal']:
        for sign in [1,-1]:
            angles=[]
            for a in range(121):
                angle=sign*a;h=rigid['bone:humerus'].transform(I.rotation_about(gh,axes[axis],angle))
                if any(I.intersection_volume(h,rigid['bone:'+b])>.001 for b in ['scapula','clavicle']):break
                angles.append(angle)
            bare.append(dict(axis=axis,sign=sign,angles_deg=angles,last_sample_clear_deg=angles[-1]))
    route=next(r for r in json.loads((D.DEST/'routes_outlet_candidate.json').read_text(encoding='utf8'))['routes'] if r['head']=='FJ1467')
    n=D.unit(route['ports']['origin']['normal_source']);trials=[]
    for reduction in np.arange(0,3.01,.5):
        displacement=-n*reduction;cap=original.translate(displacement)
        hits={k:I.intersection_volume(cap,v) for k,v in fixed.items()}
        trial=dict(boss_reduction_mm=float(reduction),displacement_source_mm=displacement.tolist(),
            fixed_rigid_intersections_mm3=hits,cap_only_paths=[],score=0.)
        trials.append(trial)
        if any(v>.001 for v in hits.values()):
            print('reduction',reduction,'fixed hits',{k:v for k,v in hits.items() if v>.001},flush=True);continue
        scores=[]
        for path in bare:
            last=0;first=None;limiter=[]
            for angle in path['angles_deg']:
                transform=I.rotation_about(gh,axes[path['axis']],angle)
                limiter=[dict(moving=k,volume_mm3=I.intersection_volume(cap,s.transform(transform))) for k,s in moving.items()]
                limiter=[h for h in limiter if h['volume_mm3']>.001]
                if limiter:first=angle;break
                last=angle
            trial['cap_only_paths'].append(dict(axis=path['axis'],sign=path['sign'],last_sample_clear_deg=last,
                first_hit_deg=first,limiter=limiter))
            scores.append(abs(last)/max(1,abs(path['last_sample_clear_deg'])))
        trial['score']=min(scores)
        print('reduction',reduction,'cap-only limits',[(p['axis'],p['sign'],p['last_sample_clear_deg']) for p in trial['cap_only_paths']],flush=True)
    selected=max(trials,key=lambda t:(t['score'],-t['boss_reduction_mm']))
    report=dict(scope='cap-only stack screening, not full rigid or drive acceptance',
        source_assembly_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),bare_paths=bare,trials=trials,
        selected=selected,coarse_step_deg=1.,stack_search_step_mm=.5,intersection_threshold_mm3=.001,
        max_boss_reduction_mm=3.,stage_pass=False,
        missing=['actual shortened pad and remapped steel/TPU geometry','complete rigid assembly and boundary refinement',
                 'mounting, material deformation, retention and physical tests'])
    folder=D.DEST/'origin_stack_revision';folder.mkdir(exist_ok=True)
    (folder/'search.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Selected cap-only reduction',selected['boss_reduction_mm'],'score',selected['score'],flush=True)


if __name__=='__main__':main()
