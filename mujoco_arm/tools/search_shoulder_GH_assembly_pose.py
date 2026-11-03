"""Find a bare-humerus assembly pose for GH caps; check return to neutral.

Only humerus is rotated about the recorded GH centre. No muscles, collar or
tendon constraints are in this stage; this is not an actuation range result.
"""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from open_shoulder_GH_assembly_entries import DEST
from check_shoulder_front_liner_keepers import dense_slide
from audit_shoulder_GH_assembly import scan


def run():
    source=DEST/'assembly_check.json';data=json.loads(source.read_text(encoding='utf8'))
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in data['sources'].items()}
    joint=json.loads((A.G1/'assembly.json').read_text(encoding='utf8'));G=np.array(joint['centre_source_mm']);n=np.array(joint['neck_axis'])
    bare={k for k in S if k.startswith(('fixture:','bench:','clamp:','liner_keeper:'))
          or k in ['bone:scapula','clavicle:fixed_support'] or k.startswith('clavicle:M4_')}
    gh={k for k in S if k.startswith('GH_')};fixed_all={k:S[k] for k in bare|gh}
    cases=[];selected=None
    for axis_label in ['flexion','abduction']:
        axis=A.frames.gh_probe_axes()[axis_label]
        for angle in [10.,20.,30.,45.,-10.,-20.,-30.]:
            hum=S['bone:humerus'].transform(A.I.rotation_about(G,axis,angle))
            # First check cap 1 including all earlier nominal nuts/cap 0.
            fixed={k:S[k] for k in bare|(gh-{'GH_cap_1','GH_top_washer_1','GH_bolt_0','GH_bolt_1'})}
            fixed['bone:humerus']=hum
            body=S['GH_cap_1']+S['GH_top_washer_1']
            quick=scan(body,n,fixed,25.,.1)
            row=dict(axis=axis_label,axis_source=axis.tolist(),angle_deg=angle,
                     cap1_straight_path_clear=quick['entire_sampled_path_clear'])
            if quick['entire_sampled_path_clear']:
                motion=[]
                for phase,values in [('to_assembly_pose',np.linspace(0.,angle,int(abs(angle))+1)),
                                     ('return_to_neutral',np.linspace(angle,0.,int(abs(angle))+1))]:
                    samples=[]
                    for value in values:
                        h=S['bone:humerus'].transform(A.I.rotation_about(G,axis,float(value)))
                        hits=[dict(part=k,volume_mm3=A.I.intersection_volume(h,s)) for k,s in fixed_all.items()]
                        hits=[r for r in hits if r['volume_mm3']>.001]
                        samples.append(dict(angle_deg=float(value),hits_over_threshold=hits))
                        if hits:break
                    motion.append(dict(phase=phase,samples=samples,
                                       sampled_path_clear=len(samples)==len(values) and not any(s['hits_over_threshold'] for s in samples)))
                row['closed_cap_humerus_motion']=motion
                if all(p['sampled_path_clear'] for p in motion):
                    dense=dense_slide({k:S[k] for k in ['GH_cap_1','GH_top_washer_1']},fixed,n,25.)
                    if dense['sampled_path_clear']:
                        row['cap1_dense_path']=dense;selected=row;cases.append(row);break
            cases.append(row)
            print(axis_label,angle,'cap clear',quick['entire_sampled_path_clear'],flush=True)
        if selected:break
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),cases=cases,selected=selected,
                closed_motion_fixed_keys=sorted(fixed_all),
                not_yet_installed_keys=sorted(set(S)-set(fixed_all)-{'bone:humerus'}),
                angle_step_deg=1.,translation_step_mm=.1,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'humerus_assembly_pose_search.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Selected',None if selected is None else (selected['axis'],selected['angle_deg']),flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
