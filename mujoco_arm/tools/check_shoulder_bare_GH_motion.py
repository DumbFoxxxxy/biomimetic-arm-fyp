"""Bare-humerus GH geometry probes, excluding all muscle/tendon attachments.

This explicitly staged test does not validate motion of the integrated
free-terminal assembly. Its ranges are neither TPU predictions nor K passes.
All excluded objects are listed; fixed joint parts are checked individually.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_GH_nut_keepers import DEST as CURRENT


def run(source,out):
    source,out=source.resolve(),out.resolve();out.mkdir(parents=True,exist_ok=True)
    data=json.loads(source.read_text(encoding='utf8'))
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in data['sources'].items()}
    fixed={k:s for k,s in S.items() if k.startswith(('fixture:','bench:','clamp:','liner_keeper:','GH_'))
           or k in ['bone:scapula','clavicle:fixed_support'] or k.startswith('clavicle:M4_')}
    moving=S['bone:humerus'];centre=A.frames.joint_centres('source')['GH'];axes=A.frames.gh_probe_axes()
    segment=A.frames.joint_centres('source')['elbow']-centre;segment/=np.linalg.norm(segment)
    fixed_bounds={k:np.array(s.bounding_box()) for k,s in fixed.items()}
    def evaluate(axis,angle):
        T=A.I.rotation_about(centre,axis,angle);posed=moving.transform(T);bounds=np.array(posed.bounding_box())
        hits=[];nearest=dict(part=None,gap_capped_2mm=2.)
        for k,obstacle in fixed.items():
            volume=A.I.intersection_volume(posed,obstacle)
            if volume>.001:hits.append(dict(part=k,volume_mm3=volume))
            bb=fixed_bounds[k];delta=np.maximum(np.maximum(bounds[:3]-bb[3:],bb[:3]-bounds[3:]),0.)
            if np.linalg.norm(delta)<nearest['gap_capped_2mm']:
                gap=A.I.gap(posed,obstacle,nearest['gap_capped_2mm'])
                if gap<nearest['gap_capped_2mm']:nearest=dict(part=k,gap_capped_2mm=gap)
        direction=T[:,:3]@segment
        return dict(angle_deg=round(float(angle),6),humerus_direction_source=direction.tolist(),
                    angle_from_neutral_bone_direction_deg=float(np.degrees(np.arccos(np.clip(segment@direction,-1.,1.)))),
                    hits_over_threshold=sorted(hits,key=lambda r:-r['volume_mm3']),minimum_gap=nearest)
    results=[]
    for label,axis in axes.items():
        for sign in [1,-1]:
            coarse=[];fine=[];clear=[];first=None
            for magnitude in range(121):
                row=evaluate(axis,float(sign*magnitude));coarse.append(row)
                if row['hits_over_threshold']:first=row;break
                clear.append(row)
            if first is not None and clear:
                base=clear[-1]['angle_deg']
                for offset in range(1,11):
                    row=evaluate(axis,round(base+sign*offset*.1,1));fine.append(row)
                    if row['hits_over_threshold']:first=row;break
                    clear.append(row)
            last=clear[-1] if clear else None
            back=[evaluate(axis,r['angle_deg']) for r in reversed(clear)]
            minimum=min((r['minimum_gap']['gap_capped_2mm'] for r in clear),default=None)
            result=dict(axis=label,sign=sign,requested_limit_deg=120.,last_contiguous_sample_clear_deg=None if last is None else last['angle_deg'],
                        first_interference=first,last_clear_pose=last,minimum_gap_along_outward_samples_mm=minimum,
                        outward_coarse_samples=coarse,boundary_samples=fine,return_samples=back,
                        return_sampled_path_clear=not any(r['hits_over_threshold'] for r in back))
            results.append(result)
            print(label,sign,'last clear',None if last is None else last['angle_deg'],
                  'first hits',None if first is None else first['hits_over_threshold'][:2],flush=True)
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),
                geometry_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/audit_shoulder_retained_joint_fit.py'),
                frames_sha256=A.hashed(Path(A.frames.__file__).resolve()),
                source_object_count=len(S),moving_keys=['bone:humerus'],fixed_keys=sorted(fixed),
                not_installed_keys=sorted(set(S)-set(fixed)-{'bone:humerus'}),
                joint_centre_source_mm=centre.tolist(),axes_source={k:v.tolist() for k,v in axes.items()},
                results=results,step_deg=1.,boundary_step_deg=.1,intersection_threshold_mm3=.001,gap_cap_mm=2.,
                anatomical_axial_sign_mapping_verified=False,stage_pass=False,actual_drive_range_deg=None,
                assumptions=['GH ideal centre fixed; real fastening and preload not solved',
                             'only bare humerus moves; all muscle ends, collar and tendon modules absent',
                             'axial sign is not mapped to anatomical internal/external rotation'],
                not_verified=['integrated or driven shoulder range, flexible contact and tendons',
                              'combined/scapular motion, loads, strength, print fit and actual performance'])
    (out/'bare_joint_motion_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');parser=argparse.ArgumentParser()
    parser.add_argument('--assembly',type=Path,default=CURRENT/'assembly_check.json')
    parser.add_argument('--out',type=Path,default=CURRENT/'bare_motion_review')
    args=parser.parse_args();run(args.assembly,args.out)
