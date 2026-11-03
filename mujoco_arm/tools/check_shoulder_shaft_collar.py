"""Screen shaft collar motion and integration, not a complete muscle assembly.

All collar bodies are ideally attached to the humerus for this screen. The
clamp friction/preload is unverified. The old muscle caps are not removed to
manufacture a whole-shoulder motion pass; their baseline limits remain separate.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_shoulder_shaft_collar as B
import audit_shoulder_retained_joint_fit as A


def load_sources(data):
    return {k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0]
            for k,r in data['sources'].items()}


def check(folder):
    source=folder/'interface_check.json'
    collar=json.loads(source.read_text(encoding='utf8'))
    assembly_file=A.ROOT/collar['assembly_source']
    if A.hashed(assembly_file)!=collar['assembly_sha256']:
        raise ValueError('Changed parent assembly')
    assembly=json.loads(assembly_file.read_text(encoding='utf8'))
    baseline_file=assembly_file.parent/'rigid_motion_check.json'
    baseline=json.loads(baseline_file.read_text(encoding='utf8'))
    if baseline['assembly_sha256']!=A.hashed(assembly_file):raise ValueError('Stale baseline motion')
    current=load_sources(assembly);new=load_sources(collar)
    fixed={k:current[k] for k in baseline['fixed_keys']}
    centre=A.frames.joint_centres('source')['GH'];axes=A.frames.gh_probe_axes()
    zero=[]
    for key,part in new.items():
        for old,solid in current.items():
            volume=A.I.intersection_volume(part,solid)
            if volume>1e-9:zero.append(dict(collar=key,reference=old,volume_mm3=volume,
                category='compliant contact requires review' if old.startswith('tpu:') or key.startswith('liner_')
                         else 'same-body or fixed rigid fit requires review'))
    def evaluate(axis,angle):
        transform=A.I.rotation_about(centre,axes[axis],angle);hits=[]
        for key,part in new.items():
            posed=part.transform(transform)
            for old,solid in fixed.items():
                volume=A.I.intersection_volume(posed,solid)
                if volume>.001:hits.append(dict(collar=key,fixed=old,volume_mm3=volume))
        return sorted(hits,key=lambda x:-x['volume_mm3'])
    results=[]
    for axis in ['flexion','abduction','longitudinal']:
        for sign in [1,-1]:
            samples=[];safe=0.;first=None
            for magnitude in range(121):
                angle=float(sign*magnitude);hits=evaluate(axis,angle)
                samples.append(dict(angle_deg=angle,hits=hits))
                if hits:first=angle;break
                safe=angle
            fine=[]
            if first is not None and first!=0:
                for step in range(1,11):
                    angle=round(safe+sign*.1,1);hits=evaluate(axis,angle)
                    fine.append(dict(angle_deg=angle,hits=hits))
                    if hits:first=angle;break
                    safe=angle
            nearest=None
            if first!=0:
                transform=A.I.rotation_about(centre,axes[axis],safe)
                for key,part in new.items():
                    posed=part.transform(transform)
                    for old,solid in fixed.items():
                        gap=A.I.gap(posed,solid,2.)
                        if nearest is None or gap<nearest['gap_capped_2mm']:
                            nearest=dict(collar=key,fixed=old,gap_capped_2mm=gap)
            results.append(dict(axis=axis,sign=sign,last_sample_clear_deg=None if first==0 else safe,
                first_hit_deg=first,last_sample_gap=nearest,coarse_samples=samples,boundary_samples=fine))
            print(folder.name,axis,sign,'component only',results[-1]['last_sample_clear_deg'],
                  'first',first,'hits',(fine[-1]['hits'] if fine else samples[-1]['hits'])[:2],flush=True)
    link_rows=[]
    for head,anchor in collar['anchors'].items():
        probe=B.rod(np.array(anchor['proposed_muscle_terminal_source_mm']),np.array(anchor['eye_centre_source_mm']),.3)
        hits=[dict(part=k,volume_mm3=A.I.intersection_volume(probe,s)) for k,s in current.items()]
        link_rows.append(dict(head=head,scope='reference chord only, NOT designed tendon',
            hits_over_threshold=[h for h in hits if h['volume_mm3']>.001]))
    report=dict(scope=__doc__,interface_sha256=A.hashed(source),assembly_sha256=A.hashed(assembly_file),
        script_sha256=A.hashed(Path(__file__).resolve()),baseline_motion_sha256=A.hashed(baseline_file),
        station_mm=collar['station_from_GH_along_anatomical_axis_mm'],
        nominal_assembly_intersections=zero,reference_link_checks=link_rows,
        collar_component_motion=results,angle_step_deg=1.,boundary_step_deg=.1,
        intersection_threshold_mm3=.001,gap_note='surface gap capped at 2mm; threshold-clear may have zero gap',
        baseline_full_rigid_limits=[{k:v for k,v in r.items() if k in
            ['axis','sign','last_sample_clear_deg','first_hit_deg']} for r in baseline['results']],
        not_verified=['old distal caps and muscle terminals have NOT been replaced or connected',
            'clamp grip/preload, TPU liner deformation, true tendon routing and terminations',
            'whole modified shoulder motion, coupled paths and physical drive range'],stage_pass=False)
    (folder/'motion_integration_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    # Source-only view of the exact new interface and its bone/joint context.
    # Omitted reference muscles/caps are stated in the scope, not accepted as
    # a complete replacement assembly or a collision-free full shoulder.
    view={k:v for k,v in assembly['sources'].items() if k.startswith(('bone:','GH_'))}
    view.update({'collar:'+k:v for k,v in collar['sources'].items()})
    (folder/'interface_view.json').write_text(json.dumps(dict(
        scope='ISOLATED COLLAR AND BONE/JOINT CONTEXT; muscle/reference caps/fixture omitted; no tendons',
        sources=view,stage_pass=False),indent=2)+'\n',encoding='utf8',newline='\n')
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stations',type=float,nargs='+',default=[40.,50.,60.])
    args=parser.parse_args()
    for station in args.stations:check(B.DEST/('station_'+format(station,'g')))
