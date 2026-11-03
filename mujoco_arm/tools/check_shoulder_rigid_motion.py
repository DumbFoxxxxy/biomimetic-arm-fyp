"""Conditional fixed-GH rigid sweeps of the actual integrated candidate.

Checks bone, caps, nominal stop placement and fixture; TPU deformation and
open-joint retention are deliberately unverified, not approximated as rigid.
"""
import hashlib
import argparse
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
    ap=argparse.ArgumentParser();ap.add_argument('--assembly',type=Path,
        default=D.DEST/'whole_shoulder_fixture_revision/actual_assembly_check.json')
    ap.add_argument('--moving-prefix',action='append',default=[],
        help='Additional ideally humerus-attached component keys, e.g. collar:')
    args=ap.parse_args();source=args.assembly;folder=source.parent
    assembly=json.loads(source.read_text(encoding='utf8'))
    if assembly.get('free_terminal_constraints_required'):
        raise ValueError('Free muscle terminals require coupled muscle/tendon constraints; '
                         'this fixed rigid-group sweep cannot validate their motion')
    solids={}
    for key,row in assembly['sources'].items():
        if key.startswith('tpu:'):continue
        path=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed rigid input '+key)
        s=D.SP.load_solid(path) if path.suffix=='.obj' else C.load_stl(path)
        if row['translation_source_mm']:s=s.translate(row['translation_source_mm'])
        solids[key]=s
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    moving_keys={'bone:humerus'}|{'cap:'+s['head']+'/'+s['end'] for s in seats if s['nearest_bone']=='humerus'}
    moving_keys|={k for k in solids if k.startswith('stop:')}
    for prefix in args.moving_prefix:
        matches={k for k in solids if k.startswith(prefix)}
        if not matches:raise ValueError('Moving prefix has no rigid components: '+prefix)
        moving_keys|=matches
    moving={k:v for k,v in solids.items() if k in moving_keys}
    fixed={k:v for k,v in solids.items() if k not in moving_keys}
    centre=D.frames.joint_centres('source')['GH'];axes=D.frames.gh_probe_axes()
    threshold=.001
    def evaluate(axis,angle):
        transform=I.rotation_about(centre,axes[axis],angle)
        hits=[]
        for k,s in moving.items():
            posed=s.transform(transform)
            for name,obstacle in fixed.items():
                volume=I.intersection_volume(posed,obstacle)
                if volume>threshold:hits.append(dict(moving=k,fixed=name,volume_mm3=volume))
        return sorted(hits,key=lambda r:-r['volume_mm3'])
    def minimum_gap(axis,angle):
        transform=I.rotation_about(centre,axes[axis],angle);nearest=None
        for k,s in moving.items():
            posed=s.transform(transform)
            for name,obstacle in fixed.items():
                gap=I.gap(posed,obstacle,2.)
                if nearest is None or gap<nearest['gap_capped_2mm']:
                    nearest=dict(moving=k,fixed=name,gap_capped_2mm=gap)
        return nearest
    results=[]
    for axis in ['flexion','abduction','longitudinal']:
        # frames calls the axial probe longitudinal; verify available names
        # instead of substituting a guessed global Euler axis.
        if axis not in axes:raise ValueError('Missing calibrated axis '+axis)
        for sign in [1,-1]:
            samples=[];safe=0.;first=None
            for magnitude in range(121):
                angle=float(sign*magnitude);hits=evaluate(axis,angle)
                samples.append(dict(angle_deg=angle,hits=hits))
                if hits:
                    first=angle;break
                safe=angle
            fine=[]
            if first is not None and first!=0:
                for offset in range(1,11):
                    angle=round(safe+sign*.1,1);hits=evaluate(axis,angle)
                    fine.append(dict(angle_deg=angle,hits=hits))
                    if hits:first=angle;break
                    safe=angle
            results.append(dict(axis=axis,sign=sign,last_sample_clear_deg=None if first==0 else safe,first_hit_deg=first,
                last_sample_minimum_gap=None if first==0 else minimum_gap(axis,safe),
                endpoint_hits=(fine[-1]['hits'] if fine else samples[-1]['hits']),
                coarse_samples=samples,boundary_samples=fine))
            print(axis,sign,'last clear',safe,'first hit',first,'pairs',results[-1]['endpoint_hits'][:3],flush=True)
    report=dict(scope='actual rigid components, fixed GH single-axis contiguous paths; not shoulder motion acceptance',
        assembly_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        moving_keys=sorted(moving),fixed_keys=sorted(fixed),GH_source_mm=centre.tolist(),
        additional_moving_prefixes=args.moving_prefix,
        axes_source={k:np.asarray(v).tolist() for k,v in axes.items()},
        step_deg=1.,boundary_step_deg=.1,intersection_threshold_mm3=threshold,
        gap_note='surface gap capped at2mm, zero when surfaces contact/intersect; report volume threshold separately',
        results=results,stage_pass=False,
        assumptions=['scapula and clavicle fixed; nominal GH centre enforced',
            'humerus-side caps and stop blanks move with humerus as an ideal installed rigid group',
            'additional moving prefixes, if present, are ideally clamped; liner preload/deformation is not solved',
            'within-group contacts retained in zero-pose report, not independently moving pairs'],
        not_verified=['joint retention/translation and real cap fastening','TPU deformation or cable drive',
            'combined GH/SC/AC motion, K1-K7 acceptance','physical hysteresis and return motion'])
    (folder/'rigid_motion_check.json').write_text(json.dumps(report,indent=2),encoding='utf8')


if __name__=='__main__':main()
