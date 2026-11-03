"""Local straight-line actuation audit, not a routed TPU or joint-contact proof.

Angles are fixed-centre axis probes. Tension coefficients have no measured
upper bounds. Positive spanning proves only conditional directional authority.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0,str(Path(__file__).resolve().parent))
import frames
import interference as I

DESIGN=frames.OUTPUT_DIR.parent/'shoulder_design'


def tension_moments(origins,insertions,centre):
    """Columns: torque per unit tensile force about centre, in mm."""
    delta=np.asarray(origins)-np.asarray(insertions)
    lengths=np.linalg.norm(delta,axis=1)
    if np.any(lengths<=1e-9):raise ValueError('Coincident force endpoints')
    directions=delta/lengths[:,None]
    return np.cross(np.asarray(insertions)-centre,directions).T,lengths


def positive_spanning(matrix):
    """Find normalized strictly positive self-stress, AND require rank3.

    This is an unbounded, tension-only torque-cone test, not a force capacity
    or socket stability calculation. Translation/contact forces are omitted.
    """
    matrix=np.asarray(matrix,float);n=matrix.shape[1]
    # max t, with M f=0, sum(f)=1 and each f_i>=t.
    eq=np.zeros((4,n+1));eq[:3,:n]=matrix;eq[3,:n]=1
    ub=np.c_[-np.eye(n),np.ones(n)]
    c=np.r_[np.zeros(n),-1.]
    result=linprog(c,A_ub=ub,b_ub=np.zeros(n),A_eq=eq,b_eq=[0,0,0,1],
                   bounds=[(0,None)]*(n+1),method='highs')
    rank=int(np.linalg.matrix_rank(matrix,tol=1e-8))
    residual=None;weights=None;margin=0.
    if result.success:
        weights=result.x[:n];residual=float(np.max(np.abs(matrix@weights)))
        margin=float(min(result.x[-1],np.min(weights)))
    ok=result.success and rank==3 and margin>1e-7 and residual<1e-7
    certificate=None
    if not ok:
        # A strict separating plane is an independent, directly checkable
        # witness: all nonnegative tension combinations stay on one side.
        separate=linprog([0,0,0,-1],A_ub=np.c_[-matrix.T,np.ones(n)],
                         b_ub=np.zeros(n),bounds=[(-1,1)]*3+[(0,None)],method='highs')
        if separate.success and separate.x[3]>1e-7:
            normal=separate.x[:3]
            certificate=dict(normal=normal.tolist(),
                             all_column_projections_mm=(normal@matrix).tolist())
    return dict(rank=rank,normalized_positive_margin=margin,
                torque_residual_mm=residual,strict_positive_span=bool(ok),
                balanced_weights=weights.tolist() if weights is not None else None,
                strict_separating_plane=certificate)


def directional_probes(matrix,axes):
    rows=[]
    for name,target in axes.items():
        for sign in [-1,1]:
            rhs=sign*np.asarray(target)
            solve=linprog(np.ones(matrix.shape[1]),A_eq=matrix,b_eq=rhs,
                          bounds=(0,None),method='highs')
            residual=float(np.max(np.abs(matrix@solve.x-rhs))) if solve.success else None
            rows.append(dict(axis=name,sign=sign,unit_torque_direction_feasible=bool(
                solve.success and residual<1e-7),residual=residual))
    return rows


def analyze(step=1.,maximum=120.):
    paths=[DESIGN/'s1_modules.json',DESIGN/'s1_gh_bare_sweep.json']
    modules=json.loads(paths[0].read_text(encoding='utf8'))
    limits=json.loads(paths[1].read_text(encoding='utf8'))
    lines=[dict(module=m['module'],head=c['head'],A=c['path_source_mm'][0],
                B=c['path_source_mm'][-1]) for m in modules['modules'] for c in m['corridors']]
    A=np.array([r['A'] for r in lines]);B=np.array([r['B'] for r in lines])
    G=frames.joint_centres('source')['GH'];axes=frames.gh_probe_axes()
    bounds={(r['axis'],r['sign']):r['first_hit_deg'] for r in limits['sweeps']}
    scans=[]
    for axis,k in axes.items():
        for sign in [-1,1]:
            hit=bounds[axis,sign];samples=[]
            for q in np.arange(0,maximum+step*.1,step):
                angle=float(sign*q);T=I.rotation_about(G,k,angle)
                moved=B@T[:,:3].T+T[:,3]
                M,L=tension_moments(A,moved,G);span=positive_spanning(M)
                valid=hit is None or abs(angle)<abs(hit)
                row=dict(deg=angle,before_known_bare_bone_first_collision=valid,
                         local_signed_moment_arms_mm=(k@M).tolist(),
                         chords_mm=L.tolist(),torque_cone=span)
                if not span['strict_positive_span']:row['direction_tests']=directional_probes(M,axes)
                samples.append(row)
            bad=next((r for r in samples if not r['torque_cone']['strict_positive_span']),None)
            scans.append(dict(axis=axis,sign=sign,bare_bone_first_hit_deg=hit,samples=samples,
                              first_sample_without_full_torque_cone_deg=bad['deg'] if bad else None))
    return dict(scope='conditional fixed-centre straight-line force model; NOT S1 acceptance',
        source_sha256={str(p.relative_to(frames.PROJECT_ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        centre_source_mm=G.tolist(),axes={k:v.tolist() for k,v in axes.items()},lines=lines,
        angle_step_deg=step,requested_probe_maximum_deg=maximum,scans=scans,
        not_verified=['actual routed cable Jacobian and flexible muscle contact',
            'effective TPU length, strain and tendon compliance',
            'contact stability, liner, translations and maintained joint centre',
            'force limits, gravity, friction, strength and combined-axis paths'],
        actual_drive_range_deg=None,stage_pass=False)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--step',type=float,default=1.)
    ap.add_argument('--max',type=float,default=120.);ap.add_argument('--json',default=str(DESIGN/'s1_local_actuation.json'))
    args=ap.parse_args()
    if not 0<args.step<=1:ap.error('Use positive step <=1 degree')
    report=analyze(args.step,args.max)
    out=Path(args.json);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2),encoding='utf8')
    for s in report['scans']:
        print(s['axis'],s['sign'],'first torque-cone failure',s['first_sample_without_full_torque_cone_deg'],
              'bare bone first collision',s['bare_bone_first_hit_deg'])
    print('Conditional audit saved; stage_pass remains False.')


if __name__=='__main__':main()
