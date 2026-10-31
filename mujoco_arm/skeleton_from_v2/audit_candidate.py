"""Sample connected zero-to-target ranges; stop at first true intersection."""
import itertools,json
import numpy as np
from scipy.optimize import least_squares,minimize
from scipy.spatial.transform import Rotation
from build_candidate import build,OUT,world,NAMES,mujoco

def run():
    m,d,solids,bodies,_=build()
    joints={mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,j):j for j in range(m.njnt)}
    def pose(q):
        d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d)
    hb=bodies['humerus_shaft'];ub=bodies['ulna_shaft'];sb=bodies['scapula_mesh']
    pose(np.zeros(7));v0=d.xpos[ub]-d.xpos[hb];v0/=np.linalg.norm(v0)
    # Fixed normal of the nominal scapular plate, independent of a fitted sum.
    sup=np.array([72.76,-82.21,268.45]);inf=np.array([44.55,-106.69,112.87]);gh=np.array([159.06,-7.3,187.54])
    w=np.cross(inf-sup,gh-sup);w/=np.linalg.norm(w)
    if w[1]>0:w=-w
    def metrics():
        v=d.xpos[ub]-d.xpos[hb];v/=np.linalg.norm(v)
        r=d.xmat[sb].reshape(3,3);quat=Rotation.from_matrix(r).as_quat()
        twist=np.degrees(2*np.arctan2(quat[:3]@w,quat[3]))
        return {'world_upper_arm_elevation_deg':float(np.degrees(np.arccos(np.clip(v0@v,-1,1)))),'scapula_up_twist_deg':float(twist)}
    # Endpoint existence test only. Geometry is assessed separately below.
    bounds=(np.array([-10,-25,-35]),np.array([55,25,35]))
    desired=Rotation.from_rotvec(w*np.radians(50)).as_matrix()
    def residual(q):
        pose(np.r_[q,np.zeros(4)]);return Rotation.from_matrix(d.xmat[sb].reshape(3,3)@desired.T).as_rotvec()
    def up(q):
        pose(np.r_[q,np.zeros(4)]);return metrics()['scapula_up_twist_deg']
    fit=minimize(lambda q:np.sum((q/50)**2),[50,-20,-20],bounds=list(zip(*bounds)),constraints=[{'type':'eq','fun':lambda q:up(q)-50}],method='SLSQP',options={'ftol':1e-10,'maxiter':300})
    scap_q=fit.x;pose(np.r_[scap_q,np.zeros(4)])
    scap_endpoint={'q_deg':scap_q.tolist(),'fit_success':bool(fit.success),**metrics()}
    def raise_res(q):
        pose(np.r_[scap_q,q,0]);return np.r_[(metrics()['world_upper_arm_elevation_deg']-170)/10,q*.0001]
    raise_fit=least_squares(raise_res,[90,70,0],bounds=([-120,-120,-90],[120,120,70]),max_nfev=500)
    actions=[]
    for name,joint,end,target,floor in [('K1','elbow',145,145,116),('K2','gh_flex',120,120,96),('K3','gh_abd',120,120,96),('K4','gh_rot',70,70,56),('K5','gh_rot',-90,90,72)]:
        q=np.zeros(7);q[joints[joint]]=end;actions.append((name,q,target,floor))
    actions.extend([('K6',np.r_[scap_q,np.zeros(4)],50,36),('K7',np.r_[scap_q,raise_fit.x,0],170,132)])
    old=json.loads((OUT.parents[1]/'cable_design_v2/action_cable_matrix.json').read_text())['actions']
    for a in old:
        if a['id'].startswith('K'):continue
        q=np.zeros(7)
        for j,v in a['q_target_deg'].items():q[joints[j]]=v
        actions.append((a['id'],q,a.get('target_deg',a.get('target_mm')),None))
    def collision(q):
        pose(q);worlds={n:world(d,s,bodies[n]) for n,s in solids.items()};hits=[];gap=2.
        for a,b in itertools.combinations(NAMES,2):
            aa=worlds[a];bb=worlds[b];ba=np.array(aa.bounding_box());bc=np.array(bb.bounding_box())
            if np.all(ba[:3]<bc[3:]) and np.all(bc[:3]<ba[3:]):
                hit=aa^bb;v=hit.volume()
                if v>.001:hits.append({'bones':[a,b],'volume_mm3':v})
            gap=min(gap,aa.min_gap(bb,2.))
        return hits,gap
    rows=[]
    for name,q,target,floor in actions:
        steps=max(1,int(np.ceil(np.max(np.abs(q)))));last=0.;first=None;mingap=2.;tested=0
        for f in np.linspace(0,1,steps+1):
            hits,gap=collision(q*f);tested+=1
            if hits:first=(float(f),hits);break
            last=float(f);mingap=min(mingap,gap)
        if first:
            lo=last;hi=first[0]
            while (hi-lo)*max(np.max(np.abs(q)),1)>.1:
                mid=(lo+hi)/2;hits,gap=collision(q*mid);tested+=1
                if hits:hi=mid;first=(mid,hits)
                else:lo=mid;mingap=min(mingap,gap)
            last=lo
        pose(q*last);actual=metrics()
        value=(abs(q[joints['elbow']]*last) if name=='K1' else abs(q[joints['gh_flex']]*last) if name=='K2' else abs(q[joints['gh_abd']]*last) if name=='K3' else abs(q[joints['gh_rot']]*last) if name in ['K4','K5'] else abs(actual['scapula_up_twist_deg']) if name=='K6' else actual['world_upper_arm_elevation_deg'] if name=='K7' else None)
        row={'action':name,'target':target,'floor':floor,'q_endpoint_deg':q.tolist(),'last_clear_fraction':last,'actual_range_deg':value,'actual_pose':actual,'first_collision':first,'minimum_tested_clear_gap_mm_capped_at_2':mingap,'tested_poses':tested,'path_reaches_endpoint':first is None}
        rows.append(row);print(json.dumps(row),flush=True)
    report={'frozen':False,'scope':'candidate bone solids only; 1 degree max joint increment, first boundary refined to 0.1 degree; no continuum proof, no force/TPU validation','scapula_endpoint_fit':scap_endpoint,'actions':rows}
    (OUT/'candidate_rom.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':run()
