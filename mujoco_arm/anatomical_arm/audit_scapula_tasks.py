"""Check additional scapular tasks without inventing extra joints or sliding bones."""
from build import *
import mujoco
from scipy.spatial.transform import Rotation

def run():
    m=mujoco.MjModel.from_xml_path(str(OUT/'arm.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);sid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'scapula');p0=d.xpos[sid].copy()*1000;R0=d.xmat[sid].reshape(3,3).copy();results=[]
    # Surface-normal axis from original scapular PCA, sign towards anterior +Y.
    v,f=read_stl(str(OUT/'assets/scapula.stl'));_,_,vt=np.linalg.svd(v-v.mean(0),full_matrices=False);normal=vt[-1];normal*=np.sign(normal[1])
    def fk(q):
        d.qpos[:]=0;d.qpos[:3]=q;mujoco.mj_forward(m,d);return d.xpos[sid].copy()*1000,d.xmat[sid].reshape(3,3).copy()
    for aid,kind,target in [('SLIDE_FORWARD','slide',15),('SLIDE_BACK','slide',-15),('SCAP_INTERNAL','rotate',20),('SCAP_EXTERNAL','rotate',-20)]:
        samples=[];q=np.zeros(3);count=int(abs(target));first=None
        def solve(value,q0):
            pt=p0+np.array([0,value,0]) if kind=='slide' else p0;rt=R0 if kind=='slide' else Rotation.from_rotvec(normal*np.radians(value)).as_matrix()@R0
            def error(x):
                p,r=fk(x);return np.r_[p-pt,30*Rotation.from_matrix(r@rt.T).as_rotvec()]
            sol=least_squares(error,q0,bounds=(np.radians([-30,-30,-20]),np.radians([60,30,60])),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=100);p,r=fk(sol.x);pe=np.linalg.norm(p-pt);re=np.degrees(Rotation.from_matrix(r@rt.T).magnitude());return sol.x,{'command':value,'q_sc_ac_deg':np.degrees(sol.x).tolist(),'position_residual_mm':float(pe),'orientation_residual_deg':float(re),'actual_translation_mm':(p-p0).tolist(),'pose_tracking_pass':bool(pe<=.1 and re<=.1)}
        for i in range(count+1):
            q,row=solve(target*i/count,q);samples.append(row)
            if not row['pose_tracking_pass'] and first is None:first=i
        prefix=abs(target) if first is None else max(0,first-1)
        if first and first>0:
            lo=first-1;hi=first
            while hi-lo>.1:
                mid=(lo+hi)/2;_,r=solve(np.sign(target)*mid,np.zeros(3))
                if r['pose_tracking_pass']:lo=mid
                else:hi=mid
            prefix=lo
        results.append({'action':aid,'target':target,'units':'mm' if kind=='slide' else 'deg','continuous_tracking_prefix':prefix,'samples':samples,'return_replayed':[solve(row['command'],np.radians(row['q_sc_ac_deg']))[1] for row in reversed(samples)],'collision_acceptance':False,'reason':'pose tracking only; paths failing position/orientation tolerance are not counted as realized movement. Bone collision audit not transferred from other actions.'})
    dump(OUT/'scapula_task_audit.json',{'tolerance_position_mm':.1,'tolerance_orientation_deg':.1,'boundary_resolution':.1,'rotation_axis_scapular_PCA_normal':normal.tolist(),'warning':'prescribed scapula pose tests, not independent SC protraction angles; no extra joint added','actions':results});print([(r['action'],r['continuous_tracking_prefix']) for r in results])

if __name__=='__main__':run()
