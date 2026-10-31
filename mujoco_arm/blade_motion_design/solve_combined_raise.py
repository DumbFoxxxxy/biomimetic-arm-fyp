"""Search the full shoulder/GH redundancy for actual 170-degree arm raising.

Keeps the original 50-degree scapular target and GH component bounds; does not
pretend their Euler angles add arithmetically. Endpoint only, not a ROM proof.
"""
import json,sys,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import minimize
from audit_arm_rom import build_arm,set_arm,world_direction,mesh_samples
from audit_geometry import mesh_local,capsule_local
from shoulder_model import inputs,EXPORT,v8
from search_gh_anchors import A,B,C
sys.stdout.reconfigure(encoding='utf-8')

def solve():
    cfg,fit,anchors=inputs();model=mujoco.MjModel.from_xml_string(ET.tostring(build_arm().getroot(),encoding='unicode'));data=mujoco.MjData(model)
    zero=np.array(fit['poses'][0]['q_rad']);up=np.array(next(p['q_rad'] for p in fit['poses'] if p['id']=='UP_ROT_1.00'))
    set_arm(model,data,zero);neutral=world_direction(model,data)
    axis=v8.direct_kinematics(cfg,zero)['body_rotation']@C
    desired=v8.rotation_matrix(axis,np.radians(170.))@neutral
    perp1=np.cross(desired,axis);perp1/=np.linalg.norm(perp1);perp2=np.cross(desired,perp1)
    body_samples=[]
    for name in ['scapula_mesh','clavicle_mesh','humerus_shaft','ulna_shaft']:
        gid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,name);p,t=mesh_local(model,data,gid);body_samples.append((int(model.geom_bodyid[gid]),p))
    qref=v8.neutral_q(cfg);fkref=v8.direct_kinematics(cfg,qref);pref=fkref['contact_mm'];centre=np.array(cfg['surface']['center_xz_mm'])
    theta_ref=np.arctan2(pref[2]-centre[1],pref[0]-centre[0]);uref=fkref['body_rotation']@v8.source_basis(cfg)[0]
    def state(z):
        set_arm(model,data,z[:5],np.degrees([z[5],z[6],0.]))
        direction=world_direction(model,data);fk=v8.direct_kinematics(cfg,z[:5]);p=fk['contact_mm'];theta=np.arctan2(p[2]-centre[1],p[0]-centre[0])
        normal=np.array([np.cos(theta),0.,np.sin(theta)]);ur=v8.rotation_matrix(np.array([0.,1.,0.]),-(theta-theta_ref))@uref
        ur=v8.normalize(ur-normal*np.dot(ur,normal));ua=fk['body_rotation']@v8.source_basis(cfg)[0];ua=v8.normalize(ua-normal*np.dot(ua,normal))
        phi=np.arctan2(normal@np.cross(ur,ua),ur@ua)
        gaps=[]
        for body,points in body_samples:
            world=points@data.xmat[body].reshape(3,3).T+data.xpos[body]
            gaps.append(float(v8.thorax_clearance(world@v8.DESIGN_TO_MJ*1000,cfg).min()))
        return direction,phi,np.array(gaps),np.linalg.norm(p[[0,2]]-centre)-400.
    def eq(z):
        d,phi,gap,offset=state(z);return np.array([(phi-np.radians(50.))*30,d@perp1*30,d@perp2*30])
    def inequality(z):
        d,phi,gap,offset=state(z);return np.r_[gap-.6,offset-2.,15.-offset,d@desired-.99,np.radians(80.)**2-np.sum(z[2:5]**2)]
    lower=np.radians([-10,-25,-65,-65,-65,-40,-20]);upper=np.radians([55,25,65,65,65,120,120])
    baseline=np.r_[up,np.radians([-34.,120.])];best=None
    for seed in [baseline,np.r_[up,np.radians([0.,120.])],np.r_[up+np.radians([5,0,0,0,0]),np.radians([-20.,120.])]]:
        result=minimize(lambda z:np.sum((z-baseline)**2),np.clip(seed,lower,upper),method='SLSQP',bounds=list(zip(lower,upper)),
            constraints=[{'type':'eq','fun':eq},{'type':'ineq','fun':inequality}],options={'maxiter':200,'ftol':1e-9})
        direction,phi,gaps,offset=state(result.x);angle=float(np.degrees(np.arccos(np.clip(direction@neutral,-1,1))))
        residual=float(np.max(np.abs(eq(result.x))));okay=bool(residual<1e-4 and np.min(inequality(result.x))>-1e-5)
        row={'pass':okay,'optimizer_success':bool(result.success),'message':str(result.message),'q_deg':np.degrees(result.x).tolist(),
            'world_raise_deg':angle,'scapular_up_rotation_deg':float(np.degrees(phi)),'body_chest_clearances_mm':gaps.tolist(),
            'marker_gap_mm':float(offset),'equality_residual':residual,'path_verified':False,'cable_feasibility_verified':False}
        print(json.dumps(row),flush=True)
        if best is None or (okay,-residual)>(best['pass'],-best['equality_residual']):best=row
        if okay:break
    (EXPORT/'combined_raise_endpoint_search.json').write_text(json.dumps(best,indent=2))
if __name__=='__main__':solve()
