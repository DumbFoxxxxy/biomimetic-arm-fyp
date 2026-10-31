"""Fit the complete 170-degree raising path, retaining the 50-degree scapula target."""
import json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import minimize
from gimbal_arm_model import build_gimbal
from audit_arm_rom import set_arm,world_direction
from audit_geometry import mesh_local,capsule_local
from shoulder_model import inputs,EXPORT,v8
from search_gh_anchors import C

def solve():
    cfg,fit,_=inputs();tree,_=build_gimbal();m=mujoco.MjModel.from_xml_string(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m)
    zero=np.array(fit['poses'][0]['q_rad']);old=np.radians(json.loads((EXPORT/'combined_raise_endpoint_search.json').read_text())['q_deg'])
    set_arm(m,d,zero);neutral=world_direction(m,d);axis=v8.direct_kinematics(cfg,zero)['body_rotation']@C
    desired=v8.rotation_matrix(axis,np.radians(170.))@neutral;p1=np.cross(desired,axis);p1/=np.linalg.norm(p1);p2=np.cross(desired,p1)
    samples=[]
    for name in ['scapula_mesh','clavicle_mesh','eye_web_S03','ac_socket','humerus_shaft','ulna_shaft']:
        g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,name)
        p=mesh_local(m,d,g)[0] if m.geom_type[g]==mujoco.mjtGeom.mjGEOM_MESH else capsule_local(m,d,g)
        samples.append((name,int(m.geom_bodyid[g]),p[::3]))
    fk0=v8.direct_kinematics(cfg,zero);centre=np.array(cfg['surface']['center_xz_mm']);pref=fk0['contact_mm'];theta0=np.arctan2(pref[2]-centre[1],pref[0]-centre[0]);u0=fk0['body_rotation']@v8.source_basis(cfg)[0]
    def endpoint(z):
        set_arm(m,d,z[:5],np.degrees([z[5],z[6],0.]));fk=v8.direct_kinematics(cfg,z[:5]);p=fk['contact_mm'];theta=np.arctan2(p[2]-centre[1],p[0]-centre[0]);normal=np.array([np.cos(theta),0.,np.sin(theta)])
        ur=v8.rotation_matrix(np.array([0.,1.,0.]),-(theta-theta0))@u0;ur=v8.normalize(ur-normal*np.dot(ur,normal));ua=fk['body_rotation']@v8.source_basis(cfg)[0];ua=v8.normalize(ua-normal*np.dot(ua,normal))
        phi=np.arctan2(normal@np.cross(ur,ua),ur@ua);direction=world_direction(m,d)
        return direction,phi,np.linalg.norm(p[[0,2]]-centre)-400.
    def eq(z):
        direction,phi,gap=endpoint(z);return np.array([(phi-np.radians(50.))*30,direction@p1*30,direction@p2*30])
    def inequality(z):
        mins=np.ones(len(samples))*1e6
        for f in np.linspace(0,1,21):
            q=zero*(1-f)+z[:5]*f;set_arm(m,d,q,np.degrees([z[5]*f,z[6]*f,0.]))
            for i,(name,b,p) in enumerate(samples):
                world=p@d.xmat[b].reshape(3,3).T+d.xpos[b]
                mins[i]=min(mins[i],float(v8.thorax_clearance(world@v8.DESIGN_TO_MJ*1000,cfg).min()))
        direction,phi,gap=endpoint(z)
        return np.r_[mins-np.array([4.5,.6,4.5,.6,.6,.6]),gap-2.,35.-gap,direction@desired-.99]
    lower=np.radians([-10,-25,-65,-65,-65,-40,-20]);upper=np.radians([55,25,65,65,65,120,120])
    result=minimize(lambda z:np.sum((z-old)**2),old,method='SLSQP',bounds=list(zip(lower,upper)),constraints=[{'type':'eq','fun':eq},{'type':'ineq','fun':inequality}],options={'maxiter':180,'ftol':1e-8})
    direction,phi,gap=endpoint(result.x);valid=bool(np.max(np.abs(eq(result.x)))<1e-4 and inequality(result.x).min()>-1e-5)
    points=np.array([np.r_[zero*(1-f)+result.x[:5]*f,result.x[5:7]*f,0.,0.] for f in np.linspace(0,1,21)])
    report={'optimizer_success':bool(result.success),'pass':valid,'message':str(result.message),'endpoint_deg':np.degrees(result.x).tolist(),'world_raise_deg':float(np.degrees(np.arccos(np.clip(direction@neutral,-1,1)))),'scapular_up_rotation_deg':float(np.degrees(phi)),'marker_gap_mm':gap,'minimum_constraints':float(inequality(result.x).min()),'points_rad':points.tolist(),'actual_mesh_and_cable_audit_pass':False}
    (EXPORT/'combined_raise_path_candidate.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='points_rad'}),flush=True)
if __name__=='__main__':solve()
