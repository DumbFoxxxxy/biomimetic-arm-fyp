"""Test all 16 spatial tendons of the full 9-DOF no-payload arm."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,sys,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation
from gimbal_arm_model import build_gimbal
from arm_targets import target_actions
from audit_arm_rom import set_arm
from shoulder_model import inputs,EXPORT,v8
sys.stdout.reconfigure(encoding='utf-8')

def native_jacobian(model,data):
    J=np.zeros((model.ntendon,model.nv))
    for i in range(model.ntendon):
        start=int(model.ten_J_rowadr[i]);count=int(model.ten_J_rownnz[i]);cols=model.ten_J_colind[start:start+count]
        J[i,cols]=data.ten_J[start:start+count]
    return J

def apply_target(model,data,z):
    set_arm(model,data,z[:5],np.degrees(z[5:8]),np.degrees(z[8]))

def audit():
    tree,specs=build_gimbal();model=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));desired=mujoco.MjData(model)
    cfg,fit,anchors=inputs();actions=target_actions();zero=next(iter(actions.values()))[0]
    humerus=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'humerus');ulna=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'ulna')
    wrist_local=v8.design_to_mj(np.array([.087155743,-.996194698,0.])*237.)
    rows=[];tables={};KP=10000.;dt=model.opt.timestep;duration=12.;hold=2.
    for action,points in actions.items():
        dense=[]
        for a,b in zip(points[:-1],points[1:]):
            for f in np.linspace(0.,1.,5,endpoint=False):dense.append(a*(1-f)+b*f)
        dense.append(points[-1]);dense=np.array(dense);control=[];preloads=[];reference=[];failure=None
        for z in dense:
            apply_target(model,desired,z);J=native_jacobian(model,desired)
            lp=linprog(np.ones(model.ntendon),A_eq=J.T,b_eq=np.zeros(model.nv),bounds=[(.5,None)]*model.ntendon,method='highs')
            if not lp.success:failure='no_positive_equilibrium_at_target';break
            if np.max(lp.x)>200.:failure='preload_exceeds_diagnostic_force_range';break
            control.append(desired.ten_length-lp.x/KP);preloads.append(lp.x)
            wrist=desired.xpos[ulna]+desired.xmat[ulna].reshape(3,3)@wrist_local
            reference.append({'humerus_R':desired.xmat[humerus].reshape(3,3).copy(),'wrist':wrist.copy(),'qpos':desired.qpos.copy()})
        if failure:
            row={'action':action,'pass':False,'prerequisite_failure':failure,'completed_targets':len(control)};rows.append(row);print(json.dumps(row),flush=True);continue
        control=np.array(control);preloads=np.array(preloads);data=mujoco.MjData(model);apply_target(model,data,zero)
        maxpos=maxangle=maxelbow=0.;finite=True;push=False
        for step in range(int((2*duration+hold)/dt)):
            t=step*dt;raw=t/duration if t<duration else (1. if t<duration+hold else max(0.,1-(t-duration-hold)/duration))
            f=raw*raw*(3-2*raw);x=f*(len(dense)-1);k=min(int(x),len(dense)-2);blend=x-k
            data.ctrl[:]=control[k]*(1-blend)+control[k+1]*blend;mujoco.mj_step(model,data)
            if step%20:continue
            z=dense[k]*(1-blend)+dense[k+1]*blend;apply_target(model,desired,z)
            wrist=data.xpos[ulna]+data.xmat[ulna].reshape(3,3)@wrist_local
            wrist_goal=desired.xpos[ulna]+desired.xmat[ulna].reshape(3,3)@wrist_local
            error=float(np.linalg.norm(wrist-wrist_goal)*1000.)
            angle=float(np.degrees(Rotation.from_matrix(data.xmat[humerus].reshape(3,3)@desired.xmat[humerus].reshape(3,3).T).magnitude()))
            eid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'elbow');elbow_error=abs(float(np.degrees(data.qpos[model.jnt_qposadr[eid]]-desired.qpos[model.jnt_qposadr[eid]])))
            maxpos=max(maxpos,error);maxangle=max(maxangle,angle);maxelbow=max(maxelbow,elbow_error)
            finite=finite and bool(np.isfinite(data.qpos).all());push=push or bool(np.any(data.actuator_force>1e-8))
        row={'action':action,'pass':bool(finite and not push and maxpos<=5. and maxangle<=2. and maxelbow<=2.),
            'maximum_wrist_tracking_error_mm':maxpos,'maximum_humerus_orientation_error_deg':maxangle,'maximum_elbow_error_deg':maxelbow,
            'maximum_preload_N':float(preloads.max()),'finite':finite,'positive_push':push}
        rows.append(row);tables[f'{action}_q']=dense;tables[f'{action}_ctrl']=control;print(json.dumps(row),flush=True)
    result={'scope':'All 16 spatial tendons, 9 scalar DOFs; no gravity or payload; qpos initialized once per independent action',
            'gh_mechanism':'user-authorized 3-axis concentric yokes with separate antagonistic sector spools',
            'nq':model.nq,'nv':model.nv,'nu':model.nu,'pass':all(r['pass'] for r in rows),'actions':rows,'full_geometry_verified':False}
    EXPORT.mkdir(parents=True,exist_ok=True);(EXPORT/'gimbal_drive_report.json').write_text(json.dumps(result,indent=2));np.savez_compressed(EXPORT/'gimbal_control_tables.npz',**tables)
    tree.write(EXPORT/'gimbal_arm.xml',encoding='utf-8',xml_declaration=True)
if __name__=='__main__':audit()
