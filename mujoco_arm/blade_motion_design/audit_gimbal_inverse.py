"""Full-arm tendon tracking with inverse-dynamics and velocity feedforward."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,sys,argparse,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation
from gimbal_arm_model import build_gimbal
from arm_targets import target_actions
from audit_gimbal_drive import native_jacobian,apply_target
from shoulder_model import EXPORT,v8
from verification_snapshot import fingerprint
sys.stdout.reconfigure(encoding='utf-8')

def audit(gravity=False,force_limit=200.):
    tree,spec=build_gimbal()
    if gravity:tree.getroot().find('option').set('gravity','0 0 -9.81')
    for actuator in tree.getroot().find('actuator'):actuator.set('forcerange',f'-{force_limit:g} 0')
    model=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));reference=mujoco.MjData(model)
    targets=target_actions();zero=next(iter(targets.values()))[0];duration=12.;hold=2.;total=26.;timegrid=np.linspace(0,total,401)
    hid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'humerus');uid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'ulna');eid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'elbow')
    wrist_local=v8.design_to_mj(np.array([.087155743,-.996194698,0.])*237.);rows=[];export={}
    for action,points in targets.items():
        ztable=[];positions=[]
        for t in timegrid:
            raw=t/duration if t<duration else (1. if t<duration+hold else max(0.,1-(t-duration-hold)/duration));f=raw*raw*(3-2*raw)
            x=f*(len(points)-1);k=min(int(x),len(points)-2);blend=x-k;z=points[k]*(1-blend)+points[k+1]*blend
            apply_target(model,reference,z);positions.append(reference.qpos.copy());ztable.append(z)
        positions=np.array(positions);ztable=np.array(ztable);velocity=np.zeros((len(timegrid),model.nv))
        for i in range(1,len(timegrid)-1):mujoco.mj_differentiatePos(model,velocity[i],timegrid[i+1]-timegrid[i-1],positions[i-1],positions[i+1])
        acceleration=np.gradient(velocity,timegrid,axis=0);controls=[];maxT=0.;failure=None
        for i in range(len(timegrid)):
            reference.qpos[:]=positions[i];reference.qvel[:]=velocity[i];reference.qacc[:]=acceleration[i]
            mujoco.mj_inverse(model,reference);J=native_jacobian(model,reference)
            if gravity:
                # Minimize the peak line force rather than sum-force. A
                # minimum-sum solution can overload one line unnecessarily.
                lp=linprog(np.r_[np.zeros(model.ntendon),1.],A_ub=np.column_stack([np.eye(model.ntendon),-np.ones(model.ntendon)]),b_ub=np.zeros(model.ntendon),A_eq=np.column_stack([J.T,np.zeros(model.nv)]),b_eq=-reference.qfrc_inverse,bounds=[(.5,None)]*model.ntendon+[(.5,force_limit)],method='highs')
                tensions=lp.x[:-1] if lp.success else None
            else:
                lp=linprog(np.ones(model.ntendon),A_eq=J.T,b_eq=-reference.qfrc_inverse,bounds=[(.5,force_limit)]*model.ntendon,method='highs');tensions=lp.x if lp.success else None
            if not lp.success:failure='positive_inverse_force_infeasible';break
            if tensions.max()>force_limit+1e-6:failure='diagnostic_actuator_saturation';break
            maxT=max(maxT,float(tensions.max()))
            controls.append(reference.ten_length-tensions/10000.+150./10000.*(J@velocity[i]))
        if failure:
            row={'action':action,'pass':False,'failure':failure,'time_grid_index':len(controls)};rows.append(row);print(json.dumps(row),flush=True);continue
        controls=np.array(controls);data=mujoco.MjData(model);apply_target(model,data,zero)
        maxpos=maxori=maxelbow=0.;finite=True;push=False;dt=model.opt.timestep
        for step in range(int(total/dt)):
            t=step*dt;x=t/total*(len(timegrid)-1);k=min(int(x),len(timegrid)-2);blend=x-k
            data.ctrl[:]=controls[k]*(1-blend)+controls[k+1]*blend;mujoco.mj_step(model,data)
            if step%20:continue
            z=ztable[k]*(1-blend)+ztable[k+1]*blend;apply_target(model,reference,z)
            wrist=data.xpos[uid]+data.xmat[uid].reshape(3,3)@wrist_local;goal=reference.xpos[uid]+reference.xmat[uid].reshape(3,3)@wrist_local
            maxpos=max(maxpos,float(np.linalg.norm(wrist-goal)*1000))
            maxori=max(maxori,float(np.degrees(Rotation.from_matrix(data.xmat[hid].reshape(3,3)@reference.xmat[hid].reshape(3,3).T).magnitude())))
            maxelbow=max(maxelbow,abs(float(np.degrees(data.qpos[model.jnt_qposadr[eid]]-reference.qpos[model.jnt_qposadr[eid]]))))
            finite=finite and bool(np.isfinite(data.qpos).all());push=push or bool(np.any(data.actuator_force>1e-8))
        row={'action':action,'pass':bool(finite and not push and maxpos<=5 and maxori<=2 and maxelbow<=2),
            'maximum_wrist_error_mm':maxpos,'maximum_humeral_orientation_error_deg':maxori,'maximum_elbow_error_deg':maxelbow,
            'maximum_feedforward_tension_N':maxT,'finite':finite,'positive_push':push}
        rows.append(row);export[f'{action}_q']=ztable;export[f'{action}_ctrl']=controls;print(json.dumps(row),flush=True)
    result={'scope':f'{"Gravity with estimated printed-part masses" if gravity else "No gravity"}/no payload; {model.ntendon} true spatial tendons; initialized once; inverse-dynamics and tendon-velocity feedforward',
        'scalar_dofs':model.nv,'actuators':model.nu,'one_way_duration_s':duration,'pass':all(r['pass'] for r in rows),
        'geometry_gate_pass':False,'actions':rows,'diagnostic_force_limit_N':force_limit,'tension_allocation':'minimax peak' if gravity else 'minimum sum','input_fingerprint':fingerprint(tree.getroot())}
    label='gimbal_gravity' if gravity else 'gimbal_inverse'
    if force_limit!=200.:label+=f'_{force_limit:g}N'
    (EXPORT/f'{label}_drive_report.json').write_text(json.dumps(result,indent=2));np.savez_compressed(EXPORT/f'{label}_control_tables.npz',time_s=timegrid,**export)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--gravity',action='store_true');parser.add_argument('--force-limit',type=float,default=200.);args=parser.parse_args();audit(args.gravity,args.force_limit)
