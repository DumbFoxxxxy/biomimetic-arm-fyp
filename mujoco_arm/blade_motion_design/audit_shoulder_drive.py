"""Actual pull-only tendon-driven simulation, never resetting tracked qpos."""
from pathlib import Path
import json,sys,math
import numpy as np
import mujoco
from scipy.spatial.transform import Rotation
from shoulder_model import compile_model,inputs,equilibrium_tensions,EXPORT,KP,v8,input_fingerprints
sys.stdout.reconfigure(encoding='utf-8')

def audit():
    config,fit,anchors=inputs(); model=compile_model(); target=mujoco.MjData(model)
    zero=np.array(fit['poses'][0]['q_rad']); samples={}; max_tension=0.
    for action in ['UP_ROT','DOWN_ROT','PROTRACT','RETRACT','ELEVATE','DEPRESS']:
        qs=np.array([zero]+[p['q_rad'] for p in fit['poses'] if p.get('action')==action]); controls=[];tensions=[]
        # Five subdivisions per original interval, with physical target lengths.
        dense=[]
        for a,b in zip(qs[:-1],qs[1:]):
            for f in np.linspace(0.,1.,5,endpoint=False):dense.append(a*(1-f)+b*f)
        dense.append(qs[-1])
        for q in dense:
            v8.set_joint_state(model,target,q);mujoco.mj_forward(model,target)
            force,J=equilibrium_tensions(model,target)
            controls.append(target.ten_length-force/KP);tensions.append(force)
        samples[action]={'q':np.array(dense),'ctrl':np.array(controls),'tensions':np.array(tensions)}
        max_tension=max(max_tension,float(np.max(tensions)))
    print('FEEDFORWARD_MAX_TENSION_N',max_tension,flush=True)
    EXPORT.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(EXPORT/'shoulder_control_tables.npz',**{f'{action}_{key}':value for action,table in samples.items() for key,value in table.items()})
    if max_tension>200:raise RuntimeError('Feedforward exceeds diagnostic actuator force range')
    rows=[]; dt=model.opt.timestep; scap=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'scapula')
    ac=np.array(config['geometry_mm']['AC']); p=np.array(config['geometry_mm']['P'])
    for action,table in samples.items():
        data=mujoco.MjData(model);v8.set_joint_state(model,data,zero);mujoco.mj_forward(model,data)
        max_pos=max_angle=0.; endpoint={}; finite=True; push=False; max_force=0.
        duration=12.; hold=2.; total=2*duration+hold
        for step in range(int(total/dt)):
            t=step*dt
            if t<duration:raw=t/duration
            elif t<duration+hold:raw=1.
            else:raw=max(0.,1.-(t-duration-hold)/duration)
            f=raw*raw*(3.-2.*raw); index=f*(len(table['q'])-1); k=min(int(index),len(table['q'])-2); blend=index-k
            q=table['q'][k]*(1-blend)+table['q'][k+1]*blend
            data.ctrl[:]=table['ctrl'][k]*(1-blend)+table['ctrl'][k+1]*blend
            mujoco.mj_step(model,data)
            if step%20:continue
            desired=v8.direct_kinematics(config,q)
            R=v8.DESIGN_TO_MJ.T@data.xmat[scap].reshape(3,3)@v8.DESIGN_TO_MJ
            actual_ac=v8.mj_to_design(data.xpos[scap]); actual_p=actual_ac+R@(p-ac)
            pos=float(np.linalg.norm(actual_p-desired['contact_mm']))
            angle=float(np.degrees(Rotation.from_matrix(R@desired['body_rotation'].T).magnitude()))
            max_pos=max(max_pos,pos);max_angle=max(max_angle,angle)
            finite=finite and bool(np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all())
            push=push or bool(np.any(data.actuator_force>1e-8));max_force=max(max_force,float(np.max(np.abs(data.actuator_force))))
            if duration+hold-.03<t<duration+hold:
                endpoint={'contact_error_mm':pos,'orientation_error_deg':angle,'qpos':data.qpos.tolist()}
        row={'action':action,'finite':finite,'positive_push_force':push,'maximum_contact_tracking_error_mm':max_pos,
             'maximum_orientation_tracking_error_deg':max_angle,'maximum_tension_N':max_force,'endpoint':endpoint,
             'pass':bool(finite and not push and max_pos<=3. and max_angle<=2.)}
        rows.append(row);print(json.dumps(row),flush=True)
    result={'scope':'No gravity/no payload; genuine tendon position actuation and prestress feedforward; qpos set only at run initialization',
            'nq':model.nq,'nv':model.nv,'nu':model.nu,'no_contact_reaction':True,
            'duration_one_way_s':duration,'max_feedforward_tension_N':max_tension,
            'pass':all(r['pass'] for r in rows),'actions':rows,'geometry_verified_here':False,'input_sha256':input_fingerprints()}
    EXPORT.mkdir(parents=True,exist_ok=True);(EXPORT/'shoulder_drive_report.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':audit()
