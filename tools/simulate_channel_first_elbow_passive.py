"""Free cap-angle response to unilateral cable length commands (proxy routing)."""
import sys
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from simulate_channel_first_scapular_pads import pad_forces
from simulate_channel_first_spool_drive import cable_tensions


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_passive_model'
    check_file=folder/'model_check.json';check=read(check_file)
    model_file=ROOT/check['model']['file'];assert provenance(model_file)==check['model']
    pad=read(ROOT/check['parent_pads']['file']);support=read(ROOT/pad['support_direction_screen']['file'])
    nominal=np.array(support['compact_force_solution']['cable_tensions_N'])
    endpoint=read(ROOT/check['endpoint_report']['file'])
    flexor=next(i for i,c in enumerate(endpoint['channels']) if c['id']=='FJ1486')
    cases=[]
    for torsion in (.02,.2,2.):
        model=load_mujoco_model(model_file);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
        initial=data.qpos.copy();available0=data.ten_length-nominal/10000.
        ids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,p['head']+'_swivel') for p in check['pivots']]
        for i in ids:model.jnt_stiffness[i]=torsion
        trace=[];saturations=0;bottomed=False;over_edge=False
        for step in range(round(2./model.opt.timestep)):
            mujoco.mj_step1(model,data)
            available=available0.copy()
            available[flexor]-=np.interp(data.time,[.3,.8,1.1,1.6],[0.,.006,.006,0.])
            requested,slack=cable_tensions(data.ten_length,data.ten_velocity,available,10000.,.5)
            saturations+=int(np.count_nonzero(requested>60.));data.ctrl[:]=np.minimum(requested,60.)
            data.qfrc_applied[:]=0.;loads=pad_forces(model,data,check['pads'],1.,3.)
            bottomed|=any(p['pad_bottomed'] for p in loads);over_edge|=any(p['over_pad_edge'] for p in loads)
            mujoco.mj_step2(model,data)
            assert np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
            if step%25==0:
                delta=np.zeros(model.nv);mujoco.mj_differentiatePos(model,delta,1.,initial,data.qpos)
                trace.append(dict(time_s=float(data.time),joint_displacement_rad=delta.tolist(),
                    cap_angles_deg={model.joint(i).name:float(np.degrees(data.qpos[model.jnt_qposadr[i]])) for i in ids},
                    tensions_N=data.ctrl.tolist(),slack=slack.tolist(),pad_loads=loads))
        peak={model.joint(i).name:max(abs(r['cap_angles_deg'][model.joint(i).name]) for r in trace) for i in ids}
        warnings={str(i):int(w.number) for i,w in enumerate(data.warning) if w.number}
        cases.append(dict(torsion_Nm_per_rad=torsion,peak_abs_cap_angles_deg=peak,trace=trace,
            saturated_channel_steps=saturations,pad_bottomed=bottomed,left_pad_footprint=over_edge,warning_counts=warnings))
        print('torsion',torsion,'cap peaks',peak,'saturated',saturations,'pad edge',over_edge,flush=True)
    save(folder/'free_swivel_response.json',dict(model_check=provenance(check_file),model=provenance(model_file),
        generator=provenance(ROOT/'tools/simulate_channel_first_elbow_passive.py'),cases=cases,
        scope='2s free response to available-length input. Joint positions never prescribed. Rotational stiffness is unmeasured equivalent TPU/end compliance, not demonstrated cap performance. Rigid collisions and all real channel/steel friction remain unchecked; no motion range acceptance.'))


if __name__=='__main__':main()
