"""Short force-driven trial of physical channel-node support, no pose commands."""
import sys
import argparse
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from simulate_channel_first_scapular_pads import pad_forces
from simulate_channel_first_spool_drive import cable_tensions


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_channel_beams'
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder')
    parser.add_argument('--duration',type=float,default=.2)
    parser.add_argument('--stroke-mm',type=float,default=2.)
    parser.add_argument('--cycle',action='store_true')
    parser.add_argument('--contact-stiffness-scale',type=float,default=1.)
    parser.add_argument('--contact-margin-mm',type=float)
    parser.add_argument('--result-name',default='short_response')
    args=parser.parse_args()
    assert args.duration>0 and args.stroke_mm>=0 and args.contact_stiffness_scale>0
    assert args.contact_margin_mm is None or args.contact_margin_mm>=0
    assert args.result_name and all(c.isalnum() or c in '_-' for c in args.result_name)
    if args.folder:folder=ROOT/args.folder
    check_file=folder/'model_check.json';check=read(check_file)
    parent=read(ROOT/check['parent']['file']);pad=read(ROOT/parent['parent_pads']['file'])
    support=read(ROOT/pad['support_direction_screen']['file'])
    model_file=ROOT/check['model']['file'];assert provenance(model_file)==check['model']
    model=load_mujoco_model(model_file)
    if args.contact_stiffness_scale!=1.:
        assert model.npair>0 and (model.pair_solref[:,0]<0).all()
        model.pair_solref[:,0]*=args.contact_stiffness_scale
        model.pair_solref[:,1]*=np.sqrt(args.contact_stiffness_scale)
    if args.contact_margin_mm is not None:
        assert model.npair>0
        model.pair_margin[:]=args.contact_margin_mm*.001
    data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    nominal=np.array(support['compact_force_solution']['cable_tensions_N'])
    available0=data.ten_length-nominal/10000.;initial=data.qpos.copy()
    endpoint=read(ROOT/parent['endpoint_report']['file'])
    flex=next(i for i,c in enumerate(endpoint['channels']) if c['id']=='FJ1486')
    beam_ids={row['head']:[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,f"{row['head']}__node_{i}")
                     for i in range(row['specification']['segments']+1)] for row in check['channels']}
    ids={row['id']:[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,name) for name in row['site_names']]
         for row in check['wire_paths']} if check.get('wire_paths') else beam_ids
    axial_q=[];rest=[];beam_slices={};strain_extrema={}
    for row in check['channels']:
        start=len(axial_q)
        for i,length in enumerate(row['specification']['segment_rest_lengths_m']):
            j=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,f"{row['head']}__axial_{i}")
            axial_q.append(model.jnt_qposadr[j]);rest.append(length)
        beam_slices[row['head']]=slice(start,len(axial_q));strain_extrema[row['head']]=[0.,0.]
    rest=np.array(rest);maximum_strain=0.;limit_hit=False
    ball_dofs=np.array([[model.jnt_dofadr[j]+k for k in range(3)] for j in range(model.njnt) if model.jnt_type[j]==mujoco.mjtJoint.mjJNT_BALL])
    hinge_dofs=np.array([model.jnt_dofadr[j] for j in range(model.njnt) if model.jnt_type[j]==mujoco.mjtJoint.mjJNT_HINGE])
    last_saved=initial.copy();maximum_saved_rotation=0.
    elbow=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'elbow');elbow_q=model.jnt_qposadr[elbow]
    trace=[];max_constraint=0.;saturated=0;max_contact_penetration=0.;contact_steps=0
    maximum_tension=np.zeros(model.nu);minimum_elbow=0.;maximum_elbow=0.;pad_edge=False;pad_bottom=False
    command_times=np.array([.1,.3,.45,.7])*args.duration if args.cycle else np.array([.25,.5,.75,1.])*args.duration
    command_values=np.array([0.,args.stroke_mm,args.stroke_mm,0.])*.001
    steps=round(args.duration/model.opt.timestep)
    clamps=[]
    for row in check.get('external_clamps',[]):
        clamps.append(dict(id=row['id'],equalities=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_EQUALITY,n) for n in row['equalities']],
            sheath_sites=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,n) for n in row['sheath_sites']],
            housing_sites=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,n) for n in row['housing_sites']]))
    def clamp_diagnostics():
        result=[]
        for clamp in clamps:
            forces=[];expected=np.zeros(model.nv);multipliers=np.zeros(data.nefc)
            for eq,sheath_site,housing_site in zip(clamp['equalities'],clamp['sheath_sites'],clamp['housing_sites']):
                mask=(data.efc_type==mujoco.mjtConstraint.mjCNSTR_EQUALITY)&(data.efc_id==eq)
                force=np.array(data.efc_force[mask]);assert len(force)==3
                forces.append(force);multipliers[mask]=force
                j1=np.zeros((3,model.nv));j2=np.zeros_like(j1);jr=np.zeros_like(j1)
                mujoco.mj_jacSite(model,data,j1,jr,sheath_site);mujoco.mj_jacSite(model,data,j2,jr,housing_site)
                expected+=(j1-j2).T@force
            actual=np.zeros(model.nv);mujoco.mj_mulJacTVec(model,data,actual,multipliers)
            reconstruction=float(np.max(abs(actual-expected),initial=0));assert reconstruction<1e-7
            outlet,seat=data.site_xpos[clamp['housing_sites']]
            result.append(dict(id=clamp['id'],force_on_sheath_at_outlet_N=forces[0].tolist(),
                force_on_sheath_at_inboard_engagement_N=forces[1].tolist(),
                net_force_on_humerus_N=(-forces[0]-forces[1]).tolist(),
                moment_on_humerus_about_seat_Nm=(-np.cross(outlet-seat,forces[0])).tolist(),
                jacobian_force_reconstruction_max_abs_mixed_units=reconstruction))
        return result
    for step in range(steps):
        mujoco.mj_step1(model,data);available=available0.copy()
        pull=float(np.interp(data.time,command_times,command_values));available[flex]-=pull
        tension,slack=cable_tensions(data.ten_length,data.ten_velocity,available,10000.,.5)
        saturated+=int(np.count_nonzero(tension>60.));data.ctrl[:]=np.minimum(tension,60.)
        data.qfrc_applied[:]=0.;loads=pad_forces(model,data,parent['pads'],1.,3.)
        maximum_tension=np.maximum(maximum_tension,data.ctrl)
        pad_edge|=any(p['over_pad_edge'] for p in loads);pad_bottom|=any(p['pad_bottomed'] for p in loads)
        equality=data.efc_type==mujoco.mjtConstraint.mjCNSTR_EQUALITY
        max_constraint=max(max_constraint,float(np.max(abs(data.efc_pos[equality]),initial=0)))
        contact_steps+=int(data.ncon>0)
        for i in range(data.ncon):max_contact_penetration=max(max_contact_penetration,-float(data.contact[i].dist))
        mujoco.mj_step2(model,data)
        assert np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
        minimum_elbow=min(minimum_elbow,float(data.qpos[elbow_q]));maximum_elbow=max(maximum_elbow,float(data.qpos[elbow_q]))
        strain=float(np.max(abs(data.qpos[axial_q]/rest)))
        all_strains=data.qpos[axial_q]/rest
        for name,sl in beam_slices.items():
            strain_extrema[name][0]=min(strain_extrema[name][0],float(np.min(all_strains[sl])))
            strain_extrema[name][1]=max(strain_extrema[name][1],float(np.max(all_strains[sl])))
        maximum_strain=max(maximum_strain,strain);limit_hit|=strain>=.49
        increment=np.zeros(model.nv);mujoco.mj_differentiatePos(model,increment,1.,last_saved,data.qpos)
        angular_increment=max(float(np.linalg.norm(increment[ball_dofs],axis=1).max()),float(abs(increment[hinge_dofs]).max()))
        if step%100==0 or np.degrees(angular_increment)>=.5 or step==steps-1:
            reactions=clamp_diagnostics()
            maximum_saved_rotation=max(maximum_saved_rotation,float(np.degrees(angular_increment)))
            last_saved=data.qpos.copy()
            mujoco.mj_forward(model,data)
            delta=np.zeros(model.nv);mujoco.mj_differentiatePos(model,delta,1.,initial,data.qpos)
            trace.append(dict(time_s=float(data.time),generalized_displacement=delta.tolist(),
                tensions_N=data.ctrl.tolist(),slack=slack.tolist(),pad_loads=loads,
                commanded_pull_mm=pull*1000.,tendon_lengths_mm=(data.ten_length*1000.).tolist(),
                clamp_reactions=reactions,clamp_reaction_evaluation_time_s=float(data.time-model.opt.timestep),
                beam_nodes_source_mm={head:(data.site_xpos[sites]*1000).tolist() for head,sites in beam_ids.items()},
                wire_nodes_source_mm={head:(data.site_xpos[sites]*1000).tolist() for head,sites in ids.items()}))
    warnings={str(i):int(w.number) for i,w in enumerate(data.warning) if w.number}
    output=folder/(args.result_name+'.json')
    save(output,dict(model_check=provenance(check_file),model=provenance(model_file),
        generator=provenance(ROOT/'tools/simulate_channel_first_elbow_channel_beams.py'),duration_s=args.duration,
        command=dict(cycle=args.cycle,times_s=command_times.tolist(),pull_m=command_values.tolist()),
        runtime_contact_parameters=dict(stiffness_scale=args.contact_stiffness_scale,pair_solref=model.pair_solref.tolist(),
            pair_margin_m=model.pair_margin.tolist()),
        elbow_range_deg=np.degrees([minimum_elbow,maximum_elbow]).tolist(),
        final_elbow_deg=float(np.degrees(data.qpos[elbow_q])),final_maximum_speed=float(np.max(abs(data.qvel))),
        peak_tensions_N=maximum_tension.tolist(),pad_edge_reached=pad_edge,pad_bottomed=pad_bottom,
        timestep_s=model.opt.timestep,saved_interval_s=.002,trace=trace,
        adaptive_rotation_trigger_deg=.5,maximum_recorded_rotation_increment_deg=maximum_saved_rotation,
        maximum_endpoint_constraint_error_m=max_constraint,saturated_channel_steps=saturated,warning_counts=warnings,
        contact_steps=contact_steps,maximum_contact_penetration_m=max_contact_penetration,
        maximum_absolute_segment_axial_strain=maximum_strain,numerical_axial_guard_hit=limit_hit,
        beam_axial_strain_extrema=strain_extrema,
        clamp_reaction_note='Pre-integration constraint forces at evaluation_time, paired with saved pose one timestep later. Cartesian signs checked by exact J.T multiplication. Only humerus distal external clamps reported; no thorax-root wrench yet.',
        state_units='MuJoCo tangent vector: rotational components radians, slide components metres',
        scope='Short numerical response only, not full cycle or accepted motion. Curve nodes move through beam axial/bending forces. Original lattice not simulated; effective tube stiffness assumed. Any contact-solved parts are specified by model_check; no general bone/self contact, lumen collapse, friction or supplied conduit reaction. Constraint accuracy and actual wire geometry need auditing.'))
    print(output.name,'complete; elbow deg',np.degrees([minimum_elbow,maximum_elbow]),'warnings',warnings,'max end constraint mm',max_constraint*1000,'saturation steps',saturated,flush=True)


if __name__=='__main__':main()
