"""Test actual activation of a proposed three-pad law along geometric branches."""
import sys
import numpy as np
from scipy.optimize import lsq_linear
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from simulate_channel_first_coordinated_spools import endpoint_moment_matrix


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/unified_SC_AC'
    support_file=folder/'scapulothoracic_support_direction_screen.json';support=read(support_file)
    endpoint_file=ROOT/support['endpoint_model']['file'];endpoint=read(endpoint_file)
    assert provenance(endpoint_file)==support['endpoint_model']
    force_file=folder/'endpoint_force_balance_paths.json';force=read(force_file)
    assert force['endpoint_model']==provenance(endpoint_file)
    chain=read(ROOT/endpoint['chain_check']['file'])
    model_file=ROOT/endpoint['generated_model']['file'];assert provenance(model_file)==endpoint['generated_model']
    model=load_mujoco_model(model_file);data=mujoco.MjData(model)
    body=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'scapula')
    origin=np.array(chain['neutral_body_origins_m']['scapula'])
    compact=support['compact_force_solution'];pads=[]
    for index in compact['candidate_indices']:
        candidate=support['candidates'][index];nominal=compact['support_forces_N'][index]
        pads.append(dict(candidate_index=index,
            local_point_m=(np.array(candidate['scapular_surface_source_mm'])*.001-origin).tolist(),
            rib_point_m=(np.array(candidate['rib_surface_source_mm'])*.001).tolist(),
            normal=candidate['normal_force_on_scapula_source'],radius_m=.012,stiffness_N_per_m=10000.,
            rest_separation_m=candidate['gap_mm']*.001+nominal/10000.,neutral_force_N=nominal))
    blocks={k:v['velocity_dof_indices'] for k,v in endpoint['neutral_joint_blocks'].items()}
    cases=[]
    for case in force['cases']:
        samples=[]
        for state in case['samples']:
            mujoco.mj_resetData(model,data);delta=np.zeros(model.nv)
            delta[blocks['SC'][case['axis']]]=np.radians(state['SC_angle_deg'])
            delta[blocks['AC'][case['axis']]]=np.radians(state['AC_angle_deg'])
            mujoco.mj_integratePos(model,data.qpos,delta,1.);mujoco.mj_forward(model,data)
            reaction=np.zeros(model.nv);contacts=[]
            for pad in pads:
                point=data.xpos[body]+data.xmat[body].reshape(3,3)@np.array(pad['local_point_m'])
                offset=point-np.array(pad['rib_point_m']);normal=np.array(pad['normal'])
                height=float(offset@normal);lateral=float(np.linalg.norm(offset-height*normal))
                active=height<pad['rest_separation_m'] and lateral<=pad['radius_m']
                compression=max(0.,pad['rest_separation_m']-height) if active else 0.
                load=pad['stiffness_N_per_m']*compression
                jac=np.zeros((3,model.nv));mujoco.mj_jac(model,data,jac,np.zeros_like(jac),point,body)
                reaction+=jac.T@normal*load
                contacts.append(dict(candidate_index=pad['candidate_index'],active=active,
                    height_mm=height*1000.,lateral_mm=lateral*1000.,compression_mm=compression*1000.,
                    force_N=load,bottomed=height<=0.))
            matrix=endpoint_moment_matrix(model,data,endpoint['channels'])
            fit=lsq_linear(matrix,data.qfrc_bias-data.qfrc_passive-reaction,bounds=(0.,60.),tol=1e-11,max_iter=1000)
            samples.append(dict(SC_angle_deg=state['SC_angle_deg'],AC_angle_deg=state['AC_angle_deg'],
                contacts=contacts,tensions_N=fit.x.tolist(),residual_norm_Nm=float(np.linalg.norm(fit.fun)),
                residual_by_joint_Nm={k:float(np.linalg.norm(fit.fun[v])) for k,v in blocks.items()}))
        cases.append(dict(axis=case['axis'],direction=case['direction'],samples=samples))
        print(case['axis'],case['direction'],'neutral',samples[0]['residual_norm_Nm'],
              'endpoint',samples[-1]['residual_norm_Nm'],'active',sum(p['active'] for p in samples[-1]['contacts']),flush=True)
    save(folder/'passive_pad_engagement_paths.json',dict(
        support=provenance(support_file),endpoint_model=provenance(endpoint_file),force_paths=provenance(force_file),
        generator=provenance(ROOT/'tools/study_channel_first_unified_pad_engagement.py'),pads=pads,cases=cases,
        scope='Static conditional three-pad law, not installed hardware. Rib-fixed finite-radius planes with compression only, no friction or tension. Forces determined by separation, not freely allocated. New reference projection performed in support study; footprint shape, mount, strain and moving contact not verified. Prescribed geometric states are not a force-driven trajectory.'))


if __name__=='__main__':main()
