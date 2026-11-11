"""Inverse geometric contact-following search, followed by solid and force checks."""
from itertools import combinations
import sys
import numpy as np
from scipy.optimize import least_squares,lsq_linear
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance
from simulate_channel_first_coordinated_spools import endpoint_moment_matrix


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/unified_SC_AC'
    pad_file=folder/'passive_pad_engagement_paths.json';pad_report=read(pad_file);pads=pad_report['pads']
    endpoint_file=ROOT/pad_report['endpoint_model']['file'];endpoint=read(endpoint_file)
    assert provenance(endpoint_file)==pad_report['endpoint_model']
    chain_file=ROOT/endpoint['chain_check']['file'];chain=read(chain_file)
    assert provenance(chain_file)==endpoint['chain_check']
    ref_file=ROOT/chain['reference']['file'];ref=read(ref_file)
    assert provenance(ref_file)==chain['reference']
    model_file=ROOT/endpoint['generated_model']['file'];assert provenance(model_file)==endpoint['generated_model']
    model=load_mujoco_model(model_file);data=mujoco.MjData(model)
    scap=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'scapula')
    hum=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'humerus')
    fore=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'forearm')
    blocks={k:v['velocity_dof_indices'] for k,v in endpoint['neutral_joint_blocks'].items()}
    origins={k:np.array(v)*1000 for k,v in chain['neutral_body_origins_m'].items()};origins['thorax']=np.zeros(3)
    bodies={k:mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,k) for k in origins}
    owners={r['source_key']:r['owner'] for r in chain['geom_owners']}
    solids={k:load_component(v) for k,v in {**ref['sources'],**ref['joint_context_sources']}.items()}
    pairs=[(a,b) for a,b in combinations(solids,2) if owners[a]!=owners[b]]
    def pose(axis,angle,ac):
        mujoco.mj_resetData(model,data);delta=np.zeros(model.nv)
        delta[blocks['SC'][axis]]=np.radians(angle);delta[blocks['AC']]=ac
        mujoco.mj_integratePos(model,data.qpos,delta,1.);mujoco.mj_forward(model,data)
    def contacts():
        out=[];reaction=np.zeros(model.nv)
        for pad in pads:
            p=data.xpos[scap]+data.xmat[scap].reshape(3,3)@pad['local_point_m']
            n=np.array(pad['normal']);offset=p-pad['rib_point_m'];h=float(offset@n)
            lateral=float(np.linalg.norm(offset-h*n));active=h<pad['rest_separation_m'] and lateral<=pad['radius_m']
            f=max(0.,pad['rest_separation_m']-h)*pad['stiffness_N_per_m'] if active else 0.
            jac=np.zeros((3,model.nv));mujoco.mj_jac(model,data,jac,np.zeros_like(jac),p,scap)
            reaction+=jac.T@n*f
            out.append(dict(height_mm=h*1000,lateral_mm=lateral*1000,force_N=f,active=active,bottomed=h<=0))
        return out,reaction
    def solid_hits():
        transforms={k:np.c_[data.xmat[i].reshape(3,3),data.xpos[i]*1000-data.xmat[i].reshape(3,3)@origins[k]] for k,i in bodies.items()}
        placed={k:s.transform(transforms[owners[k]]) for k,s in solids.items()};hits=[]
        for a,b in pairs:
            volume=float(A.I.intersection_volume(placed[a],placed[b]))
            if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
        return hits
    pose(0,0,np.zeros(3));neutral,_=contacts();target=np.array([p['height_mm'] for p in neutral])
    direction0=data.xpos[fore]-data.xpos[hum];direction0/=np.linalg.norm(direction0)
    assert not solid_hits()
    paths=[]
    for axis in range(3):
        for sign in (-1,1):
            last=0.;ac=np.zeros(3);step=.5;states=[];failure=None
            while last<10.-1e-8:
                angle=min(10.,last+step)
                def objective(candidate):
                    pose(axis,sign*angle,candidate)
                    c,_=contacts();return np.array([p['height_mm'] for p in c])-target
                fit=least_squares(objective,ac,bounds=(-np.pi/6,np.pi/6),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=300)
                angular_step=float(np.degrees(np.linalg.norm(fit.x-ac)))
                pose(axis,sign*angle,fit.x);c,reaction=contacts()
                error=float(np.max(abs(fit.fun)));collision=solid_hits() if error<=.01 else []
                reason=('contact_constraint_unsolved_in_AC_30deg_search' if error>.01 else
                        'AC_step_above_1deg' if angular_step>1. else
                        'pad_detached_or_bottomed' if any(not p['active'] or p['bottomed'] for p in c) else
                        'rigid_interference' if collision else None)
                if reason:
                    if step>.0125+1e-8:step/=2;continue
                    failure=dict(reason=reason,SC_angle_deg=sign*angle,AC_rotvec_deg=np.degrees(fit.x).tolist(),
                                 max_contact_height_error_mm=error,hits=collision)
                    break
                matrix=endpoint_moment_matrix(model,data,endpoint['channels'])
                tension=lsq_linear(matrix,data.qfrc_bias-data.qfrc_passive-reaction,bounds=(0.,60.),tol=1e-11,max_iter=1000)
                direction=data.xpos[fore]-data.xpos[hum];direction/=np.linalg.norm(direction)
                states.append(dict(SC_angle_deg=sign*angle,AC_rotvec_deg=np.degrees(fit.x).tolist(),contacts=c,
                    max_contact_height_error_mm=error,tensions_N=tension.x.tolist(),
                    residual_norm_Nm=float(np.linalg.norm(tension.fun)),
                    actual_humerus_direction_change_deg=float(np.degrees(np.arccos(np.clip(direction@direction0,-1,1)))),
                    endpoint_span_lengths_mm=(data.ten_length*1000).tolist()))
                last,ac=angle,fit.x
            for state in reversed(states):
                pose(axis,state['SC_angle_deg'],np.radians(state['AC_rotvec_deg']))
                assert not solid_hits()
            paths.append(dict(axis=axis,direction=sign,states=states,last_SC_angle_deg=sign*last,
                              first_failure=failure,return_samples=len(states)))
            print(axis,sign,'last SC',sign*last,'samples',len(states),'failure',failure['reason'] if failure else 'study cap',flush=True)
    save(folder/'contact_following_paths.json',dict(pads=provenance(pad_file),endpoint_model=provenance(endpoint_file),
        reference=provenance(ref_file),generator=provenance(ROOT/'tools/study_channel_first_contact_following_paths.py'),
        paths=paths,cross_body_pairs=len(pairs),height_tolerance_mm=.01,AC_rotvec_component_search_bound_deg=30.,
        scope='Inverse kinematic three-height constraint search only. Contact and tension evaluated after pose solution, not a driven trajectory. Entire same-body pad footprint, forces, material and mounts remain unverified. SC study cap10deg, adaptive joint increments<=1deg; failed branch bracket SC<=0.0125deg but not a global feasible-region proof. No K-range acceptance.'))


if __name__=='__main__':main()
