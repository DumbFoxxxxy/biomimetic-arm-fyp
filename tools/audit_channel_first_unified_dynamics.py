"""Independent solid audit of computed free dynamic states, not planned angles."""
from itertools import combinations
import sys
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/unified_SC_AC'
    source_file=folder/'spool_drive_dynamic_screen.json';source=read(source_file)
    endpoint_file=ROOT/source['endpoint_report']['file'];endpoint=read(endpoint_file)
    assert provenance(endpoint_file)==source['endpoint_report']
    chain_file=ROOT/endpoint['chain_check']['file'];chain=read(chain_file)
    assert provenance(chain_file)==endpoint['chain_check']
    reference_file=ROOT/chain['reference']['file'];reference=read(reference_file)
    assert provenance(reference_file)==chain['reference']
    model_file=ROOT/source['model']['file'];assert provenance(model_file)==source['model']
    model=load_mujoco_model(model_file);data=mujoco.MjData(model)
    solids={k:load_component(v) for k,v in {**reference['sources'],**reference['joint_context_sources']}.items()}
    owners={r['source_key']:r['owner'] for r in chain['geom_owners']}
    origins={k:np.array(v)*1000 for k,v in chain['neutral_body_origins_m'].items()};origins['thorax']=np.zeros(3)
    bodies={k:mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,k) for k in origins}
    pairs=[(a,b) for a,b in combinations(solids,2) if owners[a]!=owners[b]]
    rows=[]
    for trial in source['trials']:
        checked=0;failure=None
        for state in trial['trace']:
            mujoco.mj_resetData(model,data)
            mujoco.mj_integratePos(model,data.qpos,np.array(state['joint_displacement_rad']),1.)
            mujoco.mj_forward(model,data)
            transforms={k:np.c_[data.xmat[i].reshape(3,3),data.xpos[i]*1000-data.xmat[i].reshape(3,3)@origins[k]] for k,i in bodies.items()}
            placed={k:s.transform(transforms[owners[k]]) for k,s in solids.items()};hits=[]
            for a,b in pairs:
                volume=float(A.I.intersection_volume(placed[a],placed[b]))
                if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
            checked+=1
            if hits:
                failure=dict(time_s=state['time_s'],hits=hits);break
        delta=np.array([s['joint_displacement_rad'] for s in trial['trace']])
        changes=[]
        for block in endpoint['neutral_joint_blocks'].values():
            changes.append(float(np.degrees(np.linalg.norm(np.diff(delta[:,block['velocity_dof_indices']],axis=0),axis=1)).max()))
        rows.append(dict(trial=trial['trial'],checked_states=checked,available_states=len(trial['trace']),
            first_interference=failure,max_recorded_joint_rotvec_step_deg=max(changes),
            sampled_rigid_clear=failure is None))
        print(rows[-1],flush=True)
    save(folder/'spool_dynamic_solid_audit.json',dict(source=provenance(source_file),
        reference=provenance(reference_file),model=provenance(model_file),cross_body_pairs=len(pairs),
        collider_count=len(solids),threshold_mm3=.001,trials=rows,saved_sample_interval_s=.01,
        scope='Independent full rigid cross-body Boolean audit at recorded force-driven states; stops at first overlap. Visual pad envelopes excluded as deformable conditional contact law, not rigid parts. Same-body mounting fit separate. No actual cables, supply conduits, lattice or rigid collision response in simulated dynamics; an overlap invalidates that trial as a physical motion prediction.'))


if __name__=='__main__':main()
