"""Check restored rigid cap packs at independently computed dynamic poses."""
from itertools import combinations
import sys
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_passive_model'
    result_file=folder/'free_swivel_response.json';result=read(result_file)
    check_file=ROOT/result['model_check']['file'];check=read(check_file);assert provenance(check_file)==result['model_check']
    assembly_file=ROOT/check['assembly']['file'];assembly=read(assembly_file);assert provenance(assembly_file)==check['assembly']
    model_file=ROOT/check['model']['file'];assert provenance(model_file)==check['model']
    model=load_mujoco_model(model_file);data=mujoco.MjData(model)
    solids={k:load_component(v) for k,v in assembly['sources'].items()}
    owners={r['source_key']:r['owner'] for r in check['geom_owners']}
    pairs=[(a,b) for a,b in combinations(solids,2) if owners[a]!=owners[b]]
    origins={k:np.array(v)*1000 for k,v in check['neutral_origins_m'].items()}
    bodies={k:mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,k) for k in origins}
    trials=[]
    for case in result['cases']:
        delta=np.array([s['joint_displacement_rad'] for s in case['trace']])
        differences=np.diff(delta,axis=0)
        max_step=max([float(np.degrees(np.linalg.norm(differences[:,i:i+3],axis=1)).max()) for i in (0,3,6)]+
                     [float(np.degrees(abs(differences[:,9:])).max())])
        assert max_step<=1.,'Saved poses exceed one-degree sampling requirement'
        checked=0;failure=None
        for state in case['trace']:
            mujoco.mj_resetData(model,data);mujoco.mj_integratePos(model,data.qpos,np.array(state['joint_displacement_rad']),1.)
            mujoco.mj_forward(model,data)
            transforms={k:np.c_[data.xmat[i].reshape(3,3),data.xpos[i]*1000-data.xmat[i].reshape(3,3)@origins[k]] for k,i in bodies.items()}
            placed={k:s.transform(transforms[owners[k]]) for k,s in solids.items()};hits=[]
            for a,b in pairs:
                volume=float(A.I.intersection_volume(placed[a],placed[b]))
                if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
            checked+=1
            if hits:failure=dict(time_s=state['time_s'],cap_angles_deg=state['cap_angles_deg'],hits=hits);break
        trials.append(dict(torsion_Nm_per_rad=case['torsion_Nm_per_rad'],checked_states=checked,
                           available_states=len(case['trace']),first_interference=failure,max_sample_joint_step_deg=max_step))
        print(trials[-1],flush=True)
    save(folder/'free_swivel_rigid_audit.json',dict(source=provenance(result_file),model_check=provenance(check_file),
        cross_body_pairs=len(pairs),rigid_parts=len(solids),saved_interval_s=.005,trials=trials,
        scope='All cross-body rigid solids, including original passive cap packs; no whole adjacent pair exclusions. Same-body assembly fit separate. Stops at first >0.001mm3 overlap. Deforming muscle and full cable route not checked. Simulation has no rigid collision response, so overlap rejects physical motion interpretation.'))


if __name__=='__main__':main()
