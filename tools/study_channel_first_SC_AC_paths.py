"""SC/AC geometric path search; input joint coordinates are not drive evidence."""
from itertools import combinations
import sys
import numpy as np
from build_s1_teres_origin_mount import A, ROOT, read, save, load_component
from build_channel_first_inventory import provenance


def main():
    import mujoco
    sys.path.insert(0,str(ROOT))
    from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/unified_SC_AC'
    check_file=folder/'geometry_mass_body_chain_check.json';check=read(check_file)
    reference_file=ROOT/check['reference']['file'];reference=read(reference_file)
    model_file=ROOT/check['generated_model']['file']
    assert provenance(reference_file)==check['reference'] and provenance(model_file)==check['generated_model']
    model=load_mujoco_model(model_file);data=mujoco.MjData(model)
    sources={**reference['sources'],**reference['joint_context_sources']}
    solids={k:load_component(v) for k,v in sources.items()}
    owners={r['source_key']:r['owner'] for r in check['geom_owners']}
    pairs=[(a,b) for a,b in combinations(solids,2) if owners[a]!=owners[b]]
    origins={k:np.array(v)*1000 for k,v in check['neutral_body_origins_m'].items()};origins['thorax']=np.zeros(3)
    bodies={k:mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,k) for k in origins}
    dofs={k:model.jnt_dofadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,k)] for k in ('SC','AC')}
    neutral_direction=origins['forearm']-origins['humerus'];neutral_direction/=np.linalg.norm(neutral_direction)
    def probe(axis,sc,ratio,measure_gap=False):
        mujoco.mj_resetData(model,data);delta=np.zeros(model.nv)
        delta[dofs['SC']+axis]=np.radians(sc);delta[dofs['AC']+axis]=np.radians(sc*ratio)
        mujoco.mj_integratePos(model,data.qpos,delta,1.);mujoco.mj_forward(model,data)
        transforms={k:np.c_[data.xmat[i].reshape(3,3),data.xpos[i]*1000-data.xmat[i].reshape(3,3)@origins[k]] for k,i in bodies.items()}
        placed={k:s.transform(transforms[owners[k]]) for k,s in solids.items()}
        hits=[];gap=.5
        for a,b in pairs:
            volume=float(A.I.intersection_volume(placed[a],placed[b]))
            if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
            if measure_gap:gap=min(gap,A.I.gap(placed[a],placed[b],gap))
        direction=data.xpos[bodies['forearm']]-data.xpos[bodies['humerus']];direction/=np.linalg.norm(direction)
        scapula_rotation=data.xmat[bodies['scapula']].reshape(3,3)
        return dict(SC_angle_deg=float(sc),AC_angle_deg=float(sc*ratio),hits=hits,
            minimum_gap_mm_capped0_5=gap if measure_gap else None,
            humerus_direction_source=direction.tolist(),
            humerus_direction_change_deg=float(np.degrees(np.arccos(np.clip(direction@neutral_direction,-1,1)))),
            scapula_rotation_matrix=scapula_rotation.tolist(),GH_source_mm=(data.xpos[bodies['humerus']]*1000).tolist())
    assert not probe(0,0,0)['hits']
    rows=[]
    for axis in range(3):
        for sign in (-1,1):
            for ratio in (-2.,-1.,-.5,0.,.5,1.,2.):
                step=1./max(1.,abs(ratio));last=0.;first=None;out=[]
                for magnitude in np.arange(step,10.+step/2,step):
                    state=probe(axis,sign*magnitude,ratio)
                    if state['hits']:
                        first=state;break
                    last=float(magnitude);out.append(float(sign*magnitude))
                if first is not None:
                    fine_step=.1/max(1.,abs(ratio))
                    for magnitude in np.arange(last+fine_step,abs(first['SC_angle_deg'])+fine_step/2,fine_step):
                        state=probe(axis,sign*magnitude,ratio)
                        if state['hits']:first=state;break
                        last=float(magnitude);out.append(float(sign*magnitude))
                # Independent return evaluation of each previously clear sample.
                for angle in sorted(set(out+[0.]),key=abs,reverse=True):
                    assert not probe(axis,angle,ratio)['hits']
                endpoint=probe(axis,sign*last,ratio,True)
                row=dict(axis=axis,direction=sign,AC_to_SC_ratio=ratio,
                    tested_SC_limit_deg=10.,last_contiguous_clear_SC_deg=sign*last,
                    first_collision=first,endpoint=endpoint,return_samples=len(set(out+[0.])),
                    regular_max_joint_step_deg=1.,boundary_max_joint_step_deg=.1)
                rows.append(row)
            candidates=rows[-7:]
            best=max(candidates,key=lambda r:(abs(r['last_contiguous_clear_SC_deg']),-abs(r['endpoint']['AC_angle_deg'])))
            print('axis',axis,'sign',sign,'best ratio',best['AC_to_SC_ratio'],'SC',best['last_contiguous_clear_SC_deg'],flush=True)
    save(folder/'SC_AC_path_screen.json',dict(reference=provenance(reference_file),model=provenance(model_file),
        chain_check=provenance(check_file),generator=provenance(ROOT/'tools/study_channel_first_SC_AC_paths.py'),
        cross_body_pairs=len(pairs),paths=rows,
        scope='42 same-axis constant-ratio kinematic paths, outward to first collision or 10-degree SC study limit, return independently sampled. All cross-body solid pairs included. Local coordinate angles are not anatomical K targets; measured humerus direction is reported separately. No actuation, contact equilibrium, manufacturing or all-pose proof.'))


if __name__=='__main__':main()
