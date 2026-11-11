"""Check simulated curved steel centerlines against moving actual rigid parts."""
import sys
import argparse
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_channel_beams'
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder')
    parser.add_argument('--result-name',default='short_response')
    args=parser.parse_args()
    if args.folder:folder=ROOT/args.folder
    assert args.result_name and all(c.isalnum() or c in '_-' for c in args.result_name)
    output='wire_solid_screen.json' if args.result_name=='short_response' else args.result_name+'_solid_screen.json'
    result_file=folder/(args.result_name+'.json')
    archive_file=folder/(args.result_name+'_numerical_rejection.json')
    if not result_file.exists() and archive_file.exists():
        archived=read(archive_file);assert provenance(ROOT/archived['archive']['file'])==archived['archive']
        save(folder/output,dict(response_archive=provenance(archive_file),status='rejected_numerical_response',
            rejection_reasons=archived['rejection_reasons'],scope='No geometry acceptance performed on solver-reset or otherwise numerically invalid poses. Full original bytes retained in verified gzip.'))
        print('Numerically invalid archived response: geometry acceptance skipped',flush=True);return
    result=read(result_file);times=np.array([r['time_s'] for r in result['trace']]);reasons=[]
    if result['warning_counts']:reasons.append('solver_warnings')
    if result['numerical_axial_guard_hit']:reasons.append('numerical_axial_guard')
    if result['maximum_recorded_rotation_increment_deg']>1.:reasons.append('angular_sampling_exceeds_one_degree')
    if len(times)==0 or np.any(np.diff(times)<=0):reasons.append('missing_or_nonmonotonic_trace')
    if len(times) and times[-1]+2*result['timestep_s']<result['duration_s']:reasons.append('duration_not_completed')
    if result.get('termination_reason'):reasons.append(result['termination_reason'])
    if reasons:
        save(folder/output,dict(response=provenance(result_file),status='rejected_numerical_response',rejection_reasons=reasons,
            scope='Numerical validity or sampling preconditions failed. No geometry pass may be inferred.'))
        print('Numerical response rejected before geometry screen:',reasons,flush=True);return
    check_file=ROOT/result['model_check']['file'];check=read(check_file)
    assert provenance(check_file)==result['model_check']
    parent=read(ROOT/check['parent']['file']);assembly_file=ROOT/check['assembly']['file'];assembly=read(assembly_file)
    assert provenance(assembly_file)==check['assembly']
    model=load_mujoco_model(ROOT/check['model']['file']);data=mujoco.MjData(model)
    solids={k:load_component(v) for k,v in assembly['sources'].items()}
    owners={r['source_key']:r['owner'] for r in check.get('geom_owners',parent['geom_owners'])}
    origins={k:np.array(v)*1000 for k,v in parent['neutral_origins_m'].items()}
    bodies={k:mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,k) for k in origins}
    checked=({r['id']:0 for r in check['wire_paths']} if check.get('wire_paths') else {r['head']:0 for r in check['channels']})
    failures={k:None for k in checked}
    external=[r['head'] for r in check['channels'] if r.get('region')=='external_supply']
    wall_first={h:None for h in external};wall_checked={h:0 for h in external}
    rigid_pairs=[(a,b) for a,b in combinations(solids,2) if owners[a]!=owners[b]]
    keys=list(solids);indices={k:i for i,k in enumerate(keys)}
    pair_indices=np.array([(indices[a],indices[b]) for a,b in rigid_pairs],dtype=int)
    rigid_failure=None;rigid_checked=0
    states=list(result['trace'])
    if check.get('wire_paths'):
        mujoco.mj_forward(model,data)
        neutral_wire={r['id']:data.site_xpos[[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,n) for n in r['site_names']]]*1000 for r in check['wire_paths']}
        neutral_beams={r['head']:data.site_xpos[[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,f"{r['head']}__node_{i}") for i in range(r['specification']['segments']+1)]]*1000 for r in check['channels']}
        states.insert(0,dict(time_s=0.,generalized_displacement=np.zeros(model.nv).tolist(),wire_nodes_source_mm=neutral_wire,beam_nodes_source_mm=neutral_beams))
    for state in states:
        mujoco.mj_resetData(model,data);mujoco.mj_integratePos(model,data.qpos,np.array(state['generalized_displacement']),1.)
        mujoco.mj_forward(model,data)
        transforms={k:np.c_[data.xmat[i].reshape(3,3),data.xpos[i]*1000-data.xmat[i].reshape(3,3)@origins[k]] for k,i in bodies.items()}
        placed={k:s.transform(transforms[owners[k]]) for k,s in solids.items()}
        # Conservative broad phase: cache each world AABB once per pose.
        # The epsilon admits near-tangent boxes; exact material intersection
        # still decides every admitted pair. No anatomical pair exclusions.
        bounds=np.array([placed[k].bounding_box() for k in keys]);assert np.isfinite(bounds).all()
        def broadphase(left,right):
            return np.all(left[...,:3]<=right[...,3:]+1e-8,axis=-1)&np.all(right[...,:3]<=left[...,3:]+1e-8,axis=-1)
        if rigid_failure is None:
            hits=[]
            admitted=broadphase(bounds[pair_indices[:,0]],bounds[pair_indices[:,1]])
            for pair_index in np.flatnonzero(admitted):
                a,b=rigid_pairs[pair_index]
                volume=float(A.I.intersection_volume(placed[a],placed[b]))
                if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
            rigid_checked+=1
            if hits:rigid_failure=dict(time_s=state['time_s'],hits=hits)
        for head in checked:
            if failures[head] is not None:continue
            wire=A.D.outer_envelope(np.array(state['wire_nodes_source_mm'][head]),.3);hits=[]
            wire_bounds=np.array(wire.bounding_box())
            for index in np.flatnonzero(broadphase(wire_bounds,bounds)):
                key=keys[index];part=placed[key]
                volume=float(A.I.intersection_volume(wire,part))
                if volume>.001:hits.append(dict(target=key,overlap_mm3=volume))
            checked[head]+=1
            if hits:failures[head]=dict(time_s=state['time_s'],hits=hits)
        for head in external:
            if wall_first[head] is not None:continue
            from screen_channel_first_thorax_elbow_supply import flat_envelope
            envelope=flat_envelope(np.array(state['beam_nodes_source_mm'][head]),1.5)
            contacts=[];bb=np.array(envelope.bounding_box())
            for index in np.flatnonzero(broadphase(bb,bounds)):
                key=keys[index];volume=float(A.I.intersection_volume(envelope,placed[key]))
                if volume>.001:contacts.append(dict(target=key,overlap_mm3=volume))
            wall_checked[head]+=1
            if contacts:wall_first[head]=dict(time_s=state['time_s'],contacts=contacts,
                interpretation='Flexible outer envelope contact requires force/contact and lumen review. No wall reaction was solved; not automatically a physical failure or a pass.')
        if all(v is not None for v in failures.values()):break
    target=folder/output
    if target.exists():
        import hashlib
        previous=target.read_bytes()
        (folder/(target.stem+'_history_'+hashlib.sha256(previous).hexdigest()[:12]+'.json')).write_bytes(previous)
    save(folder/output,dict(response=provenance(result_file),model_check=provenance(check_file),
        assembly=provenance(assembly_file),checked_states=checked,first_interferences=failures,
        wire_radius_mm=.3,rigid_obstacles=len(solids),threshold_mm3=.001,
        rigid_cross_body_pairs=len(rigid_pairs),rigid_states_checked=rigid_checked,
        first_rigid_interference=rigid_failure,
        external_wall_states_checked=wall_checked,external_wall_first_intersections=wall_first,
        neutral_geometry_state_included=bool(check.get('wire_paths')),
        status=('geometric_interference' if rigid_failure or any(v is not None for v in failures.values()) else
                'flexible_contact_unresolved' if any(v is not None for v in wall_first.values()) else 'sampled_geometry_clear'),
        broadphase='Per-pose world AABBs, 1e-8mm conservative overlap epsilon; exact solid intersections decide all admitted pairs.',
        scope='Saved adaptive poses: centerline capsule checks against all listed moving rigid parts, including caps and fasteners. No end-zone or adjacent-body exclusions. Stops per channel at first overlap. Elastic tube wall and original lattice not contact-solved; actual terminal tails and motor supply outside beam are absent. No continuous-motion acceptance.'))
    print('Curved wire checks',checked,'first interferences',failures,flush=True)


if __name__=='__main__':main()
