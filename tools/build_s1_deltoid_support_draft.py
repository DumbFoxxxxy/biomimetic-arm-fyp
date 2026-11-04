"""Directional bone-side release relief on four metal deltoid supports.

Only supports change. Bone, flexible lattice, ports and fastener positions
stay identical. This candidate tests carrying loose supports with the muscle;
temporary holding and later module placement still need their own design.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component,export

DEST=A.DESIGN/'s1_deltoid_support_draft'


def main():
    DEST.mkdir(exist_ok=True)
    analysis_file=A.DESIGN/'s1_deltoid_installation/entry_relief_analysis.json';analysis=read(analysis_file)
    file=A.DESIGN/'s1_deltoid_distal_mount/assembly_check.json';parent=read(file)
    assert analysis['assembly_sha256']==sha(file)
    proof_file=file.parent/'whole_neutral_audit.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['neutral_pose_source_and_wire_pass']
    original={k:load_component(s) for k,s in parent['sources'].items()}
    direction=np.array(analysis['direction_source']);routes_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=read(routes_file)['routes'];assert sha(routes_file)==parent['routes_sha256']
    groups=[('origin:FJ1468/front_support','bone:clavicle','FJ1468','origin'),
            ('origin:FJ1513/front_support','bone:scapula','FJ1513','origin'),
            ('origin:FJ1467/front_support','bone:scapula','FJ1467','origin'),
            ('distal:FJ1468/front_support','bone:humerus','FJ1468','distal')]
    sources=dict(parent['sources']);solids=dict(original);trials=[];changed=set()
    offsets=np.unique(np.r_[np.arange(0.,80.0001,1.),np.arange(0.,3.0001,.1)])
    for key,bone_key,head,end in groups:
        support=original[key];bone=original[bone_key];revised=support
        for offset in offsets:
            obstacle=bone.translate(-direction*offset)
            if A.I.intersection_volume(revised,obstacle)>.0000001:revised=revised-obstacle
        n=A.D.unit(next(r for r in routes if r['head']==head)['ports'][end]['normal_source'])
        before=original[key];raw_removed=float((support-revised).volume())
        needs_change=raw_removed>1e-6
        prepared=revised.as_original().simplify(1e-6) if needs_change else support
        trial=dict(key=key,source_geometry_changed=needs_change,removed_mm3=float((support-prepared).volume()),original_volume_mm3=float(support.volume()),
            new_bone_contact_probe_mm3=float(A.I.intersection_volume(prepared.translate(-n*.02),bone)),
            original_bone_contact_probe_mm3=float(A.I.intersection_volume(before.translate(-n*.02),bone)),
            new_TPU_contact_probe_mm3=float(A.I.intersection_volume(prepared.translate(n*.02),original['tpu:deltoid'])),
            original_TPU_contact_probe_mm3=float(A.I.intersection_volume(before.translate(n*.02),original['tpu:deltoid'])),
            positive_components_over_1e_8=sum(x.volume()>1e-8 for x in prepared.decompose()),
            cleanup_added_mm3=float((prepared-revised).volume()),cleanup_removed_mm3=float((revised-prepared).volume()))
        trials.append(trial)
        save(DEST/'support_trials.json',dict(scope=__doc__,parent_sha256=sha(file),analysis_sha256=sha(analysis_file),
            script_sha256=sha(Path(__file__)),direction_source=direction.tolist(),trials=trials,complete=False,stage_pass=False))
        print(trial,flush=True)
        if trial['positive_components_over_1e_8']!=1 or trial['new_bone_contact_probe_mm3']<=.001 or trial['new_TPU_contact_probe_mm3']<=.001:
            raise ValueError('Directional relief loses connectivity or contact; keep as failed trial')
        if needs_change:
            changed.add(key)
            origin=np.array(sources[key]['translation_source_mm'] or [0.,0.,0.])
            solids[key],sources[key]=export(prepared,DEST/(key.replace(':','_').replace('/','_')+'.obj'),origin)
    hits=[];count=0;keys=list(solids)
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            if a not in changed and b not in changed:continue
            count+=1;v=float(A.I.intersection_volume(solids[a],solids[b]))
            if v>.001:hits.append(dict(a=a,b=b,volume_mm3=v))
    save(DEST/'assembly_check.json',dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),parent_neutral_proof_sha256=sha(proof_file),analysis_sha256=sha(analysis_file),script_sha256=sha(Path(__file__)),
        routes_sha256=sha(routes_file),new_keys=[],replaced_keys=sorted(changed),unchanged_parent_rows=len(sources)-len(changed),
        direction_source=direction.tolist(),support_trials=trials,sampling=dict(coarse_step_mm=1.,fine_step_mm=.1,fine_interval_mm=[0.,3.],span_mm=80.),
        affected_pairs=count,material_intersections=hits,static_geometry_pass=not hits,
        all_bone_TPU_ports_fasteners_unchanged=True,intended_body_owners=parent['intended_body_owners'],
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],adopted_as_current=False,stage_pass=False,
        actual_drive_range_deg=None,not_verified=['complete carrying/holding fixture','full staged bundle and later modules',
            'contact pressure, strength and preload','whole shoulder assembly and actuation']))
    save(DEST/'support_trials.json',dict(scope=__doc__,parent_sha256=sha(file),analysis_sha256=sha(analysis_file),
        script_sha256=sha(Path(__file__)),direction_source=direction.tolist(),trials=trials,complete=True,stage_pass=False))
    print('static hits',hits,flush=True)


if __name__=='__main__':main()
