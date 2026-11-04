"""Verify unchanged muscle plus released supports in a staged assembly context.

Only the five-part bundle approach is checked. The listed later parts remain
pending, including upper collar hardware and infraspinatus. This is not a
complete installation or a temporary holding fixture verification.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component
import mesh_geom

DEST=A.DESIGN/'s1_deltoid_support_draft'


def main():
    file=DEST/'assembly_check.json';r=read(file);parent_file=ROOT/r['parent_report'];parent=read(parent_file)
    assert sha(parent_file)==r['parent_sha256']
    assert sha(ROOT/'tools/build_s1_deltoid_support_draft.py')==r['script_sha256']
    proof_file=parent_file.parent/'whole_neutral_audit.json';proof=read(proof_file)
    assert sha(proof_file)==r['parent_neutral_proof_sha256'] and proof['neutral_pose_source_and_wire_pass']
    changed=set(r['replaced_keys']);assert len(changed)==2
    for key in set(r['sources'])-changed:assert r['sources'][key]==parent['sources'][key]
    solids={k:load_component(s) for k,s in r['sources'].items()};topology=[]
    for key in changed:
        v,f=mesh_geom.load_obj(str(ROOT/r['sources'][key]['file']));_,inv=np.unique(v,axis=0,return_inverse=True);faces=inv[f]
        _,boundary,nm=mesh_geom.edge_census(faces);xyz=v[f]
        zero=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        repeated=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
        count=sum(p.volume()>1e-8 for p in solids[key].decompose())
        topology.append(dict(key=key,boundary=boundary,nonmanifold=nm,repeated=repeated,zero_area=zero,positive_components=count))
    topology_pass=all(x['boundary']==x['nonmanifold']==x['zero_area']==x['repeated']==0 and x['positive_components']==1 for x in topology)
    base_file=A.DESIGN/'s1_deltoid_installation/support_bundle_entry_study.json';base=read(base_file)
    assert base['assembly_sha256']==r['parent_sha256']
    moving=base['moving_keys'];later=set(base['later_installation_keys'])
    added_late={'collar:jaw_upper','collar:liner_upper','tpu:infraspinatus'}
    added_late|={k for k in solids if k.startswith(('collar:bolt_','collar:top_washer_','collar:bottom_washer_','collar:nut_'))}
    assert added_late<=set(solids)
    later|=added_late;fixed=set(solids)-later-set(moving)
    assert fixed|later|set(moving)==set(solids) and not fixed&later and not set(moving)&later
    direction=np.array(r['direction_source']);offsets=np.unique(np.r_[np.arange(0.,80.0001,1.),np.arange(0.,3.0001,.1)])
    hits=[];poses=0
    for way,values in [('in',offsets[::-1]),('out',offsets)]:
        for distance in values:
            poses+=1
            for key in moving:
                posed=solids[key].translate(direction*distance)
                for target in sorted(fixed):
                    v=float(A.I.intersection_volume(posed,solids[target]))
                    if v>.001:hits.append(dict(way=way,offset_mm=float(distance),moving=key,target=target,volume_mm3=v))
            if abs(distance-round(distance))<1e-8 and int(round(distance))%10==0:
                print(way,distance,'accumulated hits',len(hits),flush=True)
    contact_pass=all(x['new_bone_contact_probe_mm3']>.001 and x['new_TPU_contact_probe_mm3']>.001 for x in r['support_trials'])
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),parent_neutral_proof_sha256=sha(proof_file),
        base_bundle_study_sha256=sha(base_file),topology=topology,topology_pass=topology_pass,unchanged_parent_rows=len(solids)-len(changed),
        all_bones_TPU_ports_fasteners_identical=True,changed_source_pairs_checked=r['affected_pairs'],
        neutral_pose_pass_inherited_and_rechecked=proof['neutral_pose_source_and_wire_pass'] and r['static_geometry_pass'],
        moving_keys=moving,fixed_keys=sorted(fixed),later_installation_keys=sorted(later),
        additional_later_keys=sorted(added_late),all_inventory_keys_accounted_for=True,
        direction_source=direction.tolist(),corridor_length_mm=80.,fine_interval_mm=[0.,3.],fine_step_mm=.1,coarse_step_mm=1.,
        poses=poses,hits=hits,source_intersection_threshold_mm3=.001,
        nominal_support_bone_and_TPU_contact_retained=contact_pass,
        conditional_bundle_placement_geometry_pass=not hits and topology_pass and r['static_geometry_pass'] and contact_pass,
        temporary_holding_fixture_designed=False,later_upper_collar_installation_verified=False,later_infraspinatus_installation_verified=False,
        subsequent_deltoid_fastening_sequence_replayed=False,whole_shoulder_assembly_verified=False,physical_wire_feeding_verified=False,
        driven_motion_verified=False,physical_strength_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',result)
    print('bundle pass',result['conditional_bundle_placement_geometry_pass'],'hits',hits[:10],flush=True)


if __name__=='__main__':main()
