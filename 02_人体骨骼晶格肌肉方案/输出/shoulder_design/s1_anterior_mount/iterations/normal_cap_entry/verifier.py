"""Read back anterior clamp: topology, contact and conditional local paths.

The deltoid placement remains unverified. Later stages are conditional on
earlier placement; they cannot turn an earlier failed approach into a pass.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,sha,save,load_component,rod
import mesh_geom

DEST=A.DESIGN/'s1_anterior_mount'


def volume(a,b):return float(A.I.intersection_volume(a,b))


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/'tools/build_s1_anterior_mount.py')==r['script_sha256']
    changed=set(r['new_keys'])|set(r['replaced_keys']);unchanged=set(parent['sources'])-changed
    for key in unchanged:assert r['sources'][key]==parent['sources'][key]
    solids={k:load_component(s) for k,s in r['sources'].items()};topology=[]
    for key in sorted(changed):
        v,f=mesh_geom.load_obj(str(ROOT/r['sources'][key]['file']))
        _,inv=np.unique(v,axis=0,return_inverse=True);faces=inv[f]
        _,boundary,nm=mesh_geom.edge_census(faces);xyz=v[f]
        repeated=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
        zero=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        count=sum(p.volume()>1e-8 for p in solids[key].decompose())
        topology.append(dict(key=key,boundary=boundary,nonmanifold=nm,repeated=repeated,zero_area=zero,positive_components=count))
    topology_pass=all(x['boundary']==x['nonmanifold']==x['repeated']==x['zero_area']==0 and x['positive_components']==1 for x in topology)
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    assert sha(route_file)==r['routes_sha256'];routes=read(route_file)['routes']
    route=next(x for x in routes if x['head']=='FJ1468');port=route['ports']['origin']
    n=np.array(port['normal_source']);I=np.array(port['inner_face_source_mm'])
    cap_key='cap:FJ1468/origin';muscle_key='tpu:deltoid';bone_key='bone:clavicle'
    front_key='origin:FJ1468/front_support';rear_key='origin:FJ1468/backing_plate'
    contacts=dict(front_bone_mm3=volume(solids[front_key].translate(-n*.02),solids[bone_key]),
        front_TPU_mm3=volume(solids[front_key].translate(n*.02),solids[muscle_key]),
        rear_bone_mm3=volume(solids[rear_key].translate(n*.02),solids[bone_key]),
        cap_TPU_mm3=volume(solids[cap_key].translate(-n*.1),solids[muscle_key]),
        cap_TPU_gap_capped_2mm=float(A.I.gap(solids[cap_key],solids[muscle_key],2.)))
    contacts_pass=all(contacts[k]>.001 for k in ['front_bone_mm3','front_TPU_mm3','rear_bone_mm3','cap_TPU_mm3'])
    previous=load_component(parent['sources'][muscle_key]);removed=previous-solids[muscle_key]
    wire_guard=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),2.)
    lost_wire_wall=volume(removed,wire_guard)
    fixed=set(solids)-set(r['new_keys'])-{muscle_key,cap_key};stages=[]
    def replay(name,keys,direction,start):
        assert not set(keys)&fixed
        hits=[];offsets=np.arange(start,-.001,-.5)
        for way,values in [('in',offsets),('out',offsets[::-1])]:
            for offset in values:
                for key in keys:
                    posed=solids[key].translate(direction*offset)
                    for target in fixed:
                        v=volume(posed,solids[target])
                        if v>.001:hits.append(dict(way=way,offset_mm=float(offset),moving=key,target=target,volume_mm3=v))
        stages.append(dict(name=name,moving=keys,already_installed=sorted(fixed),step_mm=.5,
            direction_source=direction.tolist(),poses=2*len(offsets),hits=hits,
            conditional_on_prior_stages=True))
        fixed.update(keys);print(name,'hits',len(hits),flush=True)
    replay('front_support_before_deltoid',[front_key],n,20.)
    fixed.add(muscle_key)  # Explicitly unverified transition, not a path pass.
    replay('anterior_cap',[cap_key],n,20.)
    replay('rear_backing',[rear_key],-n,35.)
    for i in range(2):
        replay(f'front_screw_{i}',[f'origin:FJ1468/bolt_{i}',f'origin:FJ1468/front_washer_{i}'],n,50.)
    for i in range(2):
        replay(f'rear_nut_{i}',[f'origin:FJ1468/back_washer_{i}',f'origin:FJ1468/nut_{i}'],-n,12.)
    assert fixed==set(solids)
    tool_hits=[];th=r['dimensions']['cap_thickness_mm']
    for i,hole in enumerate(r['selected_bone_crossings']):
        p=np.array(hole['axis_source_mm'])
        # Shaft envelope ends outside the screw head. Real socket engagement
        # and physical operator access remain unverified.
        tool=rod(p+n*(th+3.9),p+n*(th+25.),1.5)
        for key,solid in solids.items():
            v=volume(tool,solid)
            if v>.001:tool_hits.append(dict(index=i,target=key,volume_mm3=v))
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),
        topology=topology,topology_pass=topology_pass,unchanged_parent_rows=len(unchanged),
        contact_probe_shifts_mm=dict(front=.02,rear=.02,cap=.1),contact_probes=contacts,
        nominal_contact_chain_present=contacts_pass,working_wire_R2_wall_removed_mm3=lost_wire_wall,
        removed_outside_declared_terminal_mm3=r['TPU_removed_outside_terminal_mm3'],
        local_stages=stages,poses=sum(x['poses'] for x in stages),all_objects_accounted_for=True,
        local_approaches_geometry_pass=not any(x['hits'] for x in stages),
        tool_shaft_envelopes=tool_hits,real_tool_engagement_verified=False,
        static_fit_and_contact_pass=r['static_geometry_pass'] and topology_pass and contacts_pass and lost_wire_wall<.001,
        deltoid_placement_transition_verified=False,physical_preload_verified=False,
        whole_module_installation_verified=False,driven_motion_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',result)
    print('contacts',contacts,'static',result['static_fit_and_contact_pass'],'paths',result['local_approaches_geometry_pass'],flush=True)


if __name__=='__main__':main()
