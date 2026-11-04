"""Independent posterior clamp contact, topology and local staged approaches.

Front support is placed before the deltoid muscle. Muscle placement itself
is a declared unverified transition, not hidden by the local checks.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,sha,save,load_component
import mesh_geom

DEST=A.DESIGN/'s1_deltoid_posterior_mount'


def volume(a,b):return float(A.I.intersection_volume(a,b))


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/'tools/build_s1_deltoid_posterior_mount.py')==r['script_sha256']
    changed=set(r['new_keys'])|set(r['replaced_keys'])
    unchanged=set(r['sources'])-changed
    assert len(r['sources'])==218 and len(unchanged)==205
    for key in unchanged:assert r['sources'][key]==parent['sources'][key]
    for row in parent['sources'].values():assert sha(ROOT/row['file'])==row['sha256']
    solids={k:load_component(row) for k,row in r['sources'].items()}
    topology=[]
    for key in sorted(changed):
        v,f=mesh_geom.load_obj(str(ROOT/r['sources'][key]['file']))
        _,inv=np.unique(v,axis=0,return_inverse=True);faces=inv[f]
        _,boundary,nm=mesh_geom.edge_census(faces);xyz=v[f]
        repeated=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
        zero=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        count=sum(piece.volume()>1e-8 for piece in solids[key].decompose())
        topology.append(dict(key=key,boundary=boundary,nonmanifold=nm,repeated=repeated,zero_area=zero,
            positive_components_over_1e_8_mm3=count))
        assert boundary==nm==repeated==zero==0 and count==1,key
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    assert sha(route_file)==r['routes_sha256']
    route=next(row for row in read(route_file)['routes'] if row['head']=='FJ1513')
    n=np.array(route['ports']['origin']['normal_source'])
    cap_key='cap:FJ1513/origin';muscle_key='tpu:deltoid'
    front_key='origin:FJ1513/front_support';rear_key='origin:FJ1513/backing_plate'
    bone=solids['bone:scapula'];muscle=solids[muscle_key]
    contacts=dict(front_bone_shift_mm=.02,front_bone_mm3=volume(solids[front_key].translate(-n*.02),bone),
        front_TPU_shift_mm=.02,front_TPU_mm3=volume(solids[front_key].translate(n*.02),muscle),
        rear_bone_shift_mm=.02,rear_bone_mm3=volume(solids[rear_key].translate(n*.02),bone),
        cap_TPU_shift_mm=.1,cap_TPU_mm3=volume(solids[cap_key].translate(-n*.1),muscle),
        cap_TPU_gap_capped_2mm=float(A.I.gap(solids[cap_key],muscle,2.)))
    contacts_pass=all(contacts[key]>.001 for key in ['front_bone_mm3','front_TPU_mm3','rear_bone_mm3','cap_TPU_mm3'])
    # Early state, no deltoid muscle or posterior cap yet. Every omitted
    # source is named; the muscle subsequently needs its own placement check.
    fixed=set(solids)-set(r['new_keys'])-{muscle_key,cap_key}
    stages=[]
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
            direction_source=direction.tolist(),poses=2*len(offsets),hits=hits))
        fixed.update(keys);print(name,'hits',len(hits),flush=True)
    replay('front_support_before_deltoid',[front_key],n,20.)
    fixed.add(muscle_key)  # Explicitly unverified muscle placement transition.
    replay('posterior_cap',[cap_key],n,20.)
    replay('rear_backing',[rear_key],-n,30.)
    for index in range(2):
        replay(f'front_screw_{index}',[f'origin:FJ1513/bolt_{index}',f'origin:FJ1513/front_washer_{index}'],n,35.)
    for index in range(2):
        replay(f'rear_nut_{index}',[f'origin:FJ1513/back_washer_{index}',f'origin:FJ1513/nut_{index}'],-n,35.)
    assert fixed==set(solids)
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),topology=topology,
        unchanged_parent_rows=205,contact_probes=contacts,nominal_contact_chain_present=contacts_pass,
        local_stages=stages,all_218_objects_accounted_for=True,
        muscle_placement_transition_verified=False,
        local_approaches_geometry_pass=not any(row['hits'] for row in stages),
        static_fit_and_contact_pass=r['no_new_or_increased_overlap'] and r['no_active_wire_overlap'] and contacts_pass,
        physical_preload_verified=False,actual_TPU_grip_verified=False,
        whole_module_installation_verified=False,driven_motion_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',result)
    print('Contacts',contacts,'local paths',result['local_approaches_geometry_pass'],flush=True)
    if not result['static_fit_and_contact_pass'] or not result['local_approaches_geometry_pass']:
        raise ValueError('Posterior contact or approach remains unresolved; see saved evidence')


if __name__=='__main__':main()
