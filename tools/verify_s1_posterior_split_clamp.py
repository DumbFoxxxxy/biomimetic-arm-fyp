"""Review posterior split caps, pocketed nuts and rear-fed screw approaches.

Two muscle-placement transitions remain unverified. Nuts load after the empty
caps and are held until screws engage; this is not full shoulder assembly.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component,rod,shoulder_pads
import mesh_geom

DEST=A.DESIGN/'s1_posterior_split_clamp'


def volume(a,b):return float(A.I.intersection_volume(a,b))


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/r['parent_report'].replace('assembly_check.json','verification_checks.json'))==r['parent_verification_sha256']
    assert sha(ROOT/'tools/build_s1_posterior_split_clamp.py')==r['script_sha256']
    changed=set(r['new_keys'])|set(r['replaced_keys'])
    assert len(r['sources'])==217 and r['unchanged_parent_rows']==211
    for key,row in parent['sources'].items():
        assert sha(ROOT/row['file'])==row['sha256']
        if key not in r['retired_keys']+r['replaced_keys']:assert r['sources'][key]==row
    solids={key:load_component(row) for key,row in r['sources'].items()}
    topology=[]
    for key in sorted(changed):
        v,f=mesh_geom.load_obj(str(ROOT/r['sources'][key]['file']));_,inv=np.unique(v,axis=0,return_inverse=True);fi=inv[f]
        _,boundary,nm=mesh_geom.edge_census(fi);xyz=v[f]
        deg=int(np.any(np.diff(np.sort(fi,axis=1),axis=1)==0,axis=1).sum())
        area0=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        components=sum(piece.volume()>1e-8 for piece in solids[key].decompose())
        assert boundary==nm==deg==area0==0 and components==1,key
        topology.append(dict(key=key,boundary=boundary,nonmanifold=nm,repeated=deg,zero_area=area0,
            positive_components_over_1e_8_mm3=components))
    I=np.array(r['local_origin_source_mm']);F=np.array(r['local_basis_columns']);side,tangent,n=F.T
    caps={suffix:solids['cap:FJ1513/origin_'+suffix] for suffix in ['A','B']}
    cap_all=caps['A']+caps['B'];nut_checks=[];tools=[]
    for index,hole in enumerate(parent['selected_bone_crossings']):
        point=np.array(hole['axis_source_mm']);suffix=r['nut_to_half'][str(index)]
        cap=caps[suffix];nut=solids[f'origin:FJ1513/nut_{index}']
        shell=rod(point+n*.3,point+n*1.7,2.6)-rod(point+n*.2,point+n*1.8,1.75)
        missing=float((shell-cap).volume())
        rotations=[volume(nut.transform(A.I.rotation_about(point,n,sign*5.)),cap) for sign in [-1,1]]
        seat=volume(nut.translate(-n*.05),cap)
        nut_checks.append(dict(index=index,half=suffix,floor_witness_missing_mm3=missing,
            floor_probe_shift_mm=.05,floor_probe_mm3=seat,rotation_probe_deg=5.,rotation_probe_mm3=rotations,
            nominal_seat_and_antirotation_pass=missing<=.001 and seat>.001 and min(rotations)>.001))
        head_face=r['dimensions']['rear_head_base_normal_mm']-3.
        driver=rod(point+n*(head_face-.1),point+n*(head_face-25.1),1.5)
        hits=[]
        for target,solid in solids.items():
            v=volume(driver,solid)
            if v>.001:hits.append(dict(target=target,volume_mm3=v))
        tools.append(dict(index=index,radius_mm=1.5,head_engagement_included=False,hits=hits))
    contact=volume(cap_all.translate(-n*.1),solids['tpu:deltoid'])
    registration=[]
    for axis_name,axis in [('tangent',tangent),('normal',n)]:
        for sign in [-1,1]:
            registration.append(dict(axis=axis_name,sign=sign,shift_mm=.2,
                cap_intersection_mm3=volume(caps['A'].translate(axis*sign*.2),caps['B'])))
    rear_bearing=[]
    for index in range(2):
        bolt=solids[f'origin:FJ1513/bolt_{index}'];washer=solids[f'origin:FJ1513/back_washer_{index}']
        backing=solids['origin:FJ1513/backing_plate']
        rear_bearing.append(dict(index=index,head_to_washer_shift_mm=.12,
            head_to_washer_probe_mm3=volume(bolt.translate(n*.12),washer),washer_to_backing_shift_mm=.2,
            washer_to_backing_probe_mm3=volume(washer.translate(n*.2),backing)))
    install={key for key in solids if key.startswith('origin:FJ1513/')}
    install.update(r['new_keys'])
    muscles={'tpu:deltoid','tpu:infraspinatus'}
    fixed=set(solids)-install-muscles
    initial=sorted(fixed);assert len(initial)==205
    stages=[]
    def replay(name,keys,direction,start,fine=False,step=.5):
        assert not fixed&set(keys)
        samples=np.arange(0.,start+.001,step)
        if fine:samples=np.unique(np.r_[samples,np.arange(0.,3.0001,.1)])
        hits=[]
        for way,values in [('in',samples[::-1]),('out',samples)]:
            for offset in values:
                for key in keys:
                    posed=solids[key].translate(direction*offset)
                    for target in fixed:
                        v=volume(posed,solids[target])
                        if v>.001:hits.append(dict(way=way,offset_mm=float(offset),moving=key,target=target,volume_mm3=v))
        stages.append(dict(name=name,moving=keys,fixed_keys=sorted(fixed),direction_source=direction.tolist(),
            step_mm=step,fine_interval_mm=[0.,3.] if fine else None,fine_step_mm=.1 if fine else None,
            poses=2*len(samples),hits=hits))
        fixed.update(keys);print(name,'poses',2*len(samples),'hits',len(hits),flush=True)
    replay('front_support_before_two_muscles',['origin:FJ1513/front_support'],n,20.)
    fixed.update(muscles)  # Explicitly unverified muscle placement transitions.
    for suffix in r['cap_assembly_order']:
        replay('empty_cap_'+suffix,['cap:FJ1513/origin_'+suffix],np.array(r['cap_entry_directions_source'][suffix]),12.,fine=True)
    for index in range(2):
        replay(f'nut_local_axial_drop_{index}',[f'origin:FJ1513/nut_{index}'],n,3.,step=.1)
    replay('rear_backing',['origin:FJ1513/backing_plate'],-n,30.)
    for index in range(2):
        replay(f'rear_screw_and_washer_{index}',[f'origin:FJ1513/bolt_{index}',f'origin:FJ1513/back_washer_{index}'],-n,25.)
    assert fixed==set(solids)
    passed=(r['static_geometry_pass'] and not any(row['hits'] for row in stages+tools)
        and all(row['nominal_seat_and_antirotation_pass'] for row in nut_checks) and contact>.001
        and all(row['cap_intersection_mm3']>.001 for row in registration)
        and all(row['head_to_washer_probe_mm3']>.001 and row['washer_to_backing_probe_mm3']>.001 for row in rear_bearing))
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),topology=topology,
        unchanged_parent_rows=211,nut_floor_and_rotation_probes=nut_checks,cap_TPU_probe_shift_mm=.1,
        cap_TPU_probe_mm3=contact,rear_driver_envelopes=tools,initial_205_fixed_keys=initial,local_stages=stages,
        registration_transverse_probes=registration,rear_head_washer_bearing_probes=rear_bearing,
        all_217_objects_accounted_for=True,manually_held_nuts_required=True,
        muscle_transitions_unverified=sorted(muscles),operator_hand_space_verified=False,
        local_rigid_approaches_and_seats_pass=passed,permanent_nut_capture_verified=False,
        physical_preload_verified=False,whole_module_installation_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',result)
    print('Local pass',passed,'nuts',nut_checks,'tools',tools,flush=True)
    if not passed:raise ValueError('Posterior split-clamp local paths or seats need revision; see saved evidence')


if __name__=='__main__':main()
