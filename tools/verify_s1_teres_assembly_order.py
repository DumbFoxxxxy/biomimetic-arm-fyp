"""Replay ordered teres stages and terminal capture against real source parts.

This is a sampled nominal geometric sequence with 184 other parts already in
place. Wire route fit is distinct from physical wire feeding. TPU capture
probes are geometric barriers, never strength or actual retention tests.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component,rod,shoulder_pads
import mesh_geom

DEST=A.DESIGN/'s1_teres_assembly_order'


def overlap(a,b):return float(A.I.intersection_volume(a,b))


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/r['parent_report'].replace('assembly_check.json','verification_checks.json'))==r['parent_verification_sha256']
    assert sha(ROOT/'tools/build_s1_teres_assembly_order.py')==r['script_sha256']
    for key,row in parent['sources'].items():
        assert sha(ROOT/row['file'])==row['sha256']
        if key not in r['replaced_keys']:assert row==r['sources'][key]
    solids={key:load_component(row) for key,row in r['sources'].items()}
    assert len(solids)==208
    row=r['sources']['tpu:teres_major'];v,f=mesh_geom.load_obj(str(ROOT/row['file']))
    _,inv=np.unique(v,axis=0,return_inverse=True);faces=inv[f]
    _,bnd,nm=mesh_geom.edge_census(faces);xyz=v[f]
    repeated=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
    zero_area=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
    positive=sum(piece.volume()>1e-8 for piece in solids['tpu:teres_major'].decompose())
    assert bnd==nm==repeated==zero_area==0 and positive==1
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    assert sha(route_file)==r['routes_sha256']
    routes=read(route_file)['routes'];route=next(row for row in routes if row['head']=='FJ1507')
    ports=route['ports'];normals={end:A.D.unit(row['normal_source']) for end,row in ports.items()}
    muscle='tpu:teres_major';wire='wire_active:FJ1507/reference';body='stop:FJ1507/hex_body';screw='stop:FJ1507/M2x3'
    stage_keys={key for key in solids if key.startswith(('origin:FJ1507/','distal:FJ1507/','stop:FJ1507/','wire_active:FJ1507/'))}
    # Only group 0 obstructs muscle entry. Group 1 can remain pre-installed;
    # it cannot be reinserted after the other muscles in this context.
    clamp_hardware={key for key in parent['later_collar_fastener_keys'] if key.endswith('_0')}
    stage_keys|={muscle,'cap:FJ1507/origin','cap:FJ1507/distal_A','cap:FJ1507/distal_B'}|clamp_hardware
    fixed=set(solids)-stage_keys
    initial_fixed=sorted(fixed);assert len(initial_fixed)==184
    stages=[]

    def replay(name,keys,direction,distances,local_offsets=None,install=True):
        assert not set(keys)&fixed
        offsets=local_offsets or {};hits=[]
        baseline={(key,target):overlap(solids[key],solids[target]) for key in keys for target in fixed}
        for way,values in [('in',distances),('out',distances[::-1])]:
            for distance in values:
                for key in keys:
                    posed=solids[key].translate(direction*distance+offsets.get(key,np.zeros(3)))
                    for target in fixed:
                        volume=overlap(posed,solids[target]);base=baseline[key,target]
                        if volume>max(.001,base+.001):
                            hits.append(dict(way=way,offset_mm=float(distance),moving=key,target=target,volume_mm3=volume,zero_pose_mm3=base))
        stages.append(dict(name=name,moving=keys,already_installed=sorted(fixed),offsets_mm=distances.tolist(),
            direction_source=direction.tolist(),initial_offsets={key:value.tolist() for key,value in offsets.items()},
            poses=2*len(distances),hits=hits,
            inherited_zero_overlaps=[dict(moving=a,target=b,volume_mm3=value) for (a,b),value in baseline.items() if value>.001]))
        if install:fixed.update(keys)
        print(name,'poses',2*len(distances),'hits',len(hits),flush=True)

    distances=np.unique(np.r_[np.arange(0.,40.0001,.5),np.arange(0.,10.0001,.1)])[::-1]
    replay('bare_muscle_entry',[muscle],np.array(r['entry_direction_source']),distances)
    replay('origin_cap',['cap:FJ1507/origin'],normals['origin'],np.arange(20.,-.001,-.5))
    cap_side=np.array(parent['local_basis_columns'])[:,0]
    cap_offsets=np.unique(np.r_[np.arange(0.,15.0001,.5),np.arange(0.,3.0001,.1)])[::-1]
    for suffix,sign in [('A',-1),('B',1)]:
        replay('distal_cap_'+suffix,['cap:FJ1507/distal_'+suffix],cap_side*sign,cap_offsets)
    replay('origin_rear_plate',['origin:FJ1507/backing_plate'],-normals['origin'],np.arange(30.,-.001,-.5))
    for index in range(2):
        replay(f'origin_front_fastener_{index}',[f'origin:FJ1507/bolt_{index}',f'origin:FJ1507/front_washer_{index}'],
               normals['origin'],np.arange(35.,-.001,-.5))
    for index in range(2):
        replay(f'origin_rear_fastener_{index}',[f'origin:FJ1507/back_washer_{index}',f'origin:FJ1507/nut_{index}'],
               -normals['origin'],np.arange(35.,-.001,-.5))
    for index in range(2):
        replay(f'distal_front_fastener_{index}',[f'distal:FJ1507/bolt_{index}',f'distal:FJ1507/front_washer_{index}'],
               normals['distal'],np.arange(25.,-.001,-.5))
    collar_frame_file=A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    collar_axis=np.array(read(collar_frame_file)['frame_source']['basis_columns'])[:,1]
    for index in range(1):
        replay(f'collar_front_fastener_{index}',[f'collar:bolt_{index}',f'collar:top_washer_{index}'],collar_axis,np.arange(20.,-.001,-.5))
    for index in range(1):
        replay(f'collar_rear_fastener_{index}',[f'collar:bottom_washer_{index}',f'collar:nut_{index}'],-collar_axis,np.arange(20.,-.001,-.5))
    wire_hits=[]
    for target in fixed:
        volume=overlap(solids[wire],solids[target])
        if volume>.001:wire_hits.append(dict(target=target,volume_mm3=volume))
    fixed.add(wire)  # Route fit only; not a simulated physical feeding step.
    stop_report_file=A.DESIGN/'s1_teres_wire_stop/assembly_check.json'
    stop_frame=np.array(read(stop_report_file)['local_basis_columns']);x=stop_frame[:,0]
    replay('stop_feed_with_backed_screw',[body,screw],normals['distal'],np.arange(12.,-.001,-.5),
           local_offsets={screw:x*.3},install=False)
    fixed.add(body)
    replay('screw_to_first_touch',[screw],x,np.array([.3,.2,.1,0.]))
    assert fixed==set(solids),'Every final part must be accounted for, without whole-pair exclusions'

    capture=[]
    for end in ['origin','distal']:
        I=np.array(ports[end]['inner_face_source_mm']);n=normals[end];F=shoulder_pads.pad_basis(n)
        low=-20. if end=='origin' else -4.
        patch=solids[muscle]^rod(I+n*low,I+n*.2,12.)
        prefix='origin' if end=='origin' else 'distal'
        bolts=[solids[f'{prefix}:FJ1507/bolt_{index}'] for index in range(2)]
        records=[]
        for axis in [0,1]:
            for sign in [-1,1]:
                moved=patch.translate(F[:,axis]*(sign*.3))
                records.append(dict(axis=axis,sign=sign,shift_mm=.3,bolt_overlap_mm3=sum(overlap(moved,bolt) for bolt in bolts)))
        cap=solids['cap:FJ1507/origin'] if end=='origin' else solids['cap:FJ1507/distal_A']+solids['cap:FJ1507/distal_B']
        outward=overlap(patch.translate(n*.1),cap)
        backing=solids['bone:scapula'] if end=='origin' else solids['collar:jaw_lower']
        inward=overlap(patch.translate(-n*.1),backing)
        capture.append(dict(end=end,lateral_probes=records,axial_probe_mm=.1,
            outward_cap_overlap_mm3=outward,inward_backing_overlap_mm3=inward,
            all_probe_barriers_present=all(row['bolt_overlap_mm3']>.001 for row in records) and outward>.001 and inward>.001))
    stop_proof=read(stop_report_file.parent/'verification_checks.json')
    assert stop_proof['assembly_sha256']==sha(stop_report_file)
    tool_solids={key:load_component(row) for key,row in stop_proof['tool_reference_sources'].items()}
    tools=[]
    for key,tool in tool_solids.items():
        hits=[]
        for target,solid in solids.items():
            volume=overlap(tool,solid)
            if volume>.001:hits.append(dict(target=target,volume_mm3=volume))
        tools.append(dict(key=key,hits=hits))
    corridors=[]
    for row in routes:
        corridor=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        volume=overlap(corridor,solids[muscle])
        corridors.append(dict(head=row['head'],TPU_overlap_mm3=volume))
    passed=(not any(row['hits'] for row in stages+tools) and not wire_hits
        and all(row['all_probe_barriers_present'] for row in capture)
        and all(row['TPU_overlap_mm3']<=.001 for row in corridors)
        and r['removed_outside_distal_terminal_mm3']<=.001
        and r['removed_within_wire_outer_2mm_guard_mm3']<=.001
        and max(r['cleanup_added_mm3'],r['cleanup_removed_mm3'])<=.001)
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),all_other_207_sources_unchanged=True,
        topology=dict(boundary=bnd,nonmanifold=nm,repeated=repeated,zero_area=zero_area,positive_components_over_1e_8_mm3=positive),
        initial_184_fixed_keys=initial_fixed,ordered_stages=stages,final_part_count=len(fixed),
        collar_frame_sha256=sha(collar_frame_file),collar_fasteners_installed_later=sorted(clamp_hardware),
        temporary_bone_and_collar_positioning_assumed=True,operator_hand_space_verified=False,
        wire_path_fit_hits=wire_hits,physical_wire_feeding_verified=False,
        terminal_capture_probes=capture,actual_TPU_retention_verified=False,
        final_two_tool_envelopes=tools,six_active_wire_corridors=corridors,
        nominal_sequence_and_capture_geometry_pass=passed,
        whole_shoulder_installation_verified=False,physical_grip_verified=False,TPU_deformation_verified=False,
        driven_motion_verified=False,actual_drive_range_deg=None,stage_pass=False)
    save(DEST/'installation_checks.json',result)
    print('Geometry pass',passed,'Capture',capture,flush=True)
    if not passed:raise ValueError('Module geometry sequence or capture still needs revision; see saved evidence')


if __name__=='__main__':main()
