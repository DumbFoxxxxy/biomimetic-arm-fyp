"""Independently inspect collar shoe readback, local fastening and bare closing.

Local front fasteners/tools are checked against all scene objects. Separate
bare-collar closing excludes as-yet-uninstalled muscles and hardware explicitly;
it cannot establish a full shoulder assembly sequence or real threaded grip.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import (
    ROOT,A,shoulder_pads,rod,sha,read,save,load_component,
)
import mesh_geom

DEST=A.DESIGN/'s1_teres_collar_adapter'


def overlap(a,b):return float(A.I.intersection_volume(a,b))


def sweep(name,moving,solids,fixed_keys,direction,start,step=.5):
    offsets=np.linspace(start,0.,int(abs(start)/step)+1)
    hits=[]
    for way,values in [('in',offsets),('out',offsets[::-1])]:
        for value in values:
            for key in moving:
                posed=solids[key].translate(direction*value)
                for target in fixed_keys:
                    volume=overlap(posed,solids[target])
                    if volume>.001:hits.append(dict(way=way,offset_mm=float(value),moving=key,target=target,volume_mm3=volume))
    print(name,'poses',2*len(offsets),'hits',len(hits),flush=True)
    return dict(name=name,moving=moving,fixed_keys=sorted(fixed_keys),step_mm=step,
                poses=2*len(offsets),direction_source=direction.tolist(),hits=hits)


def main():
    file=DEST/'assembly_check.json';report=read(file)
    parent=read(ROOT/report['parent_report'])
    assert sha(ROOT/report['parent_report'])==report['parent_sha256']
    selection=ROOT/'design/S1_CURRENT_CANDIDATE.json'
    assert sha(selection)==parent['pointer_sha256']
    assert sha(ROOT/report['parent_report'].replace('assembly_check.json','verification_checks.json'))==report['parent_verification_sha256']
    assert sha(ROOT/'tools/build_s1_teres_collar_adapter.py')==report['script_sha256']
    changed=set(report['new_keys'])|set(report['replaced_keys'])
    unchanged=set(report['sources'])-changed
    assert len(unchanged)==198 and len(report['sources'])==204
    for key in unchanged:assert report['sources'][key]==parent['sources'][key]
    for row in parent['sources'].values():assert sha(ROOT/row['file'])==row['sha256']
    solids={key:load_component(row) for key,row in report['sources'].items()}
    topology=[]
    for key in sorted(changed):
        row=report['sources'][key];path=ROOT/row['file']
        v,f=mesh_geom.load_obj(str(path)) if path.suffix=='.obj' else mesh_geom.read_stl(path)
        _,index=np.unique(v,axis=0,return_inverse=True);faces=index[f]
        _,bnd,nm=mesh_geom.edge_census(faces)
        repeated=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
        xyz=v[f];zero_area=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        positive=[float(piece.volume()) for piece in solids[key].decompose() if piece.volume()>1e-8]
        topology.append(dict(key=key,boundary_edges=bnd,nonmanifold_edges=nm,repeated_index_faces=repeated,
                             exact_zero_area_faces=zero_area,positive_components_over_1e_8_mm3=len(positive)))
        assert bnd==nm==repeated==zero_area==0 and len(positive)==1,key
    routes_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    assert sha(routes_file)==report['routes_sha256']
    routes=read(routes_file)['routes'];route=next(row for row in routes if row['head']=='FJ1507')
    I=np.array(route['ports']['distal']['inner_face_source_mm']);n=A.D.unit(route['ports']['distal']['normal_source'])
    F=shoulder_pads.pad_basis(n);point=lambda x,y,z:I+F@np.array([x,y,z],float)
    jaw_key=report['selected_jaw'];jaw=solids[jaw_key];muscle=solids['tpu:teres_major'];cap=solids['cap:FJ1507/distal']
    original_jaw=load_component(parent['sources'][jaw_key]);original_muscle=load_component(parent['sources']['tpu:teres_major'])
    zone=rod(point(0,0,-10.),point(0,0,.2),12.)
    removed=original_muscle-muscle;outside=float((removed-zone).volume())
    lost_jaw=float((original_jaw-jaw).volume())
    for key in parent['sources']:
        if key.startswith('bone:'):assert parent['sources'][key]==report['sources'][key]
    cap_contact=overlap(cap.translate(-n*.02),muscle)
    shoe_contact=overlap(jaw.translate(n*.02),muscle)
    thread_walls=[]
    dim=report['dimensions'];front=dim['shoe_front_normal_mm'];rear=dim['shoe_rear_normal_mm']
    for index,sign in enumerate([-1,1]):
        # Witness an actual surrounding material annulus along the nominal
        # thread length. This does not prove actual threads or pullout force.
        shell=rod(point(sign*7,0,rear),point(sign*7,0,front),2.8)-rod(point(sign*7,0,rear-.1),point(sign*7,0,front+.1),1.61)
        missing=float((shell-jaw).volume())
        thread_walls.append(dict(index=index,annular_witness_missing_mm3=missing,
            nominal_engagement_length_mm=front-rear,actual_threads_verified=False,
            bolt_to_humerus_min_gap_capped_2mm=float(A.I.gap(solids[f'distal:FJ1507/bolt_{index}'],solids['bone:humerus'],2.))))
    fixed=set(solids)-set(report['new_keys']);placements=[]
    for index in range(2):
        moving=[f'distal:FJ1507/bolt_{index}',f'distal:FJ1507/front_washer_{index}']
        placements.append(sweep(f'front_fastener_{index}',moving,solids,fixed,n,25.))
        fixed.update(moving)
    tool_checks=[]
    thickness=float((np.array(route['ports']['distal']['outer_face_source_mm'])-I)@n)
    for index,sign in enumerate([-1,1]):
        tool=rod(point(sign*7,0,thickness+3.85),point(sign*7,0,thickness+28.85),1.5)
        hits=[]
        for key,solid in solids.items():
            v=overlap(tool,solid)
            if v>.001:hits.append(dict(key=key,volume_mm3=v))
        tool_checks.append(dict(index=index,radius_mm=1.5,head_engagement_included=False,hits=hits))
    interface_file=A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    collar_frame=np.array(read(interface_file)['frame_source']['basis_columns'])
    opening=-collar_frame[:,1]
    # This is an explicitly earlier assembly state: rigid skeleton and liners,
    # with collar bolts/tendon pins and muscles not yet installed. Every omitted
    # key is recorded, and the later full-muscle fit remains unverified.
    bare_fixed={key for key in solids if key.startswith(('bone:','GH_','fixture:','bench:','clamp:'))}
    bare_fixed|={'collar:jaw_upper','collar:liner_upper'}
    # Nest the actual lower liner in its own jaw and hand-place the pair.
    # Leaving it pre-positioned on bone makes the jaw scrape its back surface;
    # it is moved here, never removed from the material collision checks.
    closing_moving=[jaw_key,'collar:liner_lower']
    closing=sweep('bare_lower_jaw_and_nested_liner_closing',closing_moving,solids,bare_fixed,opening,30.)
    closing['liner_hand_nested_not_permanently_bonded']=True
    closing['omitted_keys']=sorted(set(solids)-bare_fixed-set(closing_moving))
    closing['all_omitted_parts_require_later_installation']=True
    wires=[]
    for row in routes:
        wire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        hits=[]
        for key in changed:
            v=overlap(wire,solids[key])
            if v>.001:hits.append(dict(key=key,volume_mm3=v))
        wires.append(dict(head=row['head'],hits=hits))
    local_pass=(report['no_new_or_increased_overlap'] and outside<=.001 and lost_jaw<=.001
        and cap_contact>.001 and shoe_contact>.001
        and not any(row['hits'] for row in placements+tool_checks+wires)
        and all(row['annular_witness_missing_mm3']<=.001 for row in thread_walls))
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),
        interface_frame_sha256=sha(interface_file),unchanged_parent_rows=198,default_pointer_unchanged=True,
        all_bone_sources_unchanged=True,exact_topology=topology,
        original_jaw_removed_mm3=lost_jaw,TPU_removed_outside_terminal_mm3=outside,
        cap_contact_probe_mm3=cap_contact,shoe_contact_probe_mm3=shoe_contact,probe_shift_mm=.02,
        contact_probe_is_not_preload_or_strength=True,thread_material_witnesses=thread_walls,
        fastener_placements=placements,front_tool_approaches=tool_checks,bare_collar_closing=closing,
        active_wire_replay=wires,local_interface_geometry_pass=local_pass,
        bare_collar_closing_path_pass=not closing['hits'],complete_shoulder_installation_verified=False,
        metal_thread_grip_verified=False,collar_frictional_retention_verified=False,
        TPU_deformation_verified=False,actual_drive_range_deg=None,stage_pass=False)
    save(DEST/'verification_checks.json',result)
    print({k:result[k] for k in ['local_interface_geometry_pass','bare_collar_closing_path_pass',
        'cap_contact_probe_mm3','shoe_contact_probe_mm3','TPU_removed_outside_terminal_mm3']},flush=True)
    print('Thread witnesses',thread_walls,'Tools',tool_checks,flush=True)
    if not local_pass:raise ValueError('Local collar interface needs revision; saved checks remain authoritative')


if __name__=='__main__':main()
