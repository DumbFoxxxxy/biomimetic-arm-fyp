"""Shared deltoid distal offset clamp with actual humeral bores and backing.

All three wire ports and the muscle attachment centre remain fixed. Changes
are confined to declared terminal geometry, cap and the actual fastener holes.
Nominal geometry only: no physical preload, strength or actuation acceptance.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,export,load_component,shoulder_pads

DEST=A.DESIGN/'s1_deltoid_distal_mount'


def main():
    DEST.mkdir(exist_ok=True)
    study_file=A.DESIGN/'s1_deltoid_distal_study/offset_fastener_study.json';study=read(study_file)
    assert study['script_sha256']==sha(ROOT/'tools/study_s1_deltoid_distal_offset.py')
    choice=study['selected'];assert choice and study['screen_pass']
    file=ROOT/study['parent_report'];parent=read(file);assert sha(file)==study['parent_sha256']
    original={k:load_component(s) for k,s in parent['sources'].items()}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==study['routes_sha256']
    port=next(r for r in routes if r['head']=='FJ1468')['ports']['distal']
    I=np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-np.array(port['inner_face_source_mm']))@n)
    bone=original['bone:humerus'];muscle=original['tpu:deltoid'];cap=original['cap:FJ1468/distal'];mf=A.I._manifold()
    holes=choice['holes'];points=[np.array(h['axis_source_mm']) for h in holes]
    zone=rod(I-n*18.,I+n*.2,18.)
    # Fill obsolete original metal holes without altering the three ports.
    for sign in [-1,1]:
        p=I+F[:,0]*sign*7.
        cap=cap+(rod(p-n*.05,p+n*(thickness+.05),1.71)^rod(I,I+n*thickness,11.))
    front=rod(I-n*7.5,I-n*5.5,11.)
    flange=rod(I-n*5.5,I-n*.03,10.3)
    for p,h in zip(points,holes):
        cap=cap+rod(p,p+n*thickness,4.)
        flange=flange+rod(p-n*5.5,p-n*.03,4.)
        front=front+rod(p-n*7.5,p-n*5.5,4.)+rod(p+n*(h['front_normal_mm']-.2),p-n*5.5,3.8)
    front=front-bone
    for r in routes:
        path=np.array(r['full_wire_path_source_mm'])
        front=front-A.D.outer_envelope(path,2.1)
        flange=flange-A.D.outer_envelope(path,.86)
        cap=cap-A.D.outer_envelope(path,.86)
    front_guard=front
    for axis in range(3):
        for sign in [-1,1]:front_guard=front_guard+front.translate(np.eye(3)[axis]*sign*.01)
    revised=(muscle+flange)-(front_guard^zone)-(bone^zone)
    cap=cap-revised
    lower=min(h['rear_normal_mm'] for h in holes)-2.2;upper=lower+2.
    u=A.D.unit(points[1]-points[0]);rotated=np.stack([u,np.cross(n,u),n],axis=1);midpoint=(points[0]+points[1])*.5
    backing=mf.Manifold.cube([choice['pitch_mm'],8.,2.]).translate([-choice['pitch_mm']/2.,-4.,lower]).transform(np.c_[rotated,midpoint])
    for p,h in zip(points,holes):
        backing=backing+rod(p+n*lower,p+n*upper,4.)+rod(p+n*(upper-.1),p+n*(h['rear_normal_mm']+.2),3.8)
    backing=backing-bone;drilled=bone
    length=next(t for t in [35.,40.,45.,50.,55.,60.] if t>=thickness+.75-(lower-3.15)+.5)
    for p in points:
        drill=rod(p-n*60.,p+n*10.,1.7)
        cap=cap-drill;drilled=drilled-drill;front=front-drill;backing=backing-drill;revised=revised-(drill^zone)
    for r in routes:backing=backing-A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),1.2)
    parts={'distal:FJ1468/front_support':front,'distal:FJ1468/backing_plate':backing}
    for i,p in enumerate(points):
        for name,body in M.m3_parts(p,n,thickness,lower,length).items():parts[f'distal:FJ1468/{name}_{i}']=body
    parts.update({'bone:humerus':drilled,'tpu:deltoid':revised,'cap:FJ1468/distal':cap})
    sources=dict(parent['sources']);actual=dict(original);cleanup=[]
    for key,solid in parts.items():
        prepared=solid.as_original().simplify(1e-6)
        cleanup.append(dict(key=key,added_mm3=float((prepared-solid).volume()),removed_mm3=float((solid-prepared).volume())))
        origin=np.array(sources[key]['translation_source_mm'] or [0.,0.,0.]) if key in sources else I
        actual[key],sources[key]=export(prepared,DEST/(key.replace(':','_').replace('/','_')+'.obj'),origin)
        print('exported',key,flush=True)
    pairs=[];count=0;keys=list(actual)
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            if a not in parts and b not in parts:continue
            count+=1;v=float(A.I.intersection_volume(actual[a],actual[b]))
            previous=float(A.I.intersection_volume(original[a],original[b])) if a in original and b in original else 0.
            if v>.001 or previous>.001:pairs.append(dict(a=a,b=b,volume_mm3=v,parent_mm3=previous,new_or_increased=v>max(.001,previous+.001)))
    wires=[];walls=[];removed=muscle-actual['tpu:deltoid']
    for r in routes:
        corridor=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3);hits=[]
        for key in parts:
            v=float(A.I.intersection_volume(actual[key],corridor))
            if v>.001:hits.append(dict(key=key,volume_mm3=v))
        wires.append(dict(head=r['head'],hits=hits))
        if r['head'] in ['FJ1468','FJ1467','FJ1513']:
            guard=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),2.)
            walls.append(dict(head=r['head'],removed_mm3=float(A.I.intersection_volume(removed,guard))))
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),study_sha256=sha(study_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        new_keys=sorted(set(parts)-set(original)),replaced_keys=sorted(set(parts)&set(original)),unchanged_parent_rows=len(set(original)-set(parts)),
        cap_centre_and_three_ports_unchanged=True,selected_holes=holes,numerical_cleanup=cleanup,
        dimensions=dict(pitch_mm=choice['pitch_mm'],nominal_M3_length_mm=length,cap_thickness_mm=thickness,
            rear_plate_lower_mm=lower,rear_plate_upper_mm=upper,front_support_normal_mm=[-7.5,-5.5],
            flange_normal_mm=[-5.5,-.03],flange_base_radius_mm=10.3,ear_radius_mm=4.,front_TPU_nominal_allowance_mm=.01,
            terminal_zone_radius_mm=18.,terminal_zone_normal_mm=[-18.,.2]),
        bone_removed_mm3=float((bone-actual['bone:humerus']).volume()),TPU_added_mm3=float((actual['tpu:deltoid']-muscle).volume()),
        TPU_removed_mm3=float(removed.volume()),TPU_removed_outside_terminal_mm3=float((removed-zone).volume()),
        three_working_wire_R2_wall_losses=walls,affected_pairs=count,material_intersections=pairs,six_active_wire_corridors=wires,
        static_geometry_pass=not any(x['new_or_increased'] for x in pairs) and not any(x['hits'] for x in wires)
            and all(x['removed_mm3']<.001 for x in walls),
        intended_body_owners={**parent['intended_body_owners'],**{k:'humerus_ideal_attachment' for k in parts if k.startswith('distal:')}},
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],adopted_as_current=False,stage_pass=False,
        actual_drive_range_deg=None,not_verified=['post-export contact and fastening bearing','local and whole assembly',
            'real screw threads, tool engagement, preload and strength','three actual wire stops','TPU deformation and actuation'])
    save(DEST/'assembly_check.json',report)
    print('instances',len(sources),'conflicts',[x for x in pairs if x['new_or_increased']],'wall losses',walls,flush=True)


if __name__=='__main__':main()
