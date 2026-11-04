"""Anterior deltoid clamp with fixed port, rotated bores and real supports.

Independent mechanical candidate, not a frozen S1 or an actuation result.
Original bones' placement, working muscle lattice and wire routes are retained.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,shoulder_pads,sha,read,save,export,load_component

DEST=A.DESIGN/'s1_anterior_mount'


def main():
    DEST.mkdir(exist_ok=True)
    study_file=A.DESIGN/'s1_anterior_mount_study/fastener_study.json';study=read(study_file)
    assert study['script_sha256']==sha(ROOT/'tools/study_s1_anterior_mount.py')
    file=ROOT/study['parent_report'];parent=read(file);assert sha(file)==study['parent_sha256']
    choice=study['selected'];assert choice and choice['preliminary_clear'] and choice['direction']=='front_head'
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==study['routes_sha256']
    port=next(r for r in routes if r['head']=='FJ1468')['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-I)@n)
    angle=np.deg2rad(choice['angle_deg']);u=F[:,0]*np.cos(angle)+F[:,1]*np.sin(angle)
    rotated=np.stack([u,np.cross(n,u),n],axis=1);mf=A.I._manifold()
    old={k:load_component(s) for k,s in parent['sources'].items()}
    bone=old['bone:clavicle'];muscle=old['tpu:deltoid'];cap=old['cap:FJ1468/origin']
    envelope=rod(I,I+n*thickness,11.)
    for sign in [-1,1]:
        p=I+F[:,0]*sign*7.
        fill=rod(p-n*.05,p+n*(thickness+.05),1.71)^envelope
        cap=cap+fill
    # Two real bearing ears connect the rotated screw positions to the cap.
    for h in choice['holes']:
        p=np.array(h['axis_source_mm']);cap=cap+rod(p,p+n*thickness,4.)
    zone=rod(I-n*18.,I+n*.2,12.)
    flange=rod(I-n*7.2,I-n*4.,10.3)
    shim=rod(I-n*.3,I-n*.03,5.)
    for r in routes:
        guard=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.86)
        flange=flange-guard;shim=shim-guard;cap=cap-guard
    revised=muscle+flange+(shim-cap)
    # The ear must sit outside the unchanged muscle and its wire exit.
    cap=cap-revised
    front=rod(I-n*9.2,I-n*7.2,10.8)
    for h in choice['holes']:
        p=np.array(h['axis_source_mm']);z=h['bones']['bone:clavicle'][0]['normal_mm']
        front=front+rod(p+n*(z-.2),p-n*7.2,3.8)
    front=front-bone-revised
    front=front-rod(I-n*16.,I+n*.2,3.)
    rear=[h['bones']['bone:clavicle'][1]['normal_mm'] for h in choice['holes']]
    lower=choice['rear_plate_lower_mm'];upper=lower+2.
    backing=mf.Manifold.cube([14.,8.,2.]).translate([-7.,-4.,lower]).transform(np.c_[rotated,I])
    for h,z in zip(choice['holes'],rear):
        p=np.array(h['axis_source_mm'])
        backing=backing+rod(p+n*lower,p+n*upper,4.)+rod(p+n*(upper-.1),p+n*(z+.2),3.8)
    backing=backing-bone;drilled=bone
    for h in choice['holes']:
        p=np.array(h['axis_source_mm']);drill=rod(p-n*60.,p+n*10.,1.7)
        drilled=drilled-drill;cap=cap-drill;front=front-drill;backing=backing-drill
        revised=revised-(drill^zone)
    for r in routes:
        guard=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),1.2)
        backing=backing-guard;front=front-guard
    parts={'origin:FJ1468/front_support':front,'origin:FJ1468/backing_plate':backing}
    for i,h in enumerate(choice['holes']):
        for name,solid in M.m3_parts(np.array(h['axis_source_mm']),n,thickness,lower,choice['nominal_M3_length_mm']).items():
            parts[f'origin:FJ1468/{name}_{i}']=solid
    parts.update({'bone:clavicle':drilled,'tpu:deltoid':revised,'cap:FJ1468/origin':cap})
    sources=dict(parent['sources']);solids=dict(old);cleanup=[]
    for key,solid in parts.items():
        prepared=solid.as_original().simplify(1e-6)
        cleanup.append(dict(key=key,added_mm3=float((prepared-solid).volume()),removed_mm3=float((solid-prepared).volume())))
        origin=np.array(sources[key]['translation_source_mm'] or [0.,0.,0.]) if key in sources else I
        solids[key],sources[key]=export(prepared,DEST/(key.replace(':','_').replace('/','_')+'.obj'),origin)
        print('exported',key,flush=True)
    pairs=[];count=0;keys=list(solids)
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            if a not in parts and b not in parts:continue
            count+=1;volume=float(A.I.intersection_volume(solids[a],solids[b]))
            previous=float(A.I.intersection_volume(old[a],old[b])) if a in old and b in old else 0.
            if volume>.001 or previous>.001:
                pairs.append(dict(a=a,b=b,volume_mm3=volume,parent_mm3=previous,new_or_increased=volume>max(.001,previous+.001)))
    wires=[]
    for r in routes:
        corridor=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3);hits=[]
        for key in parts:
            volume=float(A.I.intersection_volume(solids[key],corridor))
            if volume>.001:hits.append(dict(key=key,volume_mm3=volume))
        wires.append(dict(head=r['head'],hits=hits))
    save(DEST/'assembly_check.json',dict(scope=__doc__,sources=sources,n_instances=len(sources),
        parent_report=file.relative_to(ROOT).as_posix(),parent_sha256=sha(file),
        study_sha256=sha(study_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        new_keys=sorted(set(parts)-set(old)),replaced_keys=sorted(set(parts)&set(old)),
        unchanged_parent_rows=len(set(old)-set(parts)),cap_centre_and_wire_unchanged=True,
        dimensions=dict(hole_rotation_deg=choice['angle_deg'],screw_pitch_mm=14.,hole_diameter_mm=3.4,
            cap_thickness_mm=thickness,nominal_M3_length_mm=choice['nominal_M3_length_mm'],
            rear_plate_lower_mm=lower,rear_plate_upper_mm=upper,TPU_flange_radius_mm=10.3,
            TPU_flange_normal_mm=[-7.2,-4.],terminal_zone_normal_mm=[-18.,.2],terminal_zone_radius_mm=12.),
        selected_bone_crossings=choice['holes'],numerical_cleanup=cleanup,
        bone_removed_mm3=float((bone-solids['bone:clavicle']).volume()),
        TPU_added_mm3=float((solids['tpu:deltoid']-muscle).volume()),
        TPU_removed_mm3=float((muscle-solids['tpu:deltoid']).volume()),
        TPU_removed_outside_terminal_mm3=float(((muscle-solids['tpu:deltoid'])-zone).volume()),
        affected_pairs=count,material_intersections=pairs,six_active_wire_corridors=wires,
        static_geometry_pass=not any(p['new_or_increased'] for p in pairs) and not any(w['hits'] for w in wires),
        intended_body_owners={**parent['intended_body_owners'],**{k:'clavicle_ideal_fastened' for k in parts if k.startswith('origin:')}},
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None,
        not_verified=['contact and topology independently audited','local and whole assembly paths','physical strength/preload',
                      'muscle loaded deformation and actuation']))
    print('instances',len(sources),'pairs',count,'conflicts',[p for p in pairs if p['new_or_increased']],flush=True)


if __name__=='__main__':main()
