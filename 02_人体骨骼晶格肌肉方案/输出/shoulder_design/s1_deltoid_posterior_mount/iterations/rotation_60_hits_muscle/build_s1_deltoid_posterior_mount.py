"""Posterior deltoid mount with rotated hole pattern and contoured backing.

Keeps cap centre, original wire port and bones' placement. Old unused TPU
bores remain explicit. No thread/preload/strength or driven-motion claim.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,shoulder_pads,sha,read,save,export,load_component

DEST=A.DESIGN/'s1_deltoid_posterior_mount'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_teres_assembly_order/assembly_check.json';parent=read(file)
    proof_file=file.parent/'installation_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['nominal_sequence_and_capture_geometry_pass']
    study_file=A.DESIGN/'s1_deltoid_origin_holes/hole_axis_study.json';study=read(study_file)
    assert study['parent_sha256']==sha(file)
    group=next(row for row in study['groups'] if row['head']=='FJ1513')
    choice=next(row for row in group['rows'] if row['angle_deg']==60)
    assert choice['both_axes_cross_target_once']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    route=next(row for row in routes if row['head']=='FJ1513');port=route['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    u=F[:,0]*.5+F[:,1]*(np.sqrt(3)/2);v=np.cross(n,u);rotated=np.stack([u,v,n],axis=1)
    thickness=float((np.array(port['outer_face_source_mm'])-I)@n)
    original={key:load_component(row) for key,row in parent['sources'].items()}
    bone=original['bone:scapula'];muscle=original['tpu:deltoid'];cap=original['cap:FJ1513/origin']
    mf=A.I._manifold();envelope=rod(I,I+n*thickness,11.)
    cap_fills=None
    for sign in [-1,1]:
        point=I+F[:,0]*sign*7.
        fill=rod(point-n*.05,point+n*(thickness+.05),1.71)^envelope
        cap_fills=fill if cap_fills is None else cap_fills+fill
    for wire in routes:cap_fills=cap_fills-A.D.outer_envelope(np.array(wire['full_wire_path_source_mm']),.86)
    revised_cap=cap+cap_fills
    terminal_zone=rod(I-n*18.,I+n*.2,12.)
    revised_muscle=muscle;drilled=bone
    rear=[hole['bones']['bone:scapula'][1]['normal_mm'] for hole in choice['holes']]
    upper=min(rear)-.2;lower=upper-2.
    backing=mf.Manifold.cube([14.,8.,2.]).translate([-7.,-4.,lower]).transform(np.c_[rotated,I])
    for hole,z in zip(choice['holes'],rear):
        point=np.array(hole['axis_source_mm'])
        backing=backing+rod(point+n*lower,point+n*upper,4.)+rod(point+n*(upper-.1),point+n*(z+.2),3.8)
    backing=backing-bone
    for hole in choice['holes']:
        point=np.array(hole['axis_source_mm']);drill=rod(point-n*45.,point+n*10.,1.7)
        revised_cap=revised_cap-drill;drilled=drilled-drill;backing=backing-drill
        revised_muscle=revised_muscle-(drill^terminal_zone)
    for wire in routes:backing=backing-A.D.outer_envelope(np.array(wire['full_wire_path_source_mm']),1.2)
    length=next(x for x in [20.,25.,30.,35.,40.] if x>=thickness+.75-(lower-3.15)+.5)
    new={'backing_plate':backing}
    for index,hole in enumerate(choice['holes']):
        for name,solid in M.m3_parts(np.array(hole['axis_source_mm']),n,thickness,lower,length).items():new[f'{name}_{index}']=solid
    sources=dict(parent['sources']);actual=dict(original);cleanup=[]
    for key,solid,name in [('bone:scapula',drilled,'scapula_posterior_pair.obj'),
                           ('cap:FJ1513/origin',revised_cap,'posterior_rotated_cap.obj'),
                           ('tpu:deltoid',revised_muscle,'deltoid_posterior_bores.obj')]:
        prepared=solid.as_original().simplify(1e-6)
        cleanup.append(dict(key=key,added_mm3=float((prepared-solid).volume()),removed_mm3=float((solid-prepared).volume())))
        origin=np.array(sources[key]['translation_source_mm'] or [0,0,0],float)
        actual[key],sources[key]=export(prepared,DEST/name,origin)
    for name,solid in new.items():
        key='origin:FJ1513/'+name
        actual[key],sources[key]=export(solid,DEST/(name+'.obj'),I)
    changed={'bone:scapula','cap:FJ1513/origin','tpu:deltoid'}|{'origin:FJ1513/'+name for name in new}
    pairs=[];count=0;keys=list(actual)
    for index,a in enumerate(keys):
        for b in keys[index+1:]:
            if a not in changed and b not in changed:continue
            count+=1;volume=float(A.I.intersection_volume(actual[a],actual[b]))
            previous=float(A.I.intersection_volume(original[a],original[b])) if a in original and b in original else 0.
            if volume>.001 or previous>.001:
                pairs.append(dict(a=a,b=b,volume_mm3=volume,parent_mm3=previous,
                    new_or_increased=volume>max(.001,previous+.001)))
    wires=[]
    for row in routes:
        corridor=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        hits=[]
        for key in changed:
            volume=float(A.I.intersection_volume(actual[key],corridor))
            if volume>.001:hits.append(dict(key=key,volume_mm3=volume))
        wires.append(dict(head=row['head'],hits=hits))
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),parent_verification_sha256=sha(proof_file),study_sha256=sha(study_file),
        routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        new_keys=sorted(changed-{'bone:scapula','cap:FJ1513/origin','tpu:deltoid'}),
        replaced_keys=['bone:scapula','cap:FJ1513/origin','tpu:deltoid'],all_other_205_parent_rows_unchanged=True,
        dimensions=dict(hole_rotation_deg=60.,screw_pitch_mm=14.,hole_diameter_mm=3.4,nominal_M3_length_mm=length,
            rear_plate_upper_mm=upper,rear_plate_lower_mm=lower,terminal_radius_mm=12.,terminal_normal_mm=[-18.,.2]),
        selected_bone_crossings=choice['holes'],cap_centre_and_wire_unchanged=True,
        original_unused_TPU_holes_retained=True,numerical_cleanup=cleanup,
        bone_removed_mm3=float((bone-actual['bone:scapula']).volume()),
        TPU_removed_mm3=float((muscle-actual['tpu:deltoid']).volume()),
        TPU_removed_outside_terminal_mm3=float(((muscle-actual['tpu:deltoid'])-terminal_zone).volume()),
        cap_added_mm3=float((actual['cap:FJ1513/origin']-cap).volume()),
        cap_removed_mm3=float((cap-actual['cap:FJ1513/origin']).volume()),
        affected_pairs=count,material_intersections=pairs,six_active_wire_corridors=wires,
        no_new_or_increased_overlap=not any(row['new_or_increased'] for row in pairs),
        no_active_wire_overlap=not any(row['hits'] for row in wires),
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners={**parent['intended_body_owners'],**{'origin:FJ1513/'+name:'scapula_ideal_fastened' for name in new}},
        adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None,
        not_verified=['complete mount placement and tool access','preload, strength and bone edge stock after bores',
                      'TPU end retention and loaded deformation','full shoulder assembly or actuation'])
    save(DEST/'assembly_check.json',report)
    print({k:report[k] for k in ['n_instances','bone_removed_mm3','TPU_removed_mm3','affected_pairs',
        'no_new_or_increased_overlap','no_active_wire_overlap']},flush=True)
    print('New conflicts',[row for row in pairs if row['new_or_increased']],flush=True)
    if not report['no_new_or_increased_overlap'] or not report['no_active_wire_overlap']:
        raise ValueError('Posterior mount is not accepted; inspect saved conflicts')


if __name__=='__main__':main()
