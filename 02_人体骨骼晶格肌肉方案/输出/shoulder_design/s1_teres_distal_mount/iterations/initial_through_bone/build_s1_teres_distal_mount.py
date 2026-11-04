"""Construct a teres distal through-bone clamp for actual assembly review.

This is a mechanical comparison, not an accepted mount: long transhumeral
fasteners are compared with future collar adaptation. No strength/drive claim.
"""
from pathlib import Path
import json
import sys
import numpy as np

from build_s1_teres_origin_mount import (
    ROOT, A, M, shoulder_pads, rod, sha, read, save, export, load_component,
)

DEST = A.DESIGN/'s1_teres_distal_mount'


def main():
    DEST.mkdir(exist_ok=True)
    parent_file=A.DESIGN/'s1_teres_origin_mount/assembly_check.json'
    parent=read(parent_file)
    proof=read(parent_file.parent/'verification_checks.json')
    assert proof['assembly_sha256']==sha(parent_file) and proof['geometry_checks_pass']
    original={key:load_component(row) for key,row in parent['sources'].items()}
    sources=dict(parent['sources'])
    routes_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=read(routes_file)['routes']
    route=next(row for row in routes if row['head']=='FJ1507')
    port=route['ports']['distal']
    inner=np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm'])
    normal=A.D.unit(port['normal_source']);basis=shoulder_pads.pad_basis(normal)
    thickness=float((np.array(port['outer_face_source_mm'])-np.array(port['inner_face_source_mm']))@normal)
    audit_file=ROOT/'reviews/S1_2026-11-03/mechanical_interfaces/unmounted_interface_audit.json'
    interface=next(row for row in read(audit_file)['interfaces'] if row['head']=='FJ1507' and row['end']=='distal')
    bone=original['bone:humerus'];muscle=original['tpu:teres_major']
    cap=original['cap:FJ1507/distal']
    mf=A.I._manifold()
    # Preserve the original conforming TPU seat; only close the same measured
    # cap-boss gap as at the origin. All working cells remain in the union.
    patch=rod(inner-normal*.3,inner,5.)
    for row in routes:
        patch=patch-A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.86)
    revised_muscle=muscle+patch
    rear_faces=[hole['bones']['bone:humerus'][1]['cap_inner_normal_coordinate_mm'] for hole in interface['holes']]
    upper=min(rear_faces)-.2;lower=upper-2.
    backing=mf.Manifold.cube([14.,8.,2.]).translate([-7.,-4.,lower]).transform(np.c_[basis,inner])
    drilled=bone
    for hole,rear in zip(interface['holes'],rear_faces):
        assert hole['target_crossing_count']==2
        point=np.array(hole['axis_source_mm'])
        backing=backing+rod(point+normal*lower,point+normal*upper,4.)
        backing=backing+rod(point+normal*(upper-.1),point+normal*(rear+.2),3.8)
    backing=backing-bone
    for hole in interface['holes']:
        point=np.array(hole['axis_source_mm'])
        drill=rod(point+normal*10.,point-normal*65.,1.7)
        drilled=drilled-drill;backing=backing-drill
    for row in routes:
        backing=backing-A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),1.2)
    required=thickness+.75-(lower-3.15)+1.
    length=next(value for value in [30.,35.,40.,45.,50.,55.,60.] if value>=required)
    new={'backing_plate':backing}
    for i,hole in enumerate(interface['holes']):
        for name,solid in M.m3_parts(np.array(hole['axis_source_mm']),normal,thickness,lower,length).items():
            new[name+'_'+str(i)]=solid
    actual=dict(original)
    for key,solid,name in [('bone:humerus',drilled,'humerus_teres_bores.stl'),
                           ('tpu:teres_major',revised_muscle,'teres_two_end_contact.obj')]:
        origin=np.array(parent['sources'][key]['translation_source_mm'] or [0.,0.,0.])
        actual[key],sources[key]=export(solid,DEST/name,origin)
    for name,solid in new.items():
        key='distal:FJ1507/'+name
        actual[key],sources[key]=export(solid,DEST/(name+'.stl'),inner)
    changed={'bone:humerus','tpu:teres_major'}|{'distal:FJ1507/'+name for name in new}
    count=0;hits=[]
    keys=list(actual)
    for i,key in enumerate(keys):
        for other in keys[i+1:]:
            if key not in changed and other not in changed:continue
            count+=1
            volume=float(A.I.intersection_volume(actual[key],actual[other]))
            old=float(A.I.intersection_volume(original[key],original[other])) if key in original and other in original else 0.
            if volume>.001 or old>.001:
                hits.append(dict(a=key,b=other,volume_mm3=volume,parent_volume_mm3=old,
                                 new_or_increased=volume>max(.001,old+.001)))
    wires=[]
    for row in routes:
        wire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        obstacles=[]
        for key in changed:
            vol=float(A.I.intersection_volume(wire,actual[key]))
            if vol>.001:obstacles.append(dict(key=key,volume_mm3=vol))
        wires.append(dict(head=row['head'],hits=obstacles))
    observations={}
    for key in ['bone:humerus','tpu:teres_major']:
        uv=[]
        for x,y in [(0,3),(0,-3),(4,0),(-4,0),(-7,2.4),(7,2.4)]:
            point=inner+basis[:,0]*x+basis[:,1]*y
            crossings=actual[key].ray_cast(point+normal*10.,point-normal*65.)
            uv.append(dict(uv_mm=[x,y],normal_crossings_mm=[float((np.array(hit.position)-inner)@normal) for hit in crossings]))
        observations[key]=uv
    no_new=not any(row['new_or_increased'] for row in hits)
    no_wire=not any(row['hits'] for row in wires)
    report=dict(scope=__doc__,n_instances=len(sources),sources=sources,
        parent_report=parent_file.relative_to(ROOT).as_posix(),parent_sha256=sha(parent_file),
        parent_verification_sha256=sha(parent_file.parent/'verification_checks.json'),
        routes_sha256=sha(routes_file),audit_sha256=sha(audit_file),script_sha256=sha(Path(__file__)),
        new_keys=sorted(changed-{'bone:humerus','tpu:teres_major'}),replaced_keys=['bone:humerus','tpu:teres_major'],
        dimensions=dict(screw_pitch_mm=14.,humerus_bore_mm=3.4,nominal_M3_length_mm=length,
                        back_upper_mm=upper,back_lower_mm=lower),
        bone_removed_mm3=float((bone-actual['bone:humerus']).volume()),
        muscle_added_mm3=float((actual['tpu:teres_major']-muscle).volume()),
        muscle_removed_mm3=float((muscle-actual['tpu:teres_major']).volume()),
        cap_TPU_gap_capped_2mm=float(A.I.gap(cap,actual['tpu:teres_major'],2.)),
        cap_TPU_0p02_probe_mm3=float(A.I.intersection_volume(cap.translate(-normal*.02),actual['tpu:teres_major'])),
        seat_rays=observations,affected_pairs_checked=count,material_intersections=hits,
        actuation_wire_checks=wires,no_new_or_increased_overlap=no_new,no_active_wire_overlap=no_wire,
        intended_body_owners={**parent['intended_body_owners'],**{key:'humerus_ideal_fastened' for key in changed if key.startswith('distal:')}},
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        not_verified=['assembly/fastening and tool access','strength of transhumeral bores and long fasteners',
                      'preload, TPU deformation, motion and actual drive','comparison with existing shaft collar adapter'],
        adopted_as_current=False,actual_drive_range_deg=None,stage_pass=False)
    save(DEST/'assembly_check.json',report)
    print(json.dumps({k:v for k,v in report.items() if k in ['n_instances','dimensions','bone_removed_mm3',
        'muscle_added_mm3','muscle_removed_mm3','cap_TPU_gap_capped_2mm','cap_TPU_0p02_probe_mm3',
        'no_new_or_increased_overlap','no_active_wire_overlap']},indent=2),flush=True)
    print('New overlap records',json.dumps([row for row in hits if row['new_or_increased']]),flush=True)


if __name__=='__main__':main()
