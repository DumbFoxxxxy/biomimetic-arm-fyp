"""Compare an integral collar jaw extension with a local threaded metal shoe.

No humeral drill or rear fastener. Thread bores are nominal clearance display
envelopes, not manufacturing drill diameters or validated thread engagement.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import (
    ROOT,A,M,shoulder_pads,rod,sha,read,save,export,load_component,
)

DEST=A.DESIGN/'s1_teres_collar_adapter'


def main():
    DEST.mkdir(exist_ok=True)
    parent_file=A.DESIGN/'s1_teres_origin_mount/assembly_check.json'
    parent=read(parent_file)
    proof=read(parent_file.parent/'verification_checks.json')
    assert proof['assembly_sha256']==sha(parent_file) and proof['geometry_checks_pass']
    original={key:load_component(row) for key,row in parent['sources'].items()}
    sources=dict(parent['sources'])
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=read(route_file)['routes'];route=next(r for r in routes if r['head']=='FJ1507')
    port=route['ports']['distal'];I=np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm'])
    n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    cap_thickness=float((np.array(port['outer_face_source_mm'])-np.array(port['inner_face_source_mm']))@n)
    mf=A.I._manifold()
    point=lambda u,v,z:I+F@np.array([u,v,z],float)
    front=-3.2;rear=front-3.
    shoe=rod(point(0,0,rear),point(0,0,front),10.8)
    for sign in [-1,1]:
        shoe=shoe-rod(point(7*sign,0,-12),point(7*sign,0,8),1.6)
    for wire in routes:
        shoe=shoe-A.D.outer_envelope(np.array(wire['full_wire_path_source_mm']),1.2)
    target=point(-9,0,(front+rear)/2)
    choices=[]
    for key in ['collar:jaw_upper','collar:jaw_lower']:
        vertices,_=__import__('mesh_geom').manifold_mesh_to_arrays(original[key].to_mesh64())
        index=np.argmin(np.linalg.norm(vertices-target,axis=1))
        choices.append((float(np.linalg.norm(vertices[index]-target)),key,vertices[index]))
    distance,key,attachment=min(choices,key=lambda row:row[0])
    # The connector must genuinely fuse to the selected jaw. The source
    # nearest vertex is embedded by a small sphere; no disconnected islands
    # are selected or silently dropped at export.
    connector=rod(attachment,target,2.5)+mf.Manifold.sphere(2.5,32).translate(attachment)
    addition=shoe+connector
    revised_jaw=original[key]+addition
    # Restrict terminal relief to a declared end-interface volume. Any
    # collision outside it remains a real, reported unresolved clash.
    terminal_zone=rod(point(0,0,-10),point(0,0,.2),12.)
    guard=addition
    for axis in range(3):
        for sign in [-1,1]:guard=guard+addition.translate(np.eye(3)[axis]*sign*.15)
    clearance=guard^terminal_zone
    muscle=original['tpu:teres_major'];cap=original['cap:FJ1507/distal']
    patch=rod(point(0,0,-.3),point(0,0,0),5.)
    for wire in routes:patch=patch-A.D.outer_envelope(np.array(wire['full_wire_path_source_mm']),.86)
    revised_muscle=(muscle-clearance)+(patch-cap)
    actual=dict(original)
    for part,solid,name in [(key,revised_jaw,'jaw_with_teres_shoe.obj'),
                            ('tpu:teres_major',revised_muscle,'teres_terminal_relief.obj')]:
        offset=np.array(sources[part]['translation_source_mm'] or [0,0,0],float)
        actual[part],sources[part]=export(solid,DEST/name,offset)
    new={}
    for i,sign in enumerate([-1,1]):
        parts=M.m3_parts(point(7*sign,0,0),n,cap_thickness,rear,12.)
        for label in ['bolt','front_washer']:
            name=f'distal:FJ1507/{label}_{i}'
            new[name]=parts[label]
            actual[name],sources[name]=export(parts[label],DEST/f'{label}_{i}.stl',I)
    changed={key,'tpu:teres_major'}|set(new)
    count=0;hits=[]
    keys=list(actual)
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            if a not in changed and b not in changed:continue
            count+=1;v=float(A.I.intersection_volume(actual[a],actual[b]))
            old=float(A.I.intersection_volume(original[a],original[b])) if a in original and b in original else 0.
            if v>.001 or old>.001:
                hits.append(dict(a=a,b=b,volume_mm3=v,parent_volume_mm3=old,new_or_increased=v>max(.001,old+.001)))
    wires=[]
    for row in routes:
        envelope=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        collisions=[]
        for part in changed:
            v=float(A.I.intersection_volume(envelope,actual[part]))
            if v>.001:collisions.append(dict(key=part,volume_mm3=v))
        wires.append(dict(head=row['head'],hits=collisions))
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),
        parent_report=parent_file.relative_to(ROOT).as_posix(),parent_sha256=sha(parent_file),
        parent_verification_sha256=sha(parent_file.parent/'verification_checks.json'),
        script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
        replaced_keys=[key,'tpu:teres_major'],new_keys=sorted(new),
        selected_jaw=key,attachment_source_mm=attachment.tolist(),bridge_nearest_distance_mm=distance,
        dimensions=dict(shoe_front_normal_mm=front,shoe_rear_normal_mm=rear,shoe_radius_mm=10.8,
            thread_nominal='M3 in metal, 3mm nominal axial length; no physical strength claim',
            thread_display_bore_diameter_mm=3.2,display_bore_is_not_manufacturing_drill=True,
            nominal_M3_length_mm=12.,terminal_relief_radius_mm=12.,terminal_relief_normal_mm=[-10.,.2],
            discrete_relief_axis_shift_mm=.15),
        original_bone_unchanged=all(sources[k]==parent['sources'][k] for k in sources if k.startswith('bone:')),
        terminal_TPU_removed_mm3=float((muscle-actual['tpu:teres_major']).volume()),
        terminal_TPU_added_mm3=float((actual['tpu:teres_major']-muscle).volume()),
        removed_outside_terminal_mm3=float(((muscle-actual['tpu:teres_major'])-terminal_zone).volume()),
        jaw_added_mm3=float((actual[key]-original[key]).volume()),
        cap_TPU_probe_mm3=float(A.I.intersection_volume(cap.translate(-n*.02),actual['tpu:teres_major'])),
        affected_pairs_checked=count,material_intersections=hits,actuation_wire_checks=wires,
        no_new_or_increased_overlap=not any(row['new_or_increased'] for row in hits),
        no_active_wire_overlap=not any(row['hits'] for row in wires),
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners={**parent['intended_body_owners'],**{part:'humerus_ideal_fastened' for part in new}},
        adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None,
        not_verified=['metal thread manufacture, axial retention and strength','assembly and tool access',
                      'TPU end grip, deforming muscle contact and actual drive',
                      'collar closing, friction and anti-slip under muscle load'])
    save(DEST/'assembly_check.json',report)
    print({k:report[k] for k in ['n_instances','selected_jaw','bridge_nearest_distance_mm',
        'terminal_TPU_removed_mm3','removed_outside_terminal_mm3','jaw_added_mm3',
        'no_new_or_increased_overlap','no_active_wire_overlap']},flush=True)
    print('New overlaps',[r for r in hits if r['new_or_increased']],flush=True)


if __name__=='__main__':main()
