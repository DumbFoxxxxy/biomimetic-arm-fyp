"""Independent readback, interface contact and staged placement checks.

This is a nominal mechanical assembly check, not preload, strength, threads
or TPU-driven movement. Original unresolved shoulder contacts stay recorded.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'mujoco_arm/tools'), str(ROOT)]
import audit_shoulder_retained_joint_fit as A
import mesh_geom
import shoulder_pads
from build_s1_teres_origin_mount import load_component

DEST = A.DESIGN/'s1_teres_distal_mount'


def read(path): return json.loads(path.read_text(encoding='utf8'))


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def overlap(a,b): return float(A.I.intersection_volume(a,b))


def main():
    report_file = DEST/'assembly_check.json'
    report = read(report_file)
    parent = read(ROOT/report['parent_report'])
    assert sha(ROOT/report['parent_report']) == report['parent_sha256']
    selection = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    assert sha(selection) == parent['pointer_sha256']
    assert sha(ROOT/report['parent_report'].replace('assembly_check.json','verification_checks.json')) == report['parent_verification_sha256']
    for row in parent['sources'].values(): assert sha(ROOT/row['file']) == row['sha256']
    keys = set(report['sources'])
    changed = set(report['new_keys'])|set(report['replaced_keys'])
    unchanged = keys-changed
    assert len(keys) == 209 and len(unchanged) == 198
    for key in unchanged: assert report['sources'][key] == parent['sources'][key]
    solids = {key:load_component(row) for key,row in report['sources'].items()}
    topology = []
    for key in sorted(changed):
        row = report['sources'][key]
        path=ROOT/row['file']
        vertices,faces = mesh_geom.load_obj(str(path)) if path.suffix=='.obj' else mesh_geom.read_stl(path)
        _,inv = np.unique(vertices,axis=0,return_inverse=True)
        indexed = inv[faces]
        _,bnd,nm = mesh_geom.edge_census(indexed)
        degenerate = int(np.any(np.diff(np.sort(indexed,axis=1),axis=1)==0,axis=1).sum())
        topology.append(dict(key=key,exact_boundary_edges=bnd,exact_nonmanifold_edges=nm,exact_degenerate_faces=degenerate,
                             positive_components=sum(piece.volume()>1e-8 for piece in solids[key].decompose())))
        assert bnd == nm == degenerate == 0 and topology[-1]['positive_components'] == 1, key
    routes_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    assert sha(routes_file) == report['routes_sha256']
    routes = read(routes_file)['routes']
    route = next(row for row in routes if row['head']=='FJ1507')
    inner = np.array(route['ports']['distal']['inner_face_source_mm'])
    n = A.D.unit(route['ports']['distal']['normal_source'])
    basis = shoulder_pads.pad_basis(n)
    plate = solids['distal:FJ1507/backing_plate']
    bone = solids['bone:humerus']
    foot_contacts = []
    for index,sign in enumerate([-1,1]):
        point = inner+basis[:,0]*(7*sign)
        zone = A.I._manifold().Manifold.cylinder(70.,4.3,4.3,64).translate([0.,0.,-60.]).transform(np.c_[basis,point])
        v = overlap((plate.translate(n*.02)^bone),zone)
        foot_contacts.append(dict(index=index,normal_shift_mm=.02,overlap_mm3=v))
        assert v > .001, 'Each rear foot must meet bone material'
    cap = solids['cap:FJ1507/distal']
    muscle = solids['tpu:teres_major']
    bearing_volume = overlap(cap.translate(-n*.02),muscle)
    assert bearing_volume > .001
    # The complete original cap and TPU are pre-positioned for these steps;
    # this does not validate full-muscle fitting or its deformation.
    fixed = {key:value for key,value in solids.items() if key not in report['new_keys']}
    stages = [('rear_plate',['distal:FJ1507/backing_plate'],-30.)]
    for index in range(2):
        stages.append((f'front_bolt_washer_{index}',[f'distal:FJ1507/bolt_{index}',f'distal:FJ1507/front_washer_{index}'],55.))
    for index in range(2):
        stages.append((f'rear_washer_nut_{index}',[f'distal:FJ1507/back_washer_{index}',f'distal:FJ1507/nut_{index}'],-35.))
    placement = []
    for name,moving,start in stages:
        distances = np.linspace(start,0.,int(abs(start)/.5)+1)
        hits = []
        for direction,samples in [('in',distances),('out',distances[::-1])]:
            for distance in samples:
                for key in moving:
                    posed = solids[key].translate(n*distance)
                    for target,obstacle in fixed.items():
                        volume = overlap(posed,obstacle)
                        if volume>.001:
                            hits.append(dict(direction=direction,displacement_mm=float(distance),moving=key,
                                             fixed=target,volume_mm3=volume))
        placement.append(dict(name=name,moving=moving,step_mm=.5,poses=2*len(distances),hits=hits))
        fixed.update({key:solids[key] for key in moving})
        print(name,2*len(distances),'poses, hits',len(hits),flush=True)
    assert set(fixed)==keys
    # Tool approaches stop before the unmodelled head/nut engagement. Rear
    # annular envelope leaves the actual protruding screw shaft in place.
    tool_checks = []
    thickness = float((np.array(route['ports']['distal']['outer_face_source_mm'])-inner)@n)
    rear_faces = report['dimensions']['rear_seat_faces_mm']
    for index,sign in enumerate([-1,1]):
        point = inner+basis[:,0]*(7*sign)
        z0 = thickness+3.85
        front = A.I._manifold().Manifold.cylinder(25.,1.5,1.5,48).translate([0.,0.,z0]).transform(np.c_[basis,point])
        back_top = rear_faces[index]-3.25
        rear_local = (A.I._manifold().Manifold.cylinder(back_top+75.,4.25,4.25,64)-
                      A.I._manifold().Manifold.cylinder(back_top+75.+2.,1.7,1.7,64).translate([0.,0.,-1.]))
        rear = rear_local.translate([0.,0.,-75.]).transform(np.c_[basis,point])
        for side,tool in [('front',front),('rear',rear)]:
            hits = [dict(key=key,volume_mm3=overlap(tool,solid)) for key,solid in solids.items()]
            hits = [row for row in hits if row['volume_mm3']>.001]
            tool_checks.append(dict(index=index,side=side,straight_approach_hits=hits))
    # Independent replay of all six active steel corridors against every
    # new or replaced object; no anatomy/cap pair is switched off.
    wires = []
    for row in routes:
        wire = A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        hits = [dict(key=key,volume_mm3=overlap(wire,solids[key])) for key in changed]
        wires.append(dict(head=row['head'],hits=[value for value in hits if value['volume_mm3']>.001]))
    old_muscle_row = parent['sources']['tpu:teres_major']
    old_muscle = load_component(old_muscle_row)
    patch_zone = A.I._manifold().Manifold.cylinder(.3,5.,5.,96).translate([0.,0.,-.3]).transform(np.c_[basis,inner])
    unexpected_addition = float(((muscle-old_muscle)-patch_zone).volume())
    lost_stock = float((old_muscle-muscle).volume())
    results = dict(scope=__doc__,assembly_sha256=sha(report_file),verifier_sha256=sha(Path(__file__)),
        current_pointer_unchanged=True,parent_200_sources_unchanged=True,unchanged_source_count=198,
        exact_topology=topology,cap_contact_probe_mm3=bearing_volume,
        cap_contact_area_proxy_mm2=bearing_volume/.02,contact_area_proxy_is_not_pressure_or_capacity=True,
        separate_rear_foot_contacts=foot_contacts,placement_stages=placement,tool_approaches=tool_checks,
        active_wire_replay=wires,unexpected_added_TPU_outside_contact_patch_signed_mm3=unexpected_addition,
        original_TPU_removed_signed_mm3=lost_stock,
        geometry_checks_pass=(report['no_new_or_increased_overlap'] and not any(row['hits'] for row in placement)
            and not any(row['straight_approach_hits'] for row in tool_checks) and not any(row['hits'] for row in wires)
            and unexpected_addition<=.001 and lost_stock<=.001),
        physical_preload_verified=False,complete_shoulder_assembly_verified=False,
        actual_drive_range_deg=None,stage_pass=False,
        installation_context='All other 200 parent objects already pre-positioned; any blocked path must remain reported')
    save_path = DEST/'verification_checks.json'
    save_path.write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(dict(geometry_checks_pass=results['geometry_checks_pass'],foot_contacts=foot_contacts,
                          tool_checks=tool_checks,unexpected_addition_mm3=unexpected_addition,lost_TPU_mm3=lost_stock)),flush=True)
    if not results['geometry_checks_pass']: raise ValueError('Nominal mounting/placement needs revision; see saved checks')


if __name__ == '__main__': main()
