"""Close the teres-major proximal clamp with actual backing and fasteners.

Keeps the original cap/anchors/wire and conforming TPU seat. Adds only a
bounded end-boss contact patch, two scapular bores and a rear backing plate.
No manufacture, preload, strength or muscle-driven motion claim.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'mujoco_arm/tools'), str(ROOT)]
import audit_shoulder_retained_joint_fit as A
import build_shoulder_origin_mounts as M
import mesh_geom
import shoulder_pads
from build_shoulder_shaft_collar import rod

DEST = A.DESIGN/'s1_teres_origin_mount'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path): return json.loads(path.read_text(encoding='utf8'))


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')


def load_component(row):
    """Respect explicit precision for newly generated mechanical geometry."""
    path=ROOT/row['file']
    if row.get('mesh_import_mode')!='exact_indexed_coordinates':
        return A.load_checked(path,row['translation_source_mm'],row['sha256'])[0]
    assert sha(path)==row['sha256']
    vertices,faces=mesh_geom.load_obj(str(path))
    unique,inverse=np.unique(vertices,axis=0,return_inverse=True)
    faces=inverse[faces]
    _,boundary,nonmanifold=mesh_geom.edge_census(faces)
    assert boundary==nonmanifold==0
    return A.I._manifold().Manifold(A.I._manifold().Mesh64(
        unique+np.array(row['translation_source_mm']),faces.astype(np.uint32)))


def export(solid, path, origin):
    local = solid.translate(-origin)
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(local.to_mesh64())
    if path.suffix=='.obj':
        with path.open('w',encoding='ascii',newline='\n') as stream:
            stream.write('# Exact indexed mechanical contact candidate; no physical pass\n')
            for point in vertices:stream.write('v '+' '.join(format(float(value),'.17g') for value in point)+'\n')
            for face in faces:stream.write('f '+' '.join(str(int(value)+1) for value in face)+'\n')
        meta=dict(file=path.relative_to(ROOT).as_posix(),sha256=sha(path),translation_source_mm=origin.tolist(),
                  mesh_import_mode='exact_indexed_coordinates')
        actual=load_component(meta)
        if str(actual.status()).split('.')[-1]!='NoError':raise ValueError('Invalid exact indexed geometry')
    else:
        mesh_geom.write_stl(str(path), vertices, faces, header='Teres origin mechanical candidate; no physical or strength pass')
        actual, _, meta = A.load_checked(path, origin.tolist(), sha(path))
    components = [float(part.volume()) for part in actual.decompose()]
    if sum(v > 1e-8 for v in components) != 1:
        raise ValueError('Disconnected material; no automatic island deletion: '+str(path))
    meta.update(component_volumes_mm3=components,
                roundtrip_relative_volume_error=abs(actual.volume()-solid.volume())/solid.volume())
    if meta['roundtrip_relative_volume_error'] > 1e-5:
        raise ValueError('Roundtrip changes material volume')
    return actual, meta


def main():
    DEST.mkdir(exist_ok=True)
    selection_file = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    parent_file = ROOT/selection['assembly_report']
    assert sha(parent_file) == selection['assembly_sha256']
    parent = read(parent_file)
    sources = dict(parent['sources'])
    original = {}
    for key, source in sources.items():
        original[key] = A.load_checked(ROOT/source['file'], source['translation_source_mm'], source['sha256'])[0]
    routes_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes = read(routes_file)['routes']
    route = next(row for row in routes if row['head'] == 'FJ1507')
    port = route['ports']['origin']
    origin = np.array(port['inner_face_source_mm'])
    normal = A.D.unit(port['normal_source'])
    basis = shoulder_pads.pad_basis(normal)
    thickness = float((np.array(port['outer_face_source_mm'])-origin)@normal)
    audit_file = ROOT/'reviews/S1_2026-11-03/mechanical_interfaces/unmounted_interface_audit.json'
    audit = read(audit_file)
    assert audit['assembly_sha256'] == sha(parent_file)
    row = next(row for row in audit['interfaces'] if row['head'] == 'FJ1507' and row['end'] == 'origin')
    if any(hole['target_crossing_count'] != 2 for hole in row['holes']):
        raise ValueError('Review the bone crossing before constructing a backing plate')
    bone, muscle = original['bone:scapula'], original['tpu:teres_major']
    cap = original['cap:FJ1507/origin']
    # The measured original TPU already conforms to bone. Only the cap
    # boss top is extended to its real contact plane; no working cells cut.
    patch = rod(origin-normal*.3, origin, 5.)
    lumen = A.D.outer_envelope(np.array(route['full_wire_path_source_mm']), .86)
    patch = patch-lumen
    revised_muscle = muscle+patch
    upper_back = min(
        hole['bones']['bone:scapula'][1]['cap_inner_normal_coordinate_mm'] for hole in row['holes'])-.2
    lower_back = upper_back-2.
    mf = A.I._manifold()
    backing = mf.Manifold.cube([14., 8., 2.]).translate([-7., -4., lower_back]).transform(np.c_[basis, origin])
    drills = []
    for hole in row['holes']:
        point = np.array(hole['axis_source_mm'])
        rear_z = hole['bones']['bone:scapula'][1]['cap_inner_normal_coordinate_mm']
        backing = backing+rod(point+normal*lower_back, point+normal*upper_back, 4.)
        # A shallow construction overlap reaches the rear bearing surface
        # while keeping the foot below the TPU wrapping the blade edge.
        backing = backing+rod(point+normal*(upper_back-.1), point+normal*(rear_z+.2), 3.8)
        drills.append(rod(point+normal*10., point-normal*50., 1.7))
    backing = backing-bone
    drilled_bone = bone
    for drill in drills:
        backing = backing-drill
        drilled_bone = drilled_bone-drill
    # Preserve all six active-line corridors, including this muscle's line.
    for wire_route in routes:
        guard = A.D.outer_envelope(np.array(wire_route['full_wire_path_source_mm']), 1.2)
        backing = backing-guard
    minimum = thickness+.75-(lower_back-3.15)+1.
    screw_length = next(length for length in [16.,20.,25.,30.,35.,40.,45.,50.] if length >= minimum)
    new = {'backing_plate': backing}
    for index, hole in enumerate(row['holes']):
        point = np.array(hole['axis_source_mm'])
        for name, body in M.m3_parts(point, normal, thickness, lower_back, screw_length).items():
            new[name+'_'+str(index)] = body
    actual = dict(original)
    muscle_origin = np.array(parent['sources']['tpu:teres_major']['translation_source_mm'] or [0.,0.,0.])
    actual['tpu:teres_major'], sources['tpu:teres_major'] = export(revised_muscle, DEST/'teres_TPU_contact_boss.obj', muscle_origin)
    bone_origin = np.array(parent['sources']['bone:scapula']['translation_source_mm'] or [0.,0.,0.])
    actual['bone:scapula'], sources['bone:scapula'] = export(drilled_bone, DEST/'scapula_teres_two_bores.stl', bone_origin)
    for name, body in new.items():
        key = 'origin:FJ1507/'+name
        actual[key], sources[key] = export(body, DEST/(name+'.stl'), origin)
    changed = {'bone:scapula','tpu:teres_major'}|{'origin:FJ1507/'+name for name in new}
    pairs = []
    pair_count=0
    keys = list(actual)
    for i, key in enumerate(keys):
        for other in keys[i+1:]:
            if key not in changed and other not in changed: continue
            pair_count+=1
            volume = float(A.I.intersection_volume(actual[key], actual[other]))
            old_volume = float(A.I.intersection_volume(original[key], original[other])) if key in original and other in original else 0.
            if volume > .001 or old_volume > .001:
                pairs.append(dict(a=key,b=other,volume_mm3=volume,parent_volume_mm3=old_volume,
                    new_or_increased_over_0_001_mm3=volume>max(.001,old_volume+.001)))
    wire_checks = []
    for wire_route in routes:
        wire = A.D.outer_envelope(np.array(wire_route['full_wire_path_source_mm']), .3)
        hits = []
        for key in changed:
            volume = float(A.I.intersection_volume(actual[key],wire))
            if volume > .001: hits.append(dict(key=key,volume_mm3=volume))
        wire_checks.append(dict(head=wire_route['head'],hits=hits))
    report = dict(scope=__doc__,n_instances=len(sources),sources=sources,
        parent_report=selection['assembly_report'],parent_sha256=sha(parent_file),
        pointer_sha256=sha(selection_file),routes_sha256=sha(routes_file),audit_sha256=sha(audit_file),
        script_sha256=sha(Path(__file__)),new_keys=sorted(changed-{'bone:scapula','tpu:teres_major'}),
        replaced_keys=['bone:scapula','tpu:teres_major'],
        cap_and_wire_anchors_unchanged=True,GH_centre_and_bone_lengths_unchanged=True,
        dimensions=dict(boss_patch_radius_mm=5.,boss_patch_normal_interval_mm=[-.3,0.],
            screw_pitch_mm=14.,scapula_bore_diameter_mm=3.4,
            backing_base_height_mm=2.,backing_upper_normal_mm=upper_back,backing_lower_normal_mm=lower_back,
            nominal_M3_length_mm=screw_length),
        muscle_added_signed_mm3=float((actual['tpu:teres_major']-muscle).volume()),
        muscle_removed_signed_mm3=float((muscle-actual['tpu:teres_major']).volume()),
        scapula_removed_signed_mm3=float((bone-actual['bone:scapula']).volume()),
        scapula_added_signed_mm3=float((actual['bone:scapula']-bone).volume()),
        backing_to_TPU_min_gap_capped_2mm=float(A.I.gap(actual['origin:FJ1507/backing_plate'],actual['tpu:teres_major'],2.)),
        backing_foot_construction_embed_mm=.2,
        backing_relief_rule='shorten rear bearing-foot construction overlap from 1.0 to 0.2mm; do not modify or subtract the wrapping muscle',
        parent_cap_TPU_min_gap_capped_2mm=float(A.I.gap(cap,muscle,2.)),
        cap_TPU_min_gap_capped_2mm=float(A.I.gap(cap,actual['tpu:teres_major'],2.)),
        cap_TPU_0p02_normal_probe_mm3=float(A.I.intersection_volume(cap.translate(-normal*.02),actual['tpu:teres_major'])),
        backing_bone_0p02_normal_probe_mm3=float(A.I.intersection_volume(actual['origin:FJ1507/backing_plate'].translate(normal*.02),actual['bone:scapula'])),
        affected_pairs_checked=pair_count,material_intersections=pairs,actuation_wire_checks=wire_checks,
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners={**parent['intended_body_owners'],**{'origin:FJ1507/'+name:'scapula_ideal_fastened' for name in new}},
        no_new_or_increased_overlap=not any(row['new_or_increased_over_0_001_mm3'] for row in pairs),
        no_active_wire_overlap=not any(row['hits'] for row in wire_checks),
        not_verified=['real screw threads, tool access, tightening/preload or strength',
                      'complete assembly sequence, movement or physical actuation',
                      'the other five unresolved muscle/bone interfaces'],
        adopted_as_current=False,actual_drive_range_deg=None,stage_pass=False)
    save(DEST/'assembly_check.json',report)
    print(json.dumps({key:report[key] for key in ['n_instances','muscle_added_signed_mm3','muscle_removed_signed_mm3',
        'scapula_removed_signed_mm3','cap_TPU_min_gap_capped_2mm','no_new_or_increased_overlap','no_active_wire_overlap']},indent=2),flush=True)
    if not report['no_new_or_increased_overlap'] or not report['no_active_wire_overlap']:
        raise ValueError('New fit rejected; keep the report and revise actual geometry')


if __name__ == '__main__': main()
