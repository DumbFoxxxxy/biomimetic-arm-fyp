"""Measure current cap holes against actual bones before designing mounts.

No geometry changes. A hole or a clear screw axis is not a fastening or
load-bearing proof. All three bones are retained in ray observations.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'mujoco_arm/tools'))
import audit_shoulder_retained_joint_fit as A
import shoulder_pads

DEST = ROOT/'reviews/S1_2026-11-03/mechanical_interfaces'
INTERFACES = [
    ('FJ1468', 'origin', 'tpu:deltoid', 'bone:clavicle'),
    ('FJ1513', 'origin', 'tpu:deltoid', 'bone:scapula'),
    ('FJ1467', 'origin', 'tpu:deltoid', 'bone:scapula'),
    ('FJ1468', 'distal', 'tpu:deltoid', 'bone:humerus'),
    ('FJ1507', 'origin', 'tpu:teres_major', 'bone:scapula'),
    ('FJ1507', 'distal', 'tpu:teres_major', 'bone:humerus'),
]


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path): return json.loads(path.read_text(encoding='utf8'))


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    selection_file = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    assembly_file = ROOT/selection['assembly_report']
    assert sha(assembly_file) == selection['assembly_sha256']
    assembly = read(assembly_file)
    # Validate all current inputs, not only the cap subset.
    for row in assembly['sources'].values():
        assert sha(ROOT/row['file']) == row['sha256']
    route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes = read(route_file)['routes']
    needed = {bone for _, _, _, bone in INTERFACES}|{muscle for _, _, muscle, _ in INTERFACES}
    needed |= {'cap:'+head+'/'+end for head, end, _, _ in INTERFACES}
    solids = {key: A.load_checked(ROOT/assembly['sources'][key]['file'],
              assembly['sources'][key]['translation_source_mm'], assembly['sources'][key]['sha256'])[0]
              for key in needed}
    rows = []
    for head, end, muscle, target_bone in INTERFACES:
        route = next(row for row in routes if row['head'] == head)
        port = route['ports'][end]
        normal = A.D.unit(port['normal_source'])
        basis = shoulder_pads.pad_basis(normal)
        # The shared deltoid plate is centered at the seat, not at the
        # separate wire ports. Preserve the existing plate axes and centres.
        inner = np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm'])
        outer = np.array(port['outer_face_source_mm'])-np.array(port['offset_from_seat_mm'])
        cap_key = 'cap:'+head+'/'+end
        holes = []
        for sign in [-1, 1]:
            axis = inner+basis[:, 0]*sign*7.
            start = axis+normal*(float((outer-inner)@normal)+10.)
            finish = start-normal*100.
            observations = {}
            for bone in ['bone:clavicle', 'bone:scapula', 'bone:humerus']:
                crossings = []
                for hit in solids[bone].ray_cast(start, finish):
                    point = np.array(hit.position)
                    depth = float((start-point)@normal)
                    assert abs(depth-hit.distance*100.) < 1e-5
                    crossings.append(dict(point_source_mm=point.tolist(),
                        ray_depth_mm=depth, cap_inner_normal_coordinate_mm=float((point-axis)@normal)))
                observations[bone] = crossings
            target = observations[target_bone]
            # Inspect only the candidate cap-to-near-bone screw span. Do
            # not turn a long ray crossing remote material into a bolt.
            bottom_z = target[1]['cap_inner_normal_coordinate_mm']-1. if len(target)>=2 else -30.
            top_z = float((outer-inner)@normal)+3.
            shaft = A.I._manifold().Manifold.cylinder(top_z-bottom_z,1.5,1.5,64).translate([0.,0.,bottom_z]).transform(np.c_[basis,axis])
            soft_overlap = shaft ^ solids[muscle]
            overlap_bounds = None
            if soft_overlap.volume() > 1e-9:
                local_overlap = soft_overlap.transform(np.c_[basis.T,-basis.T@axis])
                overlap_bounds = list(local_overlap.bounding_box())
            holes.append(dict(sign=sign, axis_source_mm=axis.tolist(), normal_source=normal.tolist(),
                bones=observations, expected_bone=target_bone,
                target_crossing_count=len(target),
                first_target_stock_mm=target[1]['ray_depth_mm']-target[0]['ray_depth_mm'] if len(target)>=2 else None,
                cap_inner_to_first_target_mm=-target[0]['cap_inner_normal_coordinate_mm'] if target else None,
                nominal_M3_shaft_span_relative_inner_mm=[bottom_z, top_z],
                cap_shaft_overlap_mm3=float(A.I.intersection_volume(shaft, solids[cap_key])),
                TPU_shaft_overlap_mm3=float(A.I.intersection_volume(shaft, solids[muscle])),
                TPU_overlap_bounds_relative_hole_axis_mm=overlap_bounds,
                shaft_is_a_diagnostic_not_an_installed_fastener=True))
        rows.append(dict(head=head, end=end, cap=cap_key, muscle=muscle, expected_bone=target_bone,
            nominal_hole_pitch_mm=14., inner_face_source_mm=inner.tolist(), holes=holes,
            associated_origin_hardware=[key for key in assembly['sources'] if key.startswith('origin:'+head+'/')],
            fastening_chain_complete=False,
            classification='existing cap and nominal holes; complete bone-side fastening absent'))
    result = dict(scope=__doc__, selected_assembly=selection['assembly_report'], assembly_sha256=sha(assembly_file),
        current_pointer_sha256=sha(selection_file), routes_sha256=sha(route_file), script_sha256=sha(Path(__file__)),
        selected_sources_checked=len(assembly['sources']), current_sources_unchanged=True,
        interfaces=rows, teres_major_active_stop_model_present=any(
            key.startswith('stop:') and 'FJ1507' in key for key in assembly['sources']),
        no_new_geometry_generated=True, actual_drive_range_deg=None, stage_pass=False)
    (DEST/'unmounted_interface_audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    for row in rows:
        print(row['cap'], 'bone crossings', [h['target_crossing_count'] for h in row['holes']],
              'stock', [h['first_target_stock_mm'] for h in row['holes']],
              'TPU shaft overlaps', [h['TPU_shaft_overlap_mm3'] for h in row['holes']], flush=True)


if __name__ == '__main__': main()
