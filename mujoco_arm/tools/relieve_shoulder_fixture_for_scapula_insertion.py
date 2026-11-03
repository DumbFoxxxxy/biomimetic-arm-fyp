"""Local receiver opening for a verified coarse-path hidden bone barrier.

Retains the 175-object parent and all bones/muscles. Removes receiver material
only around two measured clamp regions using a continuous insertion sweep.
No strength, liner holding or complete installation/actuation acceptance.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_clavicle_mount as K
import build_shoulder_scapula_clamps as C
import build_shoulder_bench_mount as B
from build_shoulder_unified_muscle import local_export

SOURCE = K.DEST/'assembly_check.json'
DEST = A.DESIGN/'s1_fixture_assembly_review/receiver_relief'
FAILURE = DEST.parent/'iterations/coarse_hidden_barrier/selected_path_dense_check.json'


def run():
    DEST.mkdir(parents=True, exist_ok=True)
    parent = json.loads(SOURCE.read_text(encoding='utf8'))
    failure = json.loads(FAILURE.read_text(encoding='utf8'))
    if failure['assembly_sha256'] != A.hashed(SOURCE) or failure['selected_partial_path_samples_clear']:
        raise ValueError('Expected original failed dense path against this exact parent')
    direction = np.array(failure['selected_direction_source'])
    maximum = failure['maximum_translation_mm']
    helper = np.eye(3)[np.argmin(abs(direction))]
    x = A.D.unit(np.cross(helper, direction))
    y = np.cross(direction, x)
    R = np.column_stack([x, y, direction])
    prism = B.box([-.02, -.02, -.02], [.02, .02, maximum+.02]).transform(np.c_[R, np.zeros(3)])
    sources = dict(parent['sources'])
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k, r in sources.items()}
    zones_file = C.DEST/'assembly_check.json'
    zones = json.loads(zones_file.read_text(encoding='utf8'))['zones']
    modifications, zone_records = {}, []
    for tag in ['medial_mid', 'lateral_upper']:
        definition = zones[tag]
        g = np.array(definition['local_origin_source_mm'])
        F = np.array(definition['basis_columns_source'])
        crop = B.box([-8, -10, -15], [18, 10, 15]).transform(np.c_[F, g])
        stock = solids['bone:scapula'] ^ crop
        sweep = stock.minkowski_sum(prism)
        frame = 'fixture:lateral_frame' if tag == 'lateral_upper' else 'fixture:medial_frame'
        old = solids[frame]
        revised = old-sweep
        positive = [p.volume() for p in revised.decompose() if p.volume() > 1e-6]
        if len(positive) != 1:
            (DEST/'rejected_connectedness.json').write_text(json.dumps(dict(zone=tag, positive_volumes_mm3=positive), indent=2), encoding='utf8')
            raise ValueError('Insertion relief disconnects useful frame material')
        origin = np.array(sources[frame]['translation_source_mm'])
        actual, meta = local_export(revised.as_original().simplify(.001), DEST/(tag+'_receiver_frame.stl'), origin)
        if not meta['mesh_checks_pass']:
            (DEST/'rejected_mesh.json').write_text(json.dumps(meta, indent=2), encoding='utf8')
            raise ValueError('Receiver relief readback fails')
        solids[frame], sources[frame] = actual, meta
        modifications[frame] = dict(zone=tag, removed_mm3=float((old-actual).volume()),
                                    added_by_export_mm3=float((actual-old).volume()),
                                    source_sha256=parent['sources'][frame]['sha256'])
        zone_records.append(dict(zone=tag, bone_crop_local_bounds_mm=[-8, -10, -15, 18, 10, 15],
                                  source_origin_mm=g.tolist(), basis_columns_source=F.tolist(),
                                  bone_stock_volume_mm3=float(stock.volume()),
                                  continuous_sweep_volume_mm3=float(sweep.volume())))
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        if a not in modifications and b not in modifications:
            continue
        volume = A.I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=volume))
    route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    wires = []
    for row in json.loads(route_file.read_text(encoding='utf8'))['routes']:
        wire = A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), .3)
        volumes = {k:A.I.intersection_volume(wire, s) for k, s in solids.items()}
        wires.append(dict(head=row['head'], volumes_mm3=volumes,
                          hits_over_threshold=[dict(part=k, volume_mm3=v) for k, v in volumes.items() if v > .001]))
    contacts = []
    for tag, definition in zones.items():
        n = np.array(definition['basis_columns_source'])[:, 2]
        liner = 'clamp:'+tag+'/front_liner'
        frame = 'fixture:lateral_frame' if tag == 'lateral_upper' else 'fixture:medial_frame'
        contacts.append(dict(liner=liner, frame=frame,
                             nominal_gap_capped_1mm=A.I.gap(solids[liner], solids[frame], 1.),
                             normal_probe_translation_source_mm=(n*.05).tolist(),
                             probe_overlap_mm3=A.I.intersection_volume(solids[liner].translate(n*.05), solids[frame])))
    report = dict(scope=__doc__, n_instances=len(solids), sources=sources,
        parent_source=SOURCE.relative_to(A.ROOT).as_posix(), parent_sha256=A.hashed(SOURCE),
        script_sha256=A.hashed(Path(__file__).resolve()), failed_dense_path_sha256=A.hashed(FAILURE),
        clamp_zones_source_sha256=A.hashed(zones_file), routes_sha256=A.hashed(route_file),
        generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
            [A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_scapula_clamps.py']},
        scapula_insertion_relief=dict(direction_source=direction.tolist(), maximum_translation_mm=maximum,
                                     prism_transverse_and_end_clearance_mm=.02,
                                     zones=zone_records, frame_modifications=modifications),
        changed_object_count=2, changed_object_material_intersections=overlaps,
        actuation_wire_checks=wires, front_liner_receiver_contact_checks=contacts,
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners=parent['intended_body_owners'],
        stage_pass=False, actual_drive_range_deg=None,
        not_verified=['dense scapula insertion path must be recomputed after export',
                      'front liner placement/holding and complete skeleton installation',
                      'receiver wall thickness, support strength, preload and physical assembly',
                      'no new flexible muscle/tendon actuation or joint range'])
    (DEST/'assembly_check.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Receiver changes', modifications, flush=True)
    print('New overlaps >.001', [r for r in overlaps if r['volume_mm3'] > .001], flush=True)
    print('Front liner contact probes', [(r['liner'], r['probe_overlap_mm3']) for r in contacts], flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
