"""Collect current S1-Y2 parts for an isolated infraspinatus fabrication trial.

Copies exact source bytes: no scaling, mesh repair or manufacturing approval.
Hardware meshes are reference envelopes, never printable fasteners.
"""
from pathlib import Path
import csv
import hashlib
import io
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'mujoco_arm/tools'))
import mesh_geom

DEST = ROOT / 'prototype/s1_infraspinatus_trial'
PRINT_ROLES = {
    'tpu:infraspinatus': ('01_TPU', 'muscle'),
    'cap:FJ1500/origin': ('02_rigid_candidates', 'origin_cap'),
    'free_terminal:FJ1500': ('02_rigid_candidates', 'free_terminal'),
    'origin:FJ1500/front_support': ('02_rigid_candidates', 'front_support'),
    'origin:FJ1500/backing_plate': ('02_rigid_candidates', 'backing_plate'),
    'link:FJ1500/muscle_groove_ring': ('02_rigid_candidates', 'muscle_groove_ring'),
    'link:FJ1500/bone_groove_ring': ('02_rigid_candidates', 'bone_groove_ring'),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def guarded_write(path, payload):
    """Do not overwrite a manually edited or stale trial part silently."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_bytes() != payload:
        raise ValueError('Existing trial artifact differs; use a new revision: ' + str(path))
    path.write_bytes(payload)


def main():
    selection_file = ROOT / 'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    assembly_file = ROOT / selection['assembly_report']
    if sha(assembly_file) != selection['assembly_sha256']:
        raise ValueError('Stale selected assembly')
    assembly = read(assembly_file)
    selected = {key: row for key, row in assembly['sources'].items()
                if key == 'tpu:infraspinatus' or 'FJ1500' in key}
    if not PRINT_ROLES.keys() <= selected.keys():
        raise ValueError('Required trial interfaces are missing')
    rows = []
    for key, source in selected.items():
        original = ROOT / source['file']
        if sha(original) != source['sha256']:
            raise ValueError('Changed source: ' + key)
        folder, stem = PRINT_ROLES.get(key, ('03_hardware_reference_NOT_PRINT',
                                             key.replace(':', '_').replace('/', '_')))
        target = DEST / folder / (stem + original.suffix.lower())
        guarded_write(target, original.read_bytes())
        vertices, faces = (mesh_geom.load_obj(str(original)) if original.suffix.lower() == '.obj'
                           else mesh_geom.read_stl(original))
        lower, upper = vertices.min(axis=0), vertices.max(axis=0)
        translation = np.array(source['translation_source_mm'] or [0., 0., 0.])
        rows.append(dict(key=key, quantity=1, group=folder,
                         trial_file=target.relative_to(ROOT).as_posix(),
                         source_file=source['file'], source_sha256=source['sha256'],
                         trial_sha256=sha(target), triangles=len(faces),
                         source_bounds_mm=[lower.tolist(), upper.tolist()],
                         dimensions_mm=(upper-lower).tolist(),
                         source_to_assembly_translation_mm=translation.tolist(),
                         assembly_bounds_mm=[(lower+translation).tolist(), (upper+translation).tolist()],
                         print_orientation_selected=False, slicer_verified=False))
    route_file = assembly_file.parent.parent / 's1_port_path_revision/routes_compact_distal_candidate.json'
    route = next(r for r in read(route_file)['routes'] if r['head'] == 'FJ1500')
    payload = dict(scope='Current-source isolated muscle fabrication trial; not a standalone fixture or shoulder motion pass',
                   selection_source=selection_file.relative_to(ROOT).as_posix(), selection_sha256=sha(selection_file),
                   assembly_source=selection['assembly_report'], assembly_sha256=sha(assembly_file),
                   parts=rows, part_count=len(rows),
                   print_candidate_count=sum(r['group'] != '03_hardware_reference_NOT_PRINT' for r in rows),
                   excluded='Bone, GH, collar, bench and other muscles are not included; full arm assembly remains the selected assembly',
                   steel_wire_nominal_diameter_mm=0.6, internal_bore_nominal_diameter_mm=1.72,
                   active_route_source=route_file.relative_to(ROOT).as_posix(), active_route_sha256=sha(route_file),
                   active_route_source_mm=route,
                   source_geometry_changed=False, print_process_verified=False,
                   actual_drive_range_deg=None, stage_pass=False)
    guarded_write(DEST / 'manifest.json', (json.dumps(payload, ensure_ascii=False, indent=2)+'\n').encode('utf8'))
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(['key', 'group', 'quantity', 'trial_file', 'X_mm', 'Y_mm', 'Z_mm', 'source_sha256'])
    for row in rows:
        writer.writerow([row['key'], row['group'], 1, row['trial_file'],
                         *[f'{v:.4f}' for v in row['dimensions_mm']], row['source_sha256']])
    guarded_write(DEST / 'parts.csv', stream.getvalue().encode('utf-8-sig'))
    subset_sources = {}
    for row in rows:
        subset_sources[row['key']] = dict(file=row['trial_file'], sha256=row['trial_sha256'],
                                         translation_source_mm=row['source_to_assembly_translation_mm'])
    guarded_write(DEST / 'assembly_subset.json', (json.dumps(dict(
        scope='Isolated infraspinatus parts; supporting bone, collar and bench hidden; static only',
        parent_report=selection['assembly_report'], parent_sha256=sha(assembly_file),
        sources=subset_sources, n_instances=len(rows), stage_pass=False,
        actual_drive_range_deg=None), ensure_ascii=False, indent=2)+'\n').encode('utf8'))
    record = DEST / 'measurements_BLANK.json'
    if not record.exists():
        record.write_text(json.dumps(dict(scope='Unmeasured physical trial. Copy to a new record before filling.',
            manifest_sha256=sha(DEST/'manifest.json'), printer=None, nozzle_mm=None, material=None,
            slicer=None, layer_height_mm=None, orientation_description=None,
            steel_wire_actual_diameter_mm=None, terminal_fixation_description=None,
            origin_reaction_fixation_description=None, installation_natural_length_mm=None,
            threading_success=None, independent_caps_fit=None,
            cycles=[dict(cycle=i, input_pull_mm=None, belly_shortening_mm=None, force_N=None,
                         residual_after_60s_mm=None, wire_jam=None, terminal_slip=None, damage=None)
                    for i in range(1, 6)],
            video_or_photo_files=[], full_arm_installed=False, observations=None,
            actual_drive_range_deg=None, stage_pass=False), ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print(json.dumps({k:payload[k] for k in ['part_count', 'print_candidate_count', 'source_geometry_changed', 'stage_pass']}))


if __name__ == '__main__':
    main()
