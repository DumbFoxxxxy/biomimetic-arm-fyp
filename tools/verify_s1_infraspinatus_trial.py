"""Read-only verification of trial copies and static assembly, not actuation."""
from pathlib import Path
import csv
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'mujoco_arm/tools')]
import mesh_geom
from project_runtime import load_mujoco_model


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    folder = ROOT/'prototype/s1_infraspinatus_trial'
    manifest = read(folder/'manifest.json')
    selection = read(ROOT/'design/S1_CURRENT_CANDIDATE.json')
    assembly = read(ROOT/selection['assembly_report'])
    assert sha(ROOT/selection['assembly_report']) == selection['assembly_sha256'] == manifest['assembly_sha256']
    for key, source in assembly['sources'].items():
        assert sha(ROOT/source['file']) == source['sha256'], key
    csvrows = list(csv.DictReader((folder/'parts.csv').open(encoding='utf-8-sig', newline='')))
    assert len(csvrows) == len(manifest['parts']) == 27
    for row in manifest['parts']:
        source = assembly['sources'][row['key']]
        original, copy = ROOT/source['file'], ROOT/row['trial_file']
        assert original.read_bytes() == copy.read_bytes()
        assert row['source_to_assembly_translation_mm'] == (source['translation_source_mm'] or [0, 0, 0])
        vertices, faces = mesh_geom.read_stl(copy)
        assert len(faces) == row['triangles']
        assert np.array_equal(vertices.max(axis=0)-vertices.min(axis=0), row['dimensions_mm'])
        matching = next(r for r in csvrows if r['key'] == row['key'])
        assert matching['source_sha256'] == sha(copy)
        assert max(abs(float(matching[name])-row['dimensions_mm'][i])
                   for i, name in enumerate(['X_mm', 'Y_mm', 'Z_mm'])) <= .000051
    subset = read(folder/'assembly_subset.json')
    assert set(subset['sources']) == {r['key'] for r in manifest['parts']}
    for key, source in subset['sources'].items():
        row = next(r for r in manifest['parts'] if r['key'] == key)
        assert source['file'] == row['trial_file']
        assert source['sha256'] == row['trial_sha256'] == sha(ROOT/source['file'])
        assert source['translation_source_mm'] == row['source_to_assembly_translation_mm']
    assert sha(ROOT/manifest['active_route_source']) == manifest['active_route_sha256']
    route = next(r for r in read(ROOT/manifest['active_route_source'])['routes'] if r['head'] == 'FJ1500')
    assert route == manifest['active_route_source_mm']
    record = read(folder/'measurements_BLANK.json')
    assert all(all(value is None for key, value in cycle.items() if key != 'cycle') for cycle in record['cycles'])
    assert record['manifest_sha256'] == sha(folder/'manifest.json')
    model = load_mujoco_model(folder/'available_parts.xml')
    assert model.ngeom == 27 and model.nq == 0
    checks = dict(current_source_count=len(assembly['sources']), current_sources_unchanged=True,
                  all_27_copies_dimensions_transforms_and_csv_match=True,
                  subset_sources_and_route_match=True,
                  fabrication_candidates=sum(r['group'] != '03_hardware_reference_NOT_PRINT' for r in manifest['parts']),
                  hardware_reference_only=sum(r['group'] == '03_hardware_reference_NOT_PRINT' for r in manifest['parts']),
                  physical_records_unmeasured=True, mujoco_static_geoms=int(model.ngeom), mujoco_DOF=int(model.nq),
                  scope='Exact copies, dimensions, transforms and static runtime only; not slicing, fixation or actuation.',
                  manifest_sha256=sha(folder/'manifest.json'), verifier_sha256=sha(Path(__file__)),
                  actual_drive_range_deg=None, stage_pass=False)
    (folder/'verification_checks.json').write_text(json.dumps(checks, indent=2)+'\n', encoding='utf8', newline='\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
