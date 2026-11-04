"""Independently re-read the bench STL, exact topology, fits and nominal assembly paths."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'mujoco_arm/tools'), str(ROOT)]
import audit_shoulder_retained_joint_fit as A
import mesh_geom
from project_runtime import load_mujoco_model

DEST = ROOT/'prototype/s1_infraspinatus_trial/bench_adapter_v1'


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def overlap(a, b):
    # Broad phase is conservative; every potentially overlapping pair uses
    # actual material intersections, including adjacent interfaces.
    first = np.array(a.bounding_box()).reshape(2, 3)
    second = np.array(b.bounding_box()).reshape(2, 3)
    if np.any(first[1] < second[0]) or np.any(second[1] < first[0]):
        return 0.
    return float(A.I.intersection_volume(a, b))


def main():
    construction = read(DEST/'adapter_design.json')
    assembly_file = DEST/'assembly_check.json'
    assembly = read(assembly_file)
    solids = {key: A.load_checked(ROOT/row['file'], row['translation_source_mm'], row['sha256'])[0]
              for key, row in assembly['sources'].items()}
    actual = solids['bench_trial:integral_saddle_and_grip']
    for field in ['new_adapter', 'print_pose_candidate']:
        path = ROOT/construction[field]['file']
        assert sha(path) == construction[field]['sha256']
        vertices, faces = mesh_geom.read_stl(path)
        unique, inverse = np.unique(vertices, axis=0, return_inverse=True)
        indices = inverse[faces]
        assert not np.any(np.diff(np.sort(indices, axis=1), axis=1) == 0)
        _, boundary, nonmanifold = mesh_geom.edge_census(indices)
        assert boundary == nonmanifold == 0, field
    all_pairs = [dict(a=a, b=b, volume_mm3=overlap(solids[a], solids[b]))
                 for a, b in combinations(solids, 2)]
    hits = [row for row in all_pairs if row['volume_mm3'] > .001]
    normal = np.array(construction['local_basis_columns_source'])[:, 2]
    # Plain nominal screw insertion and nut+washer placement only. These do
    # not represent threaded engagement, torque, complete assembly or tools.
    motions = []
    for index in range(2):
        motions.extend([
            (f'front_screw_{index}', [f'origin:FJ1500/bolt_{index}'], 30.),
            (f'rear_washer_nut_{index}', [f'origin:FJ1500/back_washer_{index}', f'origin:FJ1500/nut_{index}'], -40.),
        ])
    paths = []
    for name, moving_keys, start in motions:
        distances = np.linspace(start, 0., round(abs(start)/.5)+1)
        pose_hits = []
        for direction, samples in [('in', distances), ('out', distances[::-1])]:
            for distance in samples:
                for key in moving_keys:
                    moved = solids[key].translate(normal*distance)
                    for fixed_key, fixed in solids.items():
                        if fixed_key in moving_keys:
                            continue  # The same moving group, not a fixed collision exclusion.
                        volume = overlap(moved, fixed)
                        if volume > .001:
                            pose_hits.append(dict(direction=direction, displacement_mm=float(distance),
                                                  moving=key, fixed=fixed_key, volume_mm3=volume))
        paths.append(dict(name=name, moving_keys=moving_keys, step_mm=.5,
                          poses=2*len(distances), hits=pose_hits))
        print(name, 'poses', 2*len(distances), 'hits', len(pose_hits), flush=True)
    bearing_probes = []
    for index in range(2):
        washer = solids[f'origin:FJ1500/back_washer_{index}']
        # Move the fixture back towards the washer by more than the explicit
        # 0.02mm nominal gap, verifying that the intended roof meets material.
        volume = overlap(actual.translate(-normal*.04), washer)
        bearing_probes.append(dict(index=index, adapter_shift_mm=-.04, overlap_mm3=volume))
        assert volume > .001
    trial = read(DEST.parent/'manifest.json')
    route = trial['active_route_source_mm']
    guard = A.D.outer_envelope(np.array(route['full_wire_path_source_mm']), 1.2)
    guard_overlap = overlap(actual, guard)
    assert guard_overlap <= .001
    origin = np.array(construction['local_origin_source_mm'])
    basis = np.array(construction['local_basis_columns_source'])
    tool_probes = []
    for index, sign in enumerate([-1, 1]):
        source = assembly['sources'][f'origin:FJ1500/nut_{index}']
        vertices, _ = mesh_geom.read_stl(ROOT/source['file'])
        local = (vertices+np.array(source['translation_source_mm'])-origin)@basis
        upper = float(local[:, 2].max())-.1
        probe = A.I._manifold().Manifold.cylinder(upper+66., 4.2, 4.2, 64).translate([sign*7., 0., -66.])
        volume = overlap(actual, probe.transform(np.c_[basis, origin]))
        assert volume <= .001
        tool_probes.append(dict(index=index, nominal_probe_diameter_mm=8.4,
                                local_z_interval_mm=[-66., upper], adapter_overlap_mm3=volume,
                                interpretation='rear straight cylindrical access envelope only; not a selected socket or gripping test'))
    model = load_mujoco_model(DEST/'available_parts.xml')
    assert model.ngeom == len(solids) == 26 and model.nq == 0
    selection = read(ROOT/'design/S1_CURRENT_CANDIDATE.json')
    parent = read(ROOT/selection['assembly_report'])
    assert sha(ROOT/selection['assembly_report']) == selection['assembly_sha256'] == construction['current_assembly_sha256']
    for row in parent['sources'].values():
        assert sha(ROOT/row['file']) == row['sha256']
    report = dict(scope=__doc__, source_design_sha256=sha(DEST/'adapter_design.json'),
        source_assembly_sha256=sha(assembly_file), verifier_sha256=sha(Path(__file__)),
        exact_STL_boundary_edges=0, exact_STL_nonmanifold_edges=0, exact_STL_degenerate_triangles=0,
        pair_count=len(all_pairs), all_pair_intersections=all_pairs, pair_hits=hits,
        nominal_screw_and_nut_paths=paths, washer_bearing_material_probes=bearing_probes,
        rear_tool_access_probes=tool_probes,
        active_line_1_2mm_radius_guard_overlap_mm3=guard_overlap,
        current_191_sources_unchanged=True, mujoco_static_geoms=int(model.ngeom), mujoco_DOF=int(model.nq),
        geometric_checks_pass=not hits and not any(row['hits'] for row in paths),
        physical_fixture_verified=False, actual_drive_range_deg=None, stage_pass=False)
    (DEST/'verification_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    if not report['geometric_checks_pass']:
        raise ValueError('Actual geometry or nominal insertion paths need revision; see recorded hits')
    print('26 objects, 325 pairs, exact STL topology, insertion paths and source hashes verified.', flush=True)


if __name__ == '__main__':
    main()
