"""Independent stop fit, placement, bearing and screw first-contact checks.

Does not infer torque, friction, strand integrity, grip force or TPU motion.
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
from verify_s1_infraspinatus_bench import overlap

DEST = ROOT/'prototype/s1_infraspinatus_trial/screw_stop_v1'


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    design = read(DEST/'stop_design.json')
    full = read(DEST/'full_shoulder_check.json')
    bench = read(DEST/'bench_check.json')
    selection = read(ROOT/'design/S1_CURRENT_CANDIDATE.json')
    parent = read(ROOT/selection['assembly_report'])
    assert sha(ROOT/selection['assembly_report']) == selection['assembly_sha256'] == design['source_assembly_sha256']
    for row in parent['sources'].values():
        assert sha(ROOT/row['file']) == row['sha256']
    current = {}
    for key, row in full['sources'].items():
        current[key] = A.load_checked(ROOT/row['file'], row['translation_source_mm'], row['sha256'])[0]
    body = current['stop_screw:FJ1500/body']
    screw = current['stop_screw:FJ1500/M2x3']
    wire = current['wire_active:FJ1500/extended_reference']
    new_keys = ['stop_screw:FJ1500/body', 'stop_screw:FJ1500/M2x3', 'wire_active:FJ1500/extended_reference']
    assert len(current) == 193
    topology = []
    for key in new_keys:
        row = full['sources'][key]
        v, f = mesh_geom.read_stl(ROOT/row['file'])
        _, inverse = np.unique(v, axis=0, return_inverse=True)
        indices = inverse[f]
        assert not np.any(np.diff(np.sort(indices, axis=1), axis=1) == 0)
        _, bnd, nm = mesh_geom.edge_census(indices)
        assert bnd == nm == 0
        topology.append(dict(key=key, boundary=0, nonmanifold=0, degenerate=0))
    pairs = []
    old_keys = [key for key in current if key not in new_keys]
    for i, key in enumerate(new_keys):
        for other in old_keys+new_keys[i+1:]:
            pairs.append(dict(a=key, b=other, overlap_mm3=overlap(current[key], current[other])))
    hits = [row for row in pairs if row['overlap_mm3'] > .001]
    assert len(pairs) == 573
    basis = np.array(design['local_basis_columns_source'])
    n, radial = basis[:, 2], basis[:, 0]
    origin = np.array(design['local_origin_source_mm'])
    # Body+initial-touch screw is inserted from outside the distal cap. The
    # unthreaded wire is not yet inserted, so omit only that absent object.
    insertion = []
    min_gap = 2.
    samples = np.arange(20., -.00001, -.1)
    stationary = {key: value for key, value in current.items() if key not in new_keys}
    for direction, distances in [('in', samples), ('out', samples[::-1])]:
        for value in distances:
            for key in new_keys[:2]:
                moved = current[key].translate(n*value)
                for other, solid in stationary.items():
                    volume = overlap(moved, solid)
                    if volume > .001:
                        insertion.append(dict(direction=direction, shift_mm=float(value), moving=key,
                                              fixed=other, overlap_mm3=volume))
                min_gap = min(min_gap, float(A.I.gap(moved, stationary['free_terminal:FJ1500'], 2.)))
    # Screw backs out radially so the wire can be fed through the body.
    # First-touch placement and 0.1mm retreat are geometric states only.
    radial_hits = []
    for value in [0., .05, .1, .15, .2]:
        moved = screw.translate(radial*value)
        for key, solid in current.items():
            if key == 'stop_screw:FJ1500/M2x3':
                continue
            volume = overlap(moved, solid)
            if volume > .001:
                radial_hits.append(dict(retreat_mm=value, fixed=key, overlap_mm3=volume))
    # A small advancement meets the steel wire. The intentional probe is
    # never included as a valid rigid state or a friction/force prediction.
    compression_probe = overlap(screw.translate(-radial*.05), wire)
    assert compression_probe > .001
    seated_probe = overlap(body.translate(-n*.06), stationary['free_terminal:FJ1500'])
    assert seated_probe > .001
    # Straight tool shank envelope, 0.9mm AF catalog candidate, held away
    # from the nominal socket face: physical engagement remains untested.
    tool_local = A.I._manifold().Manifold.cylinder(20., .9/np.sqrt(3), .9/np.sqrt(3), 6)
    R = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    tool = tool_local.transform(np.c_[R, [3.35, 0., 5.05]]).transform(np.c_[basis, origin])
    tool_hits = [dict(key=key, overlap_mm3=overlap(tool, solid)) for key, solid in current.items()]
    tool_hits = [row for row in tool_hits if row['overlap_mm3'] > .001]
    # The bench uses a different proximal support. Check all new parts
    # against that support; all other bench sources are identical copies.
    bench_adapter_key = 'bench_trial:integral_saddle_and_grip'
    row = bench['sources'][bench_adapter_key]
    adapter = A.load_checked(ROOT/row['file'], row['translation_source_mm'], row['sha256'])[0]
    bench_hits = [dict(key=key, overlap_mm3=overlap(current[key], adapter)) for key in new_keys]
    assert not any(row['overlap_mm3'] > .001 for row in bench_hits)
    # Copy blank record, never infer force from first-touch geometry.
    results = dict(scope=__doc__, design_sha256=sha(DEST/'stop_design.json'),
        full_scene_sha256=sha(DEST/'full_shoulder_check.json'), verifier_sha256=sha(Path(__file__)),
        current_191_sources_unchanged=True, current_pointer_unchanged=True,
        new_parts_topology=topology, new_pair_count=len(pairs), new_pair_checks=pairs, static_hits=hits,
        inherited_190_sources_note='Old-old pairs reuse current parent evidence; this is not a new full-shoulder collision clearance.',
        body_screw_insertion=dict(step_mm=.1, outward_length_mm=20., poses=2*len(samples),
            moving_group=new_keys[:2], minimum_cap_gap_capped_2mm=min_gap, hits=insertion,
            absent_during_insertion=['wire_active:FJ1500/extended_reference']),
        screw_retreat_checks=radial_hits, screw_0p05mm_into_wire_probe_mm3=compression_probe,
        stop_0p06mm_into_cap_probe_mm3=seated_probe,
        nominal_hex_tool_approach_hits=tool_hits, bench_adapter_cross_checks=bench_hits,
        geometry_checks_pass=not (hits or insertion or radial_hits or tool_hits),
        grip_verified=False, torque_Nm=None, grip_force_N=None, actual_drive_range_deg=None, stage_pass=False)
    (DEST/'verification_checks.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf8', newline='\n')
    print(json.dumps(dict(geometry_checks_pass=results['geometry_checks_pass'], new_pairs=len(pairs),
        insertion_poses=2*len(samples), insertion_hits=len(insertion), radial_hits=radial_hits,
        tool_hits=tool_hits, first_touch_probe_mm3=compression_probe)), flush=True)
    if not results['geometry_checks_pass']:
        raise ValueError('Candidate nominal placement or tool approach rejected; inspect report')


if __name__ == '__main__':
    main()
