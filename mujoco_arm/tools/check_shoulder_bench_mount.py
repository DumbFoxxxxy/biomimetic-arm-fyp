"""Readback and staged screw-access checks for the bench-mount candidate.

Full available rigid/compliant zero-pose pairs are checked without excluding
adjacent bodies. Screw withdrawal paths are geometric, with undeformed TPU.
They are not a full assembly plan or an actuation/strength qualification.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_shoulder_bench_mount as B
import audit_shoulder_retained_joint_fit as A


def withdrawal(screw, obstacles, axis, distance):
    """Move toward the screw head; refine the first sampled collision only."""
    direction = -np.array(axis, float)
    previous = 0.
    maximum_sampled_intersection = 0.
    first = None
    for d in np.arange(1., distance+1e-9, 1.):
        moved = screw.translate(direction*d)
        hits = [(k, A.I.intersection_volume(moved, s)) for k, s in obstacles.items()]
        maximum_sampled_intersection = max(maximum_sampled_intersection, max(v for _, v in hits))
        bad = [(k, v) for k, v in hits if v > B.THRESHOLD]
        if bad:
            for fine in np.arange(previous+.1, d+.001, .1):
                exact = [(k, A.I.intersection_volume(screw.translate(direction*fine), s))
                         for k, s in obstacles.items()]
                blocked = [(k, v) for k, v in exact if v > B.THRESHOLD]
                if blocked:
                    key, volume = max(blocked, key=lambda p:p[1])
                    first = dict(displacement_mm=round(float(fine), 1), part=key,
                                 intersection_mm3=volume,
                                 compliant=key.startswith('tpu:') or key.startswith('collar:liner'))
                    break
            if first is None:
                raise ValueError('Fine scan missed coarse collision')
            break
        previous = float(d)
    return dict(direction_source=direction.tolist(), coarse_step_mm=1.,
                boundary_step_mm=.1, maximum_translation_mm=distance,
                first_interference=first, sampled_path_clear=first is None,
                maximum_coarse_sample_intersection_mm3=maximum_sampled_intersection)


def run():
    source = B.DEST/'assembly_check.json'
    data = json.loads(source.read_text(encoding='utf8'))
    if data['script_sha256'] != A.hashed(Path(B.__file__).resolve()):
        raise ValueError('Changed bench generator')
    if data['parent_sha256'] != A.hashed(B.SOURCE):
        raise ValueError('Changed parent candidate')
    dependencies = data['generation_dependency_sha256']
    if not all(A.hashed(A.ROOT/f) == h for f, h in dependencies.items()):
        raise ValueError('Changed generation dependency')
    solids, census = {}, {}
    for key, row in data['sources'].items():
        body, _, meta = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
        positive = sum(c.volume() > 1e-6 for c in body.decompose())
        if positive != 1:
            raise ValueError('Not one connected material part: '+key)
        solids[key] = body
        census[key] = dict(boundary_edges=meta['boundary_edges'], nonmanifold_edges=meta['nonmanifold_edges'],
                           positive_material_components=positive)
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        volume = A.I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=volume,
                                 compliant=any(k.startswith('tpu:') or k.startswith('collar:liner') for k in [a, b])))
    contacts = []
    for row in data['finite_contact_checks']:
        a, b = row['a'], row['b']
        contacts.append(dict(a=a, b=b, gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                             probe_intersection_mm3=A.I.intersection_volume(
                                 solids[a].translate(row['probe_translation_source_mm']), solids[b])))
    bores = []
    for row in data['shaft_bore_checks']:
        probe = B.rod(row['start_source_mm'], row['end_source_mm'], row['shaft_radius_mm'])
        bores.append(dict(name=row['name'], intersections_mm3={k:A.I.intersection_volume(probe, solids[k])
                                                              for k in row['parts']}))
    parent = json.loads(B.SOURCE.read_text(encoding='utf8'))
    preserved = {}
    for name in ['medial', 'lateral']:
        key = 'fixture:'+name+'_frame'
        row = parent['sources'][key]
        old = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])[0]
        bb = data['base_mounts'][name]['base_bounds_source_mm']
        zone = B.box([bb[0]-8, bb[1]-8, bb[2]-10], [bb[3]+8, bb[4]+8, bb[5]+10])
        delta = (solids[key]-old)+(old-solids[key])
        preserved[key] = dict(base_zone_bounds_source_mm=list(zone.bounding_box()),
                              changed_outside_base_zone_mm3=float((delta-zone).volume()),
                              gap_to_TPU_capped_2mm={k:A.I.gap(solids[key], s, 2.)
                                                    for k, s in solids.items() if k.startswith('tpu:')})
    paths = []
    for key in solids:
        if not key.startswith('bench:') or not key.endswith('/screw'):
            continue
        name = key.split(':', 1)[1].split('/')[0]
        if name.startswith('foot_'):
            spec = data['base_mounts'][name]
        else:
            group, index = name.rsplit('_', 1)
            spec = data['base_mounts'][group]['bores'][int(index)]['fastener']
        obstacles = {k:s for k, s in solids.items() if k != key}
        fixture_only = {k:s for k, s in obstacles.items() if k.startswith(('fixture:', 'bench:'))}
        maximum = spec['nominal_length_mm']+10.
        paths.append(dict(screw=key,
                          before_bones_and_muscles=withdrawal(solids[key], fixture_only, spec['axis_source'], maximum),
                          full_undeformed_context=withdrawal(solids[key], obstacles, spec['axis_source'], maximum)))
    result = dict(scope=__doc__, assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
                  n_instances=len(solids), all_source_hashes_match=True, mesh_census=census,
                  all_pair_checks=len(solids)*(len(solids)-1)//2, material_intersections=overlaps,
                  finite_contacts=contacts, bore_checks=bores, screw_withdrawal_checks=paths,
                  bone_side_preservation_checks=preserved,
                  free_terminal_constraints_required=data['free_terminal_constraints_required'],
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['bone jaws/clavicle clamp and actual workbench anchorage',
                                'nut/washer insertion order, tools, thread engagement and preload',
                                'bench strength/stability, printed tolerances and real assembly',
                                'coupled flexible muscle/tendon motion, loading and actuation range'])
    (B.DEST/'readback_access_checks.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects;', result['all_pair_checks'], 'pairs', flush=True)
    print('Rigid overlaps > threshold', [r for r in overlaps if not r['compliant'] and r['volume_mm3'] > B.THRESHOLD], flush=True)
    print('Soft overlaps > .05 mm3', [r for r in overlaps if r['compliant'] and r['volume_mm3'] > .05], flush=True)
    for row in paths:
        print(row['screw'], 'fixture-only', row['before_bones_and_muscles']['first_interference'],
              'full context', row['full_undeformed_context']['first_interference'], flush=True)
    return result


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
