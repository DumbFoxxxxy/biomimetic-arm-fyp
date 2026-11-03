"""Independent readback/contact and straight withdrawal probes of scapula clamps.

The undeformed liners are obstacles for these geometric probes; contact does
not establish preload, frictional restraint, escape resistance or strength.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_shoulder_scapula_clamps as C
import audit_shoulder_retained_joint_fit as A
from build_shoulder_shaft_collar import rod


def translate_path(body, direction, obstacles, limit=20.):
    previous = 0.
    for d in np.arange(1., limit+.01, 1.):
        hits = [(k, A.I.intersection_volume(body.translate(direction*d), s)) for k, s in obstacles.items()]
        if any(v > .001 for _, v in hits):
            for fine in np.arange(previous+.1, d+.001, .1):
                bad = [(k, A.I.intersection_volume(body.translate(direction*fine), s))
                       for k, s in obstacles.items()]
                bad = [(k, v) for k, v in bad if v > .001]
                if bad:
                    key, volume = max(bad, key=lambda x:x[1])
                    return dict(last_sample_below_threshold_mm=round(float(fine-.1), 1),
                                first_interference_mm=round(float(fine), 1), part=key, intersection_mm3=volume)
            raise ValueError('Fine path scan missed coarse collision')
        previous = float(d)
    return dict(last_sample_below_threshold_mm=limit, first_interference_mm=None, part=None)


def run():
    source = C.DEST/'assembly_check.json'
    data = json.loads(source.read_text(encoding='utf8'))
    if data['script_sha256'] != A.hashed(Path(C.__file__).resolve()):
        raise ValueError('Changed clamp generator')
    if data['parent_sha256'] != A.hashed(C.SOURCE):
        raise ValueError('Changed bench candidate')
    if not all(A.hashed(A.ROOT/p) == h for p, h in data['generation_dependency_sha256'].items()):
        raise ValueError('Changed generation dependency')
    solids, census = {}, {}
    for key, row in data['sources'].items():
        body, _, meta = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
        npositive = sum(p.volume() > 1e-6 for p in body.decompose())
        if npositive != 1:
            raise ValueError('Disconnected material: '+key)
        solids[key] = body
        census[key] = dict(boundary_edges=meta['boundary_edges'], nonmanifold_edges=meta['nonmanifold_edges'],
                           positive_material_components=npositive)
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        vol = A.I.intersection_volume(solids[a], solids[b])
        if vol > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=vol,
                                 compliant=any(k.startswith('tpu:') or 'liner' in k for k in [a, b])))
    contacts, bores, paths, screw_paths = [], [], [], []
    for tag, zone in data['zones'].items():
        frame = 'fixture:lateral_frame' if tag == 'lateral_upper' else 'fixture:medial_frame'
        prefix = 'clamp:'+tag+'/'
        n = np.array(zone['basis_columns_source'])[:, 2]
        for a, b, direction in [
            (prefix+'front_liner', 'bone:scapula', -n),
            (prefix+'back_liner', 'bone:scapula', n),
            (prefix+'front_liner', frame, n),
            (prefix+'back_liner', prefix+'moving_jaw', -n)]:
            contacts.append(dict(a=a, b=b, gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                                 nominal_overlap_mm3=A.I.intersection_volume(solids[a], solids[b]),
                                 translation_source_mm=(direction*.05).tolist(),
                                 normal_probe_overlap_mm3=A.I.intersection_volume(solids[a].translate(direction*.05), solids[b])))
        for bolt in zone['bolts']:
            i = str(bolt['index'])
            for a, b, direction in [
                (prefix+'near_washer_'+i, frame, -n),
                (prefix+'far_washer_'+i, prefix+'moving_jaw', n),
                (prefix+'screw_'+i, prefix+'near_washer_'+i, -n),
                (prefix+'nut_'+i, prefix+'far_washer_'+i, n)]:
                contacts.append(dict(a=a, b=b, gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                                     nominal_overlap_mm3=A.I.intersection_volume(solids[a], solids[b]),
                                     translation_source_mm=(direction*.05).tolist(),
                                     normal_probe_overlap_mm3=A.I.intersection_volume(solids[a].translate(direction*.05), solids[b])))
            q = np.array(bolt['axis_source_mm'])
            a, b = q-n*(abs(zone['outer_back_local_mm'])+6), q+n*(zone['outer_front_local_mm']+6)
            shaft = rod(a, b, 1.5)
            bore = rod(a, b, 1.7)
            bores.append(dict(zone=tag, index=int(i),
                              shaft_material_intersections_mm3={k:A.I.intersection_volume(shaft, solids[k])
                                                               for k in [frame, prefix+'moving_jaw', 'bone:scapula']},
                              full_nominal_bore_bone_intersection_mm3=A.I.intersection_volume(bore, solids['bone:scapula'])))
            screw = prefix+'screw_'+i
            obstacles = {k:s for k, s in solids.items() if k != screw}
            fixture = {k:s for k, s in obstacles.items() if not k.startswith(
                ('tpu:', 'cap:', 'stop:', 'link:', 'collar:', 'free_terminal:', 'origin:'))}
            required = bolt['nominal_length_mm']+.1
            path_bare = translate_path(solids[screw], n, fixture, bolt['nominal_length_mm']+10.)
            path_full = translate_path(solids[screw], n, obstacles, bolt['nominal_length_mm']+10.)
            screw_paths.append(dict(screw=screw, direction_source=n.tolist(),
                maximum_translation_mm=bolt['nominal_length_mm']+10.,
                required_translation_to_clear_near_washer_mm=required,
                required_bare_path_clear=path_bare['first_interference_mm'] is None or path_bare['first_interference_mm'] > required,
                required_full_context_path_clear=path_full['first_interference_mm'] is None or path_full['first_interference_mm'] > required,
                bare_bone_fixture=path_bare, full_undeformed_context=path_full))
        # Screws, nuts and washers of this jaw are first removed. The liner
        # is left at its nominal bone position: no fictitious adhesive or
        # rigid liner-to-jaw constraint is introduced.
        moving = prefix+'moving_jaw'
        own_hardware = {prefix+k+'_'+str(i) for k in ['screw', 'nut', 'near_washer', 'far_washer'] for i in [0, 1]}
        obstacles = {k:s for k, s in solids.items() if k != moving and k not in own_hardware}
        bare = {k:s for k, s in obstacles.items() if not k.startswith(('tpu:', 'cap:', 'stop:', 'link:', 'collar:', 'free_terminal:', 'origin:'))}
        paths.append(dict(zone=tag, moving=moving, direction_source=(-n).tolist(), step_mm=1., boundary_step_mm=.1,
                          limit_mm=20., excluded_own_hardware=sorted(own_hardware),
                          bare_bone_fixture=translate_path(solids[moving], -n, bare),
                          full_undeformed_context=translate_path(solids[moving], -n, obstacles)))
    # The frame was changed at its bone-side end; recheck the earlier bench
    # fasteners against the new frame instead of inheriting parent access.
    parent = json.loads(C.SOURCE.read_text(encoding='utf8'))
    bench_paths = []
    for key in solids:
        if not key.startswith('bench:') or not key.endswith('/screw'):
            continue
        name = key.split(':', 1)[1].split('/')[0]
        if name.startswith('foot_'):
            spec = parent['base_mounts'][name]
        else:
            group, index = name.rsplit('_', 1)
            spec = parent['base_mounts'][group]['bores'][int(index)]['fastener']
        obstacles = {k:s for k, s in solids.items() if k != key}
        frame_only = {k:s for k, s in obstacles.items() if k.startswith(('fixture:', 'bench:'))}
        direction = -np.array(spec['axis_source'])
        distance = spec['nominal_length_mm']+10.
        bench_paths.append(dict(screw=key, maximum_translation_mm=distance,
                               before_bones_and_muscles=translate_path(solids[key], direction, frame_only, distance),
                               full_undeformed_context=translate_path(solids[key], direction, obstacles, distance)))
    report = dict(scope=__doc__, assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
                  n_instances=len(solids), all_source_hashes_match=True, mesh_census=census,
                  all_pair_count=len(solids)*(len(solids)-1)//2, material_intersections=overlaps,
                  contact_checks=contacts, external_M3_bore_checks=bores, moving_jaw_withdrawal_checks=paths,
                  clamp_screw_withdrawal_checks=screw_paths, bench_screw_withdrawal_checks=bench_paths,
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['clavicle clamp and complete fixture constraint', 'no force/preload, friction or material/printing validation',
                                'hardware insertion, tools and complete multibody assembly paths',
                                'loose liners only retained by nominal compression; real assembly procedure unverified',
                                'same unresolved channel-wall differences and coupled muscle/tendon dynamics as parent'])
    (C.DEST/'readback_contact_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects;', report['all_pair_count'], 'pairs', flush=True)
    print('Rigid overlaps > .001', [r for r in overlaps if not r['compliant'] and r['volume_mm3'] > .001], flush=True)
    print('Missing contact probes', [r for r in contacts if r['normal_probe_overlap_mm3'] <= .001], flush=True)
    print('Blocked M3 bores', [r for r in bores if max(r['shaft_material_intersections_mm3'].values()) > .001 or r['full_nominal_bore_bone_intersection_mm3'] > .001], flush=True)
    for row in paths:
        print(row['zone'], 'bare', row['bare_bone_fixture'], 'full', row['full_undeformed_context'], flush=True)
    print('M3 screw path obstructions', [r for r in screw_paths if r['full_undeformed_context']['first_interference_mm'] is not None], flush=True)
    print('M3 required withdrawal before muscles clear', all(r['required_bare_path_clear'] for r in screw_paths),
          'full-context clear', all(r['required_full_context_path_clear'] for r in screw_paths), flush=True)
    print('Bench paths before bone assembly obstructed', [r for r in bench_paths if r['before_bones_and_muscles']['first_interference_mm'] is not None], flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
