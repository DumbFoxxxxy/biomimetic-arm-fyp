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
    contacts, bores, paths = [], [], []
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
    report = dict(scope=__doc__, assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
                  n_instances=len(solids), all_source_hashes_match=True, mesh_census=census,
                  all_pair_count=len(solids)*(len(solids)-1)//2, material_intersections=overlaps,
                  contact_checks=contacts, external_M3_bore_checks=bores, moving_jaw_withdrawal_checks=paths,
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
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
