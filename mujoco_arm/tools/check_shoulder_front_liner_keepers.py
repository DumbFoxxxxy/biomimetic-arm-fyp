"""Readback, staged placement and finite capture probes of front-liner keepers.

Fasteners use nominal thread-free envelopes. Slot loading, axial clearance
and tool-envelope checks are not tightening, preload or pullout tests.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_shoulder_front_liner_keepers as L
import audit_shoulder_retained_joint_fit as A
import build_shoulder_scapula_clamps as C
import build_shoulder_clavicle_mount as K
import build_shoulder_bench_mount as B


def is_compliant(key):
    # Keeper hardware contains "liner" in its prefix but is nominal metal.
    # Only the actual muscle/liner material keys are compliant candidates.
    return (key.startswith('tpu:') or key.startswith('collar:liner') or
            key.startswith('clavicle:liner_') or
            (key.startswith('clamp:') and key.endswith(('/front_liner', '/back_liner'))))


def dense_slide(moving, fixed, direction, maximum):
    """Both directions recomputed at 0.1 mm, stop each phase at first hit."""
    values = np.linspace(0., maximum, int(round(maximum/.1))+1)
    phases = []
    for phase, positions in [('withdrawal', values), ('inverse_placement', values[::-1])]:
        samples = []
        for value in positions:
            hits = []
            maximum_overlap = 0.
            for name, body in moving.items():
                posed = body.translate(direction*value)
                for other, obstacle in fixed.items():
                    volume = A.I.intersection_volume(posed, obstacle)
                    maximum_overlap = max(maximum_overlap, volume)
                    if volume > .001:
                        hits.append(dict(moving=name, fixed=other, volume_mm3=volume))
            samples.append(dict(displacement_mm=round(float(value), 6), maximum_pair_overlap_mm3=maximum_overlap,
                                hits_over_threshold=hits))
            if hits:
                break
        phases.append(dict(name=phase, samples=samples, sample_count=len(samples),
                           entire_sampled_phase_clear=len(samples) == len(values) and not any(r['hits_over_threshold'] for r in samples)))
    return dict(direction_source=direction.tolist(), maximum_translation_mm=maximum, step_mm=.1,
                moving_keys=list(moving), fixed_keys=list(fixed), phases=phases,
                sampled_path_clear=all(p['entire_sampled_phase_clear'] for p in phases))


def first_contact(body, direction, fixed, maximum=1., step=.05):
    for value in np.arange(step, maximum+step/10, step):
        hits = [dict(part=k, volume_mm3=A.I.intersection_volume(body.translate(direction*value), s)) for k, s in fixed.items()]
        hits = [r for r in hits if r['volume_mm3'] > .001]
        if hits:
            return dict(displacement_mm=round(float(value), 6), hits=sorted(hits, key=lambda r:-r['volume_mm3']))
    return None


def run():
    source = L.DEST/'assembly_check.json'
    data = json.loads(source.read_text(encoding='utf8'))
    if data['script_sha256'] != A.hashed(Path(L.__file__).resolve()):
        raise ValueError('Changed keeper generator')
    if data['parent_sha256'] != A.hashed(L.SOURCE):
        raise ValueError('Changed receiver parent')
    if not all(A.hashed(A.ROOT/p) == h for p, h in data['generation_dependency_sha256'].items()):
        raise ValueError('Changed generation dependency')
    solids, census = {}, {}
    for key, row in data['sources'].items():
        body, _, meta = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
        positive = sum(p.volume() > 1e-6 for p in body.decompose())
        if positive != 1:
            raise ValueError('Disconnected keeper material '+key)
        solids[key] = body
        census[key] = dict(boundary_edges=meta['boundary_edges'], nonmanifold_edges=meta['nonmanifold_edges'],
                           positive_material_components=positive)
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        volume = A.I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=volume,
                                 compliant=is_compliant(a) or is_compliant(b)))
    bare_keys = {k for k in solids if k.startswith(('fixture:', 'bench:', 'liner_keeper:'))
                 or k.endswith('/front_liner') or k == 'clavicle:fixed_support' or k.startswith('clavicle:M4_')}
    declarations, contacts = {}, []
    for path, field, vector, factor in [
        (B.DEST/'assembly_check.json', 'finite_contact_checks', 'probe_translation_source_mm', 1.),
        (C.DEST/'readback_contact_checks.json', 'contact_checks', 'translation_source_mm', 1.),
        (K.DEST/'readback_contact_checks.json', 'finite_contact_checks', 'direction_source', .05)]:
        original = json.loads(path.read_text(encoding='utf8'))
        declarations[path.relative_to(A.ROOT).as_posix()] = A.hashed(path)
        for row in original[field]:
            a, b = row['a'], row['b']
            d = np.array(row[vector])*factor
            contacts.append(dict(a=a, b=b, declaration_source=path.relative_to(A.ROOT).as_posix(),
                                 probe_translation_source_mm=d.tolist(),
                                 nominal_gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                                 probe_overlap_mm3=A.I.intersection_volume(solids[a].translate(d), solids[b])))
    placements, captures, tools = [], [], []
    for tag, definition in data['keeper_definitions'].items():
        F = np.array(definition['basis_columns_source'])
        g = np.array(definition['source_origin_mm'])
        n, u = F[:, 2], F[:, 0]
        liner, frame = definition['liner_key'], definition['frame_key']
        prefix = 'liner_keeper:'+tag+'/'
        hardware = {prefix+k for k in ['nut', 'screw', 'spacer', 'washer']}
        other_stage = {k:solids[k] for k in bare_keys if k != liner and k not in hardware}
        # Load the nut while the liner is outside the frame. The new post
        # must not rely on loading a side slot that the rigid bore encloses.
        nut, screw = prefix+'nut', prefix+'screw'
        shifted_liner = solids[liner].translate(-n*20)
        nut_fixed = {**other_stage, liner:shifted_liner}
        nut_load = dense_slide({nut:solids[nut].translate(-n*20)}, nut_fixed, u, 6.)
        pack = dense_slide({liner:solids[liner], nut:solids[nut]}, other_stage, -n, 20.)
        placement = dict(zone=tag, stage_not_installed_keys=sorted(set(solids)-bare_keys),
                         nut_slot_before_liner_placement=nut_load,
                         liner_and_temporarily_held_nut=pack,
                         physical_handling_not_verified=True)
        for name, removed, limit in [('spacer', {'screw', 'washer'}, 12.),
                                     ('washer', {'screw'}, 12.), ('screw', set(), 6.1)]:
            key = prefix+name
            omitted = {prefix+k for k in removed}|{key}
            fixed = {k:solids[k] for k in bare_keys if k not in omitted}
            placement[name+'_axial_clearance'] = dense_slide({key:solids[key]}, fixed, n, limit)
        placements.append(placement)
        for a, b, d in [(prefix+'spacer', frame, -n*.05),
                        (prefix+'washer', prefix+'spacer', -n*.05),
                        (screw, prefix+'washer', -n*.05),
                        (nut, liner, n*.1)]:
            contacts.append(dict(a=a, b=b, declaration_source='new keeper nominal bearing declaration',
                                 nominal_gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                                 probe_translation_source_mm=d.tolist(),
                                 probe_overlap_mm3=A.I.intersection_volume(solids[a].translate(d), solids[b])))
        finite = []
        # Move the rigidly probed lining alone against an ideally installed
        # nut/screw group. This is only geometric capture, not TPU pullout.
        fixed = {k:solids[k] for k in bare_keys if k != liner}
        for label, axis in [('outward_from_receiver', -n), ('toward_receiver', n),
                            ('u_plus', u), ('u_minus', -u), ('v_plus', F[:, 1]), ('v_minus', -F[:, 1])]:
            finite.append(dict(probe=label, direction_source=axis.tolist(), step_mm=.05, limit_mm=1.,
                               first_interference=first_contact(solids[liner], axis, fixed)))
        # Side-slot escape after mounting, with the screw omitted: the bore
        # should block the nut before it can leave the TPU carrier entirely.
        nut_obstacles = {k:solids[k] for k in bare_keys if k not in [nut, screw]}
        radial = first_contact(solids[nut], u, nut_obstacles, maximum=6., step=.1)
        captures.append(dict(zone=tag, liner_translation_probes=finite,
                             mounted_nut_slot_exit_probe=dict(direction_source=u.tolist(), step_mm=.1,
                                                              limit_mm=6., screw_assumed_removed=True,
                                                              first_interference=radial)))
        # A deliberately stated candidate straight driver entry, outside
        # the nominal head. It is not a modeled socket or selected tool.
        back = definition['back_bearing_local_n_mm']
        probe = B.rod(g+n*(back+5.4), g+n*(back+25.4), 1.5)
        hits = {k:A.I.intersection_volume(probe, solids[k]) for k in bare_keys if k != screw}
        tools.append(dict(zone=tag, candidate_entry_diameter_mm=3., candidate_entry_length_mm=20.,
                          axis_source=n.tolist(), stage='before bones/GH/muscles',
                          intersections_mm3=hits, hits_over_threshold=[dict(part=k, volume_mm3=v) for k, v in hits.items() if v > .001],
                          actual_driver_and_socket_not_selected=True))
        print(tag, 'slot', nut_load['sampled_path_clear'], 'liner', pack['sampled_path_clear'],
              'metal', {k:r['sampled_path_clear'] for k, r in placement.items() if k.endswith('_axial_clearance')}, flush=True)
    result = dict(scope=__doc__, assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
                  n_instances=len(solids), all_source_hashes_match=True, mesh_census=census,
                  all_pair_count=len(solids)*(len(solids)-1)//2, material_intersections=overlaps,
                  explicit_material_classes={k:('compliant_candidate' if is_compliant(k) else 'rigid_or_wire_candidate') for k in solids},
                  parent_contact_declaration_sha256=declarations, rechecked_contacts=contacts,
                  partial_placement_checks=placements, finite_capture_checks=captures, candidate_tool_entries=tools,
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['thread engagement, tightening/torque, printed nut slot fit and TPU pullout or fatigue',
                                'actual retaining force/position under handling and assembly loads',
                                'complete bone/GH/rear jaw/hardware and muscle installation sequence',
                                'no new joint movement or TPU muscle driving result'])
    (L.DEST/'readback_retention_checks.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects;', result['all_pair_count'], 'pairs', flush=True)
    print('Rigid overlaps >.001', [r for r in overlaps if not r['compliant'] and r['volume_mm3'] > .001], flush=True)
    print('Contact declarations', len(contacts), 'missing', [r for r in contacts if r['probe_overlap_mm3'] <= .001], flush=True)
    print('Tool entry failures', [r for r in tools if r['hits_over_threshold']], flush=True)
    return result


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
