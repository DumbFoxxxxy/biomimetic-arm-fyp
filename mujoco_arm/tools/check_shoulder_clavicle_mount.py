"""Readback, contacts and selected staged paths for the clavicle test support.

The anatomical bone withdrawal probe is only a straight candidate path;
fixed muscles are omitted in the bare-fixture case. No physical fixation,
preload, friction, strength, tool accessibility or actuation is established.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_shoulder_clavicle_mount as K
import audit_shoulder_retained_joint_fit as A
from check_shoulder_scapula_clamps import translate_path


def run():
    source = K.DEST/'assembly_check.json'
    data = json.loads(source.read_text(encoding='utf8'))
    if data['script_sha256'] != A.hashed(Path(K.__file__).resolve()):
        raise ValueError('Changed clavicle generator')
    if data['parent_sha256'] != A.hashed(K.SOURCE):
        raise ValueError('Changed scapula-clamp parent')
    if not all(A.hashed(A.ROOT/p) == h for p, h in data['generation_dependency_sha256'].items()):
        raise ValueError('Changed generation dependency')
    solids, census = {}, {}
    for key, row in data['sources'].items():
        s, _, meta = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
        positive = sum(p.volume() > 1e-6 for p in s.decompose())
        if positive != 1:
            raise ValueError('Disconnected readback material: '+key)
        solids[key] = s
        census[key] = dict(boundary_edges=meta['boundary_edges'], nonmanifold_edges=meta['nonmanifold_edges'],
                           positive_material_components=positive)
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        vol = A.I.intersection_volume(solids[a], solids[b])
        if vol > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=vol,
                                 compliant=any(k.startswith('tpu:') or 'liner' in k for k in [a, b])))
    F = np.array(data['clavicle_fixture_frame']['basis_columns_source'])
    n = F[:, 1]
    contacts, bores, hardware_paths = [], [], []
    fixed, moving = 'clavicle:fixed_support', 'clavicle:moving_jaw'
    pairs = [('clavicle:liner_posterior', 'bone:clavicle', -n),
             ('clavicle:liner_anterior', 'bone:clavicle', n),
             ('clavicle:liner_posterior', fixed, n),
             ('clavicle:liner_anterior', moving, -n),
             (fixed, 'bench:backboard', np.array([0., 1, 0]))]
    for name, spec in data['fasteners'].items():
        label, index = name.split('_')
        prefix = 'clavicle:'+label+'_'
        if label == 'M3':
            front, back = fixed, moving
            direction = np.array(spec['direction_source'])
            near = np.array(spec['near_bearing_source_mm'])
            bore_origin = np.array(spec['bore_axis_source_mm'])
            shaft = K.rod(bore_origin-direction*20, bore_origin+direction*20, 1.5)
            bore = K.rod(bore_origin-direction*20, bore_origin+direction*20, 1.7)
            bores.append(dict(fastener=name, actual_shaft_intersections_mm3={
                k:A.I.intersection_volume(shaft, solids[k]) for k in [fixed, moving, 'bone:clavicle']},
                full_bore_bone_intersection_mm3=A.I.intersection_volume(bore, solids['bone:clavicle'])))
            required = spec['required_withdrawal_mm']
        else:
            front, back = fixed, 'bench:backboard'
            direction = np.array(spec['axis_source'])
            near = np.array(spec['near_bearing_source_mm'])
            shaft = K.rod(near-direction*2, near+direction*22, 2.)
            bores.append(dict(fastener=name, actual_shaft_intersections_mm3={
                k:A.I.intersection_volume(shaft, solids[k]) for k in [fixed, back]},
                full_bore_bone_intersection_mm3=A.I.intersection_volume(
                    K.rod(near-direction*2, near+direction*22, 2.2), solids['bone:clavicle'])))
            required = spec['nominal_length_mm']+.1
        pairs.extend([(prefix+'near_washer_'+index, front, direction),
                      (prefix+'far_washer_'+index, back, -direction),
                      (prefix+'screw_'+index, prefix+'near_washer_'+index, direction),
                      (prefix+'nut_'+index, prefix+'far_washer_'+index, -direction)])
        screw = prefix+'screw_'+index
        obstacles = {k:s for k, s in solids.items() if k != screw}
        path = translate_path(solids[screw], -direction, obstacles, spec['nominal_length_mm']+10.)
        hardware_paths.append(dict(part=screw, direction_source=(-direction).tolist(),
            required_translation_to_clear_near_washer_mm=required,
            required_sampled_path_clear=path['first_interference_mm'] is None or path['first_interference_mm'] > required,
            sampled_withdrawal=path))
    for a, b, direction in pairs:
        contacts.append(dict(a=a, b=b, gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                             nominal_overlap_mm3=A.I.intersection_volume(solids[a], solids[b]),
                             direction_source=direction.tolist(),
                             normal_probe_overlap_mm3=A.I.intersection_volume(solids[a].translate(direction*.05), solids[b])))
    clamp_hardware = {k for k in solids if k.startswith('clavicle:M3_')}
    jaw_obstacles = {k:s for k, s in solids.items() if k != moving and k not in clamp_hardware}
    jaw_path = translate_path(solids[moving], -n, jaw_obstacles, 20.)
    # Open anterior half and remove its loose liner/hardware first. The
    # posterior liner remains in place; no artificial liner adhesion.
    removed = clamp_hardware | {moving, 'clavicle:liner_anterior', 'bone:clavicle'}
    bone_obstacles = {k:s for k, s in solids.items() if k not in removed}
    bare = {k:s for k, s in bone_obstacles.items() if not k.startswith(
        ('tpu:', 'cap:', 'stop:', 'link:', 'collar:', 'free_terminal:', 'origin:'))}
    bone_paths = dict(direction_source=(-n).tolist(), limit_mm=20.,
                      excluded_objects=sorted(removed-{'bone:clavicle'}),
                      bare_bone_fixture=translate_path(solids['bone:clavicle'], -n, bare, 20.),
                      full_undeformed_context=translate_path(solids['bone:clavicle'], -n, bone_obstacles, 20.))
    # Support can only be removed after the clavicle and clamp halves are
    # taken out and both M4 groups removed. This is a separate partial step.
    removed_support = {k for k in solids if k.startswith('clavicle:')} | {'bone:clavicle'}
    support_obstacles = {k:s for k, s in solids.items() if k not in removed_support}
    support_path = translate_path(solids[fixed], np.array([0., -1, 0]), support_obstacles, 30.)
    report = dict(scope=__doc__, assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
        path_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/check_shoulder_scapula_clamps.py'),
        n_instances=len(solids), all_source_hashes_match=True, mesh_census=census,
        all_pair_count=len(solids)*(len(solids)-1)//2, material_intersections=overlaps,
        finite_contact_checks=contacts, fastener_bore_checks=bores, hardware_withdrawal_checks=hardware_paths,
        moving_jaw_withdrawal=dict(direction_source=(-n).tolist(), limit_mm=20.,
                                   excluded_hardware=sorted(clamp_hardware), sampled_path=jaw_path),
        clavicle_withdrawal=bone_paths,
        support_withdrawal=dict(direction_source=[0, -1, 0], limit_mm=30.,
                                assumed_removed_objects=sorted(removed_support-{fixed}), sampled_path=support_path),
        stage_pass=False, actual_drive_range_deg=None,
        not_verified=['full multi-part assembly/tool/nut/washer sequence and lateral screw transfer',
                      'posterior liner/rigid jaw constraint compatibility during actual assembly',
                      'friction, preload, support strength, liner compression and printed tolerances',
                      'scapula whole-bone insertion and old blocked access cases still apply',
                      'no true six-DOF fixture/load proof or new muscle actuation result'])
    (K.DEST/'readback_contact_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects;', report['all_pair_count'], 'pairs', flush=True)
    print('Rigid overlaps >.001', [r for r in overlaps if not r['compliant'] and r['volume_mm3'] > .001], flush=True)
    print('Missing contacts', [r for r in contacts if r['normal_probe_overlap_mm3'] <= .001], flush=True)
    print('Blocked bores', [r for r in bores if max(r['actual_shaft_intersections_mm3'].values()) > .001 or r['full_bore_bone_intersection_mm3'] > .001], flush=True)
    print('Hardware paths blocked before clearing:', [r for r in hardware_paths if not r['required_sampled_path_clear']], flush=True)
    print('Jaw path', jaw_path, 'bone', bone_paths['bare_bone_fixture'], bone_paths['full_undeformed_context'], 'support', support_path, flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
