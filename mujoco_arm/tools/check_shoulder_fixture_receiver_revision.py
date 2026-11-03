"""Full readback and inherited contact rechecks after receiver micro-relief.

Re-evaluates parent contact declarations using current actual meshes, rather
than inheriting unchanged pass labels. Geometry is not a preload/load proof.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relieve_shoulder_fixture_for_scapula_insertion as R
import audit_shoulder_retained_joint_fit as A
import build_shoulder_scapula_clamps as C
import build_shoulder_clavicle_mount as K
import build_shoulder_bench_mount as B


def run():
    source = R.DEST/'assembly_check.json'
    data = json.loads(source.read_text(encoding='utf8'))
    if data['script_sha256'] != A.hashed(Path(R.__file__).resolve()):
        raise ValueError('Changed receiver revision generator')
    if data['parent_sha256'] != A.hashed(R.SOURCE):
        raise ValueError('Changed parent candidate')
    if not all(A.hashed(A.ROOT/p) == h for p, h in data['generation_dependency_sha256'].items()):
        raise ValueError('Changed generation dependency')
    solids, census = {}, {}
    for key, row in data['sources'].items():
        s, _, meta = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
        positive = sum(p.volume() > 1e-6 for p in s.decompose())
        if positive != 1:
            raise ValueError('Disconnected material '+key)
        solids[key] = s
        census[key] = dict(boundary_edges=meta['boundary_edges'], nonmanifold_edges=meta['nonmanifold_edges'],
                           positive_material_components=positive)
    contacts, inputs = [], {}
    for path, field, vector_field, factor in [
        (B.DEST/'assembly_check.json', 'finite_contact_checks', 'probe_translation_source_mm', 1.),
        (C.DEST/'readback_contact_checks.json', 'contact_checks', 'translation_source_mm', 1.),
        (K.DEST/'readback_contact_checks.json', 'finite_contact_checks', 'direction_source', .05)]:
        old = json.loads(path.read_text(encoding='utf8'))
        inputs[path.relative_to(A.ROOT).as_posix()] = A.hashed(path)
        for row in old[field]:
            a, b = row['a'], row['b']
            direction = np.array(row[vector_field])*factor
            contacts.append(dict(a=a, b=b, declaration_source=path.relative_to(A.ROOT).as_posix(),
                                 translation_source_mm=direction.tolist(),
                                 gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                                 normal_probe_overlap_mm3=A.I.intersection_volume(solids[a].translate(direction), solids[b])))
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        vol = A.I.intersection_volume(solids[a], solids[b])
        if vol > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=vol,
                                 compliant=any(k.startswith('tpu:') or 'liner' in k for k in [a, b])))
    zones = json.loads((C.DEST/'assembly_check.json').read_text(encoding='utf8'))['zones']
    bores = []
    for tag, zone in zones.items():
        n = np.array(zone['basis_columns_source'])[:, 2]
        frame = 'fixture:lateral_frame' if tag == 'lateral_upper' else 'fixture:medial_frame'
        for bolt in zone['bolts']:
            q = np.array(bolt['axis_source_mm'])
            probe = B.rod(q-n*(abs(zone['outer_back_local_mm'])+6), q+n*(zone['outer_front_local_mm']+6), 1.5)
            bores.append(dict(zone=tag, index=bolt['index'], intersections_mm3={
                k:A.I.intersection_volume(probe, solids[k]) for k in [frame, 'clamp:'+tag+'/moving_jaw', 'bone:scapula']}))
    dense_file = R.DEST.parent/'selected_path_dense_check.json'
    dense = json.loads(dense_file.read_text(encoding='utf8'))
    if dense['assembly_sha256'] != A.hashed(source):
        raise ValueError('Dense path does not refer to this receiver geometry')
    result = dict(scope=__doc__, assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
                  inherited_contact_declaration_sha256=inputs,
                  n_instances=len(solids), all_source_hashes_match=True, mesh_census=census,
                  all_pair_count=len(solids)*(len(solids)-1)//2, material_intersections=overlaps,
                  rechecked_contact_pairs=contacts, rechecked_scapula_M3_bores=bores,
                  dense_path_sha256=A.hashed(dense_file), dense_partial_path_samples_clear=dense['selected_partial_path_samples_clear'],
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['temporary front liner holding, tools and remaining assembly steps',
                                'minimum wall thickness, preload/friction, material and support strength',
                                'joint movement or flexible actuation; prior blocked accesses remain unverified'])
    (R.DEST/'readback_contact_checks.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects;', result['all_pair_count'], 'pairs', flush=True)
    print('Rigid overlaps >.001', [r for r in overlaps if not r['compliant'] and r['volume_mm3'] > .001], flush=True)
    print('Contact declarations rechecked', len(contacts), 'missing probes',
          [r for r in contacts if r['normal_probe_overlap_mm3'] <= .001], flush=True)
    print('Scapula M3 bore failures', [r for r in bores if max(r['intersections_mm3'].values()) > .001], flush=True)
    return result


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
