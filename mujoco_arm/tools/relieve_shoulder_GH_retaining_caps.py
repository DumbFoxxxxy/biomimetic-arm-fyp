"""Relieve only new GH lip material around unchanged subscapularis geometry.

Keeps the original cap and bearing faces, bones, joint centre and muscles.
The box Minkowski keepout is a nominal geometric allowance, not a validated
printed clearance or deformable-muscle motion envelope.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_bench_mount as B
from build_shoulder_unified_muscle import local_export
from check_shoulder_front_liner_keepers import is_compliant

SOURCE = A.DESIGN/'s1_GH_assembly_review/retained_depth_6/assembly_check.json'
DEST = A.DESIGN/'s1_GH_assembly_review/retained_depth_6_relief'


def run():
    DEST.mkdir(parents=True, exist_ok=True)
    parent = json.loads(SOURCE.read_text(encoding='utf8'))
    original_path = A.ROOT/parent['parent_report']
    original = json.loads(original_path.read_text(encoding='utf8'))
    if A.hashed(original_path) != parent['parent_sha256']:
        raise ValueError('Original cap reference changed')
    sources = dict(parent['sources'])
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k, r in sources.items()}
    mf = A.I._manifold()
    changes = {}
    muscle = solids['tpu:subscapularis']
    for key in ['GH_cap_0', 'GH_cap_1']:
        current = solids[key]
        row = original['sources'][key]
        old = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])[0]
        lo, hi = np.array(current.bounding_box()).reshape(2, 3)
        clipped_muscle = muscle ^ B.box(lo-1., hi+1.)
        keepout = clipped_muscle.minkowski_sum(mf.Manifold.cube([.5]*3, True))
        revised = old + ((current-old)-keepout)
        if sum(p.volume() > 1e-6 for p in revised.decompose()) != 1:
            raise ValueError('Relief disconnects cap '+key)
        actual, meta = local_export(revised.as_original().simplify(.001), DEST/(key+'_retained_relief.stl'),
                                    np.array(sources[key]['translation_source_mm']))
        if not meta['mesh_checks_pass']:
            raise ValueError('Readback failed '+key)
        solids[key], sources[key] = actual, meta
        changes[key] = dict(removed_from_depth6_mm3=float((current-actual).volume()),
                            added_by_readback_mm3=float((actual-current).volume()),
                            original_cap_removed_mm3=float((old-actual).volume()),
                            muscle_overlap_mm3=A.I.intersection_volume(actual, muscle),
                            muscle_gap_capped_2mm=A.I.gap(actual, muscle, 2.))
        print(key, changes[key], flush=True)
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        volume = A.I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=volume,
                                 compliant=is_compliant(a) or is_compliant(b)))
    routes_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    wire_checks = []
    for row in json.loads(routes_file.read_text(encoding='utf8'))['routes']:
        wire = A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), .3)
        hits = {k:A.I.intersection_volume(wire, s) for k,s in solids.items()}
        wire_checks.append(dict(head=row['head'], volumes_mm3=hits,
                               hits_over_threshold=[dict(part=k, volume_mm3=v) for k,v in hits.items() if v > .001]))
    data = dict(scope=__doc__, n_instances=len(solids), sources=sources,
                parent_report=SOURCE.relative_to(A.ROOT).as_posix(), parent_sha256=A.hashed(SOURCE),
                original_cap_reference=original_path.relative_to(A.ROOT).as_posix(),
                original_cap_reference_sha256=A.hashed(original_path),
                script_sha256=A.hashed(Path(__file__).resolve()), routes_sha256=A.hashed(routes_file),
                generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
                    [A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',
                     A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py']},
                GH_retention_revision=parent['GH_retention_revision'],
                GH_lip_relief=dict(keepout='clipped unchanged subscapularis Minkowski sum with centred 0.5mm cube',
                                  allowance_source_xyz_each_side_mm=.25, cap_changes=changes),
                changed_object_count=2, pair_count_checked=len(solids)*(len(solids)-1)//2,
                material_intersections=overlaps, actuation_wire_checks=wire_checks,
                free_terminal_constraints_required=parent['free_terminal_constraints_required'],
                intended_body_owners=parent['intended_body_owners'],
                scapula_insertion_relief=parent['scapula_insertion_relief'],
                keeper_definitions=parent['keeper_definitions'], stage_pass=False, actual_drive_range_deg=None,
                not_verified=['actual-axis capture must be independently rechecked after this relief',
                              'caps, all GH fasteners, tool paths and full assembly',
                              'arbitrary escape, printed clearance, strength and wear',
                              'deformable muscle motion, force and actual drive range'])
    (DEST/'assembly_check.json').write_text(json.dumps(data, indent=2)+'\n', encoding='utf8')
    print('Rigid overlaps', [r for r in overlaps if r['volume_mm3'] > .001 and not r['compliant']], flush=True)
    print('Soft overlaps', [r for r in overlaps if r['volume_mm3'] > .001 and r['compliant']], flush=True)
    print('Wire failures', sum(bool(r['hits_over_threshold']) for r in wire_checks), flush=True)
    return data


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
