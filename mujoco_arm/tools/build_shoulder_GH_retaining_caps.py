"""Independent deeper GH retaining sectors; old mounting faces are retained.

Sphere/cup/neck and bone sources stay unchanged. Extend only retained shell
material to n=6 mm, retain the original two sectors and original fastener
seats, and clear their actual hardware. No strength or arbitrary escape proof.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_front_liner_keepers as L
import build_shoulder_bench_mount as B
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export

SOURCE = L.DEST/'assembly_check.json'
DEST = A.DESIGN/'s1_GH_assembly_review/retained_depth_6'


def run():
    DEST.mkdir(parents=True, exist_ok=True)
    parent = json.loads(SOURCE.read_text(encoding='utf8'))
    sources = dict(parent['sources'])
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k, r in sources.items()}
    reference = A.G1/'assembly.json'
    joint = json.loads(reference.read_text(encoding='utf8'))
    G, n = np.array(joint['centre_source_mm']), np.array(joint['neck_axis'])
    x = A.D.unit(np.cross(n, [0., -1, 0]))
    y = np.cross(n, x)
    F = np.column_stack([x, y, n])
    mf = A.I._manifold()
    neck_cuts = []
    for axis in joint['neck_slot']['axes']:
        for angle in range(0, 161, 2):
            direction = A.I.rotation_about([0, 0, 0], axis, float(angle))[:, :3]@n
            neck_cuts.append(rod(G, G+direction*(joint['socket_outer_radius_mm']+10.), joint['neck_slot']['radius_mm']))
    neck_slot = mf.Manifold.batch_boolean(neck_cuts, mf.OpType.Add)
    shell = (mf.Manifold.sphere(joint['socket_outer_radius_mm'], 40).translate(G)-
             mf.Manifold.sphere(joint['socket_inner_radius_mm'], 40).translate(G)-neck_slot)
    extension = shell ^ B.box([-30, -30, 2.9], [30, 30, 6.]).transform(np.c_[F, G])
    assignments = {}
    for part in extension.decompose():
        if part.volume() <= 1e-6:
            continue
        owners = [k for k in ['GH_cap_0', 'GH_cap_1'] if A.I.intersection_volume(part, solids[k]) > .001]
        if len(owners) != 1:
            raise ValueError('New retained sector does not uniquely connect to one old cap')
        assignments.setdefault(owners[0], []).append(part)
    if set(assignments) != {'GH_cap_0', 'GH_cap_1'}:
        raise ValueError('Both old retained sectors must receive connected extensions')
    changes = {}
    for key, components in assignments.items():
        old = solids[key]
        addition = mf.Manifold.batch_boolean(components, mf.OpType.Add)
        raw_volume = addition.volume()
        # Notch only the new shell around actual old screws/washers; the
        # old cap and its M3 bearing surfaces do not move or get recut.
        for label, obstacle in solids.items():
            if label.startswith('GH_') and not label.startswith('GH_cap_'):
                if A.I.intersection_volume(addition, obstacle) > 1e-9:
                    keepout = mf.Manifold.batch_boolean(
                        [obstacle]+[obstacle.translate(np.eye(3)[j]*sign*.15)
                                    for j in range(3) for sign in [-1, 1]], mf.OpType.Add)
                    addition -= keepout
        revised = old+addition
        if sum(p.volume() > 1e-6 for p in revised.decompose()) != 1:
            raise ValueError('Deeper retained cap material is disconnected')
        actual, meta = local_export(revised.as_original().simplify(.001), DEST/(key+'_retained.stl'), G)
        if not meta['mesh_checks_pass']:
            (DEST/'rejected_mesh.json').write_text(json.dumps(dict(key=key, meta=meta), indent=2), encoding='utf8')
            raise ValueError('Deeper cap readback fails '+key)
        solids[key], sources[key] = actual, meta
        changes[key] = dict(old_source_sha256=parent['sources'][key]['sha256'],
                            raw_extension_volume_mm3=float(raw_volume),
                            extension_after_hardware_clearance_mm3=float(addition.volume()),
                            actual_added_mm3=float((actual-old).volume()),
                            actual_removed_by_export_mm3=float((old-actual).volume()))
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        if a not in changes and b not in changes:
            continue
        volume = A.I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=volume))
    route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    wire_checks = []
    for row in json.loads(route_file.read_text(encoding='utf8'))['routes']:
        wire = A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), .3)
        hits = {k:A.I.intersection_volume(wire, s) for k, s in solids.items()}
        wire_checks.append(dict(head=row['head'], volumes_mm3=hits,
                                hits_over_threshold=[dict(part=k, volume_mm3=v) for k, v in hits.items() if v > .001]))
    report = dict(scope=__doc__, n_instances=len(solids), sources=sources,
        parent_report=SOURCE.relative_to(A.ROOT).as_posix(), parent_sha256=A.hashed(SOURCE),
        script_sha256=A.hashed(Path(__file__).resolve()), joint_reference_sha256=A.hashed(reference),
        routes_sha256=A.hashed(route_file),
        generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
            [A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_shaft_collar.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py']},
        GH_retention_revision=dict(retaining_depth_mm=6., prior_depth_mm=joint['cap_depth_mm'],
            ball_radius_mm=joint['ball_radius_mm'], cup_inner_radius_mm=joint['socket_inner_radius_mm'],
            nominal_unslotted_lip_radius_mm=float(np.sqrt(joint['socket_inner_radius_mm']**2-6.**2)),
            nominal_radial_capture_margin_mm=float(joint['ball_radius_mm']-np.sqrt(joint['socket_inner_radius_mm']**2-6.**2)),
            original_fastener_bearing_depth_mm=3., extension_start_overlap_mm=2.9,
            actual_axis_source=n.tolist(), cap_changes=changes),
        changed_object_count=2, changed_object_material_intersections=overlaps, actuation_wire_checks=wire_checks,
        scapula_insertion_relief=parent['scapula_insertion_relief'], keeper_definitions=parent['keeper_definitions'],
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners=parent['intended_body_owners'],
        stage_pass=False, actual_drive_range_deg=None,
        not_verified=['actual axial capture and all new cap/humerus/fastener assembly paths must be checked',
                      'arbitrary escape under coupled translations/rotations and true dynamic stability',
                      'printed clearance, materials, wear, bearing and cap/fastener strength',
                      'no inherited joint range or flexible actuation pass'])
    (DEST/'assembly_check.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Cap changes', changes, flush=True)
    print('New overlaps >.001', [r for r in overlaps if r['volume_mm3'] > .001], flush=True)
    print('Active wire failures', sum(bool(w['hits_over_threshold']) for w in wire_checks), flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
