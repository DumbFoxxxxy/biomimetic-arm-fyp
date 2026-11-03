"""Anatomy-conforming scapular clamps with external M3 ears.

Independent candidate; the source bench assembly and all bone/muscle files
are retained. PCA defines a candidate clamp frame, not an anatomical joint
axis. Discrete axial sweeps form liner cavities, not constant normal offsets.
No friction, preload, strength, printability or actuation is validated here.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_bench_mount as B
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export

SOURCE = B.DEST/'assembly_check.json'
DEST = B.DEST/'scapula_clamp_candidate'
SEEDS = dict(medial_mid=[-59.29, -9.22, 1287.22],
             superior_angle=[-62.79, -18.62, 1319.93],
             lateral_upper=[-133.08, -47.19, 1275.04])


def closest(point, triangles):
    distances, points = A.mesh_geom._closest_all(point, triangles)
    index = int(np.argmin(distances))
    return points[index], float(distances[index])


def fasteners(point, direction, length, far_face):
    """Plain M3 candidate envelopes, with nominal bearing-face contacts."""
    mf = A.I._manifold()
    p = np.asarray(point, float)
    n = np.asarray(direction, float)
    n /= np.linalg.norm(n)
    u = np.cross(np.eye(3)[np.argmin(abs(n))], n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    F = np.c_[np.column_stack([u, v, n]), p]
    profile = np.array([[0, -3.5], [2.75, -3.5], [2.75, -.5],
                        [1.5, -.5], [1.5, length-.5], [0, length-.5]])
    screw = mf.Manifold.revolve(mf.CrossSection([profile]), 64).transform(F)
    near = rod(p-n*.5, p, 3.5)-rod(p-n*.6, p+n*.1, 1.7)
    far = rod(p+n*far_face, p+n*(far_face+.5), 3.5)-rod(
        p+n*(far_face-.1), p+n*(far_face+.6), 1.7)
    nut = mf.Manifold.cylinder(2.4, 5.5/np.sqrt(3), 5.5/np.sqrt(3), 6)
    nut -= mf.Manifold.cylinder(2.6, 1.6, 1.6, 48).translate([0, 0, -.1])
    nut = nut.translate([0, 0, far_face+.5]).transform(F)
    return dict(screw=screw, near_washer=near, far_washer=far, nut=nut)


def rounded_plate(ear_u, low, high):
    points = []
    for x, y, angle in [(ear_u+3, 8, 0), (-6, 8, 90),
                         (-6, -8, 180), (ear_u+3, -8, 270)]:
        for degree in np.linspace(angle, angle+90, 9):
            t = np.radians(degree)
            points.append([x+np.cos(t), y+np.sin(t)])
    return A.I._manifold().Manifold.extrude(
        A.I._manifold().CrossSection([np.array(points)]), high-low).translate([0, 0, low])


def clearance(body, delta=.15):
    shapes = [body]+[body.translate(sign*np.eye(3)[j]*delta)
                     for j in range(3) for sign in [-1, 1]]
    mf = A.I._manifold()
    return mf.Manifold.batch_boolean(shapes, mf.OpType.Add)


def clamp_geometry(tag, bone, bone_tri, frame, inward=5., ear_u=12.):
    mf = A.I._manifold()
    seed = np.array(SEEDS[tag])
    nearby = bone_tri.reshape(-1, 3)
    nearby = nearby[np.linalg.norm(nearby-seed, axis=1) < 12.]
    values, vectors = np.linalg.eigh(np.cov(nearby.T))
    n = vectors[:, 0]
    if n[1] < 0:
        n = -n
    u = np.array([-1., 0, 0]) if tag == 'lateral_upper' else np.array([1., 0, 0])
    u -= n*np.dot(u, n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    F = np.column_stack([u, v, n])
    g = seed-inward*u
    hits = bone.ray_cast(g+n*20, g-n*20)
    if len(hits) != 2:
        raise ValueError('Clamp centre needs an explicit single blade interval: '+tag)
    crossing = [float(np.dot(np.array(h.position)-g, n)) for h in hits]
    stock = abs(crossing[0]-crossing[1])
    g = g+n*np.mean(crossing)
    local_bone = bone.transform(np.c_[F.T, -F.T@g])
    patch = local_bone ^ B.box([-7, -9, -15], [6, 9, 15])
    bb = patch.bounding_box()
    outer_front, outer_back = bb[5]+4.5, bb[2]-4.5
    if outer_front > 20 or outer_back < -20:
        raise ValueError('Unexpected blade thickness/curvature: '+tag)
    # Sweep at 0.25 mm increments. It is an axial cavity envelope, not a
    # verified uniform wall thickness or material deformation operation.
    fronts = [local_bone.translate([0, 0, z]) for z in np.arange(0, 1.501, .25)]
    backs = [local_bone.translate([0, 0, -z]) for z in np.arange(0, 1.501, .25)]
    front_pocket = mf.Manifold.batch_boolean(fronts, mf.OpType.Add)
    back_pocket = mf.Manifold.batch_boolean(backs, mf.OpType.Add)
    front = rounded_plate(ear_u, .6, outer_front)-front_pocket
    back = rounded_plate(ear_u, outer_back, -.6)-back_pocket
    liners = dict(front_liner=(front_pocket-local_bone)^B.box([-6, -6, 0], [6, 6, outer_front]),
                  back_liner=(back_pocket-local_bone)^B.box([-6, -6, outer_back], [6, 6, 0]))
    # A conforming cavity alone can lock a rigid jaw behind a local bone
    # undercut. To withdraw the jaw along -n, sweep the bone/liner pocket
    # in the opposite relative direction (+n) and open that cavity. Keep
    # the original conforming soft liner: no fictitious liner adhesion or
    # rigid liner-to-jaw motion is needed to obtain a clear path.
    local_stock = back_pocket ^ B.box([-8, -10, outer_back-1], [ear_u+5, 10, outer_front+2])
    # Use a continuous solid-prism Minkowski sweep for the rigid opening.
    # Discrete translated copies left a small undercut at superior-angle
    # liner edges. The 0.02 mm transverse clearance is explicit geometry,
    # not a claim that a printer can resolve or maintain that tolerance.
    opening = local_stock.minkowski_sum(B.box([-.02, -.02, 0], [.02, .02, 20]))
    back -= opening
    # No whole-pair exclusion: rigid jaws must not penetrate bone even where
    # the sampled axial envelope is not a true swept volume.
    front -= local_bone
    back -= local_bone
    parts = {}
    bolt_rows = []
    for i, y in enumerate([-5., 5.]):
        axis = g+F@np.array([ear_u, y, 0])
        bore = rod([ear_u, y, outer_back-6], [ear_u, y, outer_front+6], 1.7)
        front -= bore
        back -= bore
        q = axis+n*outer_front
        grip = outer_front-outer_back
        length = next(x for x in [16., 20., 25., 30., 35., 40., 45.] if x >= grip+4.)
        hardware = fasteners(q, -n, length, grip)
        parts.update({f'{key}_{i}':s for key, s in hardware.items()})
        bolt_rows.append(dict(index=i, axis_source_mm=axis.tolist(),
                              near_bearing_source_mm=q.tolist(), direction_source=(-n).tolist(),
                              grip_mm=float(grip), nominal_length_mm=length,
                              nominal_tip_beyond_nut_mm=length-.5-grip-.5-2.4))
    parts.update({key:s.transform(np.c_[F, g]) for key, s in liners.items()})
    parts['moving_jaw'] = back.transform(np.c_[F, g])
    front = front.transform(np.c_[F, g])
    # Put the connecting stem on the external ear side, rather than drawing
    # a chord through the blade. It joins a real nearby surface of the frame.
    origin = g+F@np.array([ear_u-2., 0, outer_front-1.])
    fv, ff = A.mesh_geom.manifold_mesh_to_arrays(frame.to_mesh64())
    target, distance = closest(origin, fv[ff])
    direction = (target-origin)/np.linalg.norm(target-origin)
    stem = rod(origin, target+direction*1., 3.)
    fixed = front+stem+frame
    # Open the receiver for the actual jaw withdrawal, rather than only
    # clearing its zero-pose cross-section. The continuous local sweep
    # includes 0.15 mm clearance in all three local directions.
    receiver_opening = back.minkowski_sum(B.box([-.15, -.15, -20.15], [.15, .15, .15]))
    fixed -= receiver_opening.transform(np.c_[F, g])
    # The old beam occupies the open/rear side of the new clamp. Make a
    # real receiver relief, retaining the front jaw's liner-bearing face.
    for key, body in parts.items():
        if key in ['front_liner', 'moving_jaw']:
            continue
        # Preserve the nominal near bearing plane while trimming original
        # frame material that projects into the actual washer thickness.
        fixed -= body if key.startswith('near_washer') else clearance(body)
    for bolt in bolt_rows:
        q = np.array(bolt['axis_source_mm'])
        fixed -= rod(q-n*(abs(outer_back)+6), q+n*(outer_front+6), 1.7)
    obsolete_terminal = []
    positive = sorted([p for p in fixed.decompose() if p.volume() > 1e-6],
                      key=lambda p:-p.volume())
    # One identified old rear beam terminal is superseded by the new
    # front receiver/stem. This is an explicit fixture redesign, not a
    # general disconnected-material cleanup: only this original-frame
    # piece below the moving jaw, in the measured local zone, is eligible.
    if tag == 'medial_mid' and len(positive) > 1:
        zone = B.box([4, -6, -10], [12, 6, outer_back+.01]).transform(np.c_[F, g])
        if sum(p.volume() for p in positive[1:]) > 35.:
            raise ValueError('Unexpected amount of isolated old receiver material')
        for fragment in positive[1:]:
            not_old = float((fragment-frame).volume())
            outside = float((fragment-zone).volume())
            if not_old > 1e-7 or outside > 1e-7:
                raise ValueError('Receiver fragment is new structure or outside the identified old terminal')
            obsolete_terminal.append(dict(volume_mm3=float(fragment.volume()),
                not_original_frame_mm3=not_old, outside_identified_terminal_zone_mm3=outside,
                bounds_source_mm=list(fragment.bounding_box()),
                reason='old rear beam terminal superseded by connected front jaw/stem; source frame retained'))
            fixed -= fragment
    return fixed, parts, dict(seed_source_mm=seed.tolist(), local_origin_source_mm=g.tolist(),
        basis_columns_source=F.tolist(), method='unweighted triangle-vertex PCA within 12 mm; axis sign +Y; two face ray intersections',
        local_PCA_eigenvalues=values.tolist(), local_PCA_vertex_count=len(nearby), inward_offset_mm=inward,
        centre_ray_stock_mm=stock, centre_ray_intersections_source_mm=[list(h.position) for h in hits],
        liner_axial_offset_mm=1.5, envelope_translation_step_mm=.25,
        rigid_back_jaw_relative_opening_sweep_mm=20.,
        rigid_back_jaw_opening_transverse_clearance_mm=.02,
        rigid_back_jaw_opening_method='continuous Minkowski sum with 0.04x0.04x20 mm prism',
        outer_front_local_mm=float(outer_front), outer_back_local_mm=float(outer_back),
        rigid_jaw_split_gap_outside_bone_mm=1.2, external_ear_u_mm=ear_u,
        bolts=bolt_rows, stem_start_source_mm=origin.tolist(), stem_frame_contact_source_mm=target.tolist(),
        stem_contact_distance_mm=distance, stem_radius_mm=3.,
        receiver_clearance_translation_envelope_mm=.15,
        receiver_jaw_withdrawal_opening_mm=20.,
        superseded_original_frame_terminal=obsolete_terminal,
        combined_receiver_removed_from_original_frame_mm3=float((frame-fixed).volume()),
        combined_receiver_added_to_original_frame_mm3=float((fixed-frame).volume()),
        jaw_stem_union_overlap_mm3=A.I.intersection_volume(front, stem),
        stem_frame_union_overlap_mm3=A.I.intersection_volume(stem, frame))


def run(zones, dest):
    dest.mkdir(parents=True, exist_ok=True)
    parent = json.loads(SOURCE.read_text(encoding='utf8'))
    solids, sources = {}, dict(parent['sources'])
    for key, row in sources.items():
        solids[key], _, _ = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
    bone = solids['bone:scapula']
    row = sources['bone:scapula']
    _, bone_tri, _ = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
    generated, definitions, replacements = {}, {}, {}
    for tag in zones:
        frame_key = 'fixture:lateral_frame' if tag == 'lateral_upper' else 'fixture:medial_frame'
        fixed, parts, definition = clamp_geometry(tag, bone, bone_tri, solids[frame_key],
                                                  ear_u=14. if tag == 'lateral_upper' else 12.)
        solids[frame_key] = fixed
        replacements[frame_key] = fixed
        generated.update({'clamp:'+tag+'/'+k:s for k, s in parts.items()})
        definitions[tag] = definition
    (dest/'geometry_definitions.json').write_text(json.dumps(dict(
        script_sha256=A.hashed(Path(__file__).resolve()), parent_sha256=A.hashed(SOURCE),
        zones=definitions), indent=2)+'\n', encoding='utf8', newline='\n')
    for key, solid in {**replacements, **generated}.items():
        positive = [p.volume() for p in solid.decompose() if p.volume() > 1e-6]
        if len(positive) != 1:
            (dest/'rejected_components.json').write_text(json.dumps(dict(key=key, positive_volumes_mm3=positive), indent=2)+'\n', encoding='utf8')
            raise ValueError('Clamp material disconnected: '+key)
        bb = np.array(solid.bounding_box())
        origin = (bb[:3]+bb[3:])/2
        restored, meta = local_export(solid.as_original().simplify(.001),
                                     dest/(key.replace(':', '_').replace('/', '_')+'.stl'), origin)
        if not meta['mesh_checks_pass']:
            (dest/'rejected_mesh.json').write_text(json.dumps(dict(key=key, meta=meta), indent=2)+'\n', encoding='utf8')
            raise ValueError('Clamp readback fails '+key)
        solids[key], sources[key] = restored, meta
    changed = set(replacements)|set(generated)
    intersections = []
    for a, b in itertools.combinations(solids, 2):
        if a not in changed and b not in changed:
            continue
        volume = A.I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            intersections.append(dict(a=a, b=b, volume_mm3=volume))
    wire_checks = []
    routes_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    for row in json.loads(routes_file.read_text(encoding='utf8'))['routes']:
        wire = A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), .3)
        hits = [dict(part=k, volume_mm3=A.I.intersection_volume(wire, s)) for k, s in solids.items()]
        wire_checks.append(dict(head=row['head'], hits_over_threshold=[h for h in hits if h['volume_mm3'] > .001], intersections=hits))
    report = dict(scope=__doc__, n_instances=len(solids), sources=sources,
                  parent_sha256=A.hashed(SOURCE), parent_report=SOURCE.relative_to(A.ROOT).as_posix(),
                  routes_sha256=A.hashed(routes_file),
                  generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
                      [A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py',
                       A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',
                       A.ROOT/'mujoco_arm/tools/build_shoulder_shaft_collar.py',
                       A.ROOT/'mujoco_arm/tools/audit_shoulder_retained_joint_fit.py']},
                  script_sha256=A.hashed(Path(__file__).resolve()), zones=definitions,
                  changed_object_count=len(changed), changed_object_material_intersections=intersections,
                  actuation_wire_checks=wire_checks,
                  free_terminal_constraints_required=parent['free_terminal_constraints_required'],
                  intended_body_owners={**parent['intended_body_owners'],
                                        **{k:'nominal_scapula_clamp' for k in generated}},
                  unchanged_parent_material_intersections=parent['unchanged_parent_material_intersections'],
                  static_new_pairs_below_threshold=all(r['volume_mm3'] <= .001 for r in intersections),
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['clamp preload, friction, bone/liner stress, assembly and retaining loose liners',
                                'uniform normal wall thickness and physical printability',
                                'clavicle fastening and any omitted scapula clamp zones',
                                'tool, screw and jaw assembly paths in complete anatomy',
                                'same unresolved muscle contact, channel-wall residual and actuation as parent'])
    (dest/'assembly_check.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects;', len(changed), 'new/replaced', flush=True)
    print('New overlaps > .001:', [r for r in intersections if r['volume_mm3'] > .001], flush=True)
    print('Active wire failures:', sum(bool(w['hits_over_threshold']) for w in wire_checks), flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser = argparse.ArgumentParser()
    parser.add_argument('--zones', nargs='+', choices=list(SEEDS), default=list(SEEDS))
    parser.add_argument('--out', type=Path, default=DEST)
    args = parser.parse_args()
    run(args.zones, args.out)
