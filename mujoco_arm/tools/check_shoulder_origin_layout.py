"""S1 origin-side seats: the six muscular origins on the scapula and clavicle.

§5.5 laid out the four humerus-side seats.  The six origin-side seats - five on the
scapula, one on the clavicle - had never been laid out at all, and §5.7 needs them
to bolt the end pads onto.  This script closes that gap with the same envelope and
the same thresholds as §5.5.

Each seat is a R11 x 4 mm cylinder at the area-weighted bone-surface centroid of the
muscle's origin candidate region, grown outward along the nearest bone triangle's
normal.  Checks:

1. seat vs seat - can six seats coexist;
2. seat vs bone - a flat disc on a curved bone digs in at the rim when the surface
   is concave; the intrusion volume is reported per seat rather than assumed zero;
3. seat vs another muscle's origin - a seat must not swallow a neighbouring
   muscle's attachment;
4. seat vs the moving humerus - neutral-pose clearance only, the swept check is
   S1-G1.

Candidate regions are proximity samples of a body-parts model, not confirmed
attachment sites, and the envelope has no chamfer, bracket or bolt head.

Run:
    python mujoco_arm/tools/check_shoulder_origin_layout.py
    python mujoco_arm/tools/check_shoulder_origin_layout.py --json out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import attachment_regions
import frames
import interference as I
import mesh_geom
import shoulder_pads

BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'
BONES = {'scapula': '肩胛骨.obj', 'clavicle': '锁骨.obj', 'humerus': '肱骨.obj'}

MUSCLES = {
    'FJ1468': '三角肌前束',
    'FJ1513': '三角肌后束',
    'FJ1467': '三角肌中束',
    'FJ1507': '大圆肌',
    'FJ1504': '肩胛下肌',
    'FJ1500': '冈下肌',
}
DELTOID = ('FJ1468', 'FJ1513', 'FJ1467')


def load_bone(name):
    V, F = mesh_geom.load_obj(str(BONE_DIR / BONES[name]))
    W, G = mesh_geom.weld(V, F)
    vol = mesh_geom.signed_volume(W, G)
    if vol <= 0:
        raise SystemExit('%s has a negative signed volume (%.1f); its winding is '
                         'inward, so the seat normals cannot be trusted' % (name, vol))
    return W, G


def nearest_normal(W, G, point):
    """Outward normal of the bone triangle closest to `point`.

    The mesh winding is the authority: all three bone meshes are welded closed
    solids with a positive signed volume, so the triangle normal already points
    out of the bone.  An earlier version oriented the normal away from the bone's
    overall vertex centroid, which silently inverted the seat on the clavicle and
    on the acromion - for an elongated or small-process bone the vertex centroid
    is nowhere near the contact point, so that test says nothing.  Measured
    consequence: the clavicle seat was pushed into the bone at every offset, and
    the acromion seat reported 817 mm3 of intrusion instead of 115 mm3.
    """
    tri = W[G]
    best = (np.inf, None)
    for i in range(0, len(tri), 4096):
        chunk = tri[i:i + 4096]
        d = mesh_geom.point_triangle_distance(point, chunk)
        k = int(np.argmin(d))
        if d[k] < best[0]:
            best = (float(d[k]), chunk[k])
    face = best[1]
    n = np.cross(face[1] - face[0], face[2] - face[0])
    n = n / max(np.linalg.norm(n), 1e-12)
    return n, best[0]


def points_in_solid(solid, points, r=0.05):
    """How many of `points` fall inside `solid` (tiny cube probe each)."""
    mf = I._manifold()
    cube = mf.Manifold.cube([2 * r] * 3)
    n = 0
    for p in points:
        probe = cube.translate(list(np.asarray(p, float) - r))
        if I.intersection_volume(probe, solid) > I.HARD_MM3:
            n += 1
    return n


def lift_to_clear(centre, normal, bone_solid, limit=8.0):
    """Smallest outward offset along the normal with no bone intrusion [mm].

    A flat disc tangent at a curved bone either floats on the high points or digs
    in at the rim; the lift says how much floating the flat form costs.  Returns
    (lift, residual_intrusion_at_lift).
    """
    def intrusion(offset):
        solid = shoulder_pads.build_pad(np.asarray(centre, float)
                                        + np.asarray(normal, float) * offset, normal)
        return I.intersection_volume(solid, bone_solid)

    if intrusion(0.0) <= I.HARD_MM3:
        return 0.0, 0.0
    lo, hi = 0.0, limit
    if intrusion(hi) > I.HARD_MM3:
        return None, float(intrusion(hi))
    while hi - lo > 0.02:
        mid = 0.5 * (lo + hi)
        if intrusion(mid) <= I.HARD_MM3:
            hi = mid
        else:
            lo = mid
    return hi, 0.0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    gh = frames.joint_centres('source')['GH']
    bones = {name: load_bone(name) for name in BONES}
    print('bones: %s' % ', '.join('%s %d tris' % (k, len(v[1])) for k, v in bones.items()))

    seats = []
    print('\n--- origin seats (R%.1f x %.1f mm) ---'
          % (shoulder_pads.PAD_RADIUS_MM, shoulder_pads.PAD_THICKNESS_MM))
    print('%-8s %-10s %-8s %8s %8s  %s'
          % ('muscle', 'name', 'bone', 'area', 'flat', 'centre source mm'))
    for mid, name in MUSCLES.items():
        region = attachment_regions.load_region(mid, 'origin')
        if region is None:
            print('%-8s %-10s  origin region missing' % (mid, name))
            continue
        centre, area = attachment_regions.centroid(region)
        W, G = bones[region['bone']]
        normal, surf_dist = nearest_normal(W, G, centre)
        solid = shoulder_pads.build_pad(centre, normal)
        far = shoulder_pads.pad_far_radius(centre, normal, gh)
        seat = {'id': mid, 'name': name, 'bone': region['bone'],
                'centre_source_mm': np.round(centre, 3).tolist(),
                'normal_source': np.round(normal, 4).tolist(),
                'centre_to_surface_mm': round(surf_dist, 3),
                'region_area_mm2': round(area, 1),
                'region_samples': region['n_samples'],
                'far_radius_about_gh_mm': round(far, 2),
                'origin_points': region['points'],
                '_solid': solid}
        seats.append(seat)
        print('%-8s %-10s %-8s %8.1f %8.2f  %s'
              % (mid, name, region['bone'], area, surf_dist, np.round(centre, 1)))

    print('\n--- seat vs seat ---')
    print('%-9s %-9s %-9s %-9s %s' % ('pair', 'd,mm', 'gap,mm', 'overlap', 'flag'))
    pairs = []
    worst = (1e9, None)
    for i in range(len(seats)):
        for j in range(i + 1, len(seats)):
            a, b = seats[i], seats[j]
            d = float(np.linalg.norm(np.asarray(a['centre_source_mm'])
                                     - np.asarray(b['centre_source_mm'])))
            vol = I.intersection_volume(a['_solid'], b['_solid'])
            gap = d - 2 * shoulder_pads.PAD_RADIUS_MM
            flag = 'OVERLAP' if vol > I.HARD_MM3 else ('tight' if gap < 2.0 else 'ok')
            pairs.append({'a': a['id'], 'b': b['id'], 'centre_distance_mm': round(d, 2),
                          'centre_gap_mm': round(gap, 2), 'overlap_mm3': round(vol, 4),
                          'pass': bool(vol <= I.HARD_MM3)})
            if d < worst[0]:
                worst = (d, (a['id'], b['id']))
            print('%-9s %-9.1f %-9.1f %-9.3f %s'
                  % ('%s/%s' % (a['id'][2:], b['id'][2:]), d, gap, vol, flag))
    clashes = [p for p in pairs if not p['pass']]
    print('pairs overlapping: %d / %d   closest %s at %.1f mm'
          % (len(clashes), len(pairs), '/'.join(worst[1]), worst[0]))

    print('\n--- seat vs bone ---')
    print('%-8s %-9s %-14s %-10s %s'
          % ('seat', 'own bone', 'intrusion mm3', 'lift mm', 'other bones'))
    bone_rows = []
    for s in seats:
        W, G = bones[s['bone']]
        own_solid = I._manifold().Manifold(I._manifold().Mesh64(
            W, G.astype(np.uint32)))
        own = I.intersection_volume(s['_solid'], own_solid)
        lift, residual = lift_to_clear(s['centre_source_mm'], s['normal_source'],
                                       own_solid)
        others = {}
        for name, (W2, G2) in bones.items():
            if name == s['bone']:
                continue
            solid = I._manifold().Manifold(I._manifold().Mesh64(W2, G2.astype(np.uint32)))
            v = I.intersection_volume(s['_solid'], solid)
            if v > I.HARD_MM3:
                others[name] = round(v, 3)
        bone_rows.append({'id': s['id'], 'own_bone': s['bone'],
                          'own_bone_intrusion_mm3': round(own, 3),
                          'lift_to_clear_mm': None if lift is None else round(lift, 2),
                          'residual_at_lift_mm3': round(residual, 3),
                          'other_bone_intrusion_mm3': others})
        print('%-8s %-9s %-14.3f %-10s %s'
              % (s['id'], s['bone'], own,
                 'none' if lift is None else '%.2f' % lift,
                 ', '.join('%s %.3f' % kv for kv in others.items()) or 'none'))

    print('\n--- seat vs another muscle\'s origin (fraction of its samples swallowed) ---')
    print('%-8s %-28s %s' % ('seat', 'other muscle', 'swallowed'))
    swallow_rows = []
    for s in seats:
        hits = {}
        for other in seats:
            if other['id'] == s['id']:
                continue
            n = points_in_solid(s['_solid'], other['origin_points'])
            if n:
                frac = n / max(1, len(other['origin_points']))
                hits[other['id']] = {'points': int(n),
                                     'of': int(len(other['origin_points'])),
                                     'fraction': round(frac, 3)}
        swallow_rows.append({'id': s['id'], 'swallows': hits})
        if hits:
            print('%-8s %-28s %s'
                  % (s['id'], '', ', '.join('%s %d/%d' % (k, v['points'], v['of'])
                                            for k, v in hits.items())))

    print('\n--- seat vs the moving humerus (neutral pose only) ---')
    hum = I._manifold().Manifold(I._manifold().Mesh64(bones['humerus'][0],
                                                      bones['humerus'][1].astype(np.uint32)))
    hum_rows = []
    for s in seats:
        vol = I.intersection_volume(s['_solid'], hum)
        gap = I.gap(s['_solid'], hum) if vol <= I.HARD_MM3 else 0.0
        hum_rows.append({'id': s['id'], 'overlap_mm3': round(vol, 3),
                         'gap_mm': round(gap, 2),
                         'pass': bool(vol <= I.HARD_MM3)})
        print('  %-8s overlap %7.3f mm3   gap >= %.2f mm' % (s['id'], vol, gap))
    hum_bad = [r for r in hum_rows if not r['pass']]
    print('seats colliding with the humerus at neutral: %s'
          % (', '.join(r['id'] for r in hum_bad) or 'none'))

    delt = [s for s in seats if s['id'] in DELTOID]
    delt_seats_ok = all(p['pass'] for p in pairs
                        if p['a'] in DELTOID and p['b'] in DELTOID)
    print('\n--- note on D8 ---')
    print('the three deltoid origin seats keep their own pads: %s '
          '(pairwise %s)' % ('yes' if delt_seats_ok else 'NO',
                             'clear' if delt_seats_ok else 'OVERLAP'))

    report = {
        'scope': 'R11 x 4 mm seat envelopes at the six origin candidate centroids, '
                 'neutral pose; no chamfer, bracket, bolt head, pad+cap stack or '
                 'muscle belly',
        'pad': {'radius_mm': shoulder_pads.PAD_RADIUS_MM,
                'thickness_mm': shoulder_pads.PAD_THICKNESS_MM,
                'bolt_pitch_mm': shoulder_pads.BOLT_PITCH_MM,
                'hole_radius_mm': shoulder_pads.BOLT_HOLE_RADIUS_MM},
        'seats': [{k: v for k, v in s.items() if not k.startswith('_')
                   and k != 'origin_points'} for s in seats],
        'pairs': pairs,
        'overlapping_pairs': len(clashes),
        'closest_pair_mm': round(worst[0], 2),
        'seats_feasible': not clashes,
        'bone_intrusion': bone_rows,
        'origin_swallowing': swallow_rows,
        'humerus_neutral': hum_rows,
        'deltoid_origin_seats_independent': bool(delt_seats_ok),
        'not_verified': ['candidate regions are proximity samples, not confirmed '
                         'attachment sites',
                         'the pad+cap stack is 8.15 mm tall, this checks 4 mm',
                         'no swept pose (S1-G1), no muscle envelope',
                         'no clamp force, no fixture interaction (§5.9)',
                         'no check against the §5.9 clamp zones, which are named '
                         'regions rather than geometry'],
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
