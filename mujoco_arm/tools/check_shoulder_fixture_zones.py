"""S1 fixture: can the four grip zones actually be reached from a backboard?

§5.9.3 puts a rigid backboard medial/posterior of the scapula and hangs the clamps off
it.  That section was written as prose; this script tests the consequence, because the
scapula is not a flat plate and one of the four zones sits on the far side of it.

Grip zones are defined mechanically here (no hand-picked landmasses):

    medial_mid      the medial 5% vertex band, z in the middle of that band
    superior_angle  the same band, top 8 mm
    coracoid        the 6 mm patch around the most anterior scapula vertex
    lateral_upper   the most lateral surface point 45-65% of the way from the inferior
                    angle to the GH centre

For each zone this script builds a straight 4 mm post from the zone to the backboard
plane and reports what it runs into: bone, a muscle belly, or nothing.  A zone that
cannot see the board is not a clamp location, however good it looks on paper.

Run:
    python mujoco_arm/tools/check_shoulder_fixture_zones.py
    python mujoco_arm/tools/check_shoulder_fixture_zones.py --json out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frames
import interference as I
import mesh_geom

BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'
MUSCLE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / 'muscle'
MUSCLES = {'FJ1468': '三角肌前束', 'FJ1513': '三角肌后束', 'FJ1467': '三角肌中束',
           'FJ1507': '大圆肌', 'FJ1504': '肩胛下肌', 'FJ1500': '冈下肌'}

POST_RADIUS_MM = 4.0
BOARD_CLEARANCE_MM = 2.0
#: the post starts outside the clamp pad, not on the bone surface (the pad itself
#: sits on the bone; a post modelled from the surface only measures its own grazing)
POST_START_MM = 4.0
#: a straight post clipping less than this is grazing a rounded edge, not blocked
GRAZE_MM3 = 50.0


def mf():
    return I._manifold()


def load(name, folder):
    V, F = mesh_geom.load_obj(str(folder / (name + '.obj')))
    W, G = mesh_geom.weld(V, F)
    if mesh_geom.signed_volume(W, G) <= 0:
        raise SystemExit('%s has inward winding' % name)
    return W, G, mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))


def nearest_normal(W, G, point):
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
    return n / max(np.linalg.norm(n), 1e-12), best[0]


def rod(a, b, r, segments=32):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = b - a
    L = float(np.linalg.norm(d))
    if L < 1e-9:
        return None
    cyl = mf().Manifold.cylinder(L, r, r, segments)
    z = d / L
    helper = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    x = np.cross(helper, z)
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    return cyl.transform(np.c_[np.stack([x, y, z], axis=1), a])


def zones(W, GH):
    """The four grip zones, with the criterion that produced each one."""
    out = []
    x = W[:, 0]
    medial = W[x >= np.percentile(x, 95)]
    zlo, zhi = medial[:, 2].min(), medial[:, 2].max()
    mid = medial[(medial[:, 2] >= zlo + 0.45 * (zhi - zlo))
                 & (medial[:, 2] <= zlo + 0.65 * (zhi - zlo))]
    out.append(('medial_mid', '内侧缘中段', mid.mean(0),
                'medial 5%% band, z in %.0f-%.0f' % (zlo + 0.45 * (zhi - zlo),
                                                    zlo + 0.65 * (zhi - zlo))))
    top = medial[medial[:, 2] >= zhi - 8.0]
    out.append(('superior_angle', '上角', top.mean(0),
                'medial 5%% band, top 8 mm (z >= %.0f)' % (zhi - 8.0)))

    tip = W[int(np.argmin(W[:, 1]))]
    patch = W[np.linalg.norm(W - tip, axis=1) <= 6.0]
    out.append(('coracoid', '喙突', patch.mean(0),
                '6 mm patch around the most anterior vertex'))

    inferior = W[int(np.argmin(W[:, 2]))]
    a, b = inferior, GH
    t = np.clip(((W - a) @ (b - a)) / ((b - a) @ (b - a)), 0, 1)
    proj = a + t[:, None] * (b - a)
    off = np.linalg.norm(W - proj, axis=1)
    band = (t >= 0.45) & (t < 0.65) & (off < 12.0)
    q = W[band]
    lateral = q[int(np.argmin(q[:, 0]))]
    patch = W[np.linalg.norm(W - lateral, axis=1) <= 8.0]
    out.append(('lateral_upper', '外侧缘上部', patch.mean(0),
                'most lateral surface point 45-65%% along inferior angle -> GH'))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    GH = frames.joint_centres('source')['GH']
    sW, sG, scap = load('肩胛骨', BONE_DIR)
    cW, cG, clav = load('锁骨', BONE_DIR)
    hW, hG, hum = load('肱骨', BONE_DIR)
    bellies = {}
    for mid in MUSCLES:
        _, _, solid = load(mid, MUSCLE_DIR)
        bellies[mid] = solid
    print('bones and %d bellies loaded' % len(bellies))

    dorsal = max(sW[:, 1].max(), cW[:, 1].max(),
                 max(b.bounding_box()[4] for b in bellies.values()))
    board_y = float(dorsal + BOARD_CLEARANCE_MM)
    anterior = min(sW[:, 1].min(), cW[:, 1].min(),
                   min(b.bounding_box()[1] for b in bellies.values()))
    front_y = float(anterior - BOARD_CLEARANCE_MM)
    print('dorsal-most point of scapula %.2f, clavicle %.2f, bellies %.2f mm'
          % (sW[:, 1].max(), cW[:, 1].max(),
             max(b.bounding_box()[4] for b in bellies.values())))
    print('anterior-most point  scapula %.2f, clavicle %.2f, bellies %.2f mm'
          % (sW[:, 1].min(), cW[:, 1].min(),
             min(b.bounding_box()[1] for b in bellies.values())))
    print('backboard plane y = %.2f (dorsal + %.1f);  front bar plane y = %.2f '
          '(anterior - %.1f)' % (board_y, BOARD_CLEARANCE_MM, front_y,
                                 BOARD_CLEARANCE_MM))

    print('\n--- grip zones and their reach ---')
    print('%-16s %-9s %-6s %-26s %8s %8s %s'
          % ('zone', 'side', 'n.y', 'centre source mm', 'post mm', 'bone', 'belly'))
    rows = []
    for key, label, centre, criterion in zones(sW, GH):
        normal, surf = nearest_normal(sW, sG, centre)
        for side, target_y, sign in (('dorsal', board_y, +1),
                                     ('anterior', front_y, -1)):
            start = np.array([centre[0], centre[1] + sign * POST_START_MM, centre[2]])
            target = np.array([centre[0], target_y, centre[2]])
            length = float(np.linalg.norm(target - start))
            post = rod(start, target, POST_RADIUS_MM) if length > 0.5 else None
            hits = {}
            for name, solid in (('scapula', scap), ('clavicle', clav), ('humerus', hum)):
                v = I.intersection_volume(post, solid) if post is not None else 0.0
                if v > I.HARD_MM3:
                    hits[name] = round(v, 2)
            belly_hits = {}
            for mid, solid in bellies.items():
                v = I.intersection_volume(post, solid) if post is not None else 0.0
                if v > I.HARD_MM3:
                    belly_hits[mid] = round(v, 2)
            rows.append({'zone': key, 'label': label, 'criterion': criterion,
                         'side': side,
                         'centre_source_mm': np.round(centre, 2).tolist(),
                         'normal_source': np.round(normal, 4).tolist(),
                         'normal_y': round(float(normal[1]), 3),
                         'distance_to_surface_mm': round(surf, 2),
                         'post_start_mm': POST_START_MM,
                         'post_length_mm': round(length, 2),
                         'post_radius_mm': POST_RADIUS_MM,
                         'bone_hits_mm3': hits, 'belly_hits_mm3': belly_hits,
                         'graze_threshold_mm3': GRAZE_MM3,
                         'blocked_by': [k for k, v in list(hits.items())
                                        + list(belly_hits.items()) if v > GRAZE_MM3],
                         'reachable': all(v <= GRAZE_MM3 for v in
                                          list(hits.values()) + list(belly_hits.values()))})
            print('%-16s %-9s %-6.2f %-26s %8.1f %8s %s'
                  % (key, side, normal[1], np.round(centre, 1), length,
                     ', '.join('%s %.0f' % kv for kv in hits.items()) or 'clear',
                     ', '.join('%s %.0f' % kv for kv in belly_hits.items()) or 'clear'))
        print('%-16s %-8.3f %-26s %8.1f %8s %s'
              % (key, normal[1], np.round(centre, 1), length,
                 ', '.join('%s %.0f' % kv for kv in hits.items()) or 'clear',
                 ', '.join('%s %.0f' % kv for kv in belly_hits.items()) or 'clear'))

    print('\n--- summary ---')
    reach = [r for r in rows if r['reachable']]
    print('straight posts that reach a plane: %d / %d' % (len(reach), len(rows)))
    for r in rows:
        if not r['reachable']:
            print('  BLOCKED %-16s %-9s bone %s  belly %s'
                  % (r['zone'], r['side'], r['bone_hits_mm3'] or '-',
                     r['belly_hits_mm3'] or '-'))
    by_zone = {}
    for r in rows:
        by_zone.setdefault(r['zone'], []).append(r['side'] if r['reachable'] else None)
    for zone, sides in by_zone.items():
        good = [s for s in sides if s]
        print('  %-16s reachable from: %s' % (zone, ', '.join(good) or 'neither'))

    report = {
        'scope': 'straight 4 mm posts from four mechanically defined grip zones to a '
                 'dorsal backboard plane and to an anterior bar plane; no clamp body, '
                 'jaw, bolt, bent arm or rail modelled',
        'backboard': {'dorsal_plane_y_mm': round(board_y, 3),
                      'anterior_plane_y_mm': round(front_y, 3),
                      'clearance_mm': BOARD_CLEARANCE_MM,
                      'basis': 'dorsal-most / anterior-most of scapula, clavicle and '
                               'the six muscle bellies, plus the clearance'},
        'post_radius_mm': POST_RADIUS_MM,
        'post_start_mm': POST_START_MM,
        'posts': rows,
        'reachable': ['%s/%s' % (r['zone'], r['side']) for r in reach],
        'not_verified': ['clamp jaws and the two-sided grip are not modelled, only the '
                         'reach of a straight post',
                         'bent arms and rails are not modelled; a blocked straight post '
                         'does not prove the zone is unusable',
                         'no humerus sweep against the clamps (S1-G1)',
                         'no clamp force, no fixture compliance (§5.9.4)'],
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
