"""Moment and chord diagnostics for a muscle's line of action.

Convention (identical to the project's earlier endpoint studies):

    force direction   f = unit(proximal - distal)      the wire pulls distally
    moment arm        m = cross(distal - centre, f) . axis      [mm per unit force]
    positive m        rotation about +axis
    the moving end rotates as  distal(t) = centre + R(axis, t) @ (distal0 - centre)

This is a *geometric* screen.  It says nothing about achievable force, TPU
stiffness or friction.

Run this file directly for the analytic self-test, or:
    python mujoco_arm/tools/force_line.py --muscle FJ1486 --to 150 --step 5
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frames
import mesh_geom

INTERFACES_PATH = frames.OUTPUT_DIR / 'elbow_interfaces.json'
HUMERUS_OBJ = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼' / '肱骨.obj'


# ------------------------------------------------------------------ primitives
def rotate_about(point, centre, axis, degrees):
    """Rodrigues rotation of `point` about (centre, unit axis)."""
    point = np.asarray(point, float)
    centre = np.asarray(centre, float)
    axis = frames.unit(axis)
    t = np.radians(degrees)
    v = point - centre
    cos_t, sin_t = np.cos(t), np.sin(t)
    rotated = v * cos_t + np.cross(axis, v) * sin_t + axis * (axis @ v) * (1.0 - cos_t)
    return centre + rotated


def moment_arm(centre, axis, proximal, distal):
    """Signed moment arm [mm] of a wire running proximal->distal about `axis`."""
    centre = np.asarray(centre, float)
    proximal = np.asarray(proximal, float)
    distal = np.asarray(distal, float)
    force = frames.unit(proximal - distal)
    return float(np.cross(distal - centre, force) @ frames.unit(axis))


# -------------------------------------------------------------------- profiles
def moment_profile(axis, centre, proximal, distal0, angles_deg):
    """Moment arm, chord length and chord direction over a rotation sweep."""
    angles = np.asarray(angles_deg, float)
    proximal = np.asarray(proximal, float)
    moment = np.empty_like(angles)
    chord = np.empty_like(angles)
    distal_pts = np.empty((len(angles), 3))
    for i, a in enumerate(angles):
        d = rotate_about(distal0, centre, axis, a)
        distal_pts[i] = d
        moment[i] = moment_arm(centre, axis, proximal, d)
        chord[i] = float(np.linalg.norm(proximal - d))
    return {'angles_deg': angles, 'moment_mm': moment, 'chord_mm': chord,
            'distal_points': distal_pts, 'proximal': proximal,
            'centre': np.asarray(centre, float), 'axis': frames.unit(axis)}


def sign_changes(profile):
    """Indices where the moment arm changes sign, with a linear zero crossing."""
    m = profile['moment_mm']
    a = profile['angles_deg']
    out = []
    for i in range(len(m) - 1):
        if m[i] == 0.0:
            out.append({'angle_deg': float(a[i]), 'exact': True})
        elif m[i] * m[i + 1] < 0:
            t = abs(m[i]) / (abs(m[i]) + abs(m[i + 1]))
            out.append({'angle_deg': float(a[i] + t * (a[i + 1] - a[i])), 'exact': False})
    return out


def summarize(profile):
    m = profile['moment_mm']
    c = profile['chord_mm']
    i_min = int(np.argmin(np.abs(m)))
    return {
        'moment_at_start_mm': float(m[0]),
        'moment_at_end_mm': float(m[-1]),
        'sign_changes': sign_changes(profile),
        'min_abs_moment_mm': float(abs(m[i_min])),
        'min_abs_moment_at_deg': float(profile['angles_deg'][i_min]),
        'max_abs_moment_mm': float(np.abs(m).max()),
        'chord_start_mm': float(c[0]),
        'chord_end_mm': float(c[-1]),
        'chord_change_mm': float(c[-1] - c[0]),
        'chord_change_pct': float(100.0 * (c[-1] - c[0]) / c[0]),
    }


# --------------------------------------------------------------- bone clearance
class BoneObstacle:
    """Triangle soup in SOURCE coordinates; queries take DISPLAY points."""

    def __init__(self, obj_path=HUMERUS_OBJ):
        self.verts, self.faces = mesh_geom.load_obj(obj_path)
        self.tri = self.verts[self.faces]

    def clearance_display(self, a_display, b_display):
        a = frames.display_to_source(a_display)
        b = frames.display_to_source(b_display)
        return mesh_geom.segment_min_distance(a, b, self.tri)


def bone_contact_profile(profile, obstacle, clearance_mm=0.01):
    """First angle at which the straight chord comes within `clearance_mm` of the bone."""
    hits = []
    for i, a in enumerate(profile['angles_deg']):
        d = obstacle.clearance_display(profile['proximal'], profile['distal_points'][i])
        hits.append(d)
    hits = np.asarray(hits)
    first = None
    idx = np.where(hits <= clearance_mm)[0]
    if idx.size:
        first = float(profile['angles_deg'][idx[0]])
    i_min = int(np.argmin(hits))
    return {'clearance_mm': hits, 'first_contact_deg': first,
            'min_clearance_mm': float(hits[i_min]),
            'min_clearance_deg': float(profile['angles_deg'][i_min])}


# --------------------------------------------------------------------- self test
def self_test():
    ok = True

    # 1) quarter turn about +z
    p = rotate_about(np.array([1.0, 0.0, 0.0]), np.zeros(3), np.array([0, 0, 1.0]), 90.0)
    err = float(np.abs(p - np.array([0.0, 1.0, 0.0])).max())
    print('[%s] rotate 90 deg about +z: %s' % ('PASS' if err < 1e-12 else 'FAIL', np.round(p, 12)))
    ok &= err < 1e-12

    # 2) analytic moment arm: distal (1,0,0), proximal (0,1,0) -> +sqrt(2)/2 about +z
    m = moment_arm(np.zeros(3), np.array([0, 0, 1.0]), np.array([0.0, 1.0, 0.0]),
                   np.array([1.0, 0.0, 0.0]))
    expected = 0.5 * np.sqrt(2.0)
    print('[%s] analytic moment arm %.9f (expected %.9f)'
          % ('PASS' if abs(m - expected) < 1e-12 else 'FAIL', m, expected))
    ok &= abs(m - expected) < 1e-12

    # 3) reversal: swapping which end is fixed flips the sign
    m2 = moment_arm(np.zeros(3), np.array([0, 0, 1.0]), np.array([1.0, 0.0, 0.0]),
                    np.array([0.0, 1.0, 0.0]))
    print('[%s] role swap flips sign: %.6f vs %.6f' % ('PASS' if sign_eq(m2, -m) else 'FAIL', m, m2))
    ok &= sign_eq(m2, -m)

    # 4) zero crossing detection on a synthetic sweep crossing zero
    prof = {'angles_deg': np.array([0.0, 45.0, 90.0]), 'moment_mm': np.array([1.0, 0.0, -1.0]),
            'chord_mm': np.array([1.0, 1.0, 1.0])}
    zc = sign_changes(prof)
    found = any(abs(z['angle_deg'] - 45.0) < 1e-9 for z in zc)
    print('[%s] zero crossing found at 45 deg: %s' % ('PASS' if found else 'FAIL', zc))
    ok &= found

    print('\nRESULT:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


def sign_eq(a, b):
    return (a > 0) == (b > 0)


# ------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--muscle', help='FJ1486 or FJ1477')
    ap.add_argument('--from', dest='start', type=float, default=0.0)
    ap.add_argument('--to', dest='stop', type=float, default=150.0)
    ap.add_argument('--step', type=float, default=1.0)
    ap.add_argument('--bone', action='store_true', help='also check the chord against the humerus')
    ap.add_argument('--clearance', type=float, default=0.01,
                    help='chord-to-bone distance treated as contact, mm')
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    if args.self_test or not args.muscle:
        return self_test()

    centre_src, axis_src = frames.elbow_axis()
    centre = frames.source_to_display(centre_src)

    with open(INTERFACES_PATH, 'r', encoding='utf-8') as fh:
        data = {r['id']: r for r in json.load(fh)['interfaces']}
    if args.muscle not in data:
        ap.error('unknown muscle: ' + args.muscle)
    rec = data[args.muscle]
    proximal = np.asarray(rec['proximal_TPU_pad_centre_display_mm'], float)
    distal = np.asarray(rec['distal_TPU_pad_centre_display_mm'], float)

    angles = np.arange(args.start, args.stop + 1e-9, args.step)
    prof = moment_profile(axis_src, centre, proximal, distal, angles)
    rep = summarize(prof)
    rep['muscle'] = args.muscle
    rep['elbow_axis_source'] = axis_src.tolist()
    rep['centre_display_mm'] = centre.tolist()

    print('muscle %s   sweep %.0f..%.0f deg step %.1f' % (args.muscle, args.start, args.stop, args.step))
    print('  moment@start %.3f mm   moment@end %.3f mm' % (rep['moment_at_start_mm'], rep['moment_at_end_mm']))
    print('  min |moment| %.3f mm at %.1f deg' % (rep['min_abs_moment_mm'], rep['min_abs_moment_at_deg']))
    if rep['sign_changes']:
        for z in rep['sign_changes']:
            print('  SIGN CHANGE at %.2f deg%s' % (z['angle_deg'], '' if z['exact'] else ' (interpolated)'))
    else:
        print('  no sign change in range')
    print('  chord %.3f -> %.3f mm  (%+.3f mm, %+.2f%%)'
          % (rep['chord_start_mm'], rep['chord_end_mm'], rep['chord_change_mm'], rep['chord_change_pct']))

    if args.bone:
        obstacle = BoneObstacle()
        bc = bone_contact_profile(prof, obstacle, args.clearance)
        rep['bone_first_contact_deg'] = bc['first_contact_deg']
        rep['bone_min_clearance_mm'] = bc['min_clearance_mm']
        rep['bone_min_clearance_deg'] = bc['min_clearance_deg']
        if bc['first_contact_deg'] is None:
            print('  chord never within %.2f mm of the humerus (min clearance %.2f mm at %.1f deg)'
                  % (args.clearance, bc['min_clearance_mm'], bc['min_clearance_deg']))
        else:
            print('  FIRST BONE CONTACT at %.1f deg (min clearance %.3f mm at %.1f deg)'
                  % (bc['first_contact_deg'], bc['min_clearance_mm'], bc['min_clearance_deg']))

    if args.json_out:
        with open(args.json_out, 'w', encoding='utf-8') as fh:
            json.dump(rep, fh, ensure_ascii=False, indent=2)
        print('  report written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
