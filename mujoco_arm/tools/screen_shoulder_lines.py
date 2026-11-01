"""L1 force-line screen: moment-arm sign of the six intrinsic shoulder muscles.

For each muscle the origin/insertion candidate regions in
``attachment_regions/*.npz`` are reduced to area-weighted bone-surface centroids,
then the humerus is rotated about each GH probe axis and the moment arm evaluated.

Two path models are reported:

* **straight** - the chord from origin to insertion (what the elbow studies used);
* **wrapped** - if that chord passes inside the humeral head sphere, the cable is
  routed tangent to the sphere, so the force at the insertion is taken along
  ``B -> T2`` (the outgoing tangent) instead of ``B -> A``.

This is still a SIGN SCREEN, not a force model:

* candidate regions are geometric proximity samples, not confirmed attachment sites;
* no muscle force, pennation, friction, co-contraction or scapular motion is modelled;
* the scapula and clavicle are treated as fixed (consistent with stage L1).

Run:
    python mujoco_arm/tools/screen_shoulder_lines.py
    python mujoco_arm/tools/screen_shoulder_lines.py --json out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import force_line
import frames

REGION_DIR = frames.OUTPUT_DIR / 'attachment_regions'
ANGLES = np.arange(0.0, 121.0, 1.0)
SAMPLE_DEG = (0.0, 60.0, 120.0)

#: muscle -> (display name, origin bone, intended primary axis, intended sign)
#: Signs are calibrated from the anatomy, not fitted to the data:
#:   +flexion  = anterior deltoid direction, +abduction = middle deltoid direction,
#:   -longitudinal = internal rotation (primary internal rotator = subscapularis).
#: The remaining four muscles are then checked against those conventions.
MUSCLES = {
    'FJ1468': ('三角肌前束', 'clavicle', 'flexion', +1),
    'FJ1513': ('三角肌后束', 'scapula', 'flexion', -1),
    'FJ1467': ('三角肌中束', 'scapula', 'abduction', +1),
    'FJ1507': ('大圆肌', 'scapula', 'abduction', -1),
    'FJ1504': ('肩胛下肌', 'scapula', 'longitudinal', -1),
    'FJ1500': ('冈下肌', 'scapula', 'longitudinal', +1),
}

AXES = ('flexion', 'abduction', 'longitudinal')


# ------------------------------------------------------------------ geometry
def region_centroid(mid, kind, bone):
    """Area-weighted centroid of the bone-surface points of one attachment region."""
    path = REGION_DIR / ('%s_%s_%s.npz' % (mid, kind, bone))
    if not path.is_file():
        return None
    with np.load(path) as d:
        points = d['bone_points_source_mm']
        weights = d['muscle_face_areas_mm2']
        keep = d['distances_mm'] <= np.median(d['distances_mm'])
    points, weights = points[keep], weights[keep]
    return (points * weights[:, None]).sum(0) / weights.sum()


def chord_clearance(a, b, centre):
    ab = b - a
    t = float(np.clip((centre - a) @ ab / max(ab @ ab, 1e-12), 0.0, 1.0))
    return float(np.linalg.norm(a + t * ab - centre))


def wrapped_force_direction(a, b, centre, radius):
    """Direction of the cable force at `b` when it wraps a sphere at `centre`."""
    nb = b - centre
    d_b = np.linalg.norm(nb)
    if d_b <= radius * 1.02:
        return None
    nb = nb / d_b
    na = a - centre
    na = na - (na @ nb) * nb
    n = np.linalg.norm(na)
    if n < 1e-9:
        return None
    u = na / n
    alpha = np.arccos(np.clip(radius / d_b, -1.0, 1.0))
    tangent_point = centre + radius * (np.cos(alpha) * nb + np.sin(alpha) * u)
    return frames.unit(tangent_point - b)


def moment_curve(axis, centre, proximal, distal0, angles, force_fn=None):
    values = []
    for t in angles:
        d = force_line.rotate_about(distal0, centre, axis, float(t))
        direction = force_fn(d) if force_fn else frames.unit(proximal - d)
        if direction is None:
            values.append(np.nan)
            continue
        values.append(float(np.cross(d - centre, direction) @ frames.unit(axis)))
    return np.asarray(values)


def describe(values, angles):
    finite = np.isfinite(values)
    if finite.sum() < 2:
        return {'sign': '?', 'crossings': [], 'at_deg': {}, 'magnitude_mm': 0.0}
    v = values[finite]
    a = angles[finite]
    crossings = []
    for i in range(len(v) - 1):
        if v[i] * v[i + 1] < 0:
            t = abs(v[i]) / (abs(v[i]) + abs(v[i + 1]))
            crossings.append(round(float(a[i] + t * (a[i + 1] - a[i])), 2))
    lo, hi = float(v.min()), float(v.max())
    sign = '+' if lo > 1e-9 else ('-' if hi < -1e-9 else 'flip')
    at = {}
    for s in SAMPLE_DEG:
        idx = int(np.argmin(np.abs(a - s)))
        at['%g' % s] = round(float(v[idx]), 3)
    return {'sign': sign, 'crossings': crossings, 'at_deg': at,
            'min_mm': round(lo, 3), 'max_mm': round(hi, 3),
            'magnitude_mm': round(float(np.abs(v).mean()), 3)}


# ---------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    centre = frames.joint_centres('display')['GH']
    axes = frames.gh_probe_axes()
    jc = json.loads(frames.CANDIDATES_PATH.read_text(encoding='utf-8'))
    head_radius = float(jc['GH']['radius_mm'])

    rows = []
    for mid, (label, origin_bone, role, role_sign) in MUSCLES.items():
        a_src = region_centroid(mid, 'origin', origin_bone)
        b_src = region_centroid(mid, 'insertion', 'humerus')
        if a_src is None or b_src is None:
            print('SKIP', mid, 'missing region file')
            continue
        a = frames.source_to_display(a_src)
        b = frames.source_to_display(b_src)
        clearance = chord_clearance(a, b, centre)
        wraps = clearance < head_radius

        entry = {'id': mid, 'name': label, 'intended_axis': role,
                 'intended_sign': role_sign, 'origin_bone': origin_bone,
                 'origin_display_mm': a.tolist(), 'insertion_display_mm': b.tolist(),
                 'chord_mm': round(float(np.linalg.norm(a - b)), 2),
                 'gh_chord_clearance_mm': round(clearance, 2),
                 'wrapped': bool(wraps), 'axes': {}}

        for axis_name in AXES:
            axis = axes[axis_name]
            straight = moment_curve(axis, centre, a, b, ANGLES)
            rec = {'straight': describe(straight, ANGLES)}
            # A sphere wrap is only meaningful if the insertion itself lies outside
            # the head; muscles inserting on the tubercles sit inside that sphere.
            if wraps and np.linalg.norm(b - centre) > head_radius * 1.05:
                fn = lambda d, _a=a: wrapped_force_direction(_a, d, centre, head_radius)
                wrapped = moment_curve(axis, centre, a, b, ANGLES, force_fn=fn)
                rec['wrapped'] = describe(wrapped, ANGLES)
                rec['sign_agreement'] = (rec['wrapped']['sign'] == rec['straight']['sign'])
            else:
                rec['sign_agreement'] = None
            rec['used'] = 'straight'
            entry['axes'][axis_name] = rec

        primary = entry['axes'][role]['straight']
        sign0 = primary['at_deg'].get('0')
        entry['sign_at_0'] = sign0
        entry['intended_met'] = (sign0 is not None and sign0 != 0.0 and
                                 (1 if sign0 > 0 else -1) == role_sign)
        crossings = primary['crossings']
        entry['usable_range_deg'] = [0.0, crossings[0] if crossings else 120.0]
        rows.append(entry)

    print('GH centre (display) %s    humeral head radius %.2f mm\n'
          % (np.round(centre, 2), head_radius))
    print('%-8s %-10s %-11s %-6s %-7s %-22s %-5s %s' %
          ('muscle', 'name', 'intended', 'path', 'sign0', 'moment@0/60/120 mm', 'meets', 'usable'))
    for r in rows:
        p = r['axes'][r['intended_axis']]['straight']
        at = p['at_deg']
        print('%-8s %-10s %-11s %-6s %-7s %-22s %-5s %s' % (
            r['id'], r['name'],
            '%s%s' % (r['intended_axis'][:4], '+' if r['intended_sign'] > 0 else '-'),
            'wrap' if r['wrapped'] else 'str', p['sign'],
            '%.2f / %.2f / %.2f' % (at.get('0', float('nan')), at.get('60', float('nan')),
                                    at.get('120', float('nan'))),
            'YES' if r['intended_met'] else 'NO',
            '0-%g deg' % r['usable_range_deg'][1]))

    print('\n%-8s %-9s %-9s %-9s %-10s' % ('muscle', 'flexion', 'abduct', 'axial', 'path'))
    for r in rows:
        cells = []
        for ax in AXES:
            rec = r['axes'][ax]
            used = rec[rec['used']]
            cells.append('%s%s' % (used['sign'], '' if not used['crossings'] else '[%g]' % used['crossings'][0]))
        print('%-8s %-9s %-9s %-9s %-10s' % (r['id'], cells[0], cells[1], cells[2],
                                             'wrap' if r['wrapped'] else 'straight'))

    print('\n--- antagonism on the intended axis ---')
    pairs = [('flexion', 'FJ1468', 'FJ1513'), ('abduction', 'FJ1467', 'FJ1507'),
             ('longitudinal', 'FJ1504', 'FJ1500')]
    pair_report = {}
    for axis, ma, mb in pairs:
        def sign_of(mid):
            r = next(x for x in rows if x['id'] == mid)
            v = r['axes'][axis]['straight']['at_deg'].get('0', 0.0)
            return '+' if v > 0 else ('-' if v < 0 else '0')
        sa, sb = sign_of(ma), sign_of(mb)
        ok = sa in ('+', '-') and sb in ('+', '-') and sa != sb
        pair_report[axis] = {'a': ma, 'a_sign_at_0': sa, 'b': mb, 'b_sign_at_0': sb,
                             'opposed_at_neutral': ok}
        print('  %-12s %s(%s)  vs  %s(%s)   %s' % (axis, ma, sa, mb, sb,
                                                    'OPPOSED' if ok else 'NOT OPPOSED'))

    met = sum(1 for r in rows if r['intended_met'])
    print('\nlines meeting their intended sign: %d / %d' % (met, len(rows)))

    report = {'scope': 'moment-arm sign screen only; candidate regions are proximity '
                       'samples, not confirmed attachment sites; scapula and clavicle '
                       'fixed; no muscle force, pennation, friction or scapular motion',
              'gh_centre_display_mm': centre.tolist(), 'head_radius_mm': head_radius,
              'angles_deg': [0.0, 120.0, 1.0], 'muscles': rows,
              'pairs': pair_report, 'intended_met': met, 'total': len(rows)}
    if args.json_out:
        with open(args.json_out, 'w', encoding='utf-8') as fh:
            json.dump(report, fh, ensure_ascii=False, indent=2)
        print('report written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
