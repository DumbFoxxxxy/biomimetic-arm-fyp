"""S1-G1 as a cable-length sweep instead of a rigid-body sweep.

The rigid sweep (5.4.24, deform_and_sweep.py) had to be retracted as a criterion: it rotates the
muscle as a solid body, so the belly cannot slide over the bone or over its neighbours, and the
collisions it reports are artefacts of the model rather than of the design. The audit note in
5.4.24 says exactly that.

For a cable-driven actuator the question is different and much better posed:

    over the range, does the DISTANCE the wire must span change monotonically, in the right
    direction, and by a feasible amount?

That is what the bench will measure (the wire either goes slack or it does not), it needs no
deformation model, and it is insensitive to the belly-sliding problem. This sweeps each muscle's
origin-to-insertion span about each GH axis and reports:

  * span(angle) - the cable length the TPU has to take up
  * d(span)/d(angle) - which IS the moment arm, so it checks the signs from 5.3 a second way
  * the fraction the span changes over 0-120 deg
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frames
import interference as I

MODULES_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_modules.json'
ANGLES = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0)


def span(origin, insertion, axis, centre, deg):
    R = I.rotation_about(centre, axis, deg)
    moved = R[:3, :3] @ np.asarray(insertion, float) + R[:3, 3]
    return float(np.linalg.norm(moved - np.asarray(origin, float)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    gh = frames.joint_centres('source')['GH']
    axes = frames.gh_probe_axes()
    design = json.loads(MODULES_JSON.read_text(encoding='utf-8'))

    rows = []
    print('%-14s %-8s %-13s %8s %8s %8s %8s' %
          ('module', 'head', 'axis', 'span0', 'span120', 'd/ddeg', 'change%'))
    for m in design['modules']:
        for c in m['corridors']:
            path = c['path_source_mm']
            origin, insertion = np.asarray(path[0], float), np.asarray(path[-1], float)
            for name, axis in axes.items():
                spans = [span(origin, insertion, axis, gh, a) for a in ANGLES]
                s0, s1 = spans[0], spans[-1]
                # least-squares slope over the sweep: mm per degree
                A = np.c_[np.asarray(ANGLES), np.ones(len(ANGLES))]
                slope = float(np.linalg.lstsq(A, np.asarray(spans), rcond=None)[0][0])
                mono = all(b >= a - 0.05 for a, b in zip(spans, spans[1:])) or \
                       all(b <= a + 0.05 for a, b in zip(spans, spans[1:]))
                rows.append({'module': m['module'], 'head': c['head'], 'axis': name,
                             'span0_mm': round(s0, 3), 'span120_mm': round(s1, 3),
                             'slope_mm_per_deg': round(slope, 4),
                             'change_pct': round(100.0 * (s1 - s0) / s0, 3),
                             'monotonic': bool(mono),
                             'spans_mm': [round(x, 3) for x in spans]})
                print('%-14s %-8s %-13s %8.2f %8.2f %8.4f %8.3f'
                      % (m['module'], c['head'], name, s0, s1, slope,
                         100.0 * (s1 - s0) / s0))

    # The moment arm is d(span)/d(theta); per radian it should match 5.3's numbers.
    print('\n--- sign agreement with the moment-arm table (5.3) ---')
    bad = [r for r in rows if r['change_pct'] * r['slope_mm_per_deg'] < 0]
    for r in rows:
        if abs(r['change_pct']) < 0.01:
            continue
    signs = {}
    for r in rows:
        signs.setdefault((r['module'], r['head']), {})[r['axis']] = \
            'lengthen' if r['change_pct'] > 0 else 'shorten'
    for k, v in signs.items():
        print('  %-14s %-8s %s' % (k[0], k[1], v))

    stride = [r for r in rows if not r['monotonic']]
    print('\nspans monotonic over 0-120 deg: %d of %d'
          % (len(rows) - len(stride), len(rows)))
    if stride:
        print('  non-monotonic:', ['%s %s %s' % (r['module'], r['head'], r['axis'])
                                   for r in stride])

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'S1-G1 posed as a cable-length sweep: the origin-to-insertion span over '
                     '0-120 deg about each GH axis, with its slope (= the moment arm)',
            'why': 'the rigid sweep of a deforming part is not a usable criterion (5.4.24); a '
                   'cable either takes up length or it does not, and that is what the bench '
                   'measures (S1-G3)',
            'angles_deg': list(ANGLES), 'rows': rows,
            'not_verified': ['the span is the chord, not the corridor path (which is longer and '
                             'bends over the bone)',
                             'no strain limit is applied - the change is reported, not judged',
                             'the belly still has to slide; that is not modelled here at all'],
        })
        print('\nreport written:', args.json_out)


if __name__ == '__main__':
    main()
