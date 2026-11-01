"""S1-G1 with the muscles present: three formulations, all invalid - see 5.4.23.

The bare-bone sweep (5.2.2) rotates the humerus against the scapula and clavicle.  Extending it
to the muscle modules is not a matter of adding them to the fixed set, and this tool exists to
record why the three obvious ways of doing it all fail at 0 degrees - where the design is known
valid - rather than to produce a number to quote.

1. Hold the modules fixed and rotate the humerus.  Every module inserts on the humerus, so the
   rotating bone slides through its own insertion surface and reports a hit near 1 degree.  That
   is an artefact of holding a deforming body still.

2. Test the straight chord from the fixed origin seat to the moved insertion seat against the
   bones.  It fails at 0 deg for all four modules, because the chord is not the design route -
   5.4.8 already showed relaxed corridors differ from chords (FJ1468: 150.4 mm chord vs 167.3 mm
   relaxed).

3. Relax the route at each angle (--relaxed).  Still fails at 0 deg: with a 3.2 mm criterion
   because relaxation is a soft target, and with the 1.6 mm channel-wall criterion because the
   corridor's ends sit ON the bone by construction (the seat is the bone-side attachment), so a
   bone-distance criterion cannot pass at the ends no matter what the joint angle is.

What the gate actually needs is the *deformed* belly at each angle - regenerate the envelope from
the two attachment poses - which is a bigger piece of work and the correct next step.  Until then
S1-G1's muscle half is open, and none of the numbers below should be read as a range.

Run:
    python mujoco_arm/tools/sweep_gh_muscle_paths.py                  # formulation 2
    python mujoco_arm/tools/sweep_gh_muscle_paths.py --relaxed --step 30   # formulation 3
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import attachment_regions
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP

BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'
MUSCLE_MIN_CLEAR_MM = 3.2      # relaxation target: tube radius
#: the pass criterion is the wall the printed channel needs (5.4.14), not the relaxation
#: target - requiring 3.2 mm of material everywhere fails at 0 deg even though the design
#: routes are valid at 0 deg
CRITERION_MM = 1.6

MODULES = [('teres_major', ('FJ1507',)),
           ('subscapularis', ('FJ1504',)),
           ('infraspinatus', ('FJ1500',)),
           ('deltoid', ('FJ1468', 'FJ1513', 'FJ1467'))]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--step', type=float, default=10.0)
    ap.add_argument('--max', dest='max_deg', type=float, default=120.0)
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--module', default=None, help='limit to one module key')
    ap.add_argument('--criterion-mm', type=float, default=CRITERION_MM,
                    help='material the route must keep from bone [mm]; the channel wall')
    ap.add_argument('--relaxed', action='store_true',
                    help='relax the route at each angle (the design route, ~7.5 s/angle) '
                         'instead of testing the straight chord')
    args = ap.parse_args(argv)

    gh = frames.joint_centres('source')['GH']
    axes = frames.gh_probe_axes()
    bone_tris = []
    for v in SP.BONES.values():
        BV, BF = mesh_geom.load_obj(str(BONE_DIR / v))
        BW, BG = mesh_geom.weld(BV, BF)
        bone_tris.append(BW[BG])

    def seat_of(heads, kind):
        pts, ws = [], []
        for mid in heads:
            c, a = attachment_regions.centroid(attachment_regions.load_region(mid, kind))
            pts.append(c)
            ws.append(a)
        pts = np.array(pts)
        ws = np.array(ws, float)
        return (pts * ws[:, None]).sum(0) / ws.sum()

    modules = [m for m in MODULES if args.module in (None, m[0])]
    all_tris = np.vstack(bone_tris)
    print('GH centre (source mm) %s   step %.0f deg, up to %.0f deg'
          % (np.round(gh, 2), args.step, args.max_deg))
    print('\n%-15s %-8s %-11s %9s %10s'
          % ('module', 'head', 'axis', 'first_fail', 'last_ok'))

    rows = []
    for key, heads in modules:
        origins = {mid: attachment_regions.centroid(
            attachment_regions.load_region(mid, 'origin'))[0] for mid in heads}
        for mid in heads:
            o = origins[mid]
            d0 = seat_of(heads, 'insertion')
            chord0 = float(np.linalg.norm(d0 - o))
            for name in ('flexion', 'abduction', 'longitudinal'):
                axis = axes[name]
                for sign in (+1, -1):
                    last_ok, first_fail = None, None
                    profile = []
                    for a in np.arange(0.0, args.max_deg + 1e-9, args.step):
                        ang = sign * float(a)
                        T = I.rotation_about(gh, axis, ang)
                        moved = T[:, :3] @ d0 + T[:, 3]
                        ts = np.linspace(0.0, 1.0, SP.INTERIOR_POINTS)
                        band = [o + (moved - o) * t for t in ts if 0.2 <= t <= 0.8]
                        if args.relaxed:
                            # The chord is not the design route: it cuts the bone at 0 deg
                            # even though a valid relaxed corridor exists (5.4.8), so a
                            # chord-only sweep reports failure everywhere and means nothing.
                            # Relaxing per angle costs ~7.5 s, hence the coarse step.
                            path, _ = SP.relax([o + (moved - o) * t for t in ts], all_tris,
                                               iterations=SP.ITERATIONS, step=SP.STEP_MM,
                                               lam=SP.SMOOTHING, bone_tris=bone_tris,
                                               bone_clear_mm=MUSCLE_MIN_CLEAR_MM)
                            path[-1] = moved
                            band = [path[k] for k, t in enumerate(ts) if 0.2 <= t <= 0.8]
                        worst, _, _ = SP.clearance(band, all_tris)
                        profile.append({'deg': round(ang, 1), 'chord_worst_mm': round(worst, 3)})
                        if worst <= -args.criterion_mm:
                            last_ok = ang
                        elif first_fail is None:
                            first_fail = ang
                    rows.append({'module': key, 'head': mid, 'axis': name,
                                 'sign': int(sign), 'chord0_mm': round(chord0, 2),
                                 'last_ok_deg': last_ok, 'first_fail_deg': first_fail,
                                 'profile': profile})
                    print('%-15s %-8s %-4s %5s %9s %10s'
                          % (key, mid, name[:4], '+' if sign > 0 else '-',
                             'none' if first_fail is None else '%.0f' % first_fail,
                             'none' if last_ok is None else '%.0f' % last_ok))

    clean = [r for r in rows if r['first_fail_deg'] is None]
    print('\n%-15s %s' % ('summary', 'directions fully clear to %.0f deg' % args.max_deg))
    for key, _ in modules:
        sub = [r for r in rows if r['module'] == key]
        ok = [r for r in sub if r['first_fail_deg'] is None]
        worst = min((r['first_fail_deg'] for r in sub
                     if r['first_fail_deg'] is not None), default=None)
        print('  %-15s %d/%d clear%s'
              % (key, len(ok), len(sub),
                 '' if worst is None else ', first fail at %.0f deg' % worst))

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'routability of the wire corridor at each GH angle: the insertion seat is '
                     'rotated with the humerus and the corridor rebuilt under the design bone '
                     'clearance; a rigid-body module sweep is NOT used (each module inserts on '
                     'the humerus, so that would report an artefact at 1 deg)',
            'step_deg': args.step, 'max_deg': args.max_deg,
            'relax_target_mm': MUSCLE_MIN_CLEAR_MM, 'criterion_mm': args.criterion_mm,
            'directions': rows,
            'n_directions': len(rows),
            'n_fully_clear': len(clean),
            'not_verified': ['a real belly deforms; a clear route is not a working muscle',
                             'the module envelope is not swept, only the corridor',
                             'no joint parts, no fixture pads, no contact'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
