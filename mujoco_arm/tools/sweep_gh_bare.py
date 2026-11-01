"""Bare-bone GH sweep for the open-socket L0 (D15).

D15 says the L1 joint is bone in bone: the humeral head rests in the bony glenoid
(a 2 mm liner stands in for cartilage), with no ball, no socket and no
anti-dislocation cap on the scapula side.  So the kinematic bound of that joint is
simply *bone against bone*, and that is what this script measures.

It rotates the welded humerus about the GH centre on each of the three probe axes
(`frames.gh_probe_axes`) in both directions and reports the first angle at which
the humerus penetrates the fixed scapula or clavicle above the shared threshold.

What this does NOT do:
* no liner geometry - the cartilage stand-in is only discussed, not modelled;
* no muscle envelope, no cable, no caps (they do not exist in D15);
* the GH centre is held fixed, so it says nothing about the +-2 mm centre drift
  the open socket allows (R13);
* a sampled sweep plus a boundary bisection, not a proof over all poses.

Run:
    python mujoco_arm/tools/sweep_gh_bare.py
    python mujoco_arm/tools/sweep_gh_bare.py --step 2 --max 150 --json out.json
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frames
import interference as I
import mesh_geom

BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'
HUMERUS = BONE_DIR / '肱骨.obj'
SCAPULA = BONE_DIR / '肩胛骨.obj'
CLAVICLE = BONE_DIR / '锁骨.obj'

AXES = ('flexion', 'abduction', 'longitudinal')


def bone_solid(path):
    """Welded manifold solid, with the edge census printed (R14)."""
    V, F = mesh_geom.load_obj(str(path))
    e0, b0, n0 = mesh_geom.edge_census(F)
    W, G = mesh_geom.weld(V, F)
    e1, b1, n1 = mesh_geom.edge_census(G)
    mf = I._manifold()
    solid = mf.Manifold(mf.Mesh64(W, G.astype(np.uint32)))
    print('%-9s v=%5d f=%5d  edges/boundary/non-manifold %d/%d/%d -> %d/%d/%d  '
          'volume %.1f mm3'
          % (path.stem, len(V), len(F), e0, b0, n0, e1, b1, n1, solid.volume()))
    if b1 or n1:
        raise SystemExit('welded %s is still not a closed manifold; do not trust '
                         'the sweep' % path.name)
    return solid


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--step', type=float, default=1.0)
    ap.add_argument('--max', dest='max_deg', type=float, default=150.0)
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    gh = frames.joint_centres('source')['GH']
    axes = frames.gh_probe_axes()
    print('GH centre (source mm) %s' % np.round(gh, 2))

    moving = bone_solid(HUMERUS)
    fixed = [I.Solid('scapula', bone_solid(SCAPULA)),
             I.Solid('clavicle', bone_solid(CLAVICLE))]

    print('\nangles are rotation-vector degrees about the named GH probe axis'
          '\n(+ / - are both swept; step %.1f deg, bisection %.1f deg, threshold %.3f mm3)'
          % (args.step, I.BOUNDARY_DEG, I.HARD_MM3))
    print('\n%-11s %7s %8s %9s %10s %s'
          % ('axis', 'sign', 'first', 'volume', 'worst_deg', 'first limit'))

    rows = []
    for name in AXES:
        axis = axes[name]
        for sign in (+1, -1):
            angles = [sign * float(a) for a in np.arange(0.0, args.max_deg + 1e-9, args.step)]
            t0 = time.time()
            rep = I.sweep(axis, gh, moving, fixed, angles, moving_key='humerus',
                          refine=True)
            hit = rep['first_hit_deg']
            limit = ''
            if hit is not None:
                ang = np.radians(hit)
                T = I.rotation_about(gh, axis, hit)
                moved = moving.transform(T)
                best = max(((I.intersection_volume(moved, s.geom), s.key) for s in fixed),
                           key=lambda kv: kv[0])
                limit = '%s (%.3f mm3)' % (best[1], best[0])
            worst = rep['worst']
            print('%-11s %7s %8s %9s %10s %s'
                  % (name, '+' if sign > 0 else '-',
                     '%.2f' % hit if hit is not None else 'none',
                     '%.3f' % worst['volume_mm3'] if worst else '-',
                     '%.0f' % worst['angle_deg'] if worst else '-',
                     limit if hit is not None else 'clear to %.0f deg'
                     % abs(angles[-1])))
            rows.append({'axis': name, 'sign': sign, 'axis_vector_source': axis.tolist(),
                         'sampled_deg': [angles[0], angles[-1]],
                         'step_deg': args.step,
                         'first_hit_deg': None if hit is None else round(hit, 2),
                         'worst': worst,
                         'first_limit': limit or None,
                         'elapsed_s': round(time.time() - t0, 1)})

    report = {
        'scope': 'bare bone vs bare bone at the GH centre; no liner, no ball, no socket, '
                 'no cap, no muscle envelope; GH centre held fixed',
        'l0': 'open shallow socket (D15)',
        'gh_centre_source_mm': gh.tolist(),
        'threshold_mm3': I.HARD_MM3,
        'boundary_deg': I.BOUNDARY_DEG,
        'bones': {'moving': '肱骨 FJ3368', 'fixed': ['肩胛骨 FJ3384', '锁骨 FJ3362']},
        'sweeps': rows,
        'not_verified': ['cartilage/liner thickness and contact stress',
                         'GH centre drift (R13)',
                         'muscle envelope, cable path, drive',
                         'all poses (sampled + bisection only)'],
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
