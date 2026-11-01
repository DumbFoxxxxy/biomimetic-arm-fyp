"""S1 channel path search: how much material each line is short of.

Rounds 9-10 established that three heuristics (chord, section-centroid curve,
inscribed-circle curve) all leave the belly for every one of the six lines, and that
the excursion depth differs a lot.  What is still unknown is whether a *different*
route would fit, and if not, how much material is missing.

This relaxes a path inside the bone-trimmed belly instead of guessing it: every interior
point walks along -grad(signed distance), i.e. deeper into the material, with a Laplacian
smoothing term so it does not zig-zag.  The endpoints stay pinned to the two bone-side
seats, because that is where the wire has to arrive.

Two conventions carried over, both easy to get wrong:

* ``signed_distance`` is NEGATIVE inside, so a point gets better as its signed distance
  gets more negative, and the binding number along a path is the MAXIMUM (the least
  deep point);
* the last few millimetres at each end sit on the bone surface, where the end pad and
  cap carry the wire, so the wall requirement applies to the belly only.  The search
  therefore reports the worst clearance over the interior (t in 0.15-0.85) - the same
  band §5.4.3 used.

Result per line: the best clearance reachable by local relaxation, and
``bulge_needed_mm = max(0, -1.6 - worst)``, i.e. how much the envelope has to grow at
its tightest point before the channel fits.

Run:
    python mujoco_arm/tools/search_shoulder_channel_paths.py
    python mujoco_arm/tools/search_shoulder_channel_paths.py --json out.json
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
import muscle_sections

MUSCLE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / 'muscle'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

MUSCLES = {'FJ1468': '三角肌前束', 'FJ1513': '三角肌后束', 'FJ1467': '三角肌中束',
           'FJ1507': '大圆肌', 'FJ1504': '肩胛下肌', 'FJ1500': '冈下肌'}
BONES = {'scapula': '肩胛骨.obj', 'clavicle': '锁骨.obj', 'humerus': '肱骨.obj'}

CHANNEL_OUTER_R_MM = 1.6
BAND = (0.15, 0.85)
ITERATIONS = 80
STEP_MM = 1.2
SMOOTHING = 0.25
INTERIOR_POINTS = 12


def mf():
    return I._manifold()


def load_solid(path):
    V, F = mesh_geom.load_obj(str(path))
    W, G = mesh_geom.weld(V, F)
    return mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))


def clearance(points, tri):
    """(worst, deepest) signed clearance over a set of points; negative = inside."""
    d = [mesh_geom.signed_distance(p, tri) for p in points]
    return float(max(d)), float(min(d)), d


def relax(path, tri, iterations=ITERATIONS, step=STEP_MM, lam=SMOOTHING):
    """Push interior points into the material, keep the ends pinned."""
    p = np.array(path, dtype=float)
    n = len(p)
    history = []
    for _ in range(iterations):
        moved = p.copy()
        for i in range(1, n - 1):
            d, nrm = mesh_geom.signed_distance_and_normal(p[i], tri)
            # -normal walks deeper into the material whichever side we start on
            moved[i] = p[i] - nrm * step
        # Laplacian smoothing on the interior, ends fixed
        for i in range(1, n - 1):
            moved[i] = moved[i] + lam * (moved[i - 1] + moved[i + 1] - 2.0 * moved[i])
        p = moved
        worst = max(mesh_geom.signed_distance(q, tri) for q in p[1:-1])
        history.append(round(worst, 3))
    return p, history


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    bones = {k: load_solid(BONE_DIR / v) for k, v in BONES.items()}
    rows = []
    print('channel needs %.1f mm of wall; wall judged over the belly band t %.2f-%.2f\n'
          % (CHANNEL_OUTER_R_MM, BAND[0], BAND[1]))
    print('%-8s %-10s %8s %8s %9s %9s %8s %9s %s'
          % ('muscle', 'name', 'start w', 'final w', 'improve', 'bulge', 'surplus',
             'length', 'verdict'))
    for mid, label in MUSCLES.items():
        origin = attachment_regions.load_region(mid, 'origin')
        insertion = attachment_regions.load_region(mid, 'insertion')
        o, _ = attachment_regions.centroid(origin)
        i, _ = attachment_regions.centroid(insertion)
        solid = load_solid(MUSCLE_DIR / (mid + '.obj'))
        trimmed = solid
        for b in bones.values():
            trimmed = trimmed - b
        verts, tris = mesh_geom.manifold_mesh_to_arrays(trimmed.to_mesh())
        W, G = mesh_geom.weld(verts, tris)
        tri = W[G]

        # initial path: evenly spaced points on the chord, so the relaxation starts from
        # the simplest possible guess and cannot inherit a heuristic's bias
        ts = np.linspace(0.0, 1.0, INTERIOR_POINTS)
        path = [o + (i - o) * t for t in ts]
        band = [path[k] for k, t in enumerate(ts) if BAND[0] <= t <= BAND[1]]
        start_worst, _, _ = clearance(band, tri)
        start_len = float(np.linalg.norm(i - o))

        relaxed, history = relax(path, tri)
        band = [relaxed[k] for k, t in enumerate(ts) if BAND[0] <= t <= BAND[1]]
        final_worst, final_deep, _ = clearance(band, tri)
        final_len = float(sum(np.linalg.norm(np.diff(relaxed, axis=0), axis=1)))
        # requirement is worst <= -wall, so a positive `bulge` is a deficit and a
        # negative one is spare material
        bulge = max(0.0, final_worst + CHANNEL_OUTER_R_MM)
        surplus = max(0.0, -CHANNEL_OUTER_R_MM - final_worst)
        fits = final_worst <= -CHANNEL_OUTER_R_MM
        rows.append({'id': mid, 'name': label,
                     'band': list(BAND),
                     'required_wall_mm': CHANNEL_OUTER_R_MM,
                     'start_worst_mm': round(start_worst, 3),
                     'final_worst_mm': round(final_worst, 3),
                     'final_deepest_mm': round(final_deep, 3),
                     'improvement_mm': round(start_worst - final_worst, 3),
                     'bulge_needed_mm': round(bulge, 3),
                     'surplus_mm': round(surplus, 3),
                     'chord_length_mm': round(start_len, 2),
                     'relaxed_length_mm': round(final_len, 2),
                     'length_growth_mm': round(final_len - start_len, 2),
                     'fits': bool(fits),
                     'iterations': ITERATIONS,
                     'history': history,
                     'path_source_mm': [np.round(p, 2).tolist() for p in relaxed]})
        print('%-8s %-10s %8.2f %8.2f %9.2f %9.2f %8.2f %9.1f %s'
              % (mid, label, start_worst, final_worst, start_worst - final_worst,
                 bulge, surplus, final_len,
                 'FITS' if fits else 'needs %.2f mm' % bulge))

    fits = [r['id'] for r in rows if r['fits']]
    print('\nlines whose relaxed path fits the channel: %d / %d (%s)'
          % (len(fits), len(rows), ', '.join(fits) or 'none'))
    worst = sorted(rows, key=lambda r: -r['bulge_needed_mm'])
    print('largest bulge needed: %s'
          % ', '.join('%s %.2f mm' % (r['id'], r['bulge_needed_mm'])
                      for r in worst[:3] if r['bulge_needed_mm'] > 0) or 'none')
    print('tightest surplus: %s'
          % ', '.join('%s %.2f mm' % (r['id'], r['surplus_mm'])
                      for r in sorted(rows, key=lambda r: r['surplus_mm'])[:3]))

    report = {
        'scope': 'gradient relaxation of a 12-point path inside each bone-trimmed belly; '
                 'the wall is judged over t 0.15-0.85; ends pinned to the two bone-side '
                 'seats; no lattice, no channel solid, no caps',
        'required_wall_mm': CHANNEL_OUTER_R_MM,
        'relaxation': {'iterations': ITERATIONS, 'step_mm': STEP_MM,
                       'smoothing': SMOOTHING, 'interior_points': INTERIOR_POINTS},
        'lines': rows,
        'lines_that_fit': fits,
        'not_verified': ['local relaxation, not a global search: a better route may exist',
                         'the ends are pinned on the bone surface, so the wall there is '
                         'not judged - the end pad and cap carry the wire',
                         'bulge_needed is measured against the trimmed belly; the '
                         'reconstructed belly still has to clear bone, neighbour bellies '
                         'and the fixture pads',
                         'no lattice, caps, crimp or bolt geometry',
                         'no muscle force, pennation or contact'],
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
