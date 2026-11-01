"""S1 channel routes: which of the six lines can keep a channel inside its own belly?

Round 9 showed that all three routing heuristics fail for the deltoid, because the
deltoid is a shell rather than a tube.  The obvious next question is whether that is a
deltoid problem or a shoulder-wide problem - if five of six lines are fine, only the
deltoid module has to be reconstructed; if all six fail, the whole "copy the anatomical
outer shape" premise needs restating.

So this runs the same three routes for every line, against that muscle's own belly
after trimming it by the bones:

    chord3          origin pad -> mid-section -> insertion pad
    centroid-curve  section area centroids
    medial-axis     centres of the largest inscribed circle of each section

Verdict per route: does it stay inside the belly with at least the channel wall (R1.6)?

Run:
    python mujoco_arm/tools/check_shoulder_channel_routes.py
    python mujoco_arm/tools/check_shoulder_channel_routes.py --json out.json
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
STATIONS = np.arange(0.12, 0.93, 0.10)


def mf():
    return I._manifold()


def load_solid(path):
    V, F = mesh_geom.load_obj(str(path))
    W, G = mesh_geom.weld(V, F)
    return mf().Manifold(mf().Mesh64(W, G.astype(np.uint32))), W, G


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    bones = {n: load_solid(BONE_DIR / f)[0] for n, f in BONES.items()}
    rows = []
    print('channel needs %.1f mm of wall (elbow R1.6)\n' % CHANNEL_OUTER_R_MM)
    print('%-8s %-10s %-16s %7s %7s %9s %9s %s'
          % ('muscle', 'name', 'route', 'points', 'length', 'worst', 'deepest', 'verdict'))
    for mid, label in MUSCLES.items():
        origin = attachment_regions.load_region(mid, 'origin')
        insertion = attachment_regions.load_region(mid, 'insertion')
        o, _ = attachment_regions.centroid(origin)
        i, _ = attachment_regions.centroid(insertion)
        solid, _, _ = load_solid(MUSCLE_DIR / (mid + '.obj'))
        trimmed = solid
        for n in BONES:
            trimmed = trimmed - bones[n]
        comps = len([c for c in trimmed.decompose() if c.volume() > 0.5])
        verts, tris = mesh_geom.manifold_mesh_to_arrays(trimmed.to_mesh())
        W, G = mesh_geom.weld(verts, tris)
        tri = W[G]

        for route in ('chord3', 'centroid-curve', 'medial-axis'):
            if route == 'chord3':
                chunk = muscle_sections.section_at(solid, o, i, 0.5)
                path = [o]
                if chunk and chunk.get('centroid_source_mm'):
                    path.append(np.asarray(chunk['centroid_source_mm'], float))
                path.append(i)
            else:
                key = ('centroid_source_mm' if route == 'centroid-curve'
                       else 'inradius_source_mm')
                path = [o]
                for t in STATIONS:
                    chunk = muscle_sections.section_at(solid, o, i, float(t))
                    if chunk and chunk.get(key):
                        path.append(np.asarray(chunk[key], float))
                path.append(i)
            wall = muscle_sections.path_clearance(path, tri)
            ok = wall['fits_channel']
            length = float(sum(np.linalg.norm(np.asarray(b, float) - np.asarray(a, float))
                               for a, b in zip(path, path[1:])))
            rows.append({'id': mid, 'name': label, 'route': route,
                         'trimmed_components': comps,
                         'points': len(path), 'length_mm': round(length, 2),
                         'worst_mm': wall['worst_mm'],
                         'deepest_mm': wall['deepest_mm'],
                         'inside_fraction': wall['inside_fraction'],
                         'channel_outer_r_mm': CHANNEL_OUTER_R_MM,
                         'fits': bool(ok)})
            print('%-8s %-10s %-16s %7d %7.1f %9.2f %9.2f %s'
                  % (mid, label, route, len(path), length, wall['worst_mm'],
                     wall['deepest_mm'], 'OK' if ok else 'OUT/TOO THIN'))

    best = {}
    for r in rows:
        cur = best.get(r['id'])
        if cur is None or r['worst_mm'] < cur['worst_mm']:
            best[r['id']] = r
    print('\n--- best route per line (smallest worst-point gap) ---')
    for mid, label in MUSCLES.items():
        b = best[mid]
        print('  %-8s %-10s best %-16s worst %6.2f mm  %s'
              % (mid, label, b['route'], b['worst_mm'], 'FITS' if b['fits'] else 'FAILS'))
    fits = [m for m in MUSCLES if best[m]['fits']]
    print('\nlines whose best route can carry the channel: %d / %d (%s)'
          % (len(fits), len(MUSCLES), ', '.join(fits) or 'none'))
    print('lines needing a reconstructed belly: %s'
          % ', '.join(m for m in MUSCLES if m not in fits))

    report = {
        'scope': 'three path heuristics per line, tested against that muscle\'s own '
                 'bone-trimmed belly; no lattice, no end caps, no channel geometry',
        'channel_outer_r_mm': CHANNEL_OUTER_R_MM,
        'stations': [round(float(t), 2) for t in STATIONS],
        'routes': rows,
        'best_per_line': {m: {k: best[m][k] for k in ('route', 'worst_mm', 'fits')}
                          for m in MUSCLES},
        'lines_that_fit': fits,
        'lines_needing_reconstruction': [m for m in MUSCLES if m not in fits],
        'not_verified': ['these are three heuristics, not a path search; a line that '
                         'fails here might still have a reachable route',
                         'deltoid paths are tested per head, not through a built Y',
                         'no lattice, caps, crimp or bolt geometry',
                         'no muscle force, pennation or contact'],
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
