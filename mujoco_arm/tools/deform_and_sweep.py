"""S1-G1 with muscles, done the way 5.4.23 said it must be: deform the belly, then sweep.

The three shortcuts in 5.4.23 all failed at 0 degrees because none of them moved the muscle.
This deforms it instead, with the simplest kinematic model that is consistent at both ends:

    every vertex is parameterised by its position t along the design corridor (0 = origin
    seat, 1 = insertion seat) and rotated about the GH axis by angle * s(t)

with s(t) a smoothstep, so t = 0 is untouched and stays on the scapula, t = 1 receives the full
rotation and stays glued to the humerus, and the belly interpolates between them.  That is a
first-order model: no volume preservation, no fibre reorientation, no tension.  It is enough to
ask the question the gate asks - does the belly collide with bone as the joint moves - which the
zero-pose rigid sweep could not ask at all.

Run:
    python mujoco_arm/tools/deform_and_sweep.py
    python mujoco_arm/tools/deform_and_sweep.py --step 15 --json out.json
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
import search_shoulder_channel_paths as SP

DESIGN_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_modules.json'
MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

HUMERUS_KEY = 'humerus'
FIXED_KEYS = ('scapula', 'clavicle')
#: The module STLs are re-read and re-booleaned here, and the roundtrip alone leaves slivers
#: of 0.0006-0.0219 mm3 (5.4.18).  Judging against HARD_MM3 = 0.001 reports a "hit" at 0 deg,
#: where the design is known valid - the same trap as before, so the noise floor is the gate.
THRESHOLD_MM3 = 0.05

def mf():
    return I._manifold()


def polyline_param(verts, poly):
    """t in [0, 1] along `poly` for each vertex (nearest point on the polyline)."""
    P = np.asarray(poly, float)
    seg = P[1:] - P[:-1]
    seg_len = np.linalg.norm(seg, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_len)])
    total = cum[-1] if cum[-1] > 1e-9 else 1.0
    best_d = np.full(len(verts), np.inf)
    best_t = np.zeros(len(verts))
    for k in range(len(seg)):
        a, b = P[k], P[k + 1]
        d = b - a
        L2 = float(d @ d)
        if L2 < 1e-12:
            continue
        w = np.clip(((verts - a) @ d) / L2, 0.0, 1.0)
        closest = a + w[:, None] * d
        dist = np.linalg.norm(verts - closest, axis=1)
        better = dist < best_d
        best_d[better] = dist[better]
        best_t[better] = (cum[k] + w[better] * seg_len[k]) / total
    return best_t, best_d


def smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def weight(t):
    """How much of the joint motion a point at parameter t follows.

    Zero over the proximal half, ramping to one at the insertion.  Two earlier shapes failed:

    * s = smoothstep(t): the whole belly follows the joint progressively, so the proximal belly
      swings about the GH centre and is driven into the scapula (teres major: 95 to 3854 mm3 of
      scapula overlap at 30 deg).  A muscle shortens between its attachments; it does not orbit
      the joint.
    * a translation blend, p + s(t) * delta: the distal belly is *translated* while the bone
      *rotates*, so the bone swings away from it and cuts in (217 mm3 of humerus at 30 deg).

    Following the rotation rather than the translation, and holding the proximal half still,
    avoids both.
    """
    t = np.clip(t, 0.0, 1.0)
    return smoothstep(np.clip((t - 0.5) / 0.5, 0.0, 1.0))


def deform(verts, paths, gh, axis, angle_deg):
    """Rotate each vertex about the GH axis by angle * weight(t)."""
    t_all = np.zeros(len(verts))
    d_all = np.full(len(verts), np.inf)
    for poly in paths:
        t, d = polyline_param(verts, poly)
        better = d < d_all
        t_all[better] = t[better]
        d_all[better] = d[better]
    w = weight(t_all)
    out = verts.copy()
    for k in range(len(verts)):
        if w[k] == 0.0:
            continue
        T = I.rotation_about(gh, axis, angle_deg * w[k])
        out[k] = T[:, :3] @ verts[k] + T[:, 3]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--step', type=float, default=15.0)
    ap.add_argument('--max', dest='max_deg', type=float, default=120.0)
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--module', default=None)
    args = ap.parse_args(argv)

    gh = frames.joint_centres('source')['GH']
    axes = frames.gh_probe_axes()
    design = json.loads(DESIGN_JSON.read_text(encoding='utf-8'))

    humerus = SP.load_solid(BONE_DIR / SP.BONES[HUMERUS_KEY])
    fixed = {k: SP.load_solid(BONE_DIR / SP.BONES[k]) for k in FIXED_KEYS}

    print('GH centre %s   step %.0f deg, up to %.0f deg'
          % (np.round(gh, 2), args.step, args.max_deg))
    print('\n%-15s %-8s %-6s %7s %10s %9s %9s'
          % ('module', 'head', 'axis', 'angle', 'humerus', 'scapula', 'clavicle'))

    rows = []
    for m in design['modules']:
        key = m['module']
        if args.module and args.module != key:
            continue
        stl = MODULES_DIR / ('%s.stl' % key)
        V, F = mesh_geom.read_stl(str(stl))
        W, G = mesh_geom.weld(V, F)
        paths = [c['path_source_mm'] for c in m['corridors']]
        seat0 = np.asarray(m['seat_source_mm'], float)
        for name in ('flexion', 'abduction', 'longitudinal'):
            axis = axes[name]
            for sign in (+1, -1):
                first_hit = None
                first_hit_vol = None
                worst_vol = 0.0
                for a in np.arange(0.0, args.max_deg + 1e-9, args.step):
                    ang = sign * float(a)
                    T = I.rotation_about(gh, axis, ang)
                    moved = deform(W, paths, gh, axis, ang)
                    geom = mf().Manifold(mf().Mesh64(moved, G.astype(np.uint32)))
                    hum = humerus.transform(T)
                    hit_h = I.intersection_volume(geom, hum)
                    hit_s = I.intersection_volume(geom, fixed['scapula'])
                    hit_c = I.intersection_volume(geom, fixed['clavicle'])
                    worst_vol = max(worst_vol, hit_h, hit_s, hit_c)
                    if max(hit_h, hit_s, hit_c) > THRESHOLD_MM3 and first_hit is None:
                        first_hit = ang
                        first_hit_vol = max(hit_h, hit_s, hit_c)
                        print('%-15s %-8s %-6s %7.0f %10.3f %9.3f %9.3f'
                              % (key, '-', name[:6], ang, hit_h, hit_s, hit_c))
                rows.append({'module': key, 'axis': name, 'sign': int(sign),
                             'first_hit_deg': first_hit,
                             'first_hit_mm3': (None if first_hit_vol is None
                                               else round(first_hit_vol, 3)),
                             'worst_mm3': round(worst_vol, 3)})
                print('%-15s %-8s %-6s %7s  worst %.3f mm3'
                      % (key, '-', name[:6],
                         'clean' if first_hit is None else 'hit %.0f' % first_hit,
                         worst_vol))

    clean = [r for r in rows if r['first_hit_deg'] is None]
    print('\n%-15s %s' % ('summary', 'directions clean to %.0f deg' % args.max_deg))
    for key in sorted({r['module'] for r in rows}):
        sub = [r for r in rows if r['module'] == key]
        print('  %-15s %d/%d clean, first hit %s'
              % (key, len([r for r in sub if r['first_hit_deg'] is None]), len(sub),
                 min((r['first_hit_deg'] for r in sub
                      if r['first_hit_deg'] is not None), default='none')))

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'GH sweep with the muscle bellies DEFORMED, not held rigid: each vertex '
                     'rotated by angle * smoothstep(t) about the GH axis, t measured along the '
                     'design corridor, so the origin stays on the scapula and the insertion '
                     'stays glued to the humerus',
            'model_caveat': 'first-order kinematics: no volume preservation, no fibre '
                            'reorientation, no tension',
            'step_deg': args.step, 'max_deg': args.max_deg,
            'directions': rows,
            'n_directions': len(rows), 'n_clean': len(clean),
            'not_verified': ['no joint parts, no fixture pads, no contact',
                             'a TPU belly also buckles and twists; this only moves it'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
