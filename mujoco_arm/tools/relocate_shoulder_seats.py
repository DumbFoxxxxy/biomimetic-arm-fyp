"""Relocate the two seats where no pad or cap footprint fits, onto thicker bone nearby.

5.9.14 showed that two seats cannot be fitted at any (thickness x radius): the infraspinatus origin
and the deltoid posterior origin both sit on a 2.8 mm scapular blade with the opposite-side muscle
pressed against the far face, so an R11-class clamp has nowhere to go there.

Of the three ways out (move the seat, clamp from one side, relieve the muscle), moving the seat is
the only one that keeps both verified attachments and does not add material.  This searches the
bone surface within `--radius` of each failing seat, scores candidates by whether a pad AND a cap
fit, and writes the best move as an override file that build_shoulder_ends.py consumes.

Run:
    python mujoco_arm/tools/relocate_shoulder_seats.py
    python mujoco_arm/tools/relocate_shoulder_seats.py --radius 25 --json out.json
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_shoulder_ends as BE
import check_shoulder_fixture_zones as FZ
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP

ENDPADS_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_end_pads.json'
OVERRIDES_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_seat_overrides.json'
MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'


def mf():
    return I._manifold()


def fit(p, n, bsolid, bones, modules, own, blade):
    """Best (thickness, radius) for the pad and the cap at this position; None if none fits."""
    t_max = BE.PAD_T if blade is None else float(max(1.0, min(BE.PAD_T, blade - 0.5)))
    best = None
    for cand_t in sorted({t_max, round(t_max * 0.7, 2), 1.5, 1.0}, reverse=True):
        for cand_r in (BE.PAD_R, 10.0, 9.0, 8.0, 7.0, 6.0):
            lift, _ = BE.shoulder_pads.lift_to_clear(p, n, bsolid, radius=cand_r,
                                                      thickness=cand_t)
            lift = 0.0 if lift is None else lift
            pad = BE.disc_rounded(p, n, lift, cand_t, cand_r)
            pb = max(float(I.intersection_volume(pad, b)) for b in bones.values())
            pm = max((float(I.intersection_volume(pad, s)) for k, s in modules.items()
                      if k != own), default=0.0)
            if pb > BE.SURFACE_TOL_MM3 or pm > BE.SURFACE_TOL_MM3:
                continue
            for cap_t in sorted({cand_t, round(cand_t * 0.6, 2), 1.5, 1.0}, reverse=True):
                for cap_r in (BE.PAD_R, 10.0, 9.0, 8.0, 7.0, 6.0):
                    cap = BE.disc_rounded(p, n, lift + cand_t + 0.2, cap_t, cap_r)
                    cb = max(float(I.intersection_volume(cap, b)) for b in bones.values())
                    cm = max((float(I.intersection_volume(cap, s)) for k, s in modules.items()
                              if k != own), default=0.0)
                    if cb <= BE.SURFACE_TOL_MM3 and cm <= BE.SURFACE_TOL_MM3:
                        return {'pad_thickness_mm': round(cand_t, 2),
                                'pad_radius_mm': round(cand_r, 1),
                                'cap_thickness_mm': round(cap_t, 2),
                                'cap_radius_mm': round(cap_r, 1),
                                'lift_mm': round(lift, 2),
                                'blade_thickness_mm': None if blade is None else round(blade, 2)}
            if best is None:
                best = None
    return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--radius', type=float, default=22.0)
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    design = json.loads(ENDPADS_JSON.read_text(encoding='utf-8'))
    bones = {k: SP.load_solid(BE.BONE_DIR / v) for k, v in SP.BONES.items()}
    bone_tris = {}
    for k, v in SP.BONES.items():
        V, F = mesh_geom.load_obj(str(BE.BONE_DIR / v))
        W, G = mesh_geom.weld(V, F)
        bone_tris[k] = W[G]
    modules = {}
    for key in ('deltoid', 'teres_major', 'subscapularis', 'infraspinatus'):
        p = MODULES_DIR / ('%s.stl' % key)
        if p.is_file():
            V, F = mesh_geom.read_stl(str(p))
            W, G = mesh_geom.weld(V, F)
            modules[key] = mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))

    failing = [s for s in design['seats']
               if not s['pad_footprint_cleared'] or not s['cap_footprint_cleared']]
    print('seats to relocate: %d' % len(failing))

    overrides = {}
    rows = []
    for s in failing:
        p0 = np.asarray(s['centre_source_mm'], float)
        bone_key = s['nearest_bone']
        W = SW = None
        bk_tris = bone_tris[bone_key]
        # candidate surface points: the bone's own vertices near the seat
        V, F = mesh_geom.load_obj(str(BE.BONE_DIR / SP.BONES[bone_key]))
        VW, VG = mesh_geom.weld(V, F)
        verts = VW
        d = np.linalg.norm(verts - p0, axis=1)
        near = verts[d <= args.radius]
        if len(near) == 0:
            print('  %-22s no bone vertices within %.0f mm' % (','.join(s['heads']), args.radius))
            continue
        # score by blade thickness first - a thicker blade is what makes a clamp possible
        scored = []
        for q in near:
            n = BE.bone_normal(q, [bk_tris])
            blade = BE.blade_thickness(q, n, bone_key, [None if k != bone_key else bk_tris
                                                       for k in SP.BONES])
            if blade is None:
                continue
            scored.append((blade, float(np.linalg.norm(q - p0)), q, n))
        scored.sort(key=lambda t: (-t[0], t[1]))
        print('  %-22s %d candidates, thickest %.1f mm at %.1f mm away'
              % (','.join(s['heads']), len(scored), scored[0][0], scored[0][1]))

        found = None
        for blade, dist, q, n in scored[:40]:
            got = fit(q, n, bones[bone_key], bones, modules, s['module'], blade)
            if got is not None:
                found = (dist, q, n, blade, got)
                break
        if found is None:
            print('    no fittable position within %.0f mm' % args.radius)
            rows.append({'seat': s['heads'], 'relocated': False})
            continue
        dist, q, n, blade, got = found
        key = ','.join(s['heads'])
        overrides[key] = {'centre_source_mm': np.round(q, 2).tolist(),
                          'normal_source': np.round(n, 3).tolist(), **got,
                          'moved_mm': round(dist, 2),
                          'was_centre_source_mm': np.round(p0, 2).tolist(),
                          'was_blade_mm': s['blade_thickness_mm']}
        rows.append({'seat': key, 'relocated': True, 'moved_mm': round(dist, 2),
                     'blade_mm': round(blade, 2), **got})
        print('    -> moved %.1f mm, blade %.1f mm, pad r%.1f t%.2f, cap r%.1f t%.2f'
              % (dist, blade, got['pad_radius_mm'], got['pad_thickness_mm'],
                 got['cap_radius_mm'], got['cap_thickness_mm']))

    print('\nrelocated %d of %d' % (len(overrides), len(failing)))
    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'relocation of seats where no pad/cap footprint fits at the bone-side '
                     'position (5.9.14); candidates are bone surface vertices within radius, '
                     'scored by local blade thickness then distance, and accepted only when both '
                     'the pad and the cap clear bones and the neighbouring modules',
            'radius_mm': args.radius, 'overrides': overrides, 'rows': rows,
            'not_verified': ['a moved seat changes the channel end: the corridor and the tendon '
                             'segment must be re-run for any seat that moves',
                             'no screw access check at the new position'],
        })
        print('report written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
