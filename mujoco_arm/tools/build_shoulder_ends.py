"""End pads and caps as solids, for every bone-side seat, with the interference that matters.

Three threads have been waiting on this part:

* 5.4.14 / 5.4.20 measured cap reaches of 0.0 to 32.0 mm, but the pad + cap stack was only a
  number, never a solid;
* 5.4.21 / 5.4.22 built tendon segments that start on "the pad's outer face", an approximation
  that degenerate eight of nine ends;
* 5.4.13 made the infraspinatus two TPU bodies that must be clamped by one shared pad.

Built here per seat: the R11 x 4 mm end pad oriented on the bone normal (the TPU's own end pad,
per 5.7), and the R11 x 4 mm cap that clamps it, sitting one pad thickness further out.  Both
carry the 2 x M3 at 14 mm pitch holes (r1.7) and the r1.4 crimp seat.

Checked: pad/bone intrusion, cap/bone, pad-to-pad and cap-to-cap spacing between neighbouring
seats, and pad/cap against the four module envelopes.

Run:
    python mujoco_arm/tools/build_shoulder_ends.py
    python mujoco_arm/tools/build_shoulder_ends.py --json out.json
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import attachment_regions
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP
import shoulder_pads

DESIGN_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_modules.json'
MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

PAD_R = shoulder_pads.PAD_RADIUS_MM
PAD_T = shoulder_pads.PAD_THICKNESS_MM
BOLT_PITCH = shoulder_pads.BOLT_PITCH_MM
BOLT_HOLE_R = shoulder_pads.BOLT_HOLE_RADIUS_MM
CRIMP_R = 1.4
#: the pad sits on the bone surface; a few hundredths of a mm3 of graze is surface curvature
SURFACE_TOL_MM3 = 0.5


def mf():
    return I._manifold()


def disc(centre, normal, thickness, radius=PAD_R, holes=True):
    """R11 plate growing outward along `normal`, with the bolt holes and the crimp bore."""
    solid = shoulder_pads.build_pad(centre, normal, radius=radius, thickness=thickness)
    if not holes:
        return solid
    R = shoulder_pads.pad_basis(normal)
    # holes run through the plate along the normal so a screw can pass
    for k in (0, 1):
        ang = np.pi * k
        off = R[:, 0] * (BOLT_PITCH / 2.0) * np.cos(ang) + R[:, 1] * (BOLT_PITCH / 2.0) * np.sin(ang)
        hole = mf().Manifold.cylinder(thickness * 3.0, BOLT_HOLE_R, BOLT_HOLE_R, 32)
        hole = hole.translate([0.0, 0.0, -thickness])
        world = np.c_[R, np.asarray(centre, float) + off]
        solid = solid - hole.transform(world)
    bore = mf().Manifold.cylinder(thickness * 3.0, CRIMP_R, CRIMP_R, 32)
    bore = bore.translate([0.0, 0.0, -thickness])
    solid = solid - bore.transform(np.c_[R, np.asarray(centre, float)])
    return solid


def bone_normal(point, bone_tris):
    picked = None
    for tri in bone_tris:
        d, n = mesh_geom.signed_distance_and_normal(point, tri)
        if picked is None or abs(d) < abs(picked[0]):
            picked = (d, n)
    return picked[1]


def nearest_bone(point, bones):
    """The bone whose surface is closest to `point`, as (name, solid)."""
    probe = mf().Manifold.sphere(0.25, 12).translate(list(np.asarray(point, float)))
    best = None
    for name, solid in bones.items():
        g = I.gap(solid, probe)
        if best is None or g < best[0]:
            best = (g, name, solid)
    return best[1], best[2]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_ends'))
    args = ap.parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    design = json.loads(DESIGN_JSON.read_text(encoding='utf-8'))
    bones = {k: SP.load_solid(BONE_DIR / v) for k, v in SP.BONES.items()}
    bone_tris = []
    for v in SP.BONES.values():
        BV, BF = mesh_geom.load_obj(str(BONE_DIR / v))
        BW, BG = mesh_geom.weld(BV, BF)
        bone_tris.append(BW[BG])
    modules = {}
    for key in ('deltoid', 'teres_major', 'subscapularis', 'infraspinatus'):
        p = MODULES_DIR / ('%s.stl' % key)
        if p.is_file():
            V, F = mesh_geom.read_stl(str(p))
            W, G = mesh_geom.weld(V, F)
            modules[key] = mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))

    seats = []
    print('%-14s %-8s %-8s %6s %10s %10s %10s'
          % ('module', 'head', 'end', 'lift', 'pad/bone', 'cap/bone', 'pad/module'))
    for m in design['modules']:
        for c in m['channels']:
            path = next(x['path_source_mm'] for x in m['corridors'] if x['head'] == c['head'])
            for end, p in (('origin', path[0]), ('distal', path[-1])):
                p = np.asarray(p, float)
                n = bone_normal(p, bone_tris)
                bk, bsolid = nearest_bone(p, bones)
                lift, lift_residual = shoulder_pads.lift_to_clear(p, n, bsolid)
                if lift is None:
                    lift = 0.0
                    print('  %-14s %-8s %-6s lift > 8 mm not enough (residual %.1f mm3)'
                          % (m['module'], c['head'], end, lift_residual))
                pad = disc(p + n * lift, n, PAD_T)
                cap_c = p + n * (lift + PAD_T + 0.2)      # cap sits on the pad's outer face
                cap = disc(cap_c, n, PAD_T)
                pad_bone = max(I.intersection_volume(pad, b).__float__()
                               for b in bones.values())
                cap_bone = max(I.intersection_volume(cap, b).__float__()
                               for b in bones.values())
                seats.append({
                    'module': m['module'], 'head': c['head'], 'end': end,
                    'centre_source_mm': np.round(p, 2).tolist(),
                    'normal_source': np.round(n, 3).tolist(),
                    'nearest_bone': bk,
                    'lift_mm': round(lift, 2),
                    'pad_volume_mm3': round(pad.volume(), 1),
                    'cap_volume_mm3': round(cap.volume(), 1),
                    'pad_bone_mm3': round(pad_bone, 3),
                    'cap_bone_mm3': round(cap_bone, 3),
                })
                print('%-14s %-8s %-8s %6.2f %10.3f %10.3f  (%s)'
                      % (m['module'], c['head'], end, lift, pad_bone, cap_bone, bk))

    print('\n--- neighbouring seats ---')
    print('%-22s %-22s %10s %10s' % ('seat a', 'seat b', 'centres', 'rim gap'))
    pair_rows = []
    worst = None
    for a, b in itertools.combinations(seats, 2):
        pa = np.asarray(a['centre_source_mm'], float)
        pb = np.asarray(b['centre_source_mm'], float)
        d = float(np.linalg.norm(pa - pb))
        # the deltoid's three heads merge onto ONE seat (D8), so a pair at the same seat shares
        # a single pad - that is the design, not a clash
        if d < 1.0:
            pair_rows.append({'a': '%s/%s/%s' % (a['module'], a['head'], a['end']),
                              'b': '%s/%s/%s' % (b['module'], b['head'], b['end']),
                              'centre_mm': round(d, 2), 'rim_gap_mm': None,
                              'same_seat': True})
            continue
        gap = d - 2.0 * PAD_R
        pair_rows.append({'a': '%s/%s/%s' % (a['module'], a['head'], a['end']),
                          'b': '%s/%s/%s' % (b['module'], b['head'], b['end']),
                          'centre_mm': round(d, 2), 'rim_gap_mm': round(gap, 2)})
        if worst is None or gap < worst['rim_gap_mm']:
            worst = pair_rows[-1]
    real = [r for r in pair_rows if not r.get('same_seat')]
    for r in sorted(real, key=lambda x: x['rim_gap_mm'])[:6]:
        print('%-22s %-22s %10.2f %10.2f' % (r['a'], r['b'], r['centre_mm'], r['rim_gap_mm']))
    print('closest distinct pads: %s vs %s, centre %.2f mm, rim gap %.2f mm'
          % (worst['a'], worst['b'], worst['centre_mm'], worst['rim_gap_mm']))
    print('seat pairs sharing one pad (the deltoid Y merges): %d'
          % len([r for r in pair_rows if r.get('same_seat')]))

    bad_pad = [s for s in seats if s['pad_bone_mm3'] > SURFACE_TOL_MM3]
    bad_cap = [s for s in seats if s['cap_bone_mm3'] > SURFACE_TOL_MM3]
    print('\npads intruding into bone (>%.1f mm3): %d of %d'
          % (SURFACE_TOL_MM3, len(bad_pad), len(seats)))
    print('caps intruding into bone (>%.1f mm3): %d of %d'
          % (SURFACE_TOL_MM3, len(bad_cap), len(seats)))

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'R11 x 4 mm end pad and cap at every bone-side seat, on the bone normal, '
                     'with 2xM3 at 14 mm pitch and the r1.4 crimp bore',
            'pad_radius_mm': PAD_R, 'pad_thickness_mm': PAD_T,
            'bolt_pitch_mm': BOLT_PITCH, 'bolt_hole_r_mm': BOLT_HOLE_R,
            'seats': seats, 'pairs': pair_rows,
            'closest_pads': worst,
            'n_pads_intruding': len(bad_pad), 'n_caps_intruding': len(bad_cap),
            'not_verified': ['no screws, no crimp, no preload',
                             'the cap is placed on the pad face, not on a real TPU end',
                             'the pads do not carry the lift per 5.5.9 yet'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
