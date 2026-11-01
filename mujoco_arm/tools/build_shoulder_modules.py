"""S1 modules: all four TPU muscle envelopes by one recipe.

Round 12 built the deltoid Y.  The other three modules are single-channel, so they use
the same recipe with one corridor instead of three:

    corridor   gradient relaxation inside the bone-trimmed belly, also kept at least a
               tube radius away from every bone (a corridor that hugs bone is cut away
               with the tube around it - §5.4.9)
    tube       R3.2 around the corridor: channel outer 1.6 + 1.6 mm of TPU
    envelope   anatomy trimmed by bone, unioned with the tube, trimmed by bone again
    verdict    one connected solid, bone-clear, and every corridor >= 1.6 mm inside

The deltoid is included so the generalisation can be checked against round 12's numbers:
317 261.9 mm3, one component, channels -5.57 / -6.71 / -3.05.

Run:
    python mujoco_arm/tools/build_shoulder_modules.py
    python mujoco_arm/tools/build_shoulder_modules.py --json out.json
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
import search_shoulder_channel_paths as SP

MUSCLE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / 'muscle'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

CHANNEL_OUTER_R_MM = 1.6
TUBE_R_MM = 3.2

#: module -> (label, heads).  The deltoid's heads share one distal seat; the others are
#: single muscles whose distal end is their own insertion seat.
MODULES = [
    ('deltoid', '三角肌模块（Y 形，3 通道）', ('FJ1468', 'FJ1513', 'FJ1467'), True),
    ('teres_major', '大圆肌模块', ('FJ1507',), False),
    ('subscapularis', '肩胛下肌模块', ('FJ1504',), False),
    ('infraspinatus', '冈下肌模块', ('FJ1500',), False),
]


def mf():
    return I._manifold()


def rod(a, b, r, segments=24):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = b - a
    L = float(np.linalg.norm(d))
    if L < 1e-9:
        return None
    cyl = mf().Manifold.cylinder(L, r, r, segments)
    z = d / L
    helper = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    x = np.cross(helper, z)
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    return cyl.transform(np.c_[np.stack([x, y, z], axis=1), a])


def tube(path, r=TUBE_R_MM):
    pieces = []
    for a, b in zip(path, path[1:]):
        seg = rod(a, b, r)
        if seg is not None:
            pieces.append(seg)
    for p in path:
        pieces.append(mf().Manifold.sphere(r, 20).translate(list(np.asarray(p, float))))
    return mf().Manifold.batch_boolean(pieces, mf().OpType.Add)


def shared_seat(heads):
    pts, w = [], []
    for mid in heads:
        region = attachment_regions.load_region(mid, 'insertion')
        c, a = attachment_regions.centroid(region)
        pts.append(c)
        w.append(a)
    pts = np.array(pts)
    w = np.array(w)
    return (pts * w[:, None]).sum(0) / w.sum()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_modules'))
    ap.add_argument('--bone-clear', type=float, default=TUBE_R_MM,
                    help='how far the corridor is kept from bone [mm]')
    args = ap.parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    bones = {k: SP.load_solid(BONE_DIR / v) for k, v in SP.BONES.items()}
    bone_tris = []
    for v in SP.BONES.values():
        V, F = mesh_geom.load_obj(str(BONE_DIR / v))
        W, G = mesh_geom.weld(V, F)
        bone_tris.append(W[G])

    modules = []
    for key, label, heads, merge_distal in MODULES:
        print('\n=== %s (%s) ===' % (label, '+'.join(heads)))
        # anatomy: each head trimmed by bone
        anatomy = None
        head_rows = []
        for mid in heads:
            solid = SP.load_solid(MUSCLE_DIR / (mid + '.obj'))
            trimmed = solid
            for b in bones.values():
                trimmed = trimmed - b
            comps = len([c for c in trimmed.decompose() if c.volume() > 0.5])
            anatomy = trimmed if anatomy is None else anatomy + trimmed
            head_rows.append({'head': mid, 'anatomy_mm3': round(solid.volume(), 1),
                              'trimmed_mm3': round(trimmed.volume(), 1),
                              'components': comps})
            print('  %-8s anatomy %9.1f -> trimmed %9.1f mm3 (%d piece%s)'
                  % (mid, solid.volume(), trimmed.volume(), comps,
                     '' if comps == 1 else 's'))
        seat = shared_seat(heads)
        av, af = mesh_geom.manifold_mesh_to_arrays(anatomy.to_mesh())
        AW, AG = mesh_geom.weld(av, af)
        atri = AW[AG]

        corridors = []
        trunk = None
        for mid in heads:
            origin = attachment_regions.load_region(mid, 'origin')
            o, _ = attachment_regions.centroid(origin)
            ts = np.linspace(0.0, 1.0, SP.INTERIOR_POINTS)
            start = [o + (seat - o) * t for t in ts]
            relaxed, _ = SP.relax(start, atri, iterations=SP.ITERATIONS,
                                  step=SP.STEP_MM, lam=SP.SMOOTHING,
                                  bone_tris=bone_tris, bone_clear_mm=args.bone_clear)
            relaxed[-1] = seat
            band = [relaxed[k] for k, t in enumerate(ts) if SP.BAND[0] <= t <= SP.BAND[1]]
            worst, deepest, _ = SP.clearance(band, atri)
            t_tube = tube(relaxed)
            trunk = t_tube if trunk is None else trunk + t_tube
            corridors.append({'head': mid, 'worst_belly_mm': round(worst, 3),
                              'deepest_belly_mm': round(deepest, 3),
                              'length_mm': round(float(sum(np.linalg.norm(
                                  np.diff(relaxed, axis=0), axis=1))), 2),
                              'path_source_mm': [np.round(p, 2).tolist() for p in relaxed]})
            print('  corridor %-8s worst %6.2f  deepest %6.2f  length %6.1f mm'
                  % (mid, worst, deepest, corridors[-1]['length_mm']))

        raw = trunk.volume()
        raw_comps = len([c for c in trunk.decompose() if c.volume() > 0.5])
        trimmed_trunk = trunk
        removed = {}
        for n, b in bones.items():
            removed[n] = round(I.intersection_volume(trimmed_trunk, b), 2)
            trimmed_trunk = trimmed_trunk - b
        final = (anatomy + trimmed_trunk)
        for b in bones.values():
            final = final - b
        fcomps = [c for c in final.decompose() if c.volume() > 0.5]
        added = final.volume() - anatomy.volume()

        fv, ff = mesh_geom.manifold_mesh_to_arrays(final.to_mesh())
        FW, FG = mesh_geom.weld(fv, ff)
        ftri = FW[FG]
        ch = []
        for c in corridors:
            path = c['path_source_mm']
            ts = np.linspace(0.0, 1.0, len(path))
            idx = [k for k, t in enumerate(ts) if SP.BAND[0] <= t <= SP.BAND[1]]
            band = [np.asarray(path[k]) for k in idx]
            d = [mesh_geom.signed_distance(p, ftri) for p in band]
            k_worst = int(np.argmax(d))
            ch.append({'head': c['head'], 'worst_mm': round(max(d), 3),
                       'deepest_mm': round(min(d), 3),
                       'worst_at_t': round(float(ts[idx[k_worst]]), 3),
                       'worst_to_pad_mm': round(float(
                           np.linalg.norm(np.asarray(path[-1]) - band[k_worst])), 2),
                       'fits': bool(max(d) <= -CHANNEL_OUTER_R_MM)})
        bone_hits = {n: round(I.intersection_volume(final, b), 3) for n, b in bones.items()}
        total_bone = sum(bone_hits.values())
        one = len(fcomps) == 1
        all_fit = all(x['fits'] for x in ch)
        frags = sorted((round(c.volume(), 1) for c in fcomps), reverse=True)

        target = out_dir / ('%s.stl' % key)
        mesh_geom.write_stl(str(target), fv, ff, header='S1 module %s' % key)
        RV, RF = mesh_geom.read_stl(str(target))
        RW, RG = mesh_geom.weld(RV, RF)
        rt = abs(mesh_geom.signed_volume(RW, RG) - final.volume()) / final.volume()

        print('  trunk %.1f mm3 (%d piece%s) -> trimmed at bone %s'
              % (raw, raw_comps, '' if raw_comps == 1 else 's',
                 ', '.join('%s %.1f' % kv for kv in removed.items() if kv[1] > 0)
                 or 'none'))
        print('  FINAL %.1f mm3 = anatomy %.1f + added %.1f (%.3f%%), %d piece%s, %d tris'
              % (final.volume(), anatomy.volume(), added,
                 100.0 * added / max(anatomy.volume(), 1e-9), len(fcomps),
                 '' if len(fcomps) == 1 else 's', final.num_tri()))
        if len(fcomps) > 1:
            print('    fragment volumes: %s' % ', '.join('%.1f' % v for v in frags))
            big = max(fcomps, key=lambda c: c.volume())
            for c in fcomps:
                if c is big:
                    continue
                bb = c.bounding_box()
                print('    other fragment %.1f mm3 at x %.0f..%.0f y %.0f..%.0f z %.0f..%.0f'
                      % (c.volume(), bb[0], bb[3], bb[1], bb[4], bb[2], bb[5]))
        print('  bone %s | channels %s'
              % ('clear' if total_bone <= I.HARD_MM3
                 else ', '.join('%.1f' % v for v in bone_hits.values()),
                 ', '.join('%s %s (worst at t=%.2f, %.1f mm from the pad)'
                           % (x['head'], 'fit' if x['fits'] else 'SHORT',
                              x['worst_at_t'], x['worst_to_pad_mm']) for x in ch)))

        modules.append({'module': key, 'label': label, 'heads': list(heads),
                        'merges_to_one_seat': bool(merge_distal),
                        'seat_source_mm': np.round(seat, 2).tolist(),
                        'heads_anatomy': head_rows,
                        'corridors': corridors,
                        'trunk': {'volume_untrimmed_mm3': round(raw, 1),
                                  'components_untrimmed': raw_comps,
                                  'trimmed_at_bone_mm3': removed},
                        'final': {'volume_mm3': round(final.volume(), 1),
                                  'anatomy_mm3': round(anatomy.volume(), 1),
                                  'added_mm3': round(added, 1),
                                  'added_percent_of_anatomy': round(
                                      100.0 * added / max(anatomy.volume(), 1e-9), 3),
                                  'components': len(fcomps),
                                  'fragment_volumes_mm3': frags,
                                  'triangles': int(final.num_tri()),
                                  'stl_file': target.name,
                                  'stl_roundtrip_rel_error': float('%.2e' % rt)},
                        'bone_intersections_mm3': bone_hits,
                        'channels': ch,
                        'one_piece': bool(one),
                        'bone_clear': bool(total_bone <= I.HARD_MM3),
                        'all_channels_fit': bool(all_fit)})

    print('\n=== S1 module set ===')
    print('%-15s %10s %8s %9s %6s %6s %s'
          % ('module', 'volume', 'added', 'added%', 'piece', 'bone', 'channels'))
    for m in modules:
        f = m['final']
        print('%-15s %10.1f %8.1f %9.3f %6d %6s %s'
              % (m['module'], f['volume_mm3'], f['added_mm3'],
                 f['added_percent_of_anatomy'], f['components'],
                 'clear' if m['bone_clear'] else 'HIT',
                 'all fit' if m['all_channels_fit'] else 'SHORT'))
    ok = all(m['one_piece'] and m['bone_clear'] and m['all_channels_fit'] for m in modules)
    print('\nall four modules: one piece + bone-clear + channels fit: %s'
          % ('YES' if ok else 'NO'))

    report = {
        'scope': 'one recipe for all four modules: bone-aware corridor relaxation, R3.2 '
                 'tube, anatomy trimmed by bone, union, trimmed again; no lattice, no '
                 'end caps, no origin pads',
        'channel_outer_r_mm': CHANNEL_OUTER_R_MM,
        'tube_radius_mm': TUBE_R_MM,
        'modules': modules,
        'all_modules_pass': bool(ok),
        'not_verified': ['tubes are capsule chains, not lofted bellies',
                         'no lattice, end caps, crimp, bolts or origin pads',
                         'bone interference only; neighbour bellies, fixture pads and '
                         'the 0-120 deg sweep are not checked',
                         'no muscle force, pennation, contact or printability'],
        'frozen': False,
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
