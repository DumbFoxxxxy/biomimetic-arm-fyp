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
import shoulder_pads

MUSCLE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / 'muscle'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

CHANNEL_OUTER_R_MM = 1.6
TUBE_R_MM = 3.2
#: fragments below this are trimmed slivers, not anatomical lobes [mm3]
DEBRIS_MM3 = 5.0

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


def point_in(solid, point, half=0.15):
    """Is `point` inside `solid`?  A 0.3 mm probe cube is enough for a seat centre."""
    probe = mf().Manifold.cube([2 * half] * 3).translate(
        list(np.asarray(point, float) - half))
    return float((probe ^ solid).volume()) > 0.5 * (2 * half) ** 3


def nearest_points(a, b):
    """Closest pair of points between two solids, by vertex-to-surface search."""
    va, fa = mesh_geom.manifold_mesh_to_arrays(a.to_mesh())
    vb, fb = mesh_geom.manifold_mesh_to_arrays(b.to_mesh())
    WA, GA = mesh_geom.weld(va, fa)
    WB, GB = mesh_geom.weld(vb, fb)
    triA, triB = WA[GA], WB[GB]

    def closest(verts, tri):
        best = (np.inf, None)
        for v in verts:
            d = mesh_geom.point_triangle_distance(v, tri).min()
            if d < best[0]:
                best = (float(d), np.asarray(v, float))
        return best

    da, pa = closest(WA, triB)
    db, pb = closest(WB, triA)
    return (pa, pb, da) if da <= db else (pb, pa, db)


def bridge(lobes, bone_tris, bone_solids, radius=TUBE_R_MM, offsets=(4.0, 6.0, 8.0, 12.0, 16.0),
           sink_mm=4.0):
    """Join two lobes with an arc that goes over the bone between them.

    D18 option A: the lobes are one muscle, split only because the muscle mesh
    interpenetrated the scapula and the connecting material sat inside it.  The bridge
    has to clear the ridge, so the arc apex is pushed along the bone's outward normal
    until the rods are bone-clear.

    The trick that makes it work: the rod ends must be **inside** each lobe, not spheres
    sitting on the lobe surface.  The closest-point pair is only 1.47 mm from bone, so a
    radius-3.2 sphere centred there intersects bone no matter where the apex goes - that
    is exactly why the first attempt clashed at every offset.  Each end is therefore sunk
    `sink_mm` into its own lobe, where the material is already bone-clear.
    """
    a, b, gap = nearest_points(lobes[0], lobes[1])
    mid = 0.5 * (a + b)
    picked = None
    for tri in bone_tris:
        d, n = mesh_geom.signed_distance_and_normal(mid, tri)
        if picked is None or abs(d) < abs(picked[0]):
            picked = (d, n)
    _, n_out = picked
    va, fa = mesh_geom.manifold_mesh_to_arrays(lobes[0].to_mesh())
    vb, fb = mesh_geom.manifold_mesh_to_arrays(lobes[1].to_mesh())
    WA, GA = mesh_geom.weld(va, fa)
    WB, GB = mesh_geom.weld(vb, fb)
    _, n_a = mesh_geom.signed_distance_and_normal(a, WA[GA])
    _, n_b = mesh_geom.signed_distance_and_normal(b, WB[GB])
    a_in = a - n_a * sink_mm
    b_in = b - n_b * sink_mm
    fallback = None
    for off in offsets:
        apex = mid + n_out * off
        pieces = [p for p in (rod(a_in, apex, radius), rod(apex, b_in, radius),
                              mf().Manifold.sphere(radius, 20).translate(list(a_in)),
                              mf().Manifold.sphere(radius, 20).translate(list(b_in)))
                  if p is not None]
        solid = mf().Manifold.batch_boolean(pieces, mf().OpType.Add)
        clash = sum(I.intersection_volume(solid, s) for s in bone_solids)
        rec = {'solid': solid, 'apex_source_mm': np.round(apex, 2).tolist(),
               'endpoints_sunk_mm': sink_mm, 'offset_mm': off,
               'gap_mm': round(gap, 3), 'bone_clash_mm3': round(clash, 3)}
        if clash <= I.HARD_MM3:
            return rec
        if fallback is None:
            fallback = rec
    return fallback


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
        # the origin side of the single-channel modules: that muscle's own origin pad
        seat_origin = None
        if len(heads) == 1:
            region = attachment_regions.load_region(heads[0], 'origin')
            seat_origin, _ = attachment_regions.centroid(region)
        else:
            region = attachment_regions.load_region(heads[0], 'origin')
            seat_origin, _ = attachment_regions.centroid(region)
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
        # debris: fragments far below any muscle scale are trimmed slivers, not parts.
        # Real lobes are kept and reported - dropping one that carries a seat would
        # silently delete an attachment.
        debris = [c for c in fcomps if c.volume() < DEBRIS_MM3]
        if debris:
            keep = [c for c in fcomps if c.volume() >= DEBRIS_MM3]
            if keep:
                merged_final = keep[0]
                for c in keep[1:]:
                    merged_final = merged_final + c
                final = merged_final
                fcomps = [c for c in final.decompose() if c.volume() > 0.5]
            print('  dropped %d debris fragment(s) totalling %.1f mm3 (< %.0f mm3)'
                  % (len(debris), sum(c.volume() for c in debris), DEBRIS_MM3))
        # if the lobes cannot be joined in TPU, they can still be joined by the origin
        # pad and cap that both of them sit on - that is how every muscle end is captured
        # anyway (§5.7).  Measure it rather than assume it.
        pad_rows = []
        if len(fcomps) > 1:
            probe = mf().Manifold.sphere(0.25, 12)
            for c in sorted(fcomps, key=lambda x: -x.volume()):
                d_o = I.gap(c, probe.translate(list(np.asarray(seat_origin, float))))
                pad_rows.append({'fragment_mm3': round(c.volume(), 1),
                                 'distance_to_origin_seat_mm': round(d_o, 2),
                                 'within_pad_radius': bool(
                                     d_o <= shoulder_pads.PAD_RADIUS_MM)})
            all_on_pad = all(r['within_pad_radius'] for r in pad_rows)
            print('  lobes on the shared origin pad: %s (pad R%.0f)'
                  % (', '.join('%.1f mm3 at %.2f mm' % (r['fragment_mm3'],
                                                        r['distance_to_origin_seat_mm'])
                               for r in pad_rows),
                     shoulder_pads.PAD_RADIUS_MM))
            print('  -> join through the pad and cap: %s' % ('YES' if all_on_pad else 'NO'))
        # A TPU bridge is still attempted, non-destructively: if it does not reduce the
        # piece count we keep the unbridged part and record the failure rather than
        # shipping a stub.
        lobe_rows = []
        guard = 0
        while len(fcomps) > 1 and guard < 3:
            guard += 1
            before = len(fcomps)
            lobes = sorted(fcomps, key=lambda c: -c.volume())[:2]
            br = bridge(lobes, bone_tris, list(bones.values()))
            merged_bridge = br['solid']
            for c in fcomps:
                merged_bridge = merged_bridge + c
            for b in bones.values():
                merged_bridge = merged_bridge - b
            trial = [c for c in merged_bridge.decompose() if c.volume() > 0.5]
            row = {'gap_mm': br['gap_mm'],
                   'arc_apex_source_mm': br['apex_source_mm'],
                   'arc_offset_mm': br['offset_mm'],
                   'bone_clash_mm3': br['bone_clash_mm3'],
                   'pieces_before': before, 'pieces_after': len(trial),
                   'accepted': bool(len(trial) < before)}
            lobe_rows.append(row)
            print('  bridge attempt: gap %.2f mm, apex %.0f mm out, bone clash '
                  '%.2f mm3 -> %d piece(s) %s'
                  % (br['gap_mm'], br['offset_mm'], br['bone_clash_mm3'], len(trial),
                     'ACCEPTED' if row['accepted'] else 'REJECTED'))
            if not row['accepted']:
                break
            final = merged_bridge
            fcomps = trial
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
                # a seat sits ON a surface, so point-in-solid is boundary-sensitive and
                # useless here; the fragment whose surface is nearest the seat is the one
                # that carries it
                probe = mf().Manifold.sphere(0.25, 12)
                d_o = I.gap(c, probe.translate(list(np.asarray(seat_origin, float))))
                d_s = I.gap(c, probe.translate(list(np.asarray(seat, float))))
                bb = c.bounding_box()
                print('    fragment %9.1f mm3  gap to origin seat %6.2f mm, '
                      'to distal seat %6.2f mm  x %.0f..%.0f y %.0f..%.0f z %.0f..%.0f'
                      % (c.volume(), d_o, d_s, bb[0], bb[3], bb[1], bb[4], bb[2], bb[5]))
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
                        'lobe_bridges': lobe_rows,
                        'lobe_pad_join': {
                            'pad_radius_mm': shoulder_pads.PAD_RADIUS_MM,
                            'fragments': pad_rows,
                            'joined_by_pad': bool(pad_rows) and all(
                                r['within_pad_radius'] for r in pad_rows)},
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
