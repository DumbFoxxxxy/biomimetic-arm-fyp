"""S1 deltoid Y trunk: a real reconstructed belly from the three searched corridors.

D8 says the three heads are rebuilt as one Y rather than unioned.  Round 11 supplied the
input that makes the rebuild concrete: five of six lines already have a corridor inside
their anatomical belly, and the deltoid middle head is short by 4.27 mm.  So the trunk
has to carry three channels and give the middle branch that extra material.

Construction, deliberately minimal and checkable: each head's *searched* corridor is
turned into a tube of radius r (r = bore 1.6 + wall 1.6 = 3.2 mm), the three tubes are
unioned, and the result is measured against:

1. **connectivity** - one solid or not;
2. **bone** - a tube of r = 3.2 around a corridor that was relaxed inside the belly is
   not automatically bone-clear, because the belly surface includes the bone surface;
3. **the anatomical envelope** - how much the reconstructed trunk protrudes beyond the
   trimmed union from §5.4.6.  That protrusion is the price of D8 and must be on the
   drawing, not hidden.

The bore itself is not modelled: with tube r 3.2 and bore r 1.6 the wall is exactly
1.6 mm by construction, which is the criterion, so building the bore would only add
tessellation.

Run:
    python mujoco_arm/tools/build_shoulder_y_trunk.py
    python mujoco_arm/tools/build_shoulder_y_trunk.py --json out.json
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
DELTOID = {'FJ1468': '三角肌前束', 'FJ1513': '三角肌后束', 'FJ1467': '三角肌中束'}

#: elbow channel, §5.7: the wire lumen is R0.85 and its outer wall R1.6, so the material
#: only has to clear 1.6 mm from the corridor centre.  The tube built here is R3.2, i.e.
#: the channel plus 1.6 mm of TPU all round.
CHANNEL_OUTER_R_MM = 1.6
TUBE_R_MM = 3.2


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
    """Capsule chain through the path: rods between samples plus spheres at the joints."""
    pieces = []
    for a, b in zip(path, path[1:]):
        seg = rod(a, b, r)
        if seg is not None:
            pieces.append(seg)
    for p in path:
        pieces.append(mf().Manifold.sphere(r, 20).translate(list(np.asarray(p, float))))
    return mf().Manifold.batch_boolean(pieces, mf().OpType.Add)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_deltoid_y'))
    args = ap.parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    bones = {k: SP.load_solid(BONE_DIR / v) for k, v in SP.BONES.items()}
    bone_tris = []
    for n, f in SP.BONES.items():
        V, F = mesh_geom.load_obj(str(BONE_DIR / f))
        W, G = mesh_geom.weld(V, F)
        bone_tris.append(W[G])
    print('tube radius %.1f mm around the corridor; the material must clear %.1f mm '
          '(elbow channel outer radius, §5.7)' % (TUBE_R_MM, CHANNEL_OUTER_R_MM))
    print('corridors are kept >= %.1f mm from bone as well: the belly is trimmed at the '
          'bone, so "inside the belly" alone would be cut away' % TUBE_R_MM)

    # all three corridors must end at ONE shared seat, otherwise the "Y" is three
    # separate tubes (§5.4.6 computed this seat for the merged envelope)
    ins, weights = [], []
    for mid in DELTOID:
        region = attachment_regions.load_region(mid, 'insertion')
        c, a = attachment_regions.centroid(region)
        ins.append(c)
        weights.append(a)
    pts = np.array(ins)
    w = np.array(weights)
    seat = (pts * w[:, None]).sum(0) / w.sum()
    print('shared distal seat %s' % np.round(seat, 1))

    paths, rows = {}, []
    print('\n--- corridors (same relaxation as §5.4.8, ending at the shared seat) ---')
    print('%-8s %-9s %8s %8s %7s' % ('head', 'name', 'worst', 'length', 'points'))
    for mid, label in DELTOID.items():
        origin = attachment_regions.load_region(mid, 'origin')
        o, _ = attachment_regions.centroid(origin)
        solid = SP.load_solid(MUSCLE_DIR / (mid + '.obj'))
        trimmed = solid
        for b in bones.values():
            trimmed = trimmed - b
        verts, tris = mesh_geom.manifold_mesh_to_arrays(trimmed.to_mesh())
        W, G = mesh_geom.weld(verts, tris)
        tri = W[G]
        ts = np.linspace(0.0, 1.0, SP.INTERIOR_POINTS)
        start = [o + (seat - o) * t for t in ts]
        relaxed, _ = SP.relax(start, tri, iterations=SP.ITERATIONS,
                              step=SP.STEP_MM, lam=SP.SMOOTHING,
                              bone_tris=bone_tris, bone_clear_mm=TUBE_R_MM)
        relaxed[-1] = seat
        band = [relaxed[k] for k, t in enumerate(ts) if SP.BAND[0] <= t <= SP.BAND[1]]
        worst, _, _ = SP.clearance(band, tri)
        paths[mid] = relaxed
        rows.append({'head': mid, 'name': label,
                     'worst_belly_mm': round(worst, 3),
                     'length_mm': round(float(sum(np.linalg.norm(
                         np.diff(relaxed, axis=0), axis=1))), 2),
                     'points': len(relaxed),
                     'path_source_mm': [np.round(p, 2).tolist() for p in relaxed]})
        print('%-8s %-9s %8.2f %8.1f %7d'
              % (mid, label, worst, rows[-1]['length_mm'], len(relaxed)))

    print('\n--- the Y trunk ---')
    trunk = tube(paths['FJ1468'])
    for mid in ('FJ1513', 'FJ1467'):
        trunk = trunk + tube(paths[mid])
    raw_volume = trunk.volume()
    comps_raw = len([c for c in trunk.decompose() if c.volume() > 0.5])
    print('  before trimming: volume %.1f mm3, components %d' % (raw_volume, comps_raw))

    # the tube ends sit on the bone by construction (they start at the bone-mounted
    # pads), so the part is trimmed by bone exactly like the anatomical envelope
    trimmed_trunk = trunk
    removed = {}
    for n, b in bones.items():
        v = I.intersection_volume(trimmed_trunk, b)
        removed[n] = round(v, 2)
        trimmed_trunk = trimmed_trunk - b
    trunk = trimmed_trunk
    comps = [c for c in trunk.decompose() if c.volume() > 0.5]
    print('  trimmed at bone: %s  -> volume %.1f mm3, components %d'
          % (', '.join('%s %.1f mm3' % kv for kv in removed.items()),
             trunk.volume(), len(comps)))
    print('  triangles %d' % trunk.num_tri())

    print('\n  bone interference after trimming:')
    bone_hits = {}
    for n, b in bones.items():
        v = I.intersection_volume(trunk, b)
        bone_hits[n] = round(v, 2)
        print('    %-9s %8.2f mm3' % (n, v))
    total_bone = sum(bone_hits.values())

    print('\n  protrusion beyond the anatomical envelope (§5.4.6):')
    anatomical = None
    for mid in DELTOID:
        solid = SP.load_solid(MUSCLE_DIR / (mid + '.obj'))
        trimmed = solid
        for b in bones.values():
            trimmed = trimmed - b
        anatomical = trimmed if anatomical is None else anatomical + trimmed
    outside = trunk - anatomical
    pct = 100.0 * outside.volume() / max(trunk.volume(), 1e-9)
    print('    trunk %.1f mm3, outside the anatomy %.1f mm3 (%.1f %%)'
          % (trunk.volume(), outside.volume(), pct))

    # ---- the finished envelope: anatomy kept, trunk added only where it is short
    final = (anatomical + trunk) - bones['scapula'] - bones['clavicle'] - bones['humerus']
    fcomps = [c for c in final.decompose() if c.volume() > 0.5]
    added = final.volume() - anatomical.volume()
    print('\n  finished envelope = anatomy + trunk, trimmed by bone:')
    print('    volume %.1f mm3 (anatomy %.1f + added %.1f), components %d, triangles %d'
          % (final.volume(), anatomical.volume(), added, len(fcomps), final.num_tri()))
    fverts, ftris = mesh_geom.manifold_mesh_to_arrays(final.to_mesh())
    FW, FG = mesh_geom.weld(fverts, ftris)
    ftri = FW[FG]
    print('\n  channel wall inside the finished envelope (need worst <= -%.1f):'
          % CHANNEL_OUTER_R_MM)
    ch_rows = []
    for mid, label in DELTOID.items():
        band = [paths[mid][k] for k, t in enumerate(np.linspace(0.0, 1.0, len(paths[mid])))
                if SP.BAND[0] <= t <= SP.BAND[1]]
        worst, deepest, _ = SP.clearance(band, ftri)
        ch_rows.append({'head': mid, 'name': label,
                        'worst_mm': round(worst, 3), 'deepest_mm': round(deepest, 3),
                        'fits': bool(worst <= -CHANNEL_OUTER_R_MM)})
        print('    %-8s %-9s worst %6.2f  deepest %6.2f  %s'
              % (mid, label, worst, deepest,
                 'FITS' if worst <= -CHANNEL_OUTER_R_MM else 'SHORT'))
    all_fit = all(r['fits'] for r in ch_rows)

    verts, tris = mesh_geom.manifold_mesh_to_arrays(trunk.to_mesh())
    target = out_dir / 'deltoid_Y_trunk.stl'
    mesh_geom.write_stl(str(target), verts, tris, header='S1 deltoid Y trunk')
    RV, RF = mesh_geom.read_stl(str(target))
    RW, RG = mesh_geom.weld(RV, RF)
    rt = abs(mesh_geom.signed_volume(RW, RG) - trunk.volume()) / trunk.volume()
    print('  written %s  roundtrip rel %.2e' % (target.name, rt))

    final_target = out_dir / 'deltoid_Y_final.stl'
    mesh_geom.write_stl(str(final_target), fverts, ftris, header='S1 deltoid Y envelope')
    FV, FF = mesh_geom.read_stl(str(final_target))
    FW2, FG2 = mesh_geom.weld(FV, FF)
    frt = abs(mesh_geom.signed_volume(FW2, FG2) - final.volume()) / final.volume()
    print('  written %s  roundtrip rel %.2e' % (final_target.name, frt))

    verdict_one = len(comps) == 1
    print('\n--- verdict ---')
    print('  one connected Y after trimming: %s'
          % ('YES' if verdict_one else 'NO (%d pieces)' % len(comps)))
    print('  bone-clear after trimming: %s (%.2f mm3)'
          % ('YES' if total_bone <= I.HARD_MM3 else 'NO', total_bone))

    report = {
        'scope': 'trunk built as three capsule chains of r 3.2 mm around the searched '
                 'corridors; bore not modelled (wall is 1.6 mm by construction); no '
                 'lattice, no end caps, no origin pads',
        'tube_radius_mm': TUBE_R_MM, 'channel_outer_r_mm': CHANNEL_OUTER_R_MM,
        'shared_seat_source_mm': np.round(seat, 2).tolist(),
        'corridors': rows,
        'trunk': {'volume_before_trim_mm3': round(raw_volume, 1),
                  'components_before_trim': comps_raw,
                  'volume_mm3': round(trunk.volume(), 1),
                  'components': len(comps),
                  'triangles': int(trunk.num_tri()),
                  'stl_file': target.name,
                  'stl_roundtrip_rel_error': float('%.2e' % rt)},
        'trimmed_at_bone_mm3': removed,
        'bone_intersections_mm3': bone_hits,
        'bone_intersection_total_mm3': round(total_bone, 2),
        'protrusion': {'outside_mm3': round(outside.volume(), 1),
                       'percent_of_trunk': round(pct, 1)},
        'finished_envelope': {'volume_mm3': round(final.volume(), 1),
                              'anatomy_mm3': round(anatomical.volume(), 1),
                              'added_mm3': round(added, 1),
                              'components': len(fcomps),
                              'triangles': int(final.num_tri()),
                              'stl_file': final_target.name,
                              'channels': ch_rows,
                              'all_channels_fit': bool(all_fit)},
        'one_connected_y': bool(verdict_one),
        'bone_clear': bool(total_bone <= I.HARD_MM3),
        'not_verified': ['the three corridors are relaxed per head, so the trunk is a '
                         'network of tubes, not a lofted belly',
                         'no lattice, end caps, crimp, bolt or origin pads',
                         'the bore is not cut; the 1.6 mm wall is a construction claim',
                         'no muscle force, pennation, contact or printability'],
        'frozen': False,
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
