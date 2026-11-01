"""S1 TPU modules: turn a solid envelope into a printable lattice with a wire channel.

The envelopes from §5.4 are solid.  What gets printed is the paper's X lattice with a
continuous channel: struts fill the volume, one solid tube carries the wire, and the wire
bore itself is empty.

Construction, matching the elbow generator's cell (8 struts from the cell centre to the
eight corners, r = p/2 = 1 mm, spheres at the nodes):

    struts   every cell whose centre is inside the envelope, clipped to the envelope
    channel  solid tube of R1.6 along the searched corridor, minus the R0.85 wire bore
    part     (struts  intersect  envelope)  union  channel tube  minus  bore

There is no solid skin.  The elbow's own coupons do the same - the "skin" there is extra
surface strands, not a shell - and manifold3d has no offset operator to erode a shell with.

Run one module first, then all four:
    python mujoco_arm/tools/build_shoulder_tpu.py --module teres_major
    python mujoco_arm/tools/build_shoulder_tpu.py --json out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import attachment_regions
import build_shoulder_modules as BM
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP
import shoulder_pads

MUSCLE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / 'muscle'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'
MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'

#: paper cell (x, y, z) and strut diameter p, §6.6
CELL_MM = (15.0, 10.0, 20.0)
STRUT_R_MM = 1.0
#: elbow channel: outer R1.6, wire bore R0.85, wall 0.75 (§5.7)
CHANNEL_OUTER_R_MM = 1.6
BORE_R_MM = 0.85
#: facets used for the channel tubes.  A polygonal bore is narrower than its nominal
#: radius by cos(pi/n), so the check is made against the inscribed radius, not the
#: nominal one - otherwise a correctly built channel reads as "blocked".
CHANNEL_SEGMENTS = 48
INSCRIBED = float(np.cos(np.pi / CHANNEL_SEGMENTS))
#: the wire candidate is Ø0.6 (§5.7), so a facet loss of a few hundredths on the Ø1.7 bore
#: is irrelevant; the tolerance below just keeps the check honest about tessellation
FACET_TOL_MM = 0.02
WIRE_D_MM = 0.6

MODULES = [('teres_major', ('FJ1507',), False),
           ('subscapularis', ('FJ1504',), False),
           ('infraspinatus', ('FJ1500',), False),
           ('deltoid', ('FJ1468', 'FJ1513', 'FJ1467'), True)]


def mf():
    return I._manifold()


def rod(a, b, r, segments=16):
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


def tube(path, r, segments=16):
    pieces = []
    for a, b in zip(path, path[1:]):
        seg = rod(a, b, r, segments)
        if seg is not None:
            pieces.append(seg)
    for p in path:
        # the joint sphere must be as fine as the rod, or the sphere's facets are what
        # the bore clearance measures - they are coarser for the same segment count
        pieces.append(mf().Manifold.sphere(r, segments).translate(list(np.asarray(p, float))))
    return mf().Manifold.batch_boolean(pieces, mf().OpType.Add)


def lattice(envelope, tri, cell=CELL_MM, r=STRUT_R_MM, max_cells=4000, inset_mm=0.0):
    """8-strut X cells whose centres are inside the envelope, clipped to it.

    `inset_mm` is how far inside a cell centre must be to count.  Requiring a full strut
    radius (1 mm) rejects almost every cell on a thin curved belly and leaves disconnected
    clumps (§5.4.15), so the default is 0: any centre inside counts, and the strut union is
    intersected with the envelope afterwards, which trims the straddling cells for us.
    """
    bb = envelope.bounding_box()
    lo = np.array(bb[:3]) - 1.0
    hi = np.array(bb[3:]) + 1.0
    sx, sy, sz = cell
    nx = int(np.ceil((hi[0] - lo[0]) / sx))
    ny = int(np.ceil((hi[1] - lo[1]) / sy))
    nz = int(np.ceil((hi[2] - lo[2]) / sz))
    pieces, cells, skipped = [], 0, 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                c = lo + np.array([(i + 0.5) * sx, (j + 0.5) * sy, (k + 0.5) * sz])
                if mesh_geom.signed_distance(c, tri) >= -inset_mm:
                    skipped += 1
                    continue
                if cells >= max_cells:
                    skipped += 1
                    continue
                cells += 1
                hx, hy, hz = sx / 2.0, sy / 2.0, sz / 2.0
                for dx in (-1, 1):
                    for dy in (-1, 1):
                        for dz in (-1, 1):
                            corner = c + np.array([dx * hx, dy * hy, dz * hz])
                            seg = rod(c, corner, r)
                            if seg is not None:
                                pieces.append(seg)
                pieces.append(mf().Manifold.sphere(r, 12).translate(list(c)))
    if not pieces:
        return None, cells, skipped
    struts = mf().Manifold.batch_boolean(pieces, mf().OpType.Add)
    return struts, cells, skipped


def swept_lattice(envelope, tri, path, cell=CELL_MM, r=STRUT_R_MM, max_cells=4000,
                  max_extent_mm=45.0):
    """Build the X lattice in a frame that follows the corridor.

    A world-axis grid does not follow a thin curved belly: cells straddle the surface or
    drop out, and the survivors stay unconnected (§5.4.15).  Here the grid is laid out in
    the plane perpendicular to the corridor at each station, with the paper's (x, y)
    spacing in that plane and z between stations, so the cells step *along* the belly and
    share their corner nodes.
    """
    import muscle_sections as MS
    pts = [np.asarray(p, float) for p in path]
    # stations spaced z apart along the corridor
    seg = [np.linalg.norm(b - a) for a, b in zip(pts, pts[1:])]
    total = float(sum(seg))
    n_st = max(2, int(np.ceil(total / cell[2])) + 1)
    stations = []
    for s in range(n_st):
        f = s / (n_st - 1.0)
        # position along the polyline
        acc, idx = 0.0, 0
        target = f * total
        while idx < len(seg) and acc + seg[idx] < target:
            acc += seg[idx]
            idx += 1
        idx = min(idx, len(seg) - 1)
        w = 0.0 if seg[idx] < 1e-9 else (target - acc) / seg[idx]
        p = pts[idx] + (pts[idx + 1] - pts[idx]) * np.clip(w, 0, 1)
        d = pts[min(idx + 1, len(pts) - 1)] - pts[idx]
        stations.append((p, frames_unit(d)))

    sx, sy, sz = cell
    pieces, cells, skipped = [], 0, 0

    def extent(p, axis, sign):
        """How far the envelope reaches from p along axis (bisection)."""
        lo, hi = 1.0, max_extent_mm
        if mesh_geom.signed_distance(p + axis * sign * lo, tri) >= 0:
            return 0.0
        for _ in range(9):
            mid = 0.5 * (lo + hi)
            if mesh_geom.signed_distance(p + axis * sign * mid, tri) < 0:
                lo = mid
            else:
                hi = mid
        return lo

    for p, (e1, e2, t) in stations:
        # the grid must be laid out from the section's in-centre, not from the corridor
        # point: the corridor can sit near one wall, and measuring extents from there
        # under-reports the section and rejects stations that would fit a cell
        best = (mesh_geom.signed_distance(p, tri), np.zeros(2))
        span = 24.0
        for iu in range(-4, 5):
            for iv in range(-4, 5):
                off = np.array([iu * span / 4.0, iv * span / 4.0])
                d = mesh_geom.signed_distance(p + e1 * off[0] + e2 * off[1], tri)
                if d < best[0]:
                    best = (d, off)
        base = p + e1 * best[1][0] + e2 * best[1][1]
        ex1, ex2 = extent(base, e1, +1), extent(base, e1, -1)
        ey1, ey2 = extent(base, e2, +1), extent(base, e2, -1)
        if ex1 + ex2 < sx or ey1 + ey2 < sy:
            skipped += 1
            continue
        ni = int(np.floor((ex1 + ex2) / sx))
        nj = int(np.floor((ey1 + ey2) / sy))
        for i in range(ni):
            for j in range(nj):
                du = -ex1 + (i + 0.5) * sx
                dv = -ey1 + (j + 0.5) * sy
                c = base + e1 * du + e2 * dv
                if mesh_geom.signed_distance(c, tri) >= 0:
                    skipped += 1
                    continue
                if cells >= max_cells:
                    skipped += 1
                    continue
                cells += 1
                for su in (-1, 1):
                    for sv in (-1, 1):
                        for st in (-1, 1):
                            corner = (c + e1 * (su * sx / 2) + e2 * (sv * sy / 2)
                                      + t * (st * sz / 2))
                            seg_ = rod(c, corner, r)
                            if seg_ is not None:
                                pieces.append(seg_)
                pieces.append(mf().Manifold.sphere(r, 12).translate(list(c)))
    if not pieces:
        return None, cells, skipped
    struts = mf().Manifold.batch_boolean(pieces, mf().OpType.Add)
    return struts, cells, skipped


def frames_unit(d):
    t = np.asarray(d, float)
    n = np.linalg.norm(t)
    if n < 1e-12:
        t = np.array([0.0, 0.0, 1.0])
    else:
        t = t / n
    helper = np.array([0.0, 0.0, 1.0]) if abs(t[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(helper, t)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(t, e1)
    return e1, e2, t


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--module', default='teres_major')
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_tpu'))
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--mode', choices=('world', 'swept'), default='swept',
                    help='lattice grid: world-axis (the 5.4.15 failure) or swept along the '
                         'corridor (the 5.4.15 fix)')
    ap.add_argument('--cell', type=float, nargs=3, metavar=('X', 'Y', 'Z'),
                    default=list(CELL_MM),
                    help='lattice cell in mm; default is the paper value 15/10/20')
    args = ap.parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = MODULES if args.all else [m for m in MODULES if m[0] == args.module]
    if not wanted:
        raise SystemExit('unknown module %s' % args.module)

    bones = {k: SP.load_solid(BONE_DIR / v) for k, v in SP.BONES.items()}
    rows = []
    for key, heads, merged in wanted:
        env_path = MODULES_DIR / ('%s.stl' % key)
        if not env_path.is_file():
            raise SystemExit('missing envelope %s; run build_shoulder_modules.py first'
                             % env_path)
        V, F = mesh_geom.load_obj(str(MODULES_DIR / ('%s.stl' % key))) \
            if env_path.suffix == '.obj' else mesh_geom.read_stl(str(env_path))
        W, G = mesh_geom.weld(V, F)
        tri = W[G]
        envelope = mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))
        print('\n=== %s: envelope %.1f mm3 ===' % (key, envelope.volume()))

        # channel: solid tube minus wire bore, along the searched corridor
        if merged:
            ws = np.array([attachment_regions.centroid(
                attachment_regions.load_region(m, 'insertion'))[1] for m in heads])
            seat = (np.array([attachment_regions.centroid(
                attachment_regions.load_region(m, 'insertion'))[0]
                for m in heads]) * ws[:, None]).sum(0) / ws.sum()
        else:
            seat = attachment_regions.centroid(
                attachment_regions.load_region(heads[0], 'insertion'))[0]
        bone_tris = []
        for v in SP.BONES.values():
            bv, bf = mesh_geom.load_obj(str(BONE_DIR / v))
            bw, bg = mesh_geom.weld(bv, bf)
            bone_tris.append(bw[bg])
        channels = []
        for mid in heads:
            o, _ = attachment_regions.centroid(attachment_regions.load_region(mid, 'origin'))
            ts = np.linspace(0.0, 1.0, SP.INTERIOR_POINTS)
            start = [o + (seat - o) * t for t in ts]
            relaxed, _ = SP.relax(start, tri, iterations=SP.ITERATIONS, step=SP.STEP_MM,
                                  lam=SP.SMOOTHING, bone_tris=bone_tris,
                                  bone_clear_mm=CHANNEL_OUTER_R_MM * 2.0)
            relaxed[-1] = seat
            channels.append((mid, relaxed))

        if args.mode == 'swept':
            # follow the longest corridor; for the deltoid this is the dominant branch
            main_path = max(channels, key=lambda c: float(sum(np.linalg.norm(
                np.asarray(b) - np.asarray(a)) for a, b in zip(c[1], c[1][1:]))))[1]
            struts, cells, skipped = swept_lattice(envelope, tri, main_path,
                                                   cell=tuple(args.cell))
        else:
            struts, cells, skipped = lattice(envelope, tri, cell=tuple(args.cell))
        if struts is None:
            print('  no cells landed inside the envelope - cell too coarse')
            continue
        clipped = struts ^ envelope
        print('  lattice (%s): %d cells (%d skipped), raw %.1f mm3, clipped %.1f mm3'
              % (args.mode, cells, skipped, struts.volume(), clipped.volume()))
        wall = None
        for mid, path in channels:
            t = tube(path, CHANNEL_OUTER_R_MM, segments=CHANNEL_SEGMENTS)
            wall = t if wall is None else wall + t
        bore = None
        for mid, path in channels:
            b = tube(path, BORE_R_MM, segments=CHANNEL_SEGMENTS)
            bore = b if bore is None else bore + b
        solid = (clipped + wall) - bore
        comps = [c for c in solid.decompose() if c.volume() > 0.5]
        print('  channel: %d corridor(s), wall tube %.1f mm3' % (len(channels), wall.volume()))
        print('  part: %.1f mm3, %d piece(s), material %.1f%% of the solid envelope'
              % (solid.volume(), len(comps),
                 100.0 * solid.volume() / envelope.volume()))

        # the bore must be free along every corridor
        sv, sf = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh())
        SW, SG = mesh_geom.weld(sv, sf)
        stri = SW[SG]
        # the bore may only be judged between the pads: the corridor starts and ends on the
        # bone-mounted seats, where the end cap carries the wire (§5.4.14)
        bore_rows = []
        for mid, path in channels:
            pts = []
            for a, b in zip(path, path[1:]):
                a = np.asarray(a, float)
                b = np.asarray(b, float)
                n = max(2, int(np.ceil(np.linalg.norm(b - a) / 1.0)) + 1)
                for f in np.linspace(0, 1, n):
                    pts.append(a + (b - a) * f)
            ext = BM.channel_extent([np.round(p, 4).tolist() for p in path], stri,
                                    wall_mm=CHANNEL_OUTER_R_MM)
            g0 = ext['pad_to_channel_origin_mm']
            g1 = ext['pad_to_channel_distal_mm']
            # None means "no stretch holds at all", which is NOT the same as 0.0 ("holds
            # right from the pad") - conflating them made a floating channel read as perfect
            has_run = g0 is not None and g1 is not None
            g0 = g0 or 0.0
            g1 = g1 or 0.0
            total = float(sum(np.linalg.norm(np.asarray(b) - np.asarray(a))
                              for a, b in zip(path, path[1:])))
            inner = [p for p in pts
                     if np.linalg.norm(p - np.asarray(path[0])) >= g0 - 1e-6
                     and np.linalg.norm(p - np.asarray(path[-1])) >= g1 - 1e-6] if has_run \
                else pts
            if not inner:
                inner = pts
            d = [mesh_geom.signed_distance(p, stri) for p in inner]
            clear = -min(d)          # positive = the bore point is clear of material
            effective = BORE_R_MM * INSCRIBED
            ok = clear >= BORE_R_MM - FACET_TOL_MM
            bore_rows.append({'head': mid, 'min_bore_clearance_mm': round(clear, 3),
                              'bore_r_mm': BORE_R_MM,
                              'inscribed_bore_r_mm': round(effective, 3),
                              'facet_loss_mm': round(BORE_R_MM - clear, 4),
                              'judged_between_pads_mm': [round(g0, 2), round(g1, 2)],
                              'channel_enclosed': bool(has_run),
                              'judged_length_mm': round(total - g0 - g1, 2),
                              'facets': CHANNEL_SEGMENTS,
                              'wire_d_mm': WIRE_D_MM,
                              'free': bool(ok)})
            print('    %-8s bore clearance %.3f mm over %.1f mm between the pads '
                  '(gaps %.1f/%.1f) %s'
                  % (mid, clear, total - g0 - g1, g0, g1, 'FREE' if ok else 'BLOCKED'))

        target = out_dir / ('%s_tpu.stl' % key)
        mesh_geom.write_stl(str(target), sv, sf, header='S1 TPU %s' % key)
        rows.append({'module': key, 'envelope_mm3': round(envelope.volume(), 1),
                     'cells': cells, 'cells_skipped': skipped,
                     'lattice_clipped_mm3': round(clipped.volume(), 1),
                     'wall_tube_mm3': round(wall.volume(), 1),
                     'part_mm3': round(solid.volume(), 1),
                     'components': len(comps),
                     'material_percent': round(100.0 * solid.volume()
                                               / envelope.volume(), 2),
                     'channels': bore_rows,
                     'bore_all_free': all(r['free'] for r in bore_rows),
                     'stl_file': target.name})

    if args.json_out:
        I.write_report(args.json_out, {'scope': 'X lattice (paper cell 15/10/20, p=2) plus '
                                       'continuous R1.6 channel with R0.85 bore; no end caps, '
                                       'no pads, no surface strands',
                                       'lattice_mode': args.mode,
                                       'cell_mm': list(args.cell), 'strut_r_mm': STRUT_R_MM,
                                       'channel_outer_r_mm': CHANNEL_OUTER_R_MM,
                                       'bore_r_mm': BORE_R_MM,
                                       'modules': rows})
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
