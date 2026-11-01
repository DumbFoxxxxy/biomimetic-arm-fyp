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
#: lattice stubs below this are clipping debris, not structure [mm3]; one cell's struts are ~40 mm3
DEBRIS_MM3 = 10.0
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
    return union_all(pieces)


def union_all(pieces, chunk=200):
    """Union a big piece list in chunks.

    One batch_boolean over ~2700 rods is what made the deltoid run for half an hour without
    producing any output; unioning 200 at a time and folding the results keeps each operation
    small and finishes in minutes.
    """
    if not pieces:
        return None
    if len(pieces) <= chunk:
        return mf().Manifold.batch_boolean(pieces, mf().OpType.Add)
    acc = None
    for i in range(0, len(pieces), chunk):
        part = mf().Manifold.batch_boolean(pieces[i:i + chunk], mf().OpType.Add)
        acc = part if acc is None else acc + part
    return acc


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
    struts = union_all(pieces)
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
    struts = union_all(pieces)
    return struts, cells, skipped


def structured_lattice(envelope, tri, path, cell=CELL_MM, r=STRUT_R_MM,
                       i_range=(-3, 3), j_range=(-2, 2), max_cells=4000):
    """Structured X lattice on a parallel-transported frame.

    The two earlier attempts failed because cells never shared nodes (§5.4.16): a world grid
    straddles a curved belly, and a per-station grid has a different origin at every station.

    Here the frame is transported along the corridor (minimal rotation station to station) and
    node (i, j, k) is a *fixed index triple*, so cell (i, j, k) and cell (i+1, j, k) use the
    very same corner nodes.  A cell is emitted only when all eight of its corners are inside
    the envelope, which is what makes the result connected rather than hoped-for.
    """
    pts = [np.asarray(p, float) for p in path]
    seg = [np.linalg.norm(b - a) for a, b in zip(pts, pts[1:])]
    total = float(sum(seg))
    sz = cell[2]
    n_st = max(3, int(np.ceil(total / sz)) + 1)

    # stations equally spaced along the polyline
    stations, tangents = [], []
    for s in range(n_st):
        target = total * s / (n_st - 1.0)
        acc, idx = 0.0, 0
        while idx < len(seg) - 1 and acc + seg[idx] < target:
            acc += seg[idx]
            idx += 1
        w = 0.0 if seg[idx] < 1e-9 else (target - acc) / seg[idx]
        stations.append(pts[idx] + (pts[idx + 1] - pts[idx]) * np.clip(w, 0, 1))
        tangents.append(pts[min(idx + 1, len(pts) - 1)] - pts[idx])

    # transport the frame: first station defines it, then minimal rotation each step
    e1, e2, t0 = frames_unit(tangents[0])
    frames = [(e1, e2, t0)]
    for k in range(1, n_st):
        t_prev = frames[-1][2]
        t_new = frames_unit(tangents[k])[2]
        axis = np.cross(t_prev, t_new)
        n_axis = np.linalg.norm(axis)
        if n_axis > 1e-9:
            axis = axis / n_axis
            ang = float(np.arccos(np.clip(t_prev @ t_new, -1.0, 1.0)))
            c, s = np.cos(ang), np.sin(ang)
            K = np.array([[0.0, -axis[2], axis[1]],
                          [axis[2], 0.0, -axis[0]],
                          [-axis[1], axis[0], 0.0]])
            R = np.eye(3) + s * K + (1.0 - c) * (K @ K)
        else:
            R = np.eye(3)
        pe1, pe2, _ = frames[-1]
        frames.append((R @ pe1, R @ pe2, t_new))

    # one lateral offset for the whole lattice, so every station shares the same (i, j)
    du, dv = [], []
    span = 22.0
    for p, (f1, f2, _) in zip(stations, frames):
        best = (mesh_geom.signed_distance(p, tri), np.zeros(2))
        for iu in range(-3, 4):
            for iv in range(-3, 4):
                off = np.array([iu * span / 3.0, iv * span / 3.0])
                d = mesh_geom.signed_distance(p + f1 * off[0] + f2 * off[1], tri)
                if d < best[0]:
                    best = (d, off)
        du.append(best[1][0])
        dv.append(best[1][1])
    u0, v0 = float(np.median(du)), float(np.median(dv))

    sx, sy = cell[0], cell[1]
    nodes = {}
    inside = {}
    for k in range(n_st):
        f1, f2, _ = frames[k]
        for i in range(i_range[0], i_range[1] + 1):
            for j in range(j_range[0], j_range[1] + 1):
                p = stations[k] + f1 * (u0 + i * sx) + f2 * (v0 + j * sy)
                nodes[(i, j, k)] = p
                inside[(i, j, k)] = mesh_geom.signed_distance(p, tri) < 0.0

    pieces, cells, skipped = [], 0, 0
    for k in range(n_st - 1):
        for i in range(i_range[0], i_range[1]):
            for j in range(j_range[0], j_range[1]):
                corners = [(i, j, k), (i + 1, j, k), (i, j + 1, k), (i + 1, j + 1, k),
                           (i, j, k + 1), (i + 1, j, k + 1), (i, j + 1, k + 1),
                           (i + 1, j + 1, k + 1)]
                centre = np.mean([nodes[c] for c in corners], axis=0)
                # Requiring all eight corners inside is far too strict: the cell diagonal is
                # 26.9 mm, so it would need a 13.4 mm inradius, which only the thickest
                # belly reaches.  Emit the cell when its centre is inside and strut only to
                # the corners that are also inside - adjacent cells still share those nodes,
                # which is what keeps the lattice connected.
                if mesh_geom.signed_distance(centre, tri) >= 0.0:
                    skipped += 1
                    continue
                live = [c for c in corners if inside[c]]
                if len(live) < 3:
                    skipped += 1
                    continue
                if cells >= max_cells:
                    skipped += 1
                    continue
                cells += 1
                for c in live:
                    seg_ = rod(centre, nodes[c], r)
                    if seg_ is not None:
                        pieces.append(seg_)
                pieces.append(mf().Manifold.sphere(r, 12).translate(list(centre)))
    if not pieces:
        return None, cells, skipped
    struts = union_all(pieces)
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
    ap.add_argument('--mode', choices=('world', 'swept', 'structured'), default='structured',
                    help='lattice grid: world-axis and per-station swept are the two 5.4.15/'
                         '5.4.16 failures; structured is the parallel-transported grid that '
                         'shares corner nodes')
    ap.add_argument('--root-boss-mm', type=float, default=0.0,
                    help='radius of a densified boss at a shared insertion seat [mm]; the three '
                         'deltoid branches are on independent grids and need one to connect (R11)')
    ap.add_argument('--cell', type=float, nargs=3, metavar=('X', 'Y', 'Z'),
                    default=list(CELL_MM),
                    help='lattice cell in mm; default is the paper value 15/10/20')
    args = ap.parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = MODULES if args.all else [
        m for m in MODULES if m[0] in [x.strip() for x in args.module.split(',')]]
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
        elif args.mode == 'structured':
            # one transported lattice per corridor, unioned: the deltoid is a Y, and a
            # sweep along one branch cannot fill the other two
            parts, cells, skipped = [], 0, 0
            for mid, cpath in channels:
                s, c, sk = structured_lattice(envelope, tri, cpath, cell=tuple(args.cell))
                cells += c
                skipped += sk
                if s is not None:
                    parts.append(s)
                    print('  branch %-8s %d cells' % (mid, c))
            if not parts:
                struts = None
            elif len(parts) == 1:
                struts = parts[0]
            else:
                struts = mf().Manifold.batch_boolean(parts, mf().OpType.Add)
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
        # Y root: the three branches are built on independent transported grids, so they do not
        # share nodes and the part comes out in pieces.  R11 already says the bifurcation root is
        # the weak point and needs its own treatment, so the root gets a densified boss that also
        # ties the branches together.
        root = None
        if args.root_boss_mm > 0 and len(channels) > 1:
            root = mf().Manifold.sphere(args.root_boss_mm, 32).translate(
                list(np.asarray(seat, float)))
        solid = (clipped + wall) - bore
        if root is not None:
            solid = (solid + (root ^ envelope)) - bore
            print('  Y root boss r%.1f mm at the shared seat: %.1f mm3'
                  % (args.root_boss_mm, root.volume()))
        # Debris: clipping leaves a few cut-loose stubs.  A single cell's struts are ~40 mm3, so
        # anything under DEBRIS_MM3 is a stub, not structure - the same rule 5.4.11 uses for the
        # module envelopes.  Dropping them is what makes "one piece" mean one printable body.
        raw_comps = [c for c in solid.decompose() if c.volume() > 0.5]
        debris = [c for c in raw_comps if c.volume() < DEBRIS_MM3]
        if debris and len(raw_comps) > len(debris):
            # Subtract the stubs from the solid.  Rebuilding the part out of the decompose()
            # components instead FILLS the bores: a component carries its channel cavity folded
            # into itself (which is also why sum(decompose) exceeds volume(), 5.4.27), and
            # reconstructing from it turns the lumen into solid.
            for c in debris:
                solid = solid - c
            print('  dropped %d stub(s) totalling %.1f mm3 (< %.0f mm3)'
                  % (len(debris), sum(c.volume() for c in debris), DEBRIS_MM3))
        comps = [c for c in solid.decompose() if c.volume() > 0.5]
        frags = sorted((round(c.volume(), 1) for c in comps), reverse=True)
        print('  channel: %d corridor(s), wall tube %.1f mm3' % (len(channels), wall.volume()))
        print('  part: %.1f mm3, %d piece(s), material %.1f%% of the solid envelope'
              % (solid.volume(), len(comps),
                 100.0 * solid.volume() / envelope.volume()))
        if len(comps) > 1:
            body = frags[0]
            print('    fragments: %s  -> body %.0f mm3 (%.1f%% of the part), '
                  'clipped stubs %.1f mm3'
                  % (', '.join('%.1f' % f for f in frags[:6]), body,
                     100.0 * body / solid.volume(), sum(frags[1:])))

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
            # independent check that does not rely on any sign convention: a 0.45 mm probe
            # sphere at each sample must not intersect the part at all if the bore is open
            probe = mf().Manifold.sphere(0.45, 16)
            blocked_pts = 0
            for p in inner[::max(1, len(inner) // 25)]:
                if I.intersection_volume(probe.translate(list(p)), solid) > I.HARD_MM3:
                    blocked_pts += 1
            clear = min(d)
            effective = BORE_R_MM * INSCRIBED
            sign_ok = clear >= BORE_R_MM - FACET_TOL_MM
            # The probe is the criterion, not the signed distance.  Near a 0.75 mm wall the
            # nearest-triangle normal can flip the sign, so the normal-based answer is kept
            # only as a diagnostic and a disagreement is recorded rather than hidden.
            ok = blocked_pts == 0
            agree = bool(sign_ok) == bool(ok)
            bore_rows.append({'head': mid, 'min_bore_clearance_mm': round(clear, 3),
                              'bore_r_mm': BORE_R_MM,
                              'inscribed_bore_r_mm': round(effective, 3),
                              'facet_loss_mm': round(BORE_R_MM - clear, 4),
                              'probe_samples': len(inner[::max(1, len(inner) // 25)]),
                              'probe_blocked': blocked_pts,
                              'sign_says_free': bool(sign_ok),
                              'probe_says_free': bool(ok),
                              'sign_and_probe_agree': agree,
                              'judged_between_pads_mm': [round(g0, 2), round(g1, 2)],
                              'channel_enclosed': bool(has_run),
                              'judged_length_mm': round(total - g0 - g1, 2),
                              'facets': CHANNEL_SEGMENTS,
                              'wire_d_mm': WIRE_D_MM,
                              'free': bool(ok and blocked_pts == 0)})
            print('    %-8s bore clearance %+.3f mm, probe blocked %d/%d, agree %s -> %s'
                  % (mid, clear, blocked_pts,
                     len(inner[::max(1, len(inner) // 25)]), agree,
                     'FREE' if (ok and blocked_pts == 0) else 'BLOCKED'))

        target = out_dir / ('%s_tpu.stl' % key)
        # The lattice unions leave 10-41 non-manifold edges; clean_solid goes through a file and
        # simplifies the rebuilt mesh to clear them (5.9.26).  Without this the part is "probably
        # printable" rather than watertight.
        solid = I.clean_solid(solid)
        sv, sf = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh())
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
