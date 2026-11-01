"""End-cap tendon segments: the rigid piece that spans pad-to-channel-entry.

§5.4.14 measured how far each TPU channel stops short of its bone-side seat, and §5.4.20 kept
that list current.  The distances run 0.0 to 32.0 mm, while the pad + transition ring + cap
stack is only about 13 mm - so the long ends need a separate rigid tendon segment between the
pad and the TPU, or the wire has nothing to run in.

This builds one segment per corridor end:
    outer R2.5, wire bore R0.85 (wall 1.65), running along the corridor from the seat inward by
    the measured pad-to-channel distance.

and checks the two things that could make it impossible: bone interference, and interference
with the other modules (the corridor was only kept 1.6 mm clear of neighbours, which is less
than this radius, so an overlap here would be expected rather than surprising).

Run:
    python mujoco_arm/tools/build_shoulder_tendon_segments.py
    python mujoco_arm/tools/build_shoulder_tendon_segments.py --json out.json
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

SEGMENT_R_MM = 2.5
BORE_R_MM = 0.85
SEGMENTS = 32
#: a segment shorter than this is just the pad's own boss, not a separate piece
MIN_SEGMENT_MM = 1.0


def mf():
    return I._manifold()


def rod(a, b, r, segments=SEGMENTS):
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


def build_along(path, length_mm, from_distal, bone_tris, r=SEGMENT_R_MM,
                pad_thickness_mm=4.0):
    """Tube of R around the last/first `length_mm` of the corridor, minus the wire bore.

    The corridor's end point IS the bone-side seat, so a tube of radius R starting there is
    half buried in bone - the first version of this collided by 31 to 187 mm3 on all nine
    segments.  Physically the wire leaves through the pad, so the segment starts on the pad's
    OUTER face (seat + normal * pad thickness) and the buried part of the path is dropped.
    """
    pts = [np.asarray(p, float) for p in path]
    seat = pts[-1] if from_distal else pts[0]
    picked = None
    for tri in bone_tris:
        d, n = mesh_geom.signed_distance_and_normal(seat, tri)
        if picked is None or abs(d) < abs(picked[0]):
            picked = (d, n)
    n_out = picked[1] if picked is not None else np.array([0.0, 0.0, 1.0])
    start = seat + n_out * pad_thickness_mm
    if from_distal:
        pts = pts[::-1]
    keep, acc = [pts[0]], 0.0
    for a, b in zip(pts, pts[1:]):
        seg = float(np.linalg.norm(b - a))
        if acc + seg >= length_mm:
            keep.append(a + (b - a) * ((length_mm - acc) / seg))
            break
        keep.append(b)
        acc += seg
    if len(keep) < 2:
        keep.append(pts[1] if len(pts) > 1 else pts[0])
    # Skip the leading stretch that would bury the tube, then keep everything after it.
    # Filtering every point instead (the first version) dropped the whole short segment on
    # eight of nine ends, because the point just outside the pad is only ~2.4 mm from bone -
    # under the radius but perfectly usable as the second vertex.
    dists = [min(mesh_geom.signed_distance(p, tri) for tri in bone_tris) for p in keep]
    first = next((i for i, d in enumerate(dists) if d >= r), None)
    keep = [start] + (keep[first:] if first is not None else [])
    pieces = []
    for a, b in zip(keep, keep[1:]):
        s = rod(a, b, r)
        if s is not None:
            pieces.append(s)
    for p in keep:
        pieces.append(mf().Manifold.sphere(r, SEGMENTS).translate(list(p)))
    outer = mf().Manifold.batch_boolean(pieces, mf().OpType.Add)
    bore_pieces = []
    for a, b in zip(keep, keep[1:]):
        s = rod(a, b, BORE_R_MM)
        if s is not None:
            bore_pieces.append(s)
    for p in keep:
        bore_pieces.append(mf().Manifold.sphere(BORE_R_MM, SEGMENTS).translate(list(p)))
    bore = mf().Manifold.batch_boolean(bore_pieces, mf().OpType.Add)
    length = float(sum(np.linalg.norm(b - a) for a, b in zip(keep, keep[1:])))
    return outer - bore, length, keep


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_tendons'))
    ap.add_argument('--segment-r', type=float, default=SEGMENT_R_MM,
                    help='outer radius of the tendon segment [mm]; R2.5 does not fit the '
                         'short ends, see 5.4.21')
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
    others = {}
    for key in ('deltoid', 'teres_major', 'subscapularis', 'infraspinatus'):
        p = MODULES_DIR / ('%s.stl' % key)
        if p.is_file():
            V, F = mesh_geom.read_stl(str(p))
            W, G = mesh_geom.weld(V, F)
            others[key] = mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))

    rows = []
    print('%-14s %-8s %-7s %8s %9s %8s %8s'
          % ('module', 'head', 'end', 'length', 'volume', 'bone', 'neighbour'))
    for m in design['modules']:
        paths = {c['head']: c['path_source_mm'] for c in m['corridors']}
        for c in m['channels']:
            for end, key, dist in (('origin', 'pad_to_channel_origin_mm', False),
                                   ('distal', 'pad_to_channel_distal_mm', True)):
                L = c.get(key)
                if L is None or L < MIN_SEGMENT_MM:
                    rows.append({'module': m['module'], 'head': c['head'], 'end': end,
                                 'length_mm': 0.0, 'built': False,
                                 'note': 'no segment needed (channel reaches the pad)'})
                    print('%-14s %-8s %-7s %8.1f %9s %8s %8s'
                          % (m['module'], c['head'], end, L or 0.0, '-', '-', '-'))
                    continue
                solid, length, keep = build_along(paths[c['head']], L, from_distal=dist,
                                                  bone_tris=bone_tris, r=args.segment_r)
                bone_hit = max(I.intersection_volume(solid, b).__float__()
                               for b in bones.values())
                near = {k: round(I.intersection_volume(solid, s), 3)
                        for k, s in others.items() if k != m['module']}
                worst_near = max(near.values()) if near else 0.0
                rows.append({'module': m['module'], 'head': c['head'], 'end': end,
                             'length_mm': round(length, 2), 'built': True,
                             'volume_mm3': round(solid.volume(), 1),
                             'bone_intersection_mm3': round(bone_hit, 4),
                             'neighbour_intersections_mm3': near,
                             'worst_neighbour_mm3': round(worst_near, 3)})
                print('%-14s %-8s %-7s %8.1f %9.1f %8.3f %8.3f'
                      % (m['module'], c['head'], end, length, solid.volume(),
                         bone_hit, worst_near))

    need = [r for r in rows if r.get('built')]
    print('\nsegments needed: %d of %d ends' % (len(need), len(rows)))
    if need:
        print('longest segment: %.1f mm (%s %s %s)'
              % (max(r['length_mm'] for r in need),
                 *[r[k] for k in ('module', 'head', 'end')
                   for r in [max(need, key=lambda x: x['length_mm'])]]))
        print('max bone intersection: %.4f mm3'
              % max(r['bone_intersection_mm3'] for r in need))
        worst = max(need, key=lambda r: r['worst_neighbour_mm3'])
        print('worst neighbour overlap: %.3f mm3 (%s %s %s)'
              % (worst['worst_neighbour_mm3'], worst['module'], worst['head'], worst['end']))
        bad = [r for r in need if r['bone_intersection_mm3'] > I.HARD_MM3]
        print('segments colliding with bone: %d' % len(bad))
        badn = [r for r in need if r['worst_neighbour_mm3'] > 1.0]
        print('segments overlapping a neighbour by >1 mm3: %d' % len(badn))

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'rigid tendon segments from each bone-side seat to the TPU channel entry',
            'segment_outer_r_mm': args.segment_r, 'bore_r_mm': BORE_R_MM,
            'segments': rows,
            'n_needed': len(need),
            'max_length_mm': round(max([r['length_mm'] for r in need], default=0.0), 2),
            'not_verified': ['pads, transition rings and caps are not modelled here',
                             'no wire, no crimp, no preload',
                             'no 0-120 deg sweep'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
