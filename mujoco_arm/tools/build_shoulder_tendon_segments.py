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
#: Per-end radius overrides.  D21 (=A, 5.4.32) accepts the teres major's 32.0 mm cap and strengthens
#: just that end: R2.5 fits there because a 32 mm long segment has the span for a bigger section,
#: while the short ends need R1.8 or they cannot be built at all (5.4.22).
SEGMENT_R_OVERRIDES = {'FJ1507/distal': 2.5}
BORE_R_MM = 0.85
# The posterior deltoid distal bend otherwise leaves a zero-thickness material
# pinch after STL quantization. A measured local radial overcut removes it;
# the larger 1.704mm bore is recorded per part, not hidden by a weaker check.
BORE_R_OVERRIDES = {'FJ1513/distal': 0.852}
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


def open_ended_tube(points, r, bore_r=0.85):
    """Swept tube with bores extended past both rounded outer end caps."""
    pts=[np.asarray(p,float) for p in points]
    if not 0 < bore_r < r or len(pts)<2:
        raise ValueError('Need a positive wall and at least two path vertices')
    pts=[pts[0]]+[b for a,b in zip(pts,pts[1:]) if np.linalg.norm(b-a)>1e-6]
    if len(pts)<2:raise ValueError('Zero-length tube')
    pieces=[rod(a,b,r) for a,b in zip(pts,pts[1:])]
    pieces += [mf().Manifold.sphere(r,SEGMENTS).translate(list(p)) for p in pts]
    outer=mf().Manifold.batch_boolean(pieces,mf().OpType.Add)
    bore_pieces=[rod(a,b,bore_r) for a,b in zip(pts,pts[1:])]
    bore_pieces += [mf().Manifold.sphere(bore_r,SEGMENTS).translate(list(p)) for p in pts]
    first=(pts[1]-pts[0])/np.linalg.norm(pts[1]-pts[0])
    last=(pts[-1]-pts[-2])/np.linalg.norm(pts[-1]-pts[-2])
    # Concentric small end spheres alone leave an enclosed cavity in the
    # larger end spheres. Drill right through the rounded external envelope.
    bore_pieces += [rod(pts[0]-first*(r+1),pts[0],bore_r),
                    rod(pts[-1],pts[-1]+last*(r+1),bore_r)]
    bore=mf().Manifold.batch_boolean(bore_pieces,mf().OpType.Add)
    return outer-bore


def build_along(path, length_mm, from_distal, bone_tris, r=SEGMENT_R_MM,
                pad_thickness_mm=4.0, relax=False, bore_r=BORE_R_MM):
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
    if first is not None:
        keep = [start] + keep[first:]
    else:
        # No point in the capped stretch clears by the radius.  Dropping them all (the earlier
        # behaviour) left a single sphere - a "segment" with zero length that still overlaps bone
        # wherever it sits.  Build the capped path anyway and let the relaxation try: if it cannot
        # clear, the honest report is a real segment with residual contact, not a stub.
        keep = [start] + keep[1:]
    if relax:
        keep = relax_off_bone(keep, bone_tris, r + 0.3)
    solid = open_ended_tube(keep,r,bore_r)
    length = float(sum(np.linalg.norm(b - a) for a, b in zip(keep, keep[1:])))
    return solid, length, keep


def relax_off_bone(path, bone_tris, clear_mm, step_mm=2.0, iterations=40, lam=0.3):
    """Densify and push the interior of a short path off the bone.

    The grazes left in 5.4.29 are between two consecutive vertices: a straight rod of radius 1.8
    whose endpoints clear the bone can still cut a convex surface between them.  Densifying to
    2 mm and pushing every interior point to `clear_mm` is the same treatment 5.9.11 used for the
    fixture posts.
    """
    pts = [np.asarray(p, float) for p in path]
    dense = [pts[0]]
    for a, b in zip(pts, pts[1:]):
        n = max(1, int(np.ceil(np.linalg.norm(b - a) / step_mm)))
        for k in range(1, n + 1):
            dense.append(a + (b - a) * (k / float(n)))
    for _ in range(iterations):
        moved = list(dense)
        for i in range(1, len(dense) - 1):
            best = None
            for tri in bone_tris:
                d, nrm = mesh_geom.signed_distance_and_normal(dense[i], tri)
                if d < clear_mm and (best is None or d < best[0]):
                    best = (d, nrm)
            if best is not None:
                moved[i] = dense[i] + best[1] * (clear_mm - best[0]) * 0.8
        for i in range(1, len(moved) - 1):
            moved[i] = moved[i] + lam * (moved[i - 1] + moved[i + 1] - 2.0 * moved[i])
        dense = moved
        if all(min(mesh_geom.signed_distance(p, tri) for tri in bone_tris) >= clear_mm
               for p in dense[1:-1]):
            break
    return dense


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_tendons'))
    ap.add_argument('--standoff-mm', type=float, default=4.0,
                    help='how far the segment starts from the seat along the bone normal [mm]')
    ap.add_argument('--relax', action='store_true',
                    help='densify and push each segment off the bone (the 5.9.11 method)')
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
                r_end = SEGMENT_R_OVERRIDES.get('%s/%s' % (c['head'], end), args.segment_r)
                bore_end = BORE_R_OVERRIDES.get('%s/%s' % (c['head'], end), BORE_R_MM)
                solid, length, keep = build_along(paths[c['head']], L, from_distal=dist,
                                                  bone_tris=bone_tris, r=r_end,
                                                  relax=args.relax,
                                                  pad_thickness_mm=args.standoff_mm,bore_r=bore_end)
                # A calculated solid is not a manufacturing deliverable until its
                # actual mesh and centreline have been exported and read back.
                # Small parts must be tessellated near the origin. Float32 STL
                # at source Z~1300mm collapses close vertices. Store the exact
                # placement explicitly instead of rounding global geometry.
                translation = np.asarray(keep[0],float)
                solid = I.clean_solid(solid.translate(-translation))
                V, F = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh())
                filename = '%s_%s_tendon.stl' % (c['head'], end)
                removed_micro_shell_mm3 = 0.0
                for attempt in range(3):
                    V, F = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh())
                    mesh_geom.write_stl(str(out_dir / filename), V, F, header='S1 tendon ' + filename)
                    RV, RF = mesh_geom.read_stl(str(out_dir / filename))
                    RW, RG = mesh_geom.weld(RV, RF)
                    _, boundary, nonmanifold = mesh_geom.edge_census(RG)
                    restored = mf().Manifold(mf().Mesh64(RW, RG.astype(np.uint32)))
                    relative_error = abs(restored.volume()-solid.volume()) / solid.volume()
                    pieces = restored.decompose()
                    if not boundary and not nonmanifold and len(pieces)==1 and relative_error<=1e-5:
                        break
                    if attempt<2 and (boundary or nonmanifold):
                        # Repeat the same fixed-tolerance repair after float32
                        # roundtrip; do not relax the acceptance threshold.
                        solid=I.clean_solid(restored,eps=1e-3)
                        continue
                    # Float32 STL can leave sub-micron sliver shells. Remove
                    # only explicitly recorded minute shells, never real parts
                    # or a large sealed lumen (the original endpoint defect).
                    pieces = sorted(pieces,key=lambda p:abs(p.volume()),reverse=True)
                    micro = sum(abs(p.volume()) for p in pieces[1:])
                    if (attempt<2 and not boundary and not nonmanifold and len(pieces)>1
                            and micro<1e-3 and micro/max(1,restored.volume())<1e-5):
                        removed_micro_shell_mm3 += micro
                        solid=pieces[0]
                        continue
                    raise ValueError('Invalid exported tendon mesh: ' + filename)
                else:
                    raise ValueError('Tendon mesh repair did not converge: ' + filename)
                solid = restored.translate(translation)
                bone_hit = max(I.intersection_volume(solid, b).__float__()
                               for b in bones.values())
                near = {k: round(I.intersection_volume(solid, s), 3)
                        for k, s in others.items() if k != m['module']}
                worst_near = max(near.values()) if near else 0.0
                rows.append({'module': m['module'], 'head': c['head'], 'end': end,
                             'length_mm': round(length, 2), 'built': True,
                             'file': filename, 'path_source_mm': [p.tolist() for p in keep],
                             'coordinate_frame': 'local_mm',
                             'translation_source_mm': translation.tolist(),
                             'outer_radius_mm': r_end,
                             'bore_radius_mm': bore_end,
                             'boundary_edges': int(boundary), 'nonmanifold_edges': int(nonmanifold),
                             'material_components': len(restored.decompose()),
                             'roundtrip_relative_error': relative_error,
                             'removed_micro_shell_absolute_mm3': removed_micro_shell_mm3,
                             'stl_cleanup_passes': attempt+1,
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
            'standoff_mm': args.standoff_mm, 'relaxed': bool(args.relax),
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
