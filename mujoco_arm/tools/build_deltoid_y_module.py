"""S1 deltoid module: can the three heads become one printable Y?

D8 says the three deltoid heads must be rebuilt as one Y-shaped belly - one shared
distal seat on the humerus, three origins on the clavicle and scapula - rather than a
boolean union of three anatomical meshes.  Round 3 showed the sliver joint that makes
D8's argument: the heads touch through only 24.1 / 113.1 mm3 of overlap.

This script builds the envelope that the Y has to fit inside, and answers four things:

1. after trimming each head by the bones (§5.4.1), is the union still one piece;
2. how much material is there, and does it survive trimming;
3. can three channels run from the three origin pads, through their own head, into the
   one shared distal pad - i.e. is there a continuous channel inside the union;
4. what does the module still collide with (the other three modules, and the bones).

The channel is modelled as a three-point polyline per head (origin pad -> the head's
own mid-section centroid -> shared distal pad).  That is a routing feasibility test,
not a fibre direction and not a tendon path.

Run:
    python mujoco_arm/tools/build_deltoid_y_module.py
    python mujoco_arm/tools/build_deltoid_y_module.py --out <dir>
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
import muscle_sections
import shoulder_pads

MUSCLE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / 'muscle'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

DELTOID = {'FJ1468': '三角肌前束', 'FJ1513': '三角肌后束', 'FJ1467': '三角肌中束'}
OTHER_MODULES = {'FJ1507': '大圆肌', 'FJ1504': '肩胛下肌', 'FJ1500': '冈下肌'}
BONES = {'scapula': '肩胛骨.obj', 'clavicle': '锁骨.obj', 'humerus': '肱骨.obj'}

#: elbow channel bore (R0.85) and wall (R1.6), §5.7
CHANNEL_INNER_R_MM = 0.85
CHANNEL_OUTER_R_MM = 1.6
#: the three bores sit on this radius inside the shared R11 distal pad [mm]
BORE_RING_R_MM = 4.0
CHANNEL_SAMPLE_MM = 1.5


def mf():
    return I._manifold()


def load_solid(path):
    V, F = mesh_geom.load_obj(str(path))
    W, G = mesh_geom.weld(V, F)
    return mf().Manifold(mf().Mesh64(W, G.astype(np.uint32))), W, G


def polyline_wall(points, solid, tri, step=CHANNEL_SAMPLE_MM):
    """Thin wrapper kept for the printed table; the check lives in muscle_sections."""
    return muscle_sections.path_clearance(points, tri, step=step)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_deltoid_y'))
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    bones = {}
    for name, fname in BONES.items():
        solid, W, G = load_solid(BONE_DIR / fname)
        bones[name] = (solid, W, G)

    # ---- origin pads and the shared distal seat
    origins, heads = {}, {}
    for mid, label in DELTOID.items():
        region = attachment_regions.load_region(mid, 'origin')
        o, _ = attachment_regions.centroid(region)
        origins[mid] = o
        solid, W, G = load_solid(MUSCLE_DIR / (mid + '.obj'))
        heads[mid] = (solid, W, G)
    ins = []
    for mid in DELTOID:
        region = attachment_regions.load_region(mid, 'insertion')
        c, _ = attachment_regions.centroid(region)
        ins.append((c, region['area_mm2']))
    pts = np.array([c for c, _ in ins])
    w = np.array([a for _, a in ins])
    distal = (pts * w[:, None]).sum(0) / w.sum()
    # the shared seat's normal: area-weighted mean of the three insertion normals
    hW, hG, htri = None, None, None
    hum_solid, hW, hG = bones['humerus']
    htri = hW[hG]
    normals = []
    hum_centre = hW.mean(axis=0)
    for (c, _), mid in zip(ins, DELTOID):
        region = attachment_regions.load_region(mid, 'insertion', 'humerus')
        idx = int(np.argmax(region['weights']))
        face = htri[region['tris'][idx]]
        n = np.cross(face[1] - face[0], face[2] - face[0])
        n = n / np.linalg.norm(n)
        if n @ (c - hum_centre) < 0:
            n = -n
        normals.append(n)
    seat_normal = frames.unit((np.array(normals) * w[:, None]).sum(0))
    seat = shoulder_pads.build_pad(distal, seat_normal)
    print('shared distal seat at %s, normal %s' % (np.round(distal, 1),
                                                   np.round(seat_normal, 3)))

    # ---- trim each head by the bones, then union
    trimmed = {}
    print('\n--- per head: trim by bone ---')
    print('%-8s %-10s %10s %10s %10s %6s' % ('head', 'name', 'anatomy', 'trimmed',
                                             'removed', 'comps'))
    for mid, label in DELTOID.items():
        solid, _, _ = heads[mid]
        cut = solid
        for name in BONES:
            cut = cut - bones[name][0]
        comps = len([c for c in cut.decompose() if c.volume() > 0.5])
        trimmed[mid] = cut
        print('%-8s %-10s %10.1f %10.1f %10.1f %6d'
              % (mid, label, solid.volume(), cut.volume(),
                 solid.volume() - cut.volume(), comps))

    merged = trimmed['FJ1468']
    for mid in ('FJ1513', 'FJ1467'):
        merged = merged + trimmed[mid]
    merged = merged + seat
    comps = [c for c in merged.decompose() if c.volume() > 0.5]
    print('\n--- merged envelope ---')
    print('  volume %.1f mm3, components %d, triangles %d'
          % (merged.volume(), len(comps), merged.num_tri()))
    if len(comps) != 1:
        print('  WARNING: the three trimmed heads do NOT form one connected solid;')
        print('           the Y has to be reconstructed, not unioned (D8).')

    print('\n--- residual interference after trimming ---')
    print('  bones: %s' % (', '.join(
        '%s %.2f mm3' % (n, I.intersection_volume(merged, bones[n][0]))
        for n in BONES)))
    other = {}
    for mid, label in OTHER_MODULES.items():
        solid, _, _ = load_solid(MUSCLE_DIR / (mid + '.obj'))
        v = I.intersection_volume(merged, solid)
        other[mid] = round(v, 2)
    print('  other modules: %s' % ', '.join('%s %.1f mm3' % kv for kv in other.items()))

    # ---- three channels into one shared seat
    print('\n--- channel routing: three origins -> one shared seat ---')
    print('%-8s %-9s %7s %8s %9s %9s %s'
          % ('head', 'route', 'samples', 'worst', 'deepest', 'inside', 'verdict'))
    merged_arrays = mesh_geom.manifold_mesh_to_arrays(merged.to_mesh())
    _, MG = mesh_geom.weld(*merged_arrays)
    mW = merged_arrays[0]
    tri = mW[MG]
    channels = []
    for mid, label in DELTOID.items():
        solid, _, _ = heads[mid]
        for route in ('chord3', 'centroid-curve', 'medial-axis'):
            if route == 'chord3':
                chunk = muscle_sections.section_at(solid, origins[mid], distal, 0.5)
                perp = np.array([0.0, 0.0, 0.0])
                if chunk and chunk.get('centroid_source_mm'):
                    perp = np.asarray(chunk['centroid_source_mm'], float)
                path = [origins[mid], perp, distal]
            else:
                # both curves are section-centre heuristics, kept side by side because
                # they fail differently: the area centroid can sit OUTSIDE a concave
                # section, and the inscribed-circle centre jumps between lobes.
                key = ('centroid_source_mm' if route == 'centroid-curve'
                       else 'inradius_source_mm')
                path = [origins[mid]]
                for t in np.arange(0.12, 0.93, 0.10):
                    chunk = muscle_sections.section_at(solid, origins[mid], distal,
                                                       float(t))
                    if chunk and chunk.get(key):
                        path.append(np.asarray(chunk[key], float))
                path.append(distal)
            wall = polyline_wall(path, merged, tri)
            bore_ok = wall['fits_channel']
            channels.append({'head': mid, 'name': label, 'route': route,
                             'path_points': len(path),
                             'path_source_mm': [np.round(p, 2).tolist() for p in path],
                             'path_length_mm': round(float(sum(np.linalg.norm(
                                 np.asarray(b, float) - np.asarray(a, float))
                                 for a, b in zip(path, path[1:]))), 2),
                             'worst_mm': wall['worst_mm'],
                             'deepest_mm': wall['deepest_mm'],
                             'inside_fraction': wall['inside_fraction'],
                             'channel_outer_r_mm': CHANNEL_OUTER_R_MM,
                             'fits_channel': bool(bore_ok),
                             'profile_mm': wall['profile_mm']})
            print('%-8s %-15s %7d %8.2f %9.2f %9.2f %s'
                  % (mid, route, wall['samples'], wall['worst_mm'], wall['deepest_mm'],
                     wall['inside_fraction'], 'OK' if bore_ok else 'OUT/TOO THIN'))

    # ---- can three bores share the R11 seat?
    pitch = 2.0 * BORE_RING_R_MM * np.sin(np.pi / 3.0)
    bore_gap = pitch - 2 * CHANNEL_INNER_R_MM
    edge = shoulder_pads.PAD_RADIUS_MM - BORE_RING_R_MM - CHANNEL_INNER_R_MM
    print('\n--- shared seat can carry three bores ---')
    print('  bore centres on r %.1f mm (120 deg apart): pitch %.2f mm, '
          'bore-to-bore gap %.2f mm, edge clearance %.2f mm'
          % (BORE_RING_R_MM, pitch, bore_gap, edge))
    seat_ok = min(bore_gap, edge) > 1.0
    print('  verdict: %s' % ('OK' if seat_ok else 'TOO TIGHT'))

    # ---- write the envelope
    oname = out_dir / 'deltoid_Y_envelope.stl'
    verts, tris = mesh_geom.manifold_mesh_to_arrays(merged.to_mesh())
    mesh_geom.write_stl(str(oname), verts, tris, header='S1 deltoid Y envelope')
    RV, RF = mesh_geom.read_stl(str(oname))
    RW, RG = mesh_geom.weld(RV, RF)
    rt = abs(mesh_geom.signed_volume(RW, RG) - merged.volume()) / merged.volume()

    report = {
        'scope': 'envelope only: three anatomical heads trimmed by bone, unioned, plus '
                 'the shared R11 seat; no lattice, no end caps, no channel geometry',
        'shared_seat': {'centre_source_mm': np.round(distal, 2).tolist(),
                        'normal_source': np.round(seat_normal, 4).tolist(),
                        'radius_mm': shoulder_pads.PAD_RADIUS_MM,
                        'thickness_mm': shoulder_pads.PAD_THICKNESS_MM},
        'per_head_trim_mm3': {mid: {'anatomy': round(heads[mid][0].volume(), 1),
                                    'trimmed': round(trimmed[mid].volume(), 1)}
                              for mid in DELTOID},
        'merged': {'volume_mm3': round(merged.volume(), 1),
                   'components': len(comps),
                   'triangles': int(merged.num_tri()),
                   'stl_file': oname.name,
                   'stl_roundtrip_rel_error': float('%.2e' % rt)},
        'bone_intersections_mm3': {n: round(I.intersection_volume(merged, bones[n][0]), 3)
                                   for n in BONES},
        'other_module_intersections_mm3': other,
        'channels': channels,
        'seat_bores': {'ring_radius_mm': BORE_RING_R_MM,
                       'pitch_mm': round(pitch, 3),
                       'bore_to_bore_gap_mm': round(bore_gap, 3),
                       'edge_clearance_mm': round(edge, 3),
                       'ok': bool(seat_ok)},
        'not_verified': ['the union is not a built Y: it is the trimmed anatomy plus a '
                         'seat, kept only to test whether a Y is even possible',
                         'no lattice, end caps, crimp or bolt geometry',
                         'channel paths are 3-point polylines, not tendon routes',
                         'no muscle force, pennation or contact'],
        'frozen': False,
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    print('envelope written:', oname, ' roundtrip rel %.2e' % rt)
    return 0


if __name__ == '__main__':
    sys.exit(main())
