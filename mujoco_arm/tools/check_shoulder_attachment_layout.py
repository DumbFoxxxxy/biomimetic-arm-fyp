"""L1 bone-side interface study: can six muscle lines land on the proximal humerus?

Three questions, answered with solids, in this order:

1. **pad vs pad** - six 22 mm end pads at the six insertion candidates, and the
   four-pad grouped layout (the three deltoid heads share one seat, D4/D8).
2. **pad vs bearing** - the pads are humerus-side parts, so the real obstacle is
   the *scapula*-side GH hardware.  Two bearings are tested:
   a. ``--bearing g1``: the existing G1 GH solid joint (r=16 mm ball, r=16.4/19.5
      socket, two anti-dislocation caps, M3 fasteners), loaded from
      ``GH_connection/``.  Static neutral pose.
   b. ``--bearing cap``: an open shallow cup on the *anatomical* head (bore
      r_head + 0.4 mm, wall 3 mm), swept over rim half-angle alpha, to find the
      largest alpha that touches neither the four pads nor the humerus bone.
3. **layout table** - each insertion in a GH spherical frame: radius, polar angle
   from the glenoid normal, azimuth, and the pad's far-corner radius, which is
   the number that decides whether a socket bore can clear it at all.

Scope and limits, to be repeated in any report:

* the pad is a plain cylinder - no chamfer, no bracket, no bolt head, no tendon;
* insertion centroids are proximity-sample centroids of candidate regions, not
  confirmed attachment sites;
* static neutral pose only - no sweep, no compliant routing, no assembly order;
* a pass here does not mean the muscle envelope or the drive fits.

Run:
    python mujoco_arm/tools/check_shoulder_attachment_layout.py
    python mujoco_arm/tools/check_shoulder_attachment_layout.py --json out.json
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frames
import interference as I
import mesh_geom

REGION_DIR = frames.OUTPUT_DIR / 'attachment_regions'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'
HUMERUS_OBJ = BONE_DIR / '肱骨.obj'
SCAPULA_OBJ = BONE_DIR / '肩胛骨.obj'
GH_DIR = (frames.OUTPUT_DIR / 'joint_model' / 'SC_connection' / 'AC_compact'
          / 'GH_connection')

PAD_RADIUS_MM = 11.0
PAD_THICKNESS_MM = 4.0
BOLT_PITCH_MM = 14.0
BOLT_HOLE_RADIUS_MM = 1.7

#: radial clearance between the anatomical head and the socket bore [mm]
SOCKET_BORE_CLEARANCE_MM = 0.4
#: socket wall thickness [mm]
SOCKET_WALL_MM = 3.0
#: glenoid normal cut-offs for the sensitivity band [mm]
GLENOID_CUTOFFS_MM = (26.0, 28.0, 30.0)

#: parts that move with the scapula in the G1 assembly (part_owners: AC)
G1_FIXED_PARTS = ('scapula_AC_GH_sockets', 'GH_cap_0', 'GH_cap_1',
                  'GH_bolt_0', 'GH_bolt_1', 'GH_nut_0', 'GH_nut_1',
                  'GH_top_washer_0', 'GH_top_washer_1',
                  'GH_bottom_washer_0', 'GH_bottom_washer_1')

MUSCLES = {
    'FJ1468': '三角肌前束',
    'FJ1513': '三角肌后束',
    'FJ1467': '三角肌中束',
    'FJ1507': '大圆肌',
    'FJ1504': '肩胛下肌',
    'FJ1500': '冈下肌',
}
DELTOID = ('FJ1468', 'FJ1513', 'FJ1467')


# ------------------------------------------------------------------ geometry
def insertion_region(mid):
    path = REGION_DIR / ('%s_insertion_humerus.npz' % mid)
    if not path.is_file():
        return None
    with np.load(path) as d:
        points = d['bone_points_source_mm']
        weights = d['muscle_face_areas_mm2']
        keep = d['distances_mm'] <= np.median(d['distances_mm'])
        tris = d['bone_triangle_indices']
    return points[keep], weights[keep], tris[keep]


def region_centroid(mid, kind):
    """Area-weighted bone-surface centroid of an origin/insertion candidate region."""
    for bone in ('scapula', 'clavicle', 'humerus'):
        path = REGION_DIR / ('%s_%s_%s.npz' % (mid, kind, bone))
        if not path.is_file():
            continue
        with np.load(path) as d:
            points = d['bone_points_source_mm']
            weights = d['muscle_face_areas_mm2']
            keep = d['distances_mm'] <= np.median(d['distances_mm'])
        return ((points[keep] * weights[keep][:, None]).sum(0) / weights[keep].sum(),
                bone, int(keep.sum()), float(weights[keep].sum()))
    return None, None, 0, 0.0


def pad_basis(normal):
    n = frames.unit(normal)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = frames.unit(np.cross(n, helper))
    y = np.cross(n, x)
    return np.stack([x, y, n], axis=1)


def build_pad(centre_display, normal_source):
    """Cylinder pad sitting on the bone surface, growing outward along the normal."""
    mf = I._manifold()
    centre_source = frames.display_to_source(centre_display)
    R = pad_basis(normal_source)
    solid = mf.Manifold.cylinder(PAD_THICKNESS_MM, PAD_RADIUS_MM, PAD_RADIUS_MM, 48)
    solid = solid.translate([0.0, 0.0, PAD_THICKNESS_MM / 2.0])
    T = np.c_[R, centre_source]
    return solid.transform(T)


def pad_far_radius(centre_source, normal_source, centre_ref):
    """Largest distance from `centre_ref` to any point of the pad envelope [mm].

    This is the number a socket bore has to clear: the pad lies *inside* the
    head radius at its centre, but its rim corners stick out.
    """
    c = np.asarray(centre_source, dtype=float)
    R = pad_basis(normal_source)
    worst = 0.0
    for t in (0.0, PAD_THICKNESS_MM):
        for k in range(72):
            a = 2.0 * math.pi * k / 72.0
            p = c + R[:, 2] * t + R[:, 0] * PAD_RADIUS_MM * math.cos(a) \
                + R[:, 1] * PAD_RADIUS_MM * math.sin(a)
            worst = max(worst, float(np.linalg.norm(p - centre_ref)))
    return worst


def glenoid_normal(gh, radius, scaps):
    """Unit vector GH -> centroid of scapula surface within a cutoff of the head."""
    d = np.linalg.norm(scaps - gh, axis=1)
    out = []
    for cut in GLENOID_CUTOFFS_MM:
        m = d <= cut
        c = scaps[m].mean(0)
        out.append((cut, int(m.sum()), frames.unit(c - gh), float(d.min())))
    return out


def fossa_cap(gh, scaps, band=2.0):
    """Angular radius of the bony glenoid fossa [deg], as the alpha reference.

    The fossa is taken as the scapula surface within `band` mm of its closest
    approach to the head centre; its angular radius about the fossa direction is
    the coverage a real glenoid gives, i.e. the anatomical scale for alpha.
    """
    d = np.linalg.norm(scaps - gh, axis=1)
    m = d <= d.min() + band
    c = scaps[m].mean(0)
    n = frames.unit(c - gh)
    ang = np.degrees(np.arccos(np.clip((scaps[m] - gh) @ n / d[m], -1.0, 1.0)))
    return {'band_mm': band, 'vertices': int(m.sum()), 'normal': n.tolist(),
            'closest_mm': float(d.min()),
            'angular_radius_median_deg': float(np.median(ang)),
            'angular_radius_p95_deg': float(np.percentile(ang, 95)),
            'angular_radius_max_deg': float(ang.max())}


def head_cap_error(V, gh, normal, alpha_deg, radius):
    """How spherical is the *bone* inside the contact cap?  [mm]

    A printed cup cannot run on the raw bone if this spread is large, so the
    spread decides whether a spherical sleeve is needed between bone and cup.
    Vertices farther than 1.5 x the fitted radius are the shaft and are skipped.
    """
    v = V - gh
    r = np.linalg.norm(v, axis=1)
    th = np.degrees(np.arccos(np.clip(v @ normal / np.maximum(r, 1e-9), -1.0, 1.0)))
    m = (th <= alpha_deg) & (r <= 1.5 * radius)
    dev = r[m] - radius
    return {'alpha_deg': alpha_deg, 'vertices': int(m.sum()),
            'radial_min_mm': float(dev.min()), 'radial_max_mm': float(dev.max()),
            'radial_rms_mm': float(np.sqrt(np.mean(dev ** 2))),
            'nominal_radius_mm': radius}


def spherical_coordinates(gh, normal, p):
    """(r, theta from normal, azimuth in degrees) of `p` about the GH centre."""
    v = np.asarray(p, dtype=float) - gh
    r = float(np.linalg.norm(v))
    n = np.asarray(normal, dtype=float)
    theta = math.degrees(math.acos(np.clip(v @ n / max(r, 1e-12), -1.0, 1.0)))
    helper = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = frames.unit(np.cross(n, helper))
    e2 = np.cross(n, e1)
    phi = math.degrees(math.atan2(v @ e2, v @ e1))
    return r, theta, phi, e1, e2


def build_cup(alpha_deg, gh_source, bore_radius, normal):
    """Open shallow cup: spherical shell about GH, kept only inside a cone."""
    mf = I._manifold()
    centre = np.asarray(gh_source, dtype=float)      # everything stays in source mm
    outer = bore_radius + SOCKET_WALL_MM
    shell = (mf.Manifold.sphere(outer, 96) - mf.Manifold.sphere(bore_radius, 96))
    shell = shell.translate(list(centre))
    ang = math.radians(min(alpha_deg, 88.0))
    length = outer + 6.0
    cone = mf.Manifold.cylinder(length, 0.05, length * math.tan(ang), 96)
    n = frames.unit(normal)
    helper = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    x = frames.unit(np.cross(n, helper))
    y = np.cross(n, x)
    cone = cone.transform(np.c_[np.stack([x, y, n], axis=1), centre])
    return shell ^ cone


# --------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--bearing', choices=('both', 'g1', 'cap'), default='both')
    args = ap.parse_args(argv)

    V, F = mesh_geom.load_obj(str(HUMERUS_OBJ))
    tri = V[F]
    humerus_centre = V.mean(axis=0)
    gh = frames.joint_centres('source')['GH']
    jc = json.loads(frames.CANDIDATES_PATH.read_text(encoding='utf-8'))
    head_radius = float(jc['GH']['radius_mm'])
    bore = head_radius + SOCKET_BORE_CLEARANCE_MM

    scaps, scap_f = mesh_geom.load_obj(str(SCAPULA_OBJ))
    glen = glenoid_normal(gh, head_radius, scaps)
    normal = glen[-1][2]
    spread = max(math.degrees(math.acos(np.clip(g[2] @ normal, -1, 1))) for g in glen)
    closest_fossa = min(g[3] for g in glen)
    fossa = fossa_cap(gh, scaps)

    hum_w, hum_f = mesh_geom.weld(V, F)
    sc_w, sc_ff = mesh_geom.weld(scaps, scap_f)
    print('weld check: humerus %s -> %s, scapula %s -> %s '
          '(edges, boundary, non-manifold)'
          % (mesh_geom.edge_census(F), mesh_geom.edge_census(hum_f),
             mesh_geom.edge_census(scap_f), mesh_geom.edge_census(sc_ff)))

    print('humerus %d verts   GH %s   fitted head radius %.2f mm'
          % (len(V), np.round(gh, 1), head_radius))
    print('closest scapula surface to the head centre: %.2f mm '
          '(= %.2f mm outside the fitted head sphere)' % (closest_fossa,
                                                          closest_fossa - head_radius))
    print('glenoid normal (GH -> fossa) %s   definition spread %.1f deg'
          % (np.round(normal, 4), spread))
    for cut, cnt, nv, dmin in glen:
        print('    cutoff %4.0f mm: %4d verts  n = %s  (%.1f deg from chosen)'
              % (cut, cnt, np.round(nv, 4),
                 math.degrees(math.acos(np.clip(nv @ normal, -1, 1)))))
    print('bony glenoid fossa within %.1f mm of its closest approach: %d verts, '
          'angular radius median %.1f / p95 %.1f / max %.1f deg'
          % (fossa['band_mm'], fossa['vertices'],
             fossa['angular_radius_median_deg'], fossa['angular_radius_p95_deg'],
             fossa['angular_radius_max_deg']))

    # ---------------------------------------------------------- pads
    pads = []
    for mid, name in MUSCLES.items():
        region = insertion_region(mid)
        if region is None:
            print('SKIP', mid)
            continue
        points, weights, tris = region
        centroid = (points * weights[:, None]).sum(0) / weights.sum()
        idx = int(np.argmax(weights))          # largest sample of the region
        face = tri[tris[idx]]
        nrm = np.cross(face[1] - face[0], face[2] - face[0])
        nrm = nrm / max(np.linalg.norm(nrm), 1e-12)
        if nrm @ (centroid - humerus_centre) < 0:
            nrm = -nrm
        r, theta, phi, e1, e2 = spherical_coordinates(gh, normal, centroid)
        pads.append({'id': mid, 'name': name,
                     'centre_source_mm': centroid.tolist(),
                     'centre_display_mm': frames.source_to_display(centroid).tolist(),
                     'normal_source': nrm.tolist(),
                     'distance_to_gh_mm': round(float(np.linalg.norm(centroid - gh)), 2),
                     'radius_mm': round(r, 2),
                     'polar_from_glenoid_deg': round(theta, 1),
                     'azimuth_deg': round(phi, 1),
                     'pad_far_radius_mm': round(pad_far_radius(centroid, nrm, gh), 2),
                     'region_samples': int(len(points)),
                     'region_area_mm2': round(float(weights.sum()), 1)})

    print('\n--- insertion layout in the GH spherical frame (theta from glenoid normal) ---')
    print('%-8s %-10s %7s %7s %8s %9s %10s'
          % ('muscle', 'name', 'r/mm', 'theta', 'phi', 'pad rmax', 'rmax-bore'))
    for p in pads:
        print('%-8s %-10s %7.2f %7.1f %8.1f %9.2f %10.2f'
              % (p['id'], p['name'], p['radius_mm'], p['polar_from_glenoid_deg'],
                 p['azimuth_deg'], p['pad_far_radius_mm'],
                 p['pad_far_radius_mm'] - bore))

    for p in pads:
        p['_solid'] = build_pad(np.asarray(p['centre_display_mm']),
                                np.asarray(p['normal_source']))

    # ------------------------------------------------ cap-free retention test
    print('\n--- cap-free retention: does every line pull the head into the socket? ---')
    print('    pull = unit(origin - insertion);  compression = pull . n_glenoid')
    print('%-8s %-10s %-9s %10s %8s %9s'
          % ('line', 'name', 'origin', 'compression', 'shear', 'seat angle'))
    comp_rows = []
    for p in pads:
        origin, bone, ns, area = region_centroid(p['id'], 'origin')
        if origin is None:
            print('%-8s %-10s  origin region missing' % (p['id'], p['name']))
            continue
        pull = frames.unit(origin - np.asarray(p['centre_source_mm']))
        c = float(pull @ normal)
        shear = float(np.linalg.norm(pull - c * normal))
        seat = math.degrees(math.atan2(shear, c))
        comp_rows.append({'id': p['id'], 'name': p['name'], 'origin_bone': bone,
                          'origin_source_mm': origin.tolist(),
                          'compression': round(c, 4), 'shear': round(shear, 4),
                          'seat_angle_deg': round(seat, 1),
                          'seats_ball': bool(c > 0.0),
                          'origin_region_samples': ns,
                          'origin_region_area_mm2': round(area, 1)})
        print('%-8s %-10s %-9s %10.3f %8.3f %9.1f'
              % (p['id'], p['name'], bone, c, shear, seat))
    n_seat = sum(1 for r in comp_rows if r['seats_ball'])
    print('lines whose pull has a positive socket-axis component: %d / %d'
          % (n_seat, len(comp_rows)))
    print('minimum compression %.3f  ->  cap-free retention %s'
          % (min(r['compression'] for r in comp_rows) if comp_rows else float('nan'),
             'holds on muscle preload alone' if comp_rows and n_seat == len(comp_rows)
             else 'NOT guaranteed by the six lines alone'))

    print('\n--- six separate pads: pairwise ---')
    print('%-9s %-9s %-10s %-9s %s' % ('pair', 'd,mm', 'needed,mm', 'gap,mm', 'overlap'))
    centre_pairs = []
    worst = (1e9, None)
    for i in range(len(pads)):
        for j in range(i + 1, len(pads)):
            a, b = pads[i], pads[j]
            d = float(np.linalg.norm(np.asarray(a['centre_source_mm']) -
                                     np.asarray(b['centre_source_mm'])))
            vol = I.intersection_volume(a['_solid'], b['_solid'])
            centre_pairs.append({'a': a['id'], 'b': b['id'],
                                 'centre_distance_mm': round(d, 2),
                                 'centre_gap_mm': round(d - 2 * PAD_RADIUS_MM, 2),
                                 'overlap_mm3': round(vol, 4),
                                 'pass': bool(vol <= I.HARD_MM3)})
            if d < worst[0]:
                worst = (d, (a['id'], b['id']))
            flag = 'OVERLAP' if vol > I.HARD_MM3 else ('tight' if d - 2 * PAD_RADIUS_MM < 2.0 else 'ok')
            print('%-9s %-9.1f %-10.1f %-9.1f %s'
                  % ('%s/%s' % (a['id'][2:], b['id'][2:]), d, 2 * PAD_RADIUS_MM,
                     d - 2 * PAD_RADIUS_MM, '%.3f mm3 %s' % (vol, flag)))
    clashes = [c for c in centre_pairs if not c['pass']]
    print('pairs overlapping: %d / %d   closest pair %s at %.1f mm'
          % (len(clashes), len(centre_pairs), '/'.join(worst[1]), worst[0]))

    # ------------------------------------------- grouped: deltoid = one seat
    delt = [p for p in pads if p['id'] in DELTOID]
    pts = np.array([p['centre_source_mm'] for p in delt])
    w = np.array([p['region_area_mm2'] for p in delt])
    merged = (pts * w[:, None]).sum(0) / w.sum()
    merged_normal = frames.unit((np.array([p['normal_source'] for p in delt])
                                 * w[:, None]).sum(0))
    groups = [{'id': 'DELTOID', 'name': '三角肌（三束合并）',
               'centre_source_mm': merged.tolist(),
               'centre_display_mm': frames.source_to_display(merged).tolist(),
               'normal_source': merged_normal.tolist(),
               'members': [p['id'] for p in delt],
               'radius_mm': round(float(np.linalg.norm(merged - gh)), 2),
               'polar_from_glenoid_deg': round(spherical_coordinates(
                   gh, normal, merged)[1], 1),
               'pad_far_radius_mm': round(pad_far_radius(merged, merged_normal, gh), 2)}]
    for p in pads:
        if p['id'] not in DELTOID:
            groups.append({'id': p['id'], 'name': p['name'],
                           'centre_source_mm': p['centre_source_mm'],
                           'centre_display_mm': p['centre_display_mm'],
                           'normal_source': p['normal_source'],
                           'members': [p['id']],
                           'radius_mm': p['radius_mm'],
                           'polar_from_glenoid_deg': p['polar_from_glenoid_deg'],
                           'pad_far_radius_mm': p['pad_far_radius_mm']})
    for g in groups:
        g['_solid'] = build_pad(np.asarray(g['centre_display_mm']),
                                np.asarray(g['normal_source']))

    print('\n--- grouped layout: %d seats ---' % len(groups))
    print('%-9s %-9s %-10s %-11s %s'
          % ('pair', 'd,mm', 'gap,mm', 'r_max,mm', 'overlap'))
    grouped_pairs = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            a, b = groups[i], groups[j]
            d = float(np.linalg.norm(np.asarray(a['centre_source_mm']) -
                                     np.asarray(b['centre_source_mm'])))
            vol = I.intersection_volume(a['_solid'], b['_solid'])
            grouped_pairs.append({'a': a['id'], 'b': b['id'],
                                  'centre_distance_mm': round(d, 2),
                                  'centre_gap_mm': round(d - 2 * PAD_RADIUS_MM, 2),
                                  'overlap_mm3': round(vol, 4),
                                  'pass': bool(vol <= I.HARD_MM3)})
            print('%-9s %-9.1f %-10.1f %-11.2f %.3f mm3 %s'
                  % ('%s/%s' % (a['id'][:6], b['id'][:6]), d, d - 2 * PAD_RADIUS_MM,
                     max(a['pad_far_radius_mm'], b['pad_far_radius_mm']), vol,
                     'OVERLAP' if vol > I.HARD_MM3 else 'ok'))
    grouped_clashes = [c for c in grouped_pairs if not c['pass']]
    over_bore = [g['id'] for g in groups
                 if g['pad_far_radius_mm'] > bore and g['radius_mm'] < 40.0]
    print('grouped pairs overlapping: %d / %d' % (len(grouped_clashes), len(grouped_pairs)))
    print('seats near the head whose pad corner reaches past the socket bore (r=%.2f): %s'
          % (bore, ', '.join('%s %.2f' % (g['id'], g['pad_far_radius_mm'])
                             for g in groups if g['id'] in over_bore) or 'none'))

    report = {'scope': 'static neutral pose; plain cylindrical pad envelopes at the '
                       'insertion candidate centroids; proximity samples are not '
                       'confirmed attachment sites; no bracket, bolt head, tendon, '
                       'compliant routing or assembly order',
              'pad': {'radius_mm': PAD_RADIUS_MM, 'thickness_mm': PAD_THICKNESS_MM,
                      'bolt_pitch_mm': BOLT_PITCH_MM,
                      'hole_radius_mm': BOLT_HOLE_RADIUS_MM},
              'gh': {'centre_source_mm': gh.tolist(), 'fitted_head_radius_mm': head_radius,
                     'socket_bore_mm': bore,
                     'socket_wall_mm': SOCKET_WALL_MM,
                     'glenoid_normal': normal.tolist(),
                     'glenoid_normal_definition_spread_deg': round(spread, 2),
                     'glenoid_cutoff_sensitivity': [
                         {'cutoff_mm': c, 'vertices': n_, 'normal': v.tolist(),
                          'off_from_chosen_deg': round(math.degrees(
                              math.acos(np.clip(v @ normal, -1, 1))), 2)}
                         for c, n_, v, _ in glen],
                     'closest_scapula_surface_mm': round(closest_fossa, 2),
                     'fossa_cap': fossa},
              'muscles': [{k: v for k, v in p.items() if not k.startswith('_')}
                          for p in pads],
              'pairs': centre_pairs,
              'overlapping_pairs': len(clashes),
              'closest_pair_mm': round(worst[0], 2),
              'six_separate_pads_feasible': not clashes,
              'grouped': {'pads': [{k: v for k, v in g.items() if not k.startswith('_')}
                                   for g in groups],
                          'pairs': grouped_pairs,
                          'overlapping_pairs': len(grouped_clashes),
                          'pads_past_socket_bore': over_bore,
                          'feasible': not grouped_clashes},
              'cap_free_retention': {'normal': 'glenoid normal, GH -> fossa',
                                     'lines': comp_rows,
                                     'lines_seating': n_seat,
                                     'all_lines_seat': bool(comp_rows and n_seat == len(comp_rows))}}

    # ------------------------------- pads vs the fixed scapula bone itself
    mf = I._manifold()
    scap_solid = mf.Manifold(mf.Mesh64(sc_w, sc_ff.astype(np.uint32)))
    print('\n--- pads vs the fixed scapula bone (neutral; independent of bearing) ---')
    print('%-9s %-10s %-12s %s' % ('seat', 'r/mm', 'overlap', 'note'))
    scap_rows = []
    for g in groups:
        vol = I.intersection_volume(g['_solid'], scap_solid)
        gapv = I.gap(g['_solid'], scap_solid) if vol <= I.HARD_MM3 else 0.0
        scap_rows.append({'id': g['id'], 'radius_mm': g['radius_mm'],
                          'overlap_mm3': round(vol, 4), 'gap_mm': round(gapv, 3),
                          'pass': bool(vol <= I.HARD_MM3)})
        print('%-9s %-10.2f %-12.3f %s'
              % (g['id'], g['radius_mm'], vol,
                 'CLEAR >= %.2f mm' % gapv if vol <= I.HARD_MM3 else 'COLLIDES'))
    scap_bad = [r for r in scap_rows if not r['pass']]
    print('seats colliding with the scapula bone: %s'
          % (', '.join(r['id'] for r in scap_bad) or 'none'))
    report['scapula_bone'] = {'seats': scap_rows,
                              'colliding_seats': [r['id'] for r in scap_bad]}

    # ---------------------------------------------- bearing A: existing G1
    if args.bearing in ('both', 'g1') and GH_DIR.is_dir():
        mf = I._manifold()
        shift = np.c_[np.eye(3), frames.SHIFT]      # G1 parts are saved in display mm
        g1 = {}
        for name in G1_FIXED_PARTS:
            path = GH_DIR / (name + '.obj')
            if not path.is_file():
                print('G1 part missing:', name)
                continue
            vv, ff = mesh_geom.load_obj(str(path))
            g1[name] = mf.Manifold(mf.Mesh64(vv, ff.astype(np.uint32))).transform(shift)
        fixed = g1['scapula_AC_GH_sockets']
        for name, solid in g1.items():
            if name != 'scapula_AC_GH_sockets':
                fixed = fixed + solid
        print('\n--- bearing A: existing G1 solid GH joint (static, neutral) ---')
        print('    ball r=16.0  socket 16.4/19.5  two caps depth 3.0  M3 fasteners')
        print('%-9s %-10s %-12s %s' % ('seat', 'r_max,mm', 'intrusion', 'hits'))
        g1_rows = []
        for g in groups:
            vol = I.intersection_volume(g['_solid'], fixed)
            hits = {}
            for name, solid in g1.items():
                v2 = I.intersection_volume(g['_solid'], solid)
                if v2 > I.HARD_MM3:
                    hits[name] = round(v2, 4)
            g1_rows.append({'id': g['id'], 'pad_far_radius_mm': g['pad_far_radius_mm'],
                            'intrusion_mm3': round(vol, 4), 'hits': hits,
                            'pass': bool(vol <= I.HARD_MM3)})
            print('%-9s %-10.2f %-12.3f %s'
                  % (g['id'], g['pad_far_radius_mm'], vol,
                     ', '.join('%s %.3f' % (k, v) for k, v in hits.items()) or '-'))
        g1_bad = [r for r in g1_rows if not r['pass']]
        print('G1 bearing compatible with the four seats: %s (%d/%d collide)'
              % ('YES' if not g1_bad else 'NO', len(g1_bad), len(g1_rows)))
        report['bearing_g1'] = {
            'source': str(GH_DIR.relative_to(frames.PROJECT_ROOT)).replace('\\', '/'),
            'ball_radius_mm': 16.0, 'socket_inner_mm': 16.4, 'socket_outer_mm': 19.5,
            'pose': 'neutral, static',
            'seats': g1_rows,
            'colliding_seats': [r['id'] for r in g1_bad],
            'compatible': not g1_bad}

    # ------------- bearing B: how deep could *any* socket be before it hits us?
    if args.bearing in ('both', 'cap'):
        mf = I._manifold()
        hum_solid = mf.Manifold(mf.Mesh64(hum_w, hum_f.astype(np.uint32)))

        def cup_clear(alpha, b):
            cup = build_cup(float(alpha), gh, b, normal)
            if I.intersection_volume(hum_solid, cup) > I.HARD_MM3:
                return False
            return all(I.intersection_volume(g['_solid'], cup) <= I.HARD_MM3
                       for g in groups)

        def alpha_max_for(b):
            """Largest rim half-angle with no pad and no bone penetration (0 = none)."""
            if not cup_clear(0.5, b):
                return 0.0
            lo, hi = 0.5, None
            for a in np.arange(2.5, 90.1, 2.5):
                if cup_clear(float(a), b):
                    lo = float(a)
                else:
                    hi = float(a)
                    break
            if hi is None:
                return 90.0
            while hi - lo > I.BOUNDARY_DEG:
                mid = 0.5 * (lo + hi)
                if cup_clear(mid, b):
                    lo = mid
                else:
                    hi = mid
            return lo

        # (label, bore, why) - the fitted head is NOT a usable journal surface, see
        # head_sphericity below, so a bore that actually clears the bone is also tried.
        max_bone = head_radius + max(
            head_cap_error(V, gh, normal, a, head_radius)['radial_max_mm']
            for a in (25.0, 40.0))
        bores = [('fitted-head', bore,
                  'bore = fitted head %.2f + %.1f' % (head_radius,
                                                      SOCKET_BORE_CLEARANCE_MM)),
                 ('bone-clearance', round(max_bone + 0.1, 2),
                  'bore clears the worst bone vertex inside the contact cap')]
        print('\n--- bearing B: how deep may a socket be before it touches a pad or the bone? ---')
        print('    (shell: spherical, centred at GH, wall %.1f mm, open cone, no cap)' % SOCKET_WALL_MM)
        cup_report = []
        for label, b, why in bores:
            amax = alpha_max_for(b)
            hits = []
            for g in groups:
                a = 0.5
                while a < 89.0:
                    if I.intersection_volume(g['_solid'], build_cup(a, gh, b, normal)) > I.HARD_MM3:
                        hits.append({'id': g['id'], 'first_hit_deg': round(a, 1)})
                        break
                    a += 1.0
            row = {'label': label, 'bore_mm': b, 'why': why,
                   'alpha_max_deg': round(amax, 2),
                   'seat_first_hit': hits}
            cup_report.append(row)
            print('  %-14s bore %6.2f mm: alpha_max = %5.2f deg   seat hits: %s'
                  % (label, b, amax,
                     ', '.join('%s@%.0f' % (h['id'][:6], h['first_hit_deg'])
                               for h in hits) or 'none'))

        report['bearing_cup'] = {
            'wall_mm': SOCKET_WALL_MM,
            'variants': cup_report,
            'head_sphericity': [head_cap_error(V, gh, normal, a, head_radius)
                                for a in (25.0, 35.0, 40.0, 45.0)],
            'scope': 'spherical shell vs pad envelopes and vs the welded humerus bone, '
                     'neutral pose only; no brace, rim or fixation geometry'}
        print('\nhead sphericity inside the contact cap (deviation from r = %.2f mm):'
              % head_radius)
        for row in report['bearing_cup']['head_sphericity']:
            print('  alpha %4.0f deg: %4d verts  dev %+.2f .. %+.2f mm  rms %.2f mm'
                  % (row['alpha_deg'], row['vertices'], row['radial_min_mm'],
                     row['radial_max_mm'], row['radial_rms_mm']))
        print('worst bone radius inside the contact cap: %.2f mm = fitted %.2f + %.2f'
              % (max_bone, head_radius, max_bone - head_radius))


    if args.json_out:
        with open(args.json_out, 'w', encoding='utf-8') as fh:
            json.dump(report, fh, ensure_ascii=False, indent=2)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
