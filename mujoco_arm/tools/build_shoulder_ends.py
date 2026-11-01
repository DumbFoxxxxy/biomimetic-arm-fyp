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


def disc(centre, normal, thickness, radius=PAD_R, holes=True, boss_mm=0.0, boss_r=5.0,
         taper_mm=0.0):
    """R11 plate growing outward along `normal`, with the bolt holes and the crimp bore.

    `boss_mm` adds a smaller cylinder on the plate's outer face.  5.4.31 needs the tendon segment
    to start 7 mm above the bone while the fitted plate is only 4 mm thick, so the missing 3 mm has
    to be a boss - the same feature 5.9.18's one-sided clamp needs anyway.

    `taper_mm` makes the plate a frustum, `radius` at the bone end shrinking by `taper_mm` at the
    outer end.  5.9.20 left one seat whose cap intrudes 0.611 mm3 with every plain disc; that is a
    rim contact, so the rim is what has to move.
    """
    solid = shoulder_pads.build_pad(centre, normal, radius=radius, thickness=thickness,
                                    radius_top=None if taper_mm <= 0 else radius - taper_mm)
    if boss_mm > 0.0:
        R0 = shoulder_pads.pad_basis(normal)
        boss = mf().Manifold.cylinder(boss_mm, boss_r, boss_r, 48)
        boss = boss.translate([0.0, 0.0, thickness + boss_mm / 2.0])
        solid = solid + boss.transform(np.c_[R0, np.asarray(centre, float)])
    if not holes:
        return solid
    R = shoulder_pads.pad_basis(normal)
    total = thickness + boss_mm
    # holes run through the plate (and the boss) along the normal so a screw can pass
    for k in (0, 1):
        ang = np.pi * k
        off = R[:, 0] * (BOLT_PITCH / 2.0) * np.cos(ang) + R[:, 1] * (BOLT_PITCH / 2.0) * np.sin(ang)
        hole = mf().Manifold.cylinder(total * 3.0, BOLT_HOLE_R, BOLT_HOLE_R, 32)
        hole = hole.translate([0.0, 0.0, -total])
        world = np.c_[R, np.asarray(centre, float) + off]
        solid = solid - hole.transform(world)
    bore = mf().Manifold.cylinder(total * 3.0, CRIMP_R, CRIMP_R, 32)
    bore = bore.translate([0.0, 0.0, -total])
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


def blade_thickness(point, normal, bone_key, bone_tris, limit=15.0):
    """How far the bone extends inward along -normal from the seat [mm], or None.

    A flat pad thicker than this must emerge on the far side - which is what 5.9.12 class B
    turned out to be: a 4 mm pad on a 2.8 mm scapular blade, not a placement error.
    """
    tris = bone_tris[list(SP.BONES).index(bone_key)]
    p = np.asarray(point, float)
    n = np.asarray(normal, float)
    t = 0.3
    if mesh_geom.signed_distance(p - n * t, tris) > 0.0:
        return None                     # not on this bone's surface in this direction
    while t < limit:
        if mesh_geom.signed_distance(p - n * t, tris) > 0.0:
            return float(t)
        t += 0.1
    return None


def disc_rounded(centre, normal, lift, thickness, radius, boss_mm=0.0):
    """The disc exactly as it will be recorded and re-read by the other tools.

    The fit loops used to measure unrounded values and then store rounded ones, so a candidate
    that "passed" could still be a part that overlaps 87.95 mm3 - which is precisely the 10x
    disagreement between this tool and the assembly checker.  Measuring the recorded form removes
    the gap instead of explaining it.
    """
    c = np.round(np.asarray(centre, float) + np.asarray(normal, float) * round(lift, 2), 2)
    n = np.round(np.asarray(normal, float), 3)
    return disc(c, n, round(thickness, 2), radius=round(radius, 1), boss_mm=round(boss_mm, 2))


def conforming_disc(centre, normal, height, radius, bones, own_module, modules,
                    boss_mm=0.0):
    """A pad whose underside is the bone surface instead of a flat face lifted clear of it.

    5.9.28 measured 2 988 mm3 of air under the ten flat pads - a 4.69 mm lift means only one rim
    point touches the bone.  Subtracting the bones from a prism that spans the whole lift makes the
    underside conform: the pad now fills what used to be air, no lift is needed to keep it out of
    the bone, and the load path is the whole footprint.

    `height` is measured from the seat, so the caller passes lift + pad thickness: the part of the
    prism below the bone surface is removed, and the outer face stays where the flat pad's was.
    """
    R = shoulder_pads.pad_basis(normal)
    c = np.round(np.asarray(centre, float), 2)
    prism = mf().Manifold.cylinder(round(height, 2), round(radius, 1), round(radius, 1), 64)
    solid = prism.transform(np.c_[R, c])
    for k, b in bones.items():
        solid = solid - b
    for k, s in modules.items():
        if k != own_module:
            solid = solid - s
    total = round(height, 2) + boss_mm
    if boss_mm > 0.0:
        boss = mf().Manifold.cylinder(boss_mm, 5.0, 5.0, 48)
        boss = boss.translate([0.0, 0.0, round(height, 2) + boss_mm / 2.0])
        solid = solid + boss.transform(np.c_[R, c])
    for k in (0, 1):
        ang = np.pi * k
        off = R[:, 0] * (BOLT_PITCH / 2.0) * np.cos(ang) + R[:, 1] * (BOLT_PITCH / 2.0) * np.sin(ang)
        hole = mf().Manifold.cylinder(total * 3.0, BOLT_HOLE_R, BOLT_HOLE_R, 32)
        hole = hole.translate([0.0, 0.0, -total])
        solid = solid - hole.transform(np.c_[R, c + off])
    bore = mf().Manifold.cylinder(total * 3.0, CRIMP_R, CRIMP_R, 32)
    bore = bore.translate([0.0, 0.0, -total])
    solid = solid - bore.transform(np.c_[R, c])
    return solid


def main(argv=None):
    # A seat the fitter could not solve can still take a smaller, laterally offset clamp at the
    # same attachment point: 5.9.18 found the deltoid posterior origin needs only a 1.0 mm offset,
    # which costs 2.68% of its abduction moment arm - against 45.1% for the 17.39 mm seat move D22
    # rejected.  Those solved clamps are read here and used verbatim.
    clamp_specs = {}
    _clamp_path = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_one_sided_clamp.json'
    if _clamp_path.is_file():
        _c = json.loads(_clamp_path.read_text(encoding='utf-8'))
        # two shapes: one solved seat ('seat' + 'found') or several ('solved': {seat: spec})
        if _c.get('solved'):
            clamp_specs.update({k: v for k, v in _c['solved'].items() if v})
        elif _c.get('found'):
            clamp_specs[_c['seat']] = _c['found']
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--pad-boss-mm', type=float, default=3.0,
                    help='boss height on the pad; 3 mm is the specified value - 5.4.31 needs the '
                         'tendon to start 7 mm above the bone and the plate is only 4 mm')
    ap.add_argument('--conforming', action='store_true',
                    help='pad underside = bone surface (5.9.28); removes the up-to-4.69 mm lift and '
                         'the 2 988 mm3 of air under the flat pads')
    ap.add_argument('--conforming-clear-mm', type=float, default=0.0,
                    help='nominal bedding clearance on the conforming face (not a true offset)')
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

    # --- collect the seats, then merge duplicates ---
    # The deltoid's three heads share one insertion seat (D8).  Building a pad per head put three
    # coincident pads there, which 5.9.12 measured as six pairs of 1516.19 mm3.  A seat gets one
    # pad and one cap, and the heads' bores simply converge inside it.
    raw = []
    for m in design['modules']:
        for c in m['channels']:
            path = next(x['path_source_mm'] for x in m['corridors'] if x['head'] == c['head'])
            for end, p in (('origin', path[0]), ('distal', path[-1])):
                raw.append({'module': m['module'], 'head': c['head'], 'end': end,
                            'p': np.asarray(p, float)})

    merged = []
    for r in raw:
        for g in merged:
            if np.linalg.norm(r['p'] - g['p']) < 1.0:
                g['members'].append('%s/%s' % (r['head'], r['end']))
                break
        else:
            merged.append({'p': r['p'], 'module': r['module'],
                           'members': ['%s/%s' % (r['head'], r['end'])]})
    print('seats: %d raw -> %d distinct (merged: %s)'
          % (len(raw), len(merged),
             '; '.join(', '.join(g['members']) for g in merged if len(g['members']) > 1)))

    seats = []
    for g in merged:
        p = g['p']
        seat_p = p.copy()
        n = bone_normal(p, bone_tris)
        bk, bsolid = nearest_bone(p, bones)
        g_module = g['module']
        spec = clamp_specs.get(','.join(g['members']))
        if spec is not None:
            p = seat_p + np.asarray(spec['offset_source_mm'], float)
            print('  %-22s using the solved one-sided clamp: offset %.1f mm'
                  % (','.join(g['members']), spec['offset_mm']))
        # The pad must clear the bones AND the neighbouring modules.  On a 2.8 mm blade with the
        # opposite-side muscle pressed against the far face even a 2.3 mm disc reaches it, so the
        # footprint is fitted rather than fixed.
        blade = blade_thickness(p, n, bk, bone_tris)
        t_max = PAD_T if blade is None else float(max(1.0, min(PAD_T, blade - 0.5)))
        # The pad and the cap have to be fitted TOGETHER.  The cap's height is lift + pad thickness
        # + boss, so a thinner pad lowers the cap - and whether that helps depends on the seat.  The
        # first version fixed the pad at the thickest candidate that cleared the bone and fitted the
        # cap afterwards, so a thinner pad could never be used to bring the cap down; that is what
        # left the deltoid middle origin with "no cap" (5.9.21).  Now every pad candidate that
        # clears the bone is carried through to the cap, and the first pair that clears both wins.
        pad_r, thick, pad, lift = PAD_R, t_max, None, 0.0
        cap_r, cap_t, cap, cap_hit = PAD_R, t_max, None, None
        pad_cleared = cap_cleared = False
        if spec is not None:
            # verbatim from the solved clamp: no search, no re-fitting
            pad_r = float(spec['pad_radius_mm'])
            thick = float(spec['pad_thickness_mm'])
            lift, _ = shoulder_pads.lift_to_clear(p, n, bsolid, radius=pad_r, thickness=thick)
            lift = 0.0 if lift is None else lift
            pad = disc_rounded(p, n, lift, thick, pad_r, boss_mm=args.pad_boss_mm)
            cap_r = float(spec['cap_radius_mm'])
            cap_t = float(spec['cap_thickness_mm'])
            cap = disc_rounded(p, n, lift + thick + args.pad_boss_mm + 0.2, cap_t, cap_r)
            pad_cleared = cap_cleared = pad is not None and cap is not None
        else:
            best = None
            for cand_t in sorted({t_max, round(t_max * 0.7, 2), 2.5, 2.0, 1.5, 1.0}, reverse=True):
                for cand_r in (PAD_R, 10.0, 9.0, 8.0, 7.0, 6.0):
                    cand_lift, _ = shoulder_pads.lift_to_clear(p, n, bsolid, radius=cand_r,
                                                               thickness=cand_t)
                    cand_lift = 0.0 if cand_lift is None else cand_lift
                    if args.conforming:
                        # 5.9.28: the flat pad floats above the bone (up to 4.69 mm of air), so the
                        # underside IS the bone surface instead.  The prism spans the same height
                        # the flat pad occupied (lift + thickness), so the outer face and therefore
                        # the cap stay exactly where they were - what changes is that the space
                        # below is now material instead of air.
                        geom = conforming_disc(p, n, cand_lift + cand_t, cand_r,
                                              bones, g['module'], modules,
                                              boss_mm=args.pad_boss_mm)
                    else:
                        geom = disc_rounded(p, n, cand_lift, cand_t, cand_r,
                                            boss_mm=args.pad_boss_mm)
                    hit = max(max(I.intersection_volume(geom, b).__float__()
                                  for b in bones.values()),
                              max((I.intersection_volume(geom, s).__float__()
                                   for k, s in modules.items() if k != g_module), default=0.0))
                    if hit > SURFACE_TOL_MM3:
                        continue
                    # pad clears - now the cap at this pad's height
                    found_cap = None
                    for c_t in sorted({cand_t, round(cand_t * 0.6, 2), 1.5, 1.0}, reverse=True):
                        for c_r in (PAD_R, 10.0, 9.0, 8.0, 7.0, 6.0):
                            c_geom = disc_rounded(
                                p, n, cand_lift + cand_t + args.pad_boss_mm + 0.2, c_t, c_r)
                            c_hit = max(max(I.intersection_volume(c_geom, b).__float__()
                                            for b in bones.values()),
                                        max((I.intersection_volume(c_geom, s).__float__()
                                             for k, s in modules.items() if k != g_module),
                                            default=0.0))
                            if c_hit <= SURFACE_TOL_MM3:
                                found_cap = (c_r, c_t, c_geom, c_hit)
                                break
                        if found_cap:
                            break
                    if found_cap is not None:
                        best = (cand_r, cand_t, geom, cand_lift, *found_cap)
                        break
                    if best is None:
                        best = (cand_r, cand_t, geom, cand_lift, PAD_R, cand_t, None, None)
                if best is not None and best[6] is not None:
                    break
            if best is not None:
                pad_r, thick, pad, lift, cap_r, cap_t, cap, cap_hit = best
        pad_bone = max(I.intersection_volume(pad, b).__float__() for b in bones.values())
        pad_mod = max((I.intersection_volume(pad, s).__float__() for k, s in modules.items()
                       if k != g['module']), default=0.0)
        pad_cleared = pad_bone <= SURFACE_TOL_MM3 and pad_mod <= SURFACE_TOL_MM3
        if cap is None:
            cap = disc_rounded(p, n, lift + thick + args.pad_boss_mm + 0.2, cap_t, cap_r)
        cap_bone = max(I.intersection_volume(cap, b).__float__() for b in bones.values())
        cap_mod = max((I.intersection_volume(cap, s).__float__() for k, s in modules.items()
                       if k != g['module']), default=0.0)
        cap_cleared = cap_bone <= SURFACE_TOL_MM3 and cap_mod <= SURFACE_TOL_MM3

        seats.append({
            'module': g['module'], 'heads': g['members'],
            'head': g['members'][0].split('/')[0], 'end': g['members'][0].split('/')[1],
            'centre_source_mm': np.round(p, 2).tolist(),
            'seat_source_mm': np.round(seat_p, 2).tolist(),
            'normal_source': np.round(n, 3).tolist(),
            'nearest_bone': bk,
            'lift_mm': round(lift, 2),
            'blade_thickness_mm': None if blade is None else round(blade, 2),
            'pad_thickness_mm': round(thick, 2),
            # the shape must travel with the numbers: a conforming pad is not a flat disc, and a
            # reader that rebuilds a disc measures a different part (5.9.29)
            'conforming': bool(args.conforming),
            'pad_footprint_cleared': bool(pad_cleared),
            'cap_footprint_cleared': bool(cap_cleared),
            'pad_radius_mm': round(pad_r, 1),      # the fitted radius, not the constant
            'cap_radius_mm': round(cap_r, 2),
            'cap_thickness_mm': round(cap_t, 2),
            'pad_volume_mm3': round(pad.volume(), 1),
            'cap_volume_mm3': round(cap.volume(), 1),
            'pad_bone_mm3': round(pad_bone, 3),
            'pad_other_module_mm3': round(pad_mod, 3),
            'cap_bone_mm3': round(cap_bone, 3),
            'cap_other_module_mm3': round(cap_mod, 3),
        })
        r = seats[-1]
        print('  %-24s lift %4.2f  blade %5s  pad t %4.2f  cap r %4.1f  '
              'pad/bone %.2f  cap/bone %.2f  cap/mod %.2f'
              % (','.join(g['members']), r['lift_mm'], r['blade_thickness_mm'],
                 r['pad_thickness_mm'], r['cap_radius_mm'], r['pad_bone_mm3'],
                 r['cap_bone_mm3'], r['cap_other_module_mm3']))

    print('\n--- self-check: rebuild every pad and cap from its RECORDED numbers ---')
    # The fit loops above report on the candidate they were evaluating, which is not necessarily
    # what got written down (rounded centre/normal/lift/radius/thickness).  That mismatch is what
    # made this tool record 8.33 mm3 on a seat the assembly checker measured at 87.95.  Reporting
    # is done here, from the same numbers the other tools read, so agreement is by construction.
    mismatches = 0
    for s in seats:
        c = np.asarray(s['centre_source_mm'], float)
        n = np.asarray(s['normal_source'], float)
        lift = float(s['lift_mm'])
        if args.conforming:
            # rebuild the RECORDED form: a conforming pad is not a flat disc, so measuring a disc
            # here would report a mismatch that is really a different part.
            pad = conforming_disc(c, n, lift + s['pad_thickness_mm'], s['pad_radius_mm'],
                                  bones, s['module'], modules, boss_mm=args.pad_boss_mm)
        else:
            pad = disc(c + n * lift, n, s['pad_thickness_mm'], radius=s['pad_radius_mm'])
        cap = disc(c + n * (lift + s['pad_thickness_mm'] + args.pad_boss_mm + 0.2), n,
                   s['cap_thickness_mm'], radius=s['cap_radius_mm'])
        pad_bone = max(I.intersection_volume(pad, b).__float__() for b in bones.values())
        cap_bone = max(I.intersection_volume(cap, b).__float__() for b in bones.values())
        pad_mod = max((I.intersection_volume(pad, m).__float__() for k, m in modules.items()
                       if k != s['module']), default=0.0)
        cap_mod = max((I.intersection_volume(cap, m).__float__() for k, m in modules.items()
                       if k != s['module']), default=0.0)
        if abs(pad_mod - s['pad_other_module_mm3']) > 0.5 \
                or abs(cap_mod - s['cap_other_module_mm3']) > 0.5:
            mismatches += 1
            print('  %-22s was %.2f/%.2f -> is %.2f/%.2f (pad/cap vs other module)'
                  % (','.join(s['heads']), s['pad_other_module_mm3'],
                     s['cap_other_module_mm3'], pad_mod, cap_mod))
        s['pad_bone_mm3'] = round(pad_bone, 3)
        s['pad_other_module_mm3'] = round(pad_mod, 3)
        s['cap_bone_mm3'] = round(cap_bone, 3)
        s['cap_other_module_mm3'] = round(cap_mod, 3)
    print('  seats whose recorded numbers did not match what was reported: %d' % mismatches)

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
            'conforming': bool(args.conforming),
            'pad_boss_mm': args.pad_boss_mm,
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
