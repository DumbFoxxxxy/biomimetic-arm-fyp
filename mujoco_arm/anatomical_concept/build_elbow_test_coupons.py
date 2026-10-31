"""Connected TPU coupons for empirical calibration; no measured performance claim.

X is this generator's loading axis, not an inferred axis mapping to the 2026
paper. The literal 15/10/30 reference is a comparison, not a reproduced specimen.
"""
from build_direct_muscle_yoke import ROOT, OUT, mf, np, rod, write, intersection
from build_continuous_muscle_channels import check_stl
import json, csv

DEST = OUT / 'elbow_test_coupons'


def make_coupon(cell, count, skin=False, guides=False, continuous=False):
    L, W, H = cell
    length = L * count
    pieces = []
    corners = []
    for i in range(count):
        centre = np.array([(i + .5) * L, 0., 0.])
        vertices = [np.array([x, y, z]) for x in [i*L, (i+1)*L]
                    for y in [-W/2, W/2] for z in [-H/2, H/2]]
        pieces += [rod(centre, p, 1., 20) for p in vertices]
        pieces += [mf.Manifold.sphere(1., 16).translate(p) for p in vertices+[centre]]
        corners += vertices
        if skin:
            # Four explicit face diamonds, joined to actual cell corner nodes.
            for dim, half, other, span in [(1, W/2, 2, H/2), (2, H/2, 1, W/2)]:
                for sign in [-1, 1]:
                    mid = centre.copy(); mid[dim] = sign*half
                    for x in [i*L, (i+1)*L]:
                        for s in [-1, 1]:
                            p = mid.copy(); p[0] = x; p[other] = s*span
                            pieces.append(rod(mid, p, 1., 20))
                    pieces.append(mf.Manifold.sphere(1., 16).translate(mid))
        if guides:
            pieces.append(rod(centre-[1.5,0,0], centre+[1.5,0,0], 2., 32))
    # 3mm gripping plates; cell corner spheres overlap them, not point-contact.
    for x in [-3., length]:
        pieces.append(mf.Manifold.cube([3., W+4, H+4]).translate([x,-W/2-2,-H/2-2]))
    if continuous:
        pieces.append(rod([-2,0,0], [length+2,0,0], 1.6, 32))
    s = mf.Manifold.batch_boolean(pieces, mf.OpType.Add)
    if guides:
        s = s - rod([-4,0,0], [length+4,0,0], .85, 32)
    s = s.as_original().simplify(.001)
    assert s.status() == mf.Error.NoError
    components = len([p for p in s.decompose() if p.volume() > .001])
    assert components == 1, components
    wire_overlap = intersection(s, rod([-4,0,0], [length+4,0,0], .3, 32)) if guides else None
    if guides: assert wire_overlap < 1e-6
    # Loading X lies on bed; non-loading Z is vertical. Explicit bed translation.
    lo = np.array(s.bounding_box()[:3]); translation = [-length/2, 0., -float(lo[2])]
    printed = s.translate(translation)
    return printed, dict(cell_loadingX_widthY_heightZ_mm=cell, axial_cells=count,
        nominal_cell_gauge_length_mm=length, grip_plate_thickness_mm=3,
        strut_diameter_mm=2, face_diamonds=skin, short_guides=guides,
        continuous_wall=continuous, bore_diameter_mm=1.7 if guides else None,
        continuous_wall_outer_diameter_mm=3.2 if continuous else None,
        steel_diameter_candidate_mm=.6 if guides else None,
        material_components=components, material_volume_mm3=float(s.volume()),
        steel_intersection_mm3=wire_overlap, print_bounds_mm=list(printed.bounding_box()),
        design_to_print_translation_mm=translation,
        measured_stroke_mm=None, measured_force_N=None, performance_verified=False)


def run():
    DEST.mkdir(exist_ok=True)
    recipes = [('C0_sparse', (10,15,20), 6, False, False, False),
               ('C1_short_guides', (10,15,20), 6, True, True, False),
               ('C2_continuous_channel', (10,15,20), 6, True, True, True),
               ('C3_literal_15_10_30', (15,10,30), 4, False, False, False)]
    rows = []
    for name, cell, count, skin, guides, continuous in recipes:
        folder = DEST/name; folder.mkdir(exist_ok=True)
        solid, info = make_coupon(cell, count, skin, guides, continuous)
        meta = write(folder/'coupon.stl', solid)
        info.update(id=name, mesh=meta, roundtrip=check_stl(folder/'coupon.stl',solid))
        assert info['roundtrip']['pass']
        (folder/'design.json').write_text(json.dumps(info,indent=2), encoding='utf8')
        rows.append(info)
        print(name, 'components', info['material_components'], 'volume',round(solid.volume(),2),flush=True)
    (DEST/'manifest.json').write_text(json.dumps(dict(scope='empirical test candidates, not qualified printing instructions or measured TPU response',
        material_candidate='TPU95A; printer, brand, nozzle and slicer unknown',
        source='local 2026 Yoshimura paper section 3.1; continuous tubular walls and empirical parameter mapping',
        reference_limit='C3 preserves literal older PSM cell dimensions along local X/Y/Z only; not a claim of reproducing 2026 orientation or manufacture',
        independent_specimens_per_recipe=3, coupons=rows),indent=2),encoding='utf8')
    columns=['coupon','specimen','printer','nozzle_mm','material_brand_batch','hardness','print_orientation',
             'layer_mm','temperature_C','measured_gauge_mm','cycle','mode','displacement_mm','force_N',
             'residual_after_60s_mm','wire_pull_mm','wire_jam','end_damage','notes']
    # Preserve any experimental data on generator reruns.
    target = DEST/'measurements_BLANK.csv'
    if not target.exists():
        with target.open('w',newline='',encoding='utf-8-sig') as stream:
            writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
            for name,*_ in recipes:
                for specimen in range(1,4): writer.writerow(dict(coupon=name,specimen=specimen))


if __name__ == '__main__': run()
