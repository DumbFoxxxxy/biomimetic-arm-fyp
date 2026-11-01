"""S1-P0: the z-calibration coupons (D7) of the shoulder TPU muscle.

D7 fixes the calibration scheme: four coupons with z = 20 / 25 / 30 / 40 mm, three
specimens each, and x = 15, y = 10, p = 2 held at the paper's values.  This script
builds the meshes and the paperwork; it does NOT claim any measured stiffness.

The lattice cell is the same one the elbow uses (``build_elbow_test_coupons.py``):
eight r = 1 mm struts from the cell centre to the eight corners, r = 1 mm spheres at
all nine nodes, so the strut diameter is the paper's p = 2 mm.  Rather than fork that
cell silently, this script rebuilds the two elbow recipes (C0 and C3) and compares
their volumes against the stored evidence - if the numbers disagree, the cell is not
the same cell and the calibration would not transfer.

Coupon geometry chosen here:

    x = 15 mm, y = 10 mm (paper)   4 cells along x, 3 along y, 2 along z
    z = 20 / 25 / 30 / 40 mm       -> 60 x 30 x (40..80) mm envelope
    3 mm grip plates on all four in-plane faces

One coupon therefore gives both E_x and E_y, as §6.6 requires, without a second
specimen.  Two cells through the thickness matches the S1 bellies, whose mean
cross-section is ~25-35 mm across; a one-cell plate would be surface-dominated.

Run:
    python mujoco_arm/tools/build_shoulder_z_coupons.py --verify-elbow
    python mujoco_arm/tools/build_shoulder_z_coupons.py --out <dir>
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frames
import interference as I
import mesh_geom

DEFAULT_OUT = frames.OUTPUT_DIR / 'shoulder_z_coupons'
ELBOW_COUPONS = frames.OUTPUT_DIR / 'elbow_test_coupons'

#: paper-reference parameters held fixed for the calibration (D7 / §6.6) [mm]
X_MM, Y_MM, P_MM = 15.0, 10.0, 2.0
Z_VALUES_MM = (20.0, 25.0, 30.0, 40.0)
SPECIMENS_PER_Z = 3
#: lattice coupons [mm]
GRIP_PLATE_MM = 3.0
CELLS = {'x': 4, 'y': 3, 'z': 2}

#: the two elbow recipes, rebuilt here as a cross-check of the cell [mm]
ELBOW_RECIPES = {'C0_sparse': {'cell': (10.0, 15.0, 20.0), 'count': 6,
                               'volume_mm3': 4557.204},
                 'C3_literal_15_10_30': {'cell': (15.0, 10.0, 30.0), 'count': 4,
                                         'volume_mm3': 4455.283}}


def _mf():
    return I._manifold()


def rod(a, b, r, segments=20):
    """Cylinder from a to b with radius r."""
    mf = _mf()
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    d = b - a
    length = float(np.linalg.norm(d))
    if length < 1e-12:
        return None
    cyl = mf.Manifold.cylinder(length, r, r, segments)
    z = d / length
    helper = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    x = np.cross(helper, z)
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    return cyl.transform(np.c_[np.stack([x, y, z], axis=1), a])


def cell_pieces(centre, half):
    """The project's X cell: 8 struts centre -> corners, r1 nodes, at one cell."""
    mf = _mf()
    c = np.asarray(centre, dtype=float)
    hx, hy, hz = half
    corners = [c + np.array([sx * hx, sy * hy, sz * hz])
               for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    pieces = [rod(c, p, P_MM / 2.0) for p in corners]
    pieces += [mf.Manifold.sphere(P_MM / 2.0, 16).translate(p) for p in corners + [c]]
    return pieces


def lattice_block(cell, nx, ny, nz, plates='all'):
    """A block of cells with 3 mm grip plates.

    plates='axial' reproduces the elbow coupon: plates only on the two x ends.
    plates='all'   adds plates on all four in-plane faces, so one coupon can be
                   compressed along x as well as along y.
    """
    mf = _mf()
    lx, ly, lz = cell
    pieces = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                centre = np.array([(i + 0.5) * lx, (j + 0.5) * ly, (k + 0.5) * lz])
                pieces += cell_pieces(centre, (lx / 2.0, ly / 2.0, lz / 2.0))
    sx, sy = nx * lx, ny * ly
    t = GRIP_PLATE_MM
    pieces.append(mf.Manifold.cube([t, sy + 4.0, lz + 4.0]).translate([-t, -2.0, -2.0]))
    pieces.append(mf.Manifold.cube([t, sy + 4.0, lz + 4.0]).translate([sx, -2.0, -2.0]))
    if plates == 'all':
        pieces.append(mf.Manifold.cube([sx + 4.0, t, lz + 4.0]).translate([-2.0, -t, -2.0]))
        pieces.append(mf.Manifold.cube([sx + 4.0, t, lz + 4.0]).translate([-2.0, sy, -2.0]))
    solid = mf.Manifold.batch_boolean(pieces, I._manifold().OpType.Add)
    solid = solid.as_original().simplify(0.001)
    if solid.status() != mf.Error.NoError:
        raise SystemExit('boolean failed: %s' % solid.status())
    return solid


def components(solid):
    return len([c for c in solid.decompose() if c.volume() > 0.001])


def verify_elbow():
    """Rebuild the two stored elbow coupons and compare volumes."""
    print('--- cell cross-check against the stored elbow coupons ---')
    ok = True
    for name, recipe in ELBOW_RECIPES.items():
        path = ELBOW_COUPONS / name / 'design.json'
        stored = recipe['volume_mm3']
        if path.is_file():
            with open(path, 'r', encoding='utf-8') as fh:
                stored = float(json.load(fh)['material_volume_mm3'])
        solid = lattice_block(recipe['cell'], recipe['count'], 1, 1, plates='axial')
        vol = solid.volume()
        rel = abs(vol - stored) / stored
        comp = components(solid)
        good = rel < 1e-4 and comp == 1
        ok &= good
        print('[%s] %-22s rebuilt %10.3f mm3  stored %10.3f  rel %.2e  components %d'
              % ('PASS' if good else 'FAIL', name, vol, stored, rel, comp))
    print('RESULT:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


def build(out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    print('--- S1-P0 coupons: x %.0f y %.0f p %.0f, %dx%dx%d cells ---'
          % (X_MM, Y_MM, P_MM, CELLS['x'], CELLS['y'], CELLS['z']))
    print('%-8s %6s %8s %10s %10s %10s  %s'
          % ('z,mm', 'cells', 'comps', 'volume', 'bbox x', 'bbox y', 'bbox z'))
    for z in Z_VALUES_MM:
        cell = (X_MM, Y_MM, z)
        solid = lattice_block(cell, CELLS['x'], CELLS['y'], CELLS['z'], plates='all')
        comp = components(solid)
        if comp != 1:
            raise SystemExit('z=%.0f coupon is %d components, not 1' % (z, comp))
        verts, tris = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh())
        bb = solid.bounding_box()
        name = 'z%02d' % int(z)
        folder = out_dir / name
        folder.mkdir(exist_ok=True)
        mesh_geom.write_stl(str(folder / 'coupon.stl'), verts, tris,
                            header='S1-P0 z=%d coupon, body parts 3D derived' % int(z))
        RV, RF = mesh_geom.read_stl(str(folder / 'coupon.stl'))
        RW, RG = mesh_geom.weld(RV, RF)
        stl_vol = abs(mesh_geom.signed_volume(RW, RG))
        rel = abs(stl_vol - solid.volume()) / max(solid.volume(), 1e-9)
        if rel > 1e-5:
            raise SystemExit('z=%.0f STL roundtrip relative error %.2e' % (z, rel))
        info = {'id': name, 'stage': 'S1-P0',
                'parameters_mm': {'x': X_MM, 'y': Y_MM, 'z': z, 'p': P_MM},
                'cells': dict(CELLS),
                'strut_diameter_mm': P_MM,
                'grip_plate_mm': GRIP_PLATE_MM,
                'grip_plates': 'all four in-plane faces',
                'material_components': comp,
                'material_volume_mm3': round(float(solid.volume()), 3),
                'stl_roundtrip_volume_mm3': round(stl_vol, 3),
                'stl_roundtrip_rel_error': float('%.2e' % rel),
                'stl_file': 'coupon.stl (generated; not committed from this checkout '
                            'because new LFS objects cannot be added here)',
                'triangle_count': int(len(tris)),
                'bounding_box_mm': [round(float(v), 2) for v in bb],
                'specimens_required': SPECIMENS_PER_Z,
                'load_directions': ['x (compression)', 'y (contraction)'],
                'print_orientation': 'z vertical, x horizontal (same convention as the '
                                     'product muscle, §5.8)',
                'measured': None,
                'performance_verified': False}
        (folder / 'design.json').write_text(json.dumps(info, ensure_ascii=False, indent=2),
                                            encoding='utf-8')
        rows.append(info)
        print('%-8.0f %6d %8d %10.3f %10.1f %10.1f %10.1f'
              % (z, CELLS['x'] * CELLS['y'] * CELLS['z'], comp,
                 solid.volume(), bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2]))

    manifest = {'scope': 'S1-P0 z-calibration coupons; manufacturing geometry only, no '
                         'measured stiffness, no printing qualification',
                'decision': 'D7 - four z values, monotonic interpolation inside 20-40 mm',
                'held_fixed_mm': {'x': X_MM, 'y': Y_MM, 'p': P_MM},
                'z_values_mm': list(Z_VALUES_MM),
                'specimens_per_z': SPECIMENS_PER_Z,
                'total_specimens': SPECIMENS_PER_Z * len(Z_VALUES_MM),
                'cell_source': 'same X cell as build_elbow_test_coupons.py, verified by '
                               'rebuilding C0 and C3 (see --verify-elbow)',
                'print_orientation': 'z vertical, x horizontal; same printer, nozzle, '
                                     'layer height and material batch as the product',
                'coupons': rows,
                'pending': ['actual printing', 'E_x/E_y', 'kink strain', 'hysteresis',
                            '60 s residual', 'monotonic fit and residual',
                            'z=25 vs z=40 stiffness ratio must land in 0.2-0.5'],
                'frozen': False}
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2),
                                           encoding='utf-8')
    readme = out_dir / 'measurements_BLANK.csv'
    columns = ['coupon', 'specimen', 'printer', 'nozzle_mm', 'material_brand_batch',
               'hardness', 'print_orientation', 'layer_mm', 'temperature_C',
               'load_axis', 'measured_gauge_mm', 'displacement_mm', 'force_N', 'E_MPa',
               'kink_strain', 'cycle', 'hysteresis_mm', 'residual_after_60s_mm', 'notes']
    if not readme.exists():
        with readme.open('w', newline='', encoding='utf-8-sig') as fh:
            w = csv.DictWriter(fh, fieldnames=columns)
            w.writeheader()
            for z in Z_VALUES_MM:
                for s in range(1, SPECIMENS_PER_Z + 1):
                    for axis in ('x', 'y'):
                        w.writerow({'coupon': 'z%02d' % int(z), 'specimen': s,
                                    'load_axis': axis})
    print('\nmanifest and blank measurement sheet written to', out_dir)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify-elbow', action='store_true')
    ap.add_argument('--out', default=str(DEFAULT_OUT))
    args = ap.parse_args(argv)
    if args.verify_elbow:
        return verify_elbow()
    return build(args.out)


if __name__ == '__main__':
    sys.exit(main())
