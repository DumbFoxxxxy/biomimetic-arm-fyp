"""Export the end pads and caps as STL, one file per part.

Until now the ten pads and ten caps existed only as parameters in s1_end_pads.json and were
rebuilt on demand - there was nothing a machinist or a printer could open.  This writes each part
from exactly the recorded numbers (the same disc() the fitter used, hole and boss included), so
the files and the evidence cannot diverge.

Run:
    python mujoco_arm/tools/export_shoulder_parts.py
    python mujoco_arm/tools/export_shoulder_parts.py --json out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_shoulder_ends as BE
import frames
import interference as I
import mesh_geom

ENDPADS_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_end_pads.json'


def safe(name):
    return name.replace('/', '_').replace(',', '+')


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--boss-mm', type=float, default=3.0)
    ap.add_argument('--out', default=str(frames.OUTPUT_DIR / 'shoulder_parts'))
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    design = json.loads(ENDPADS_JSON.read_text(encoding='utf-8'))
    rows = []
    print('%-34s %-10s %10s %10s' % ('seat', 'part', 'volume mm3', 'tris'))
    for s in design['seats']:
        c = np.asarray(s['centre_source_mm'], float)
        n = np.asarray(s['normal_source'], float)
        lift = float(s['lift_mm'])
        pad_t = float(s['pad_thickness_mm'])
        name = safe(','.join(s['heads']))

        pad = BE.disc_rounded(c, n, lift, pad_t, float(s['pad_radius_mm']),
                              boss_mm=args.boss_mm)
        cap = BE.disc_rounded(c, n, lift + pad_t + args.boss_mm + 0.2,
                              float(s['cap_thickness_mm']), float(s['cap_radius_mm']))
        for kind, geom in (('pad', pad), ('cap', cap)):
            V, F = mesh_geom.manifold_mesh_to_arrays(geom.to_mesh())
            target = out / ('%s_%s.stl' % (name, kind))
            mesh_geom.write_stl(str(target), V, F, header='S1 %s %s' % (kind, name))
            RV, RF = mesh_geom.read_stl(str(target))
            RW, RG = mesh_geom.weld(RV, RF)
            rt = abs(mesh_geom.signed_volume(RW, RG) - geom.volume()) / geom.volume()
            e, bnd, nm = mesh_geom.edge_census(RG)
            rows.append({'seat': ','.join(s['heads']), 'part': kind,
                         'file': target.name, 'volume_mm3': round(geom.volume(), 1),
                         'triangles': int(len(F)),
                         'boundary_edges': int(bnd), 'nonmanifold_edges': int(nm),
                         'roundtrip_rel_error': float('%.2e' % rt)})
            print('%-34s %-10s %10.1f %10d%s'
                  % (name, kind, geom.volume(), len(F),
                     '' if bnd == 0 and nm == 0 else '  ** %d boundary / %d non-manifold'
                     % (bnd, nm)))

    bad = [r for r in rows if r['boundary_edges'] or r['nonmanifold_edges']]
    print('\nfiles written: %d (%d pads, %d caps)'
          % (len(rows), sum(1 for r in rows if r['part'] == 'pad'),
             sum(1 for r in rows if r['part'] == 'cap')))
    print('watertight: %d / %d' % (len(rows) - len(bad), len(rows)))
    if bad:
        print('  not watertight:', [r['file'] for r in bad])

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'STL export of the ten end pads and ten caps, built from the recorded '
                     'definitions with the 3 mm boss, holes and crimp bore included',
            'boss_mm': args.boss_mm, 'out_dir': str(out),
            'parts': rows, 'n_files': len(rows),
            'n_not_watertight': len(bad),
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
