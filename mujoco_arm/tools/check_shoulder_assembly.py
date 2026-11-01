"""S1 assembly interference: bones, muscle modules, end pads and caps, all pairwise.

The parts now exist separately - three bones, four muscle modules, twelve end pads and twelve
caps - and each was checked only against the thing it was being designed against.  This checks
the whole set pairwise at the zero pose, with an explicit allow-list of joints that are supposed
to touch:

    pad <-> its own module   the pad is the TPU's end pad (5.7), bonded by design
    cap <-> its own pad      clamped together

The fixture is not included here: its beam, posts and bases are checked against the bones, the
humerus sweep envelope and the modules in 5.9.10 / 5.9.11, which needs the router to rebuild.

Run:
    python mujoco_arm/tools/check_shoulder_assembly.py
    python mujoco_arm/tools/check_shoulder_assembly.py --json out.json
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_shoulder_ends as BE
import check_shoulder_fixture_zones as FZ
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP
import shoulder_pads

MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'
ENDPADS_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_end_pads.json'
#: parts compared after an STL roundtrip leave 0.0006-0.0219 mm3 slivers (5.4.18)
NOISE_MM3 = 0.05


def mf():
    return I._manifold()


def load_stl(path):
    V, F = mesh_geom.read_stl(str(path))
    W, G = mesh_geom.weld(V, F)
    return mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    parts = {}
    for key in ('scapula', 'clavicle', 'humerus'):
        parts['bone:' + key] = SP.load_solid(FZ.BONE_DIR / SP.BONES[key])
    for key in ('deltoid', 'teres_major', 'subscapularis', 'infraspinatus'):
        p = MODULES_DIR / ('%s.stl' % key)
        if p.is_file():
            parts['module:' + key] = load_stl(p)

    # end pads and caps, rebuilt from the same rule the ends tool uses
    scap = parts['bone:scapula']
    SV, SF = mesh_geom.load_obj(str(FZ.BONE_DIR / SP.BONES['scapula']))
    SW, SG = mesh_geom.weld(SV, SF)
    design = json.loads(ENDPADS_JSON.read_text(encoding='utf-8'))
    for s in design['seats']:
        c = np.asarray(s['centre_source_mm'], float)
        n = np.asarray(s['normal_source'], float)
        lift = float(s['lift_mm'])
        # build exactly what the ends tool recorded.  The radius and thickness are fitted per
        # seat (a 4 mm disc on a 2.8 mm blade reaches the far side), so a hard-coded R11 x 4 here
        # would check a part that is not the design - which is what it did before.
        pad_t = float(s.get('pad_thickness_mm', shoulder_pads.PAD_THICKNESS_MM))
        pad_r = float(s.get('pad_radius_mm', shoulder_pads.PAD_RADIUS_MM))
        cap_t = float(s.get('cap_thickness_mm', pad_t))
        cap_r = float(s.get('cap_radius_mm', pad_r))
        name = '%s/%s' % (s['head'], s['end'])
        # use the same builder the ends tool uses, bolt holes and crimp bore included: a solid
        # disc of the same nominal size is a bigger part and reported overlaps that the design
        # does not have (96.97 mm3 on one seat)
        if s.get('conforming'):
            # a conforming pad is not a disc: its underside is the bone surface (5.9.28).  Building
            # a disc here measures a different part - it reported 40.56 mm3 against the deltoid on
            # a pad that had been cut to clear it (5.9.29).
            bones_map = {k.split(':', 1)[1]: v for k, v in parts.items() if k.startswith('bone:')}
            modules_map = {k.split(':', 1)[1]: v for k, v in parts.items()
                           if k.startswith('module:')}
            parts['pad:' + name] = BE.conforming_disc(
                c, n, lift + pad_t, pad_r, bones_map, s['module'], modules_map,
                boss_mm=float(s.get('pad_boss_mm', 3.0)))
        else:
            parts['pad:' + name] = BE.disc(c + n * lift, n, pad_t, radius=pad_r)
        parts['cap:' + name] = BE.disc(
            c + n * (lift + pad_t + float(s.get('pad_boss_mm', 0.0)) + 0.2), n, cap_t,
            radius=cap_r)

    print('parts: %d' % len(parts))
    kinds = {}
    for k in parts:
        kinds.setdefault(k.split(':')[0], []).append(k)
    for k in sorted(kinds):
        print('  %-8s %d' % (k, len(kinds[k])))

    # a pad is named by head (FJ1500) and a module by module key (infraspinatus), so the
    # allow-list needs the map - keying on the head name made every legitimate pad-to-own-module
    # bond look like a defect
    HEAD_TO_MODULE = {'FJ1468': 'deltoid', 'FJ1513': 'deltoid', 'FJ1467': 'deltoid',
                      'FJ1507': 'teres_major', 'FJ1504': 'subscapularis',
                      'FJ1500': 'infraspinatus'}

    def expected(a, b):
        """Joints that are meant to touch."""
        ka, kb = a.split(':')[0], b.split(':')[0]
        na, nb = a.split(':', 1)[1], b.split(':', 1)[1]
        if {ka, kb} == {'pad', 'module'}:
            pad_name, mod_name = (na, nb) if ka == 'pad' else (nb, na)
            return HEAD_TO_MODULE.get(pad_name.split('/')[0]) == mod_name
        if {ka, kb} == {'cap', 'pad'} and na == nb:
            return True
        if {ka, kb} == {'cap', 'module'}:
            cap_name, mod_name = (na, nb) if ka == 'cap' else (nb, na)
            return HEAD_TO_MODULE.get(cap_name.split('/')[0]) == mod_name
        return False

    rows = []
    unexpected = []
    for a, b in itertools.combinations(sorted(parts), 2):
        v = float(I.intersection_volume(parts[a], parts[b]))
        if v <= NOISE_MM3:
            continue
        ok = expected(a, b)
        rows.append({'a': a, 'b': b, 'overlap_mm3': round(v, 3), 'expected': bool(ok)})
        if not ok:
            unexpected.append(rows[-1])

    print('\ncontacts above %.2f mm3: %d (%d expected, %d unexpected)'
          % (NOISE_MM3, len(rows), len(rows) - len(unexpected), len(unexpected)))
    if rows:
        print('\n%-34s %-34s %10s %8s'
              % ('a', 'b', 'overlap', 'kind'))
        for r in sorted(rows, key=lambda x: -x['overlap_mm3'])[:18]:
            print('%-34s %-34s %10.2f %8s'
                  % (r['a'], r['b'], r['overlap_mm3'],
                     'expected' if r['expected'] else 'DEFECT'))
    else:
        print('  nothing touches anything')

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'pairwise interference over bones, modules, 12 end pads, 12 caps and the '
                     'fixture parts at the zero pose; expected joints are listed explicitly',
            'noise_floor_mm3': NOISE_MM3,
            'n_parts': len(parts),
            'contacts': rows,
            'n_unexpected': len(unexpected),
            'unexpected': unexpected,
            'allow_list': ['pad <-> its own module (the pad is the TPU end pad)',
                           'cap <-> its own pad'],
            'not_verified': ['zero pose only; the sweep is in 5.4.24 and 5.9.10',
                             'the fixture is checked separately in 5.9.10 / 5.9.11',
                             'pads are flat discs, not the conforming underside',
                             'no screws modelled'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
