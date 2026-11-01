"""Do the four S1 muscle modules collide with each other, and with the bones?

Every module so far was checked against the bones one at a time.  Nothing has checked them
against each other, which matters because the deltoid wraps over the infraspinatus and the
subscapularis lies on the other side of the same scapula.

Each module envelope is built by trimming its own anatomy at the bone, so two modules can
still share space wherever their bellies overlapped in the source data - the anatomical
meshes were never guaranteed to be disjoint (5.4.4 already recorded a 121.5 mm3
infraspinatus residual against the deltoid).

Run:
    python mujoco_arm/tools/check_shoulder_module_interference.py
    python mujoco_arm/tools/check_shoulder_module_interference.py --json out.json
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP

MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

MODULES = [('deltoid', '三角肌'), ('teres_major', '大圆肌'),
           ('subscapularis', '肩胛下肌'), ('infraspinatus', '冈下肌')]

#: These modules are compared after an STL roundtrip, and the roundtrip alone leaves
#: slivers of 6e-4 to 2.2e-2 mm3.  The build-time threshold HARD_MM3 = 0.001 is therefore
#: below the noise floor here; anything under this is not a real collision.
ROUNDTRIP_MM3 = 0.05
#: deepest muscle first: the more superficial module yields where they overlap
TRIM_ORDER = ['infraspinatus', 'subscapularis', 'teres_major', 'deltoid']


def mf():
    return I._manifold()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--resolve', action='store_true',
                    help='trim the more superficial module wherever two overlap, deepest '
                         'first, then re-verify')
    args = ap.parse_args(argv)

    solids = {}
    for key, label in MODULES:
        p = MODULES_DIR / ('%s.stl' % key)
        if not p.is_file():
            raise SystemExit('missing %s; run build_shoulder_modules.py first' % p)
        V, F = mesh_geom.read_stl(str(p))
        W, G = mesh_geom.weld(V, F)
        solids[key] = mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))
        print('%-15s %-10s %9.1f mm3  %d piece(s)'
              % (key, label, solids[key].volume(),
                 len([c for c in solids[key].decompose() if c.volume() > 0.5])))

    bones = {k: SP.load_solid(BONE_DIR / v) for k, v in SP.BONES.items()}

    print('\n--- module to module ---')
    print('%-16s %-16s %14s %8s' % ('a', 'b', 'overlap', 'verdict'))
    pairs = []
    for (ka, la), (kb, lb) in itertools.combinations(MODULES, 2):
        v = I.intersection_volume(solids[ka], solids[kb])
        verdict = 'CLEAR' if v <= I.HARD_MM3 else 'OVERLAP'
        pairs.append({'a': ka, 'b': kb, 'overlap_mm3': round(v, 4),
                      'clear': bool(v <= I.HARD_MM3)})
        print('%-16s %-16s %14.4f %8s' % (ka, kb, v, verdict))

    print('\n--- module to bone (must be zero by construction) ---')
    bone_rows = []
    for key, label in MODULES:
        hits = {n: I.intersection_volume(solids[key], b) for n, b in bones.items()}
        worst = max(hits.values())
        bone_rows.append({'module': key, 'intersections_mm3': 
                          {k: round(v, 4) for k, v in hits.items()},
                          'clear': bool(worst <= ROUNDTRIP_MM3)})
        print('%-15s %-10s %s' % (key, label,
                                  'clear' if worst <= ROUNDTRIP_MM3
                                  else ', '.join('%s %.4f' % (n, v)
                                                 for n, v in hits.items() if v > 0)))

    total_overlap = sum(p['overlap_mm3'] for p in pairs)
    worst_pair = max(pairs, key=lambda p: p['overlap_mm3'])
    print('\nworst module-module overlap: %s/%s %.2f mm3'
          % (worst_pair['a'], worst_pair['b'], worst_pair['overlap_mm3']))
    print('all module pairs clear: %s'
          % ('YES' if total_overlap <= ROUNDTRIP_MM3 else 'NO'))
    print('all modules bone-clear: %s'
          % ('YES' if all(b['clear'] for b in bone_rows) else 'NO'))

    resolved = None
    if args.resolve:
        print('\n--- trimming the more superficial module (deepest first) ---')
        trimmed = dict(solids)
        removed = {}
        for i, a in enumerate(TRIM_ORDER):
            for b in TRIM_ORDER[i + 1:]:
                before = trimmed[b].volume()
                trimmed[b] = trimmed[b] - trimmed[a]
                cut = before - trimmed[b].volume()
                if cut > ROUNDTRIP_MM3:
                    removed['%s<- %s' % (b, a)] = round(cut, 3)
                    print('  %-14s -= %-14s removed %8.3f mm3' % (b, a, cut))
        print('  (nothing removed)' if not removed else '')
        # re-verify the trimmed set
        rows = []
        for key, label in MODULES:
            comps = [c for c in trimmed[key].decompose() if c.volume() > 0.5]
            hits = max(I.intersection_volume(trimmed[key], bb).__float__()
                       for bb in bones.values())
            rows.append({'module': key, 'volume_mm3': round(trimmed[key].volume(), 1),
                         'components': len(comps),
                         'worst_bone_mm3': round(hits, 4)})
            print('  %-15s %9.1f mm3  %d piece(s)  worst bone %.4f mm3'
                  % (key, trimmed[key].volume(), len(comps), hits))
        pair_after = []
        for (ka, la), (kb, lb) in itertools.combinations(MODULES, 2):
            v = I.intersection_volume(trimmed[ka], trimmed[kb])
            pair_after.append({'a': ka, 'b': kb, 'overlap_mm3': round(v, 4)})
        worst_after = max(pair_after, key=lambda p: p['overlap_mm3'])
        print('  worst pair after trimming: %s/%s %.4f mm3'
              % (worst_after['a'], worst_after['b'], worst_after['overlap_mm3']))
        print('  all pairs clear after trimming: %s'
              % ('YES' if worst_after['overlap_mm3'] <= ROUNDTRIP_MM3 else 'NO'))
        resolved = {'removed_mm3': removed, 'modules_after': rows,
                    'pairs_after': pair_after,
                    'worst_pair_after_mm3': worst_after['overlap_mm3'],
                    'all_pairs_clear_after': bool(
                        worst_after['overlap_mm3'] <= ROUNDTRIP_MM3)}
        print('  NOTE: the channel walls are NOT re-checked on the trimmed solids')

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'envelope-to-envelope overlap between the four S1 modules and '
                     'module-to-bone; the printed lattice parts are not checked here',
            'pairs': pairs,
            'module_vs_bone': bone_rows,
            'worst_pair_overlap_mm3': worst_pair['overlap_mm3'],
            'roundtrip_noise_mm3': ROUNDTRIP_MM3,
            'all_pairs_clear': bool(total_overlap <= ROUNDTRIP_MM3),
            'all_bone_clear': all(b['clear'] for b in bone_rows),
            'resolved': resolved,
            'not_verified': ['the TPU lattice parts are not checked, only their envelopes',
                             'fixture pads (D17) are not included',
                             'no 0-120 deg sweep'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
