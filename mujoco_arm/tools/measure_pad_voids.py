"""Measure the void under each end pad against the corrected seat height.

The 1024 mm3 in 5.4.28 was measured before 5.9.19 fixed build_pad (pads used to float 2 mm above
their seat) and before the 3 mm boss.  The lifts are now real - up to 4.69 mm - so the voids under
flat pads are certainly larger.  This recomputes them per seat and totals what a gasket would fill.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, r'D:\XJTLU\PSP\mujoco_arm\tools')

import build_shoulder_ends as BE
import frames
import interference as I
import search_shoulder_channel_paths as SP
import shoulder_pads

ENDPADS = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_end_pads.json'
mf = I._manifold()


def main():
    d = json.loads(ENDPADS.read_text(encoding='utf-8'))
    bones = {k: SP.load_solid(BE.BONE_DIR / v) for k, v in SP.BONES.items()}
    total = 0.0
    rows = []
    print('%-40s %6s %8s %10s' % ('seat', 'lift', 'pad_t', 'void mm3'))
    for s in d['seats']:
        c = np.asarray(s['centre_source_mm'], float)
        n = np.asarray(s['normal_source'], float)
        lift = float(s['lift_mm'])
        pad_t = float(s['pad_thickness_mm'])
        r = float(s['pad_radius_mm'])
        name = ','.join(s['heads'])
        if lift <= 0.01:
            rows.append({'seat': name, 'lift_mm': lift, 'void_mm3': 0.0})
            print('%-40s %6.2f %8.2f %10s' % (name, lift, pad_t, 'none'))
            continue
        # from the seat up to the pad's inner face, minus bone, minus the pad itself
        prism = mf.Manifold.cylinder(lift, r, r, 48)
        Rb = shoulder_pads.pad_basis(n)
        prism = prism.transform(np.c_[Rb, c])
        for b in bones.values():
            prism = prism - b
        pad = shoulder_pads.build_pad(c + n * lift, n, radius=r, thickness=pad_t)
        prism = prism - pad
        v = max(0.0, prism.volume())
        total += v
        rows.append({'seat': name, 'lift_mm': lift, 'pad_t_mm': pad_t,
                     'pad_r_mm': r, 'void_mm3': round(v, 1)})
        print('%-40s %6.2f %8.2f %10.1f' % (name, lift, pad_t, v))
    print('\ntotal void under the pads: %.1f mm3' % total)
    need = [r for r in rows if r['void_mm3'] > 1.0]
    print('seats needing a filler: %d of %d' % (len(need), len(rows)))
    I.write_report(str(frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_pad_voids.json'), {
        'scope': 'void between each pad underside and the bone, recomputed AFTER the 5.9.19 seat '
                 'fix and the 3 mm boss - the 1024 mm3 in 5.4.28 predates both',
        'method': 'cylinder(radius=pad_r, height=lift) at the seat, minus the bones, minus the pad',
        'total_void_mm3': round(total, 1), 'seats': rows,
        'not_verified': ['no gasket or conforming underside is modelled',
                         'the void is measured normal to the seat, not to the bone everywhere'],
    })
    print('\nreport written: s1_pad_voids.json')


if __name__ == '__main__':
    main()
