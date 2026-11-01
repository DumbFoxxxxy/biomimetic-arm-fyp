"""Flat lifted pad vs conforming pad, on all ten seats.

5.9.28 measures ~2950 mm3 of void under the flat pads.  A pad whose underside IS the bone surface
would have no void at all - this quantifies both options per seat:

  flat       : cylinder from the seat, height lift + pad_t, lifted by `lift` as today
  conforming : the same cylinder minus the bones, so it hugs the surface with no gap

For each: volume, bone interference, and the void that remains underneath.
"""
import json
import sys

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
    rows = []
    tot_flat = tot_conf = tot_void_flat = 0.0
    print('%-40s %9s %9s %9s %8s' % ('seat', 'flat mm3', 'conf mm3', 'flat bone', 'void'))
    for s in d['seats']:
        c = np.asarray(s['centre_source_mm'], float)
        n = np.asarray(s['normal_source'], float)
        lift = max(float(s['lift_mm']), 0.05)
        t = float(s['pad_thickness_mm'])
        r = float(s['pad_radius_mm'])
        name = ','.join(s['heads'])
        R = shoulder_pads.pad_basis(n)

        prism = mf.Manifold.cylinder(lift + t, r, r, 48).transform(np.c_[R, c])
        flat = shoulder_pads.build_pad(c + n * lift, n, radius=r, thickness=t)
        flat_bone = max(I.intersection_volume(flat, b).__float__() for b in bones.values())
        # void under the flat pad = the prism from seat to the pad underside, minus bone, minus pad
        gap = mf.Manifold.cylinder(lift, r, r, 48).transform(np.c_[R, c])
        for b in bones.values():
            gap = gap - b
        gap = gap - flat
        void = max(0.0, gap.volume())

        conf = prism
        for b in bones.values():
            conf = conf - b
        conf_bone = max(I.intersection_volume(conf, b).__float__() for b in bones.values())

        tot_flat += flat.volume()
        tot_conf += conf.volume()
        tot_void_flat += void
        rows.append({'seat': name, 'lift_mm': round(lift, 2), 'pad_t_mm': t, 'pad_r_mm': r,
                     'flat_mm3': round(flat.volume(), 1),
                     'flat_bone_mm3': round(flat_bone, 3),
                     'void_under_flat_mm3': round(void, 1),
                     'conforming_mm3': round(conf.volume(), 1),
                     'conforming_bone_mm3': round(conf_bone, 3)})
        print('%-40s %9.1f %9.1f %9.3f %8.1f' % (name, flat.volume(), conf.volume(), flat_bone, void))

    print('\nflat total %.1f mm3, conforming total %.1f mm3 (-%.1f%%)'
          % (tot_flat, tot_conf, 100.0 * (tot_flat - tot_conf) / tot_flat))
    print('void under the flat pads: %.1f mm3' % tot_void_flat)
    print('void under the conforming pads: 0 by construction (the underside IS the surface)')
    I.write_report(str(frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_pad_conforming.json'), {
        'scope': 'flat lifted pad vs bone-conforming pad on all ten seats, measured after the '
                 '5.9.19 seat fix and the 3 mm boss',
        'flat_total_mm3': round(tot_flat, 1), 'conforming_total_mm3': round(tot_conf, 1),
        'void_under_flat_mm3': round(tot_void_flat, 1),
        'seats': rows,
        'not_verified': ['a print/contact clearance on the conforming face is not modelled '
                         '(the underside touches the bone with zero gap)',
                         'cap positions and the tendon 7 mm standoff would shift if this is adopted',
                         'no adhesion or bedding compound is specified'],
    })
    print('\nreport written: s1_pad_conforming.json')


if __name__ == '__main__':
    main()
