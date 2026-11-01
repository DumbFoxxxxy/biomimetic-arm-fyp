"""S1 end caps: can the elbow cap form be transplanted to the shoulder?

The elbow cap is fully specified in
``mujoco_arm/anatomical_concept/build_shaped_elbow_muscles.py``:

    line 70-72  p = end surface point + 9 mm outward
                TPU pad  = cylinder r 11.0, thickness 4.0, centred on p
                16 transition rods r 1.2 on a ring of radius 7.0, pad -> surface
    line 141-150 rigid cap = cylinder r 11.0, thickness 4.0, from p+2.15 to p+6.15
                - two M3 through holes on a 14.0 mm pitch (r 1.7)
                - central wire bore r 0.8, crimp barrel r 1.25
                - proximal cable port boss r 3.0 x 6.0 long, bore r 1.1
                - distal crimp seat recess r 1.4

So the end form is three stacked pieces (pad / rods / cap), and the only one that
depends on the belly shape is the **transition ring at radius 7.0 mm**.  The pad
and the cap are R11 regardless, because §5.5 already sized the bone-side pad at
R11 - the cap is bolted onto that pad, not onto the belly.

This script therefore asks one question per muscle end: does the anatomical end
section have a 7 mm inradius, so the transition ring can land on it?  If not, the
belly needs a local flare (a stalk), and the script reports how much radial growth
that takes.

Run:
    python mujoco_arm/tools/size_shoulder_end_caps.py
    python mujoco_arm/tools/size_shoulder_end_caps.py --json out.json
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import attachment_regions
import frames
import interference as I
import mesh_geom
import muscle_sections

MUSCLE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / 'muscle'

MUSCLES = {
    'FJ1468': '三角肌前束',
    'FJ1513': '三角肌后束',
    'FJ1467': '三角肌中束',
    'FJ1507': '大圆肌',
    'FJ1504': '肩胛下肌',
    'FJ1500': '冈下肌',
}

#: elbow cap form, measured from build_shaped_elbow_muscles.py lines 70-72 / 141-150 [mm]
ELBOW_CAP = {'radius_mm': 11.0, 'thickness_mm': 4.0,
             'bolt_pitch_mm': 14.0, 'bolt_hole_r_mm': 1.7,
             'wire_bore_r_mm': 0.8, 'cable_boss_r_mm': 3.0, 'cable_boss_len_mm': 6.0,
             'cable_bore_r_mm': 1.1, 'crimp_seat_r_mm': 1.4,
             'crimp_barrel_r_mm': 1.25, 'wire_r_mm': 0.45}
#: elbow TPU end pad: R11 x 4 mm, sitting 9 mm out from the muscle end surface [mm]
ELBOW_PAD = {'radius_mm': 11.0, 'thickness_mm': 4.0,
             'offset_from_surface_mm': 9.0,
             'transition_ring_r_mm': 7.0, 'transition_rod_r_mm': 1.2, 'rods': 16}

#: the transition ring is what has to land on the belly, so that is the threshold
RING_R_MM = ELBOW_PAD['transition_ring_r_mm']
#: two M3 on a 14 mm pitch need this much radius to keep the holes inside the cap
MIN_RADIUS_TWO_BOLT_MM = ELBOW_CAP['bolt_pitch_mm'] / 2.0 + ELBOW_CAP['bolt_hole_r_mm']


def classify(inradius):
    """Can the 16-rod transition ring land on this end section?"""
    if inradius is None:
        return 'no-section'
    if inradius >= RING_R_MM:
        return 'ring-ok'
    return 'flare-needed'


def solid_of(path):
    V, F = mesh_geom.load_obj(str(path))
    W, G = mesh_geom.weld(V, F)
    mf = I._manifold()
    return mf.Manifold(mf.Mesh64(W, G.astype(np.uint32)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    print('elbow end form: pad r%.1f x %.1f at +%.1f from the surface, %d rods r%.1f '
          'on a r%.1f ring, cap r%.1f x %.1f, 2 x M3 on %.1f pitch'
          % (ELBOW_PAD['radius_mm'], ELBOW_PAD['thickness_mm'],
             ELBOW_PAD['offset_from_surface_mm'], ELBOW_PAD['rods'],
             ELBOW_PAD['transition_rod_r_mm'], RING_R_MM,
             ELBOW_CAP['radius_mm'], ELBOW_CAP['thickness_mm'],
             ELBOW_CAP['bolt_pitch_mm']))
    print('the belly-dependent part is the r%.1f transition ring; '
          'pad and cap stay R%.1f because the bone pad is R%.1f (§5.5)'
          % (RING_R_MM, ELBOW_CAP['radius_mm'], ELBOW_PAD['radius_mm']))

    rows = []
    print('\n--- end sections (t=0.02 origin end, t=0.98 insertion end) ---')
    print('%-8s %-9s %-6s %8s %8s %8s %8s  %s'
          % ('muscle', 'end', 'bone', 'area', 'wall', 'inrad', 'growth', 'verdict'))
    for mid, name in MUSCLES.items():
        solid = solid_of(MUSCLE_DIR / (mid + '.obj'))
        origin = attachment_regions.load_region(mid, 'origin')
        insertion = attachment_regions.load_region(mid, 'insertion')
        o, _ = attachment_regions.centroid(origin)
        i, _ = attachment_regions.centroid(insertion)
        for label, t, bone in (('origin', 0.02, origin['bone']),
                               ('insertion', 0.98, insertion['bone'])):
            sec = muscle_sections.section_at(solid, o, i, t)
            if sec is None:
                continue
            verdict = classify(sec['inradius_mm'])
            growth = None if sec['inradius_mm'] is None else \
                round(max(0.0, RING_R_MM - sec['inradius_mm']), 2)
            rows.append({'id': mid, 'name': name, 'end': label, 'bone': bone, 't': t,
                         'area_mm2': sec['area_mm2'], 'wall_mm': sec['wall_mm'],
                         'inradius_mm': sec['inradius_mm'],
                         'inradius_xy_mm': sec['inradius_xy_mm'],
                         'inradius_source_mm': sec['inradius_source_mm'],
                         'centroid_source_mm': sec['centroid_source_mm'],
                         'radial_growth_to_ring_mm': growth,
                         'verdict': verdict})
            print('%-8s %-9s %-6s %8.1f %8.2f %8s %8s  %s'
                  % (mid, label, bone, sec['area_mm2'], sec['wall_mm'],
                     '-' if sec['inradius_mm'] is None else '%.2f' % sec['inradius_mm'],
                     '-' if growth is None else '%.2f' % growth, verdict))

    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    print('\nverdicts: %s' % ', '.join('%s %d' % kv for kv in sorted(counts.items())))
    ok = [r for r in rows if r['verdict'] == 'ring-ok']
    flare = [r for r in rows if r['verdict'] == 'flare-needed']
    print('ends that take the elbow transition unchanged: %d / %d'
          % (len(ok), len(rows)))
    for r in flare:
        print('  flare needed: %-8s %-9s inradius %s -> growth %s mm'
              % (r['id'], r['end'],
                 '-' if r['inradius_mm'] is None else '%.2f' % r['inradius_mm'],
                 '-' if r['radial_growth_to_ring_mm'] is None
                 else '%.2f' % r['radial_growth_to_ring_mm']))

    modules = [
        {'module': '三角肌（Y 形，3 通道）', 'members': ['FJ1468', 'FJ1513', 'FJ1467'],
         'caps': '3 个起点端（各 1 根钢丝）+ 1 个合并止点端（3 根钢丝）'},
        {'module': '大圆肌', 'members': ['FJ1507'], 'caps': '2 个端'},
        {'module': '肩胛下肌', 'members': ['FJ1504'], 'caps': '2 个端'},
        {'module': '冈下肌', 'members': ['FJ1500'], 'caps': '2 个端'},
    ]
    for m in modules:
        ends = [r for r in rows if r['id'] in m['members']]
        worst = min(ends, key=lambda r: (1e9 if r['inradius_mm'] is None
                                         else r['inradius_mm']))
        m['tightest'] = {'id': worst['id'], 'end': worst['end'],
                         'inradius_mm': worst['inradius_mm'],
                         'radial_growth_to_ring_mm': worst['radial_growth_to_ring_mm'],
                         'verdict': worst['verdict']}
        m['needs_flare'] = ['%s/%s' % (r['id'], r['end'])
                            for r in ends if r['verdict'] != 'ring-ok']
    print('\n--- per module ---')
    for m in modules:
        print('  %-18s 最紧 %s/%s inradius %s -> %s；需加座: %s'
              % (m['module'], m['tightest']['id'], m['tightest']['end'],
                 m['tightest']['inradius_mm'], m['tightest']['verdict'],
                 ', '.join(m['needs_flare']) or '无'))

    report = {
        'scope': 'end sections of the anatomical bellies, perpendicular to the '
                 'origin-insertion chord; elbow end form read from the build script; '
                 'no pad, rod, cap, bolt or flare geometry actually modelled',
        'elbow_cap_mm': ELBOW_CAP,
        'elbow_pad_mm': ELBOW_PAD,
        'thresholds_mm': {'transition_ring_r': RING_R_MM,
                          'cap_two_bolt_min_r': round(MIN_RADIUS_TWO_BOLT_MM, 2)},
        'inradius_step_mm': muscle_sections.INRADIUS_STEP_MM,
        'ends': rows,
        'verdict_counts': counts,
        'ring_ok_count': len(ok),
        'flare_needed_count': len(flare),
        'modules': modules,
        'not_verified': ['inradius is a 1 mm grid search, under-reports by up to ~0.7 mm',
                         'the flare (stalk) shape, its bond to the belly and its '
                         'interference with bone are NOT modelled',
                         'cap/pad stack height is 8.15 mm but §5.5 only cleared 4 mm',
                         'origin-side bone pads have not been laid out at all',
                         'no clamp force, creep or crimp grip'],
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
