"""Can the deltoid posterior origin take a one-sided clamp?  Search an offset cap for it.

D22 sent this seat to option B because moving it would cost 45.1% of its abduction moment arm
(5.9.16).  B means: no bolt through the 2.8 mm blade, and the TPU end is pressed from the outer
side only.

5.9.14 established that no coaxial cap fits at any (thickness, radius) - the cap's rim reaches the
infraspinatus, which wraps around that part of the scapula.  So this searches a cap that is
*offset laterally* as well: a plate whose centre is displaced in the pad plane, in eight
directions, at several distances, with the usual thickness and radius ladder.

If a combination clears, a one-sided clamp exists and its parameters are the answer.  If none does,
the seat cannot be clamped at all and the muscle end itself has to change - which is a different
decision and has to be said out loud rather than absorbed into a tolerance.

Run:
    python mujoco_arm/tools/fit_one_sided_clamp.py
    python mujoco_arm/tools/fit_one_sided_clamp.py --json out.json
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
import search_shoulder_channel_paths as SP

MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'
ENDPADS_JSON = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_end_pads.json'
SEAT = 'FJ1513/origin'


def mf():
    return I._manifold()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--seat', default=SEAT)
    args = ap.parse_args(argv)

    design = json.loads(ENDPADS_JSON.read_text(encoding='utf-8'))
    seat = next(s for s in design['seats'] if ','.join(s['heads']) == args.seat)
    p = np.asarray(seat.get('seat_source_mm', seat['centre_source_mm']), float)
    n = np.asarray(seat['normal_source'], float)

    bones = {k: SP.load_solid(BE.BONE_DIR / v) for k, v in SP.BONES.items()}
    modules = {}
    for key in ('deltoid', 'teres_major', 'subscapularis', 'infraspinatus'):
        f = MODULES_DIR / ('%s.stl' % key)
        if f.is_file():
            V, F = mesh_geom.read_stl(str(f))
            W, G = mesh_geom.weld(V, F)
            modules[key] = mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))

    own = seat['module']
    # Do NOT seed the search from the seat's current fitted pad: this tool's output is written back
    # into that same JSON by build_shoulder_ends, so reading it made the search depend on its own
    # previous result and the run stopped reproducing.  Everything below comes from the blade
    # thickness alone.
    print('seat %s at %s, n %s' % (args.seat, np.round(p, 1), np.round(n, 2)))
    print('recorded (ignored as input): pad r%s t%s, overlap %s mm3'
          % (seat['pad_radius_mm'], seat['pad_thickness_mm'],
             seat['pad_other_module_mm3']))
    lift = 0.0

    # in-plane basis for the offset
    B = BE.shoulder_pads.pad_basis(n)
    e1, e2 = B[:, 0], B[:, 1]

    def hits(geom):
        bone = max(float(I.intersection_volume(geom, b)) for b in bones.values())
        mod = max((float(I.intersection_volume(geom, s)) for k, s in modules.items()
                   if k != own), default=0.0)
        return bone, mod

    found = None
    tried = 0
    # The pad overlaps by 8.33 mm3 on its own, so an offset cap alone cannot fix this seat: the
    # offset has to move the whole clamp (pad and cap together).  Small offsets are searched first
    # because the offset is itself an attachment move - D22 rejected 17.39 mm at 45.1% of the
    # abduction moment arm, so anything comparable is out.
    blade = seat['blade_thickness_mm']
    t_max = BE.PAD_T if blade is None else float(max(1.0, min(BE.PAD_T, blade - 0.5)))
    for r_off in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
        for k in range(8):
            a = 2.0 * np.pi * k / 8.0
            off = e1 * (r_off * np.cos(a)) + e2 * (r_off * np.sin(a))
            for pad_t in sorted({t_max, round(t_max * 0.7, 2), 1.5, 1.0}, reverse=True):
                for pad_r in (BE.PAD_R, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0):
                    tried += 1
                    c = p + off
                    pad = BE.disc_rounded(c, n, lift, pad_t, pad_r)
                    pb, pm = hits(pad)
                    if pb > BE.SURFACE_TOL_MM3 or pm > BE.SURFACE_TOL_MM3:
                        continue
                    for cap_t in sorted({pad_t, round(pad_t * 0.6, 2), 1.5, 1.0}, reverse=True):
                        for cap_r in (BE.PAD_R, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0):
                            tried += 1
                            cap = BE.disc_rounded(c, n, lift + pad_t + 0.2, cap_t, cap_r)
                            cb, cm = hits(cap)
                            if cb <= BE.SURFACE_TOL_MM3 and cm <= BE.SURFACE_TOL_MM3:
                                found = {'offset_mm': round(float(r_off), 1),
                                         'offset_dir_deg': round(float(np.degrees(a)), 0),
                                         'offset_source_mm': np.round(off, 2).tolist(),
                                         'pad_radius_mm': round(pad_r, 1),
                                         'pad_thickness_mm': round(pad_t, 2),
                                         'cap_radius_mm': round(cap_r, 1),
                                         'cap_thickness_mm': round(cap_t, 2),
                                         'clamp_centre_source_mm': np.round(c, 2).tolist()}
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if found:
            break

    print('\ntried %d cap configurations' % tried)
    if found:
        print('a one-sided clamp exists: pad r%.1f t%.2f + cap r%.1f t%.2f, offset %.1f mm at %.0f deg'
              % (found['pad_radius_mm'], found['pad_thickness_mm'],
                 found['cap_radius_mm'], found['cap_thickness_mm'],
                 found['offset_mm'], found['offset_dir_deg']))
        print('  clamp centre %s' % found['clamp_centre_source_mm'])
    else:
        print('NO cap configuration clears: this seat cannot take a cap at all')

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'search for a laterally offset cap at the deltoid posterior origin, the seat '
                     'D22 sent to option B (one-sided clamp)',
            'seat': args.seat, 'tried': tried, 'found': found,
            'pad': {'radius_mm': seat['pad_radius_mm'], 'thickness_mm': pad_t,
                    'other_module_mm3': seat['pad_other_module_mm3']},
            'not_verified': ['bosses, screws and the crimp are not modelled',
                             'the TPU end thickness the clamp must trap is not modelled',
                             'no clamping force'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
