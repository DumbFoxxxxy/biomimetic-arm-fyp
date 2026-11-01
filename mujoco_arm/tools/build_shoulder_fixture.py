"""D17 v3 fixture: two clamp groups, no structure across the subscapular fossa.

5.9.8 established the rule the hard way - a frame that spans the subscapular fossa crushes the
subscapularis by 6076.6 mm3, and removing the beams leaves 6.2 mm3, so the beams were the whole
problem.  v3 therefore splits into a medial group (medial border mid + superior angle, joined by a
short beam along the border) and a lateral group (lateral border upper), each on its own base.

5.9.8 also listed what was still missing, and this tool does those three things:

* the outward routing for each group, checked against the humerus SWEEP envelope (not just the
  zero pose - the first frame version died at 37.9 deg of abduction because a vertical rod lay
  across the humerus);
* the base itself;
* the lift at each pad, reported so the wedge gap can be filled.

Run:
    python mujoco_arm/tools/build_shoulder_fixture.py
    python mujoco_arm/tools/build_shoulder_fixture.py --json out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import check_shoulder_fixture_zones as FZ
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP
import shoulder_pads

MODULES_DIR = frames.OUTPUT_DIR / 'shoulder_modules'

POST_R_MM = 4.0
BEAM_R_MM = 4.0
BASE_SIZE = (40.0, 30.0, 6.0)
SWEEP_STEP = 30.0
SWEEP_MAX = 120.0
#: a few hundredths of a mm3 is surface curvature, not a collision
TOL_MM3 = 0.5


def mf():
    return I._manifold()


def rod(a, b, r):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = b - a
    L = float(np.linalg.norm(d))
    if L < 1e-9:
        return None
    cyl = mf().Manifold.cylinder(L, r, r, 32)
    z = d / L
    helper = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    x = np.cross(helper, z)
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    return cyl.transform(np.c_[np.stack([x, y, z], axis=1), a])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    gh = frames.joint_centres('source')['GH']
    axes = frames.gh_probe_axes()
    scap = SP.load_solid(FZ.BONE_DIR / SP.BONES['scapula'])
    clav = SP.load_solid(FZ.BONE_DIR / SP.BONES['clavicle'])
    humerus = SP.load_solid(FZ.BONE_DIR / SP.BONES['humerus'])
    SV, SF = mesh_geom.load_obj(str(FZ.BONE_DIR / SP.BONES['scapula']))
    SW, SG = mesh_geom.weld(SV, SF)

    zones = FZ.zones(SW, gh)
    want = {'medial_mid', 'superior_angle', 'lateral_upper'}
    picks = {k: p for k, _, p, _ in zones if k in want}
    print('clamp points')
    for k in ('medial_mid', 'superior_angle', 'lateral_upper'):
        print('  %-15s %s' % (k, np.round(picks[k], 1)))

    # humerus sweep envelope: flexion and abduction, both signs
    poses = []
    for name in ('flexion', 'abduction'):
        for sign in (+1, -1):
            for a in np.arange(0.0, SWEEP_MAX + 1e-9, SWEEP_STEP):
                T = I.rotation_about(gh, axes[name], sign * float(a))
                poses.append(humerus.transform(T))
    envelope = mf().Manifold.batch_boolean(poses, mf().OpType.Add)
    print('\nhumerus sweep envelope: %d poses, %.0f mm3' % (len(poses), envelope.volume()))

    bones = {'scapula': scap, 'clavicle': clav}
    modules = {}
    for key in ('deltoid', 'teres_major', 'subscapularis', 'infraspinatus'):
        p = MODULES_DIR / ('%s.stl' % key)
        if p.is_file():
            V, F = mesh_geom.read_stl(str(p))
            W, G = mesh_geom.weld(V, F)
            modules[key] = mf().Manifold(mf().Mesh64(W, G.astype(np.uint32)))

    def normal_at(p):
        picked = None
        for tri in (SW[SG],):
            d, n = mesh_geom.signed_distance_and_normal(p, tri)
            if picked is None or abs(d) < abs(picked[0]):
                picked = (d, n)
        return picked[1]

    parts = {}
    rows = []
    for k in ('medial_mid', 'superior_angle', 'lateral_upper'):
        p = np.asarray(picks[k], float)
        n = normal_at(p)
        lift, _ = shoulder_pads.lift_to_clear(p, n, scap)
        if lift is None:
            lift = 0.0
        pad = shoulder_pads.build_pad(p + n * lift, n)
        parts[k + '_pad'] = pad
        rows.append({'zone': k, 'centre_source_mm': np.round(p, 2).tolist(),
                     'lift_mm': round(lift, 2),
                     'pad_bone_mm3': round(I.intersection_volume(pad, scap), 3)})
        print('  pad %-15s lift %.2f mm, bone %.3f mm3' % (k, lift, rows[-1]['pad_bone_mm3']))

    # beams: one short beam along the medial border, joining the two medial pads
    pm, ps = np.asarray(picks['medial_mid'], float), np.asarray(picks['superior_angle'], float)
    beam = rod(pm, ps, BEAM_R_MM)
    if beam is not None:
        parts['medial_beam'] = beam

    # outward routing: medial group leaves along +x, lateral group along -x and down
    base_med = pm + np.array([55.0, 0.0, -10.0])
    base_lat = np.asarray(picks['lateral_upper'], float) + np.array([-45.0, 0.0, -25.0])
    posts = {}
    for name, start, end in (('medial_post', ps, base_med), ('lateral_post',
                                                             np.asarray(picks['lateral_upper'],
                                                                        float), base_lat)):
        r = rod(start, end, POST_R_MM)
        if r is not None:
            posts[name] = r
            parts[name] = r

    for name, c in (('medial_base', base_med), ('lateral_base', base_lat)):
        box = mf().Manifold.cube(list(BASE_SIZE)).translate(list(c))
        parts[name] = box

    print('\n--- the three things 5.9.8 listed as missing ---')
    out = {}
    for name, geom in parts.items():
        hit_bone = max(I.intersection_volume(geom, b).__float__() for b in bones.values())
        hit_hum = I.intersection_volume(geom, envelope)
        hit_mod = max((I.intersection_volume(geom, m).__float__()
                       for m in modules.values()), default=0.0)
        out[name] = {'bone_mm3': round(hit_bone, 3),
                     'humerus_sweep_mm3': round(hit_hum, 3),
                     'module_mm3': round(hit_mod, 3)}
        print('  %-16s bone %8.2f   humerus sweep %10.2f   module %10.2f'
              % (name, hit_bone, hit_hum, hit_mod))

    sweep_bad = [k for k, v in out.items() if v['humerus_sweep_mm3'] > TOL_MM3]
    bone_bad = [k for k, v in out.items() if v['bone_mm3'] > TOL_MM3]
    mod_bad = [k for k, v in out.items() if v['module_mm3'] > TOL_MM3]
    print('\nparts hit by the humerus sweep: %d %s' % (len(sweep_bad), sweep_bad))
    print('parts hitting bone: %d %s' % (len(bone_bad), bone_bad))
    print('parts hitting a module: %d %s' % (len(mod_bad), mod_bad))

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'D17 v3 fixture: 3 pads (medial mid, superior angle, lateral upper), a short '
                     'medial beam, one outward post per group, and two bases; checked against the '
                     'humerus sweep envelope and the modules',
            'post_r_mm': POST_R_MM, 'beam_r_mm': BEAM_R_MM, 'base_mm': list(BASE_SIZE),
            'sweep': {'axes': ['flexion', 'abduction'], 'step_deg': SWEEP_STEP,
                      'max_deg': SWEEP_MAX, 'poses': len(poses),
                      'envelope_mm3': round(envelope.volume(), 1)},
            'pads': rows, 'parts': out,
            'sweep_hits': sweep_bad, 'bone_hits': bone_bad, 'module_hits': mod_bad,
            'not_verified': ['no jaws, screws, rails or clamp body',
                             'the wedge gap under each lifted pad is not filled',
                             'single-axis sweeps only, no combined poses',
                             'clamp force and stiffness are not modelled'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
