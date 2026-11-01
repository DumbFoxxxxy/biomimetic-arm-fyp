"""S1 muscle envelopes: can the anatomical bellies be used as the TPU outer shape?

The scheme-two premise is "copy the anatomical outer shape, engineer the ends and
the channels".  That premise only holds if the anatomical volumes actually fit, so
this script measures them against the four things that can reject them:

1. **bones** - a belly that already intersects humerus/scapula/clavicle in the
   neutral pose has to be trimmed before it can be a printed part;
2. **the joint sphere** - radius = fitted head + measured cartilage gap; anything
   a belly puts inside it is squeezed by the articulation, not by a bone;
3. **each other** - four modules (D4) means four separate solids;
4. **the steel channel** - the chord from origin to insertion must run *inside*
   the belly with enough material around it to carry a channel.

The last one is the reason this script also computes a signed distance along the
chord: it reports the minimum wall available, which is what sets the channel
radius.  It does NOT infer a fibre direction (PCA on a muscle mesh says nothing
about pennation), and a chord is not a tendon path.

Run:
    python mujoco_arm/tools/check_shoulder_muscle_envelopes.py
    python mujoco_arm/tools/check_shoulder_muscle_envelopes.py --json out.json
"""
import argparse
import json
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
BONE_DIR = frames.PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输入' / '骨骼'

MUSCLES = {
    'FJ1468': '三角肌前束',
    'FJ1513': '三角肌后束',
    'FJ1467': '三角肌中束',
    'FJ1507': '大圆肌',
    'FJ1504': '肩胛下肌',
    'FJ1500': '冈下肌',
}
DELTOID = ('FJ1468', 'FJ1513', 'FJ1467')
BONES = {'humerus': '肱骨.obj', 'scapula': '肩胛骨.obj', 'clavicle': '锁骨.obj'}

CHORD_STEP_MM = 1.0
#: wire channel bore used by the elbow design [mm] (§5.7: R0.85 inner / R1.6 outer)
CHANNEL_INNER_R_MM = 0.85
CHANNEL_OUTER_R_MM = 1.6


def solid_and_mesh(path):
    V, F = mesh_geom.load_obj(str(path))
    W, G = mesh_geom.weld(V, F)
    mf = I._manifold()
    return mf.Manifold(mf.Mesh64(W, G.astype(np.uint32))), W, G


def section_profile(solid, origin, insertion, stations=13):
    return muscle_sections.section_profile(solid, origin, insertion, stations=stations)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args(argv)

    gh = frames.joint_centres('source')['GH']
    jc = json.loads(frames.CANDIDATES_PATH.read_text(encoding='utf-8'))
    head = float(jc['GH']['radius_mm'])
    scaps, _ = mesh_geom.load_obj(str(BONE_DIR / BONES['scapula']))
    cartilage = float(np.linalg.norm(scaps - gh, axis=1).min()) - head
    joint_r = head + cartilage
    print('GH %s  fitted head %.2f mm  cartilage gap %.2f mm  joint sphere r=%.2f mm'
          % (np.round(gh, 2), head, cartilage, joint_r))

    bones = {}
    for name, fname in BONES.items():
        solid, _, _ = solid_and_mesh(BONE_DIR / fname)
        bones[name] = solid
    print('bones: %s' % ', '.join('%s %.0f mm3' % (k, v.volume())
                                  for k, v in bones.items()))

    mf = I._manifold()
    joint_sphere = mf.Manifold.sphere(joint_r, 96).translate(list(gh))

    mrows, solids = [], {}
    print('\n--- envelopes: cross-sections perpendicular to the O-I chord ---')
    print('%-8s %-10s %8s %7s %6s %7s %7s %7s %8s %8s %7s'
          % ('muscle', 'name', 'vol,cm3', 'O-I mm', 'sec', 'wallmin', 'bellymin',
             'offmax', 'areamax', 'trim%', 'joint'))
    for mid, name in MUSCLES.items():
        solid, W, G = solid_and_mesh(MUSCLE_DIR / (mid + '.obj'))
        solids[mid] = solid
        origin = attachment_regions.load_region(mid, 'origin')
        insertion = attachment_regions.load_region(mid, 'insertion')
        o, _ = attachment_regions.centroid(origin)
        i, _ = attachment_regions.centroid(insertion)
        prof = section_profile(solid, o, i)
        hit_bones, bone_vol = {}, 0.0
        for bname, bsolid in bones.items():
            v = I.intersection_volume(solid, bsolid)
            if v > I.HARD_MM3:
                hit_bones[bname] = round(v, 2)
                bone_vol += v
        joint_vol = I.intersection_volume(solid, joint_sphere)
        row = {'id': mid, 'name': name,
               'volume_mm3': round(solid.volume(), 1),
               'triangles': int(len(G)),
               'bbox_source_mm': [np.round(solid.bounding_box()[:3], 1).tolist(),
                                  np.round(solid.bounding_box()[3:], 1).tolist()],
               'origin_bone': origin['bone'], 'insertion_bone': insertion['bone'],
               'origin_source_mm': o.tolist(), 'insertion_source_mm': i.tolist(),
               'sections': prof,
               'bone_intersections_mm3': hit_bones,
               'bone_intersection_total_mm3': round(bone_vol, 2),
               'joint_sphere_intrusion_mm3': round(joint_vol, 2)}
        row['bone_intersection_pct'] = round(100.0 * bone_vol / max(solid.volume(), 1e-9), 2)
        row['net_volume_mm3'] = round(solid.volume() - bone_vol, 1)
        mrows.append(row)
        print('%-8s %-10s %8.1f %7.1f %3d/%-2d %7.2f %7s %7.2f %8.1f %8.1f %7.1f'
              % (mid, name, solid.volume() / 1000.0, prof['length_mm'],
                 prof.get('sections_inside', 0), prof.get('sections_total', 0),
                 prof.get('wall_min_mm', -1.0),
                 '%.2f' % prof['wall_min_belly_mm']
                 if prof.get('wall_min_belly_mm') is not None else '-',
                 prof.get('chord_offset_max_mm', -1.0),
                 prof.get('area_max_mm2', -1.0), row['bone_intersection_pct'],
                 joint_vol))

    print('\n--- bone intersections (neutral pose) ---')
    for r in mrows:
        print('  %-8s trimmed %5.1f%% (%7.1f of %7.1f mm3)   %s'
              % (r['id'], r['bone_intersection_pct'],
                 r['bone_intersection_total_mm3'], r['volume_mm3'],
                 ', '.join('%s %.1f mm3' % kv
                           for kv in r['bone_intersections_mm3'].items())
                 or 'none'))

    print('\n--- muscle vs muscle ---')
    pairs = []
    ids = [r['id'] for r in mrows]
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            ka, kb = ids[a], ids[b]
            vol = I.intersection_volume(solids[ka], solids[kb])
            gap = I.gap(solids[ka], solids[kb]) if vol <= I.HARD_MM3 else 0.0
            flag = 'OVERLAP' if vol > I.HARD_MM3 else ('touch' if gap <= 0.05 else 'clear')
            pairs.append({'a': ka, 'b': kb, 'overlap_mm3': round(vol, 3),
                          'gap_mm': round(gap, 2), 'flag': flag,
                          'deltoid_pair': bool(ka in DELTOID and kb in DELTOID)})
            if vol > I.HARD_MM3 or gap < 1.0:
                print('  %-8s %-8s overlap %10.3f mm3  gap %.2f mm  %s'
                      % (ka, kb, vol, gap, flag))
    n_ov = sum(1 for p in pairs if p['overlap_mm3'] > I.HARD_MM3)
    print('  pairs overlapping: %d / %d' % (n_ov, len(pairs)))

    delt_pairs = [p for p in pairs if p['deltoid_pair']]
    print('\n--- deltoid Y-merge input (D8) ---')
    for p in delt_pairs:
        print('  %-8s %-8s overlap %10.3f mm3  gap %.2f mm'
              % (p['a'], p['b'], p['overlap_mm3'], p['gap_mm']))
    n_delt_ov = sum(1 for p in delt_pairs if p['overlap_mm3'] > I.HARD_MM3)
    print('  deltoid head pairs already merged in the source data: %d / %d'
          % (n_delt_ov, len(delt_pairs)))

    report = {
        'scope': 'anatomical muscle meshes in the source neutral pose, welded; '
                 'envelope only - no TPU lattice, no end cap, no channel geometry',
        'gh_centre_source_mm': gh.tolist(),
        'fitted_head_radius_mm': head,
        'cartilage_gap_mm': round(cartilage, 3),
        'joint_sphere_radius_mm': round(joint_r, 3),
        'channel': {'inner_r_mm': CHANNEL_INNER_R_MM,
                    'outer_r_mm': CHANNEL_OUTER_R_MM},
        'muscles': mrows,
        'pairs': pairs,
        'overlapping_pairs': n_ov,
        'deltoid_pairs_overlapping': n_delt_ov,
        'not_verified': ['candidate regions are proximity samples, not confirmed '
                         'attachment sites',
                         'no muscle force, pennation, friction or contact',
                         'chord is not a tendon path and PCA is not a fibre direction',
                         'neutral pose only; swept poses belong to S1-G1'],
    }
    if args.json_out:
        I.write_report(args.json_out, report)
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
