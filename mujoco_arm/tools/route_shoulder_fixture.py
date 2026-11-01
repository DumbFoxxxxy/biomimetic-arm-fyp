"""Route the D17 v3 fixture exits with an obstacle-avoiding planner instead of straight lines.

5.9.9 measured what straight lines cost: the medial beam cuts 866.7 mm3 into the blade because the
medial border is curved and a chord crosses inside it, and the lateral base landed inside the
humerus sweep (364.0 mm3) and on the deltoid (982.9 mm3).

The channel relaxation (SP.relax) cannot be reused as-is: it pushes points *into* a belly while
staying off bone, and a fixture lives in free space with no belly.  So this uses the opposite
constraint - pure repulsion from obstacles - in two forms:

    beam   densify between the two medial pads, then push each interior point out to the beam
           radius from the scapula, which makes the beam follow the border instead of chording it
    post   a greedy step planner: from the pad's outer face, try headings in a fan around the
           current one and take the first step that keeps the post radius clear of the bones AND
           of the humerus sweep envelope; stop when far enough out, and put the base there

The base position therefore falls out of the route, which is what 5.9.9 asked for.

Run:
    python mujoco_arm/tools/route_shoulder_fixture.py
    python mujoco_arm/tools/route_shoulder_fixture.py --json out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_shoulder_fixture as FX
import check_shoulder_fixture_zones as FZ
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP
import shoulder_pads

CLEAR_MM = 4.5          # beam clearance: beam radius 4.0 plus a little
#: The posts need their own, smaller clearance.  Sharing one number made the two requirements
#: fight: at 4.5 the lateral post cannot leave its pad at all, and at 3.5 the beam sinks back
#: 17.21 mm3 into the blade it had just been routed around.  The lateral zone is simply tighter.
POST_CLEAR_MM = 3.5
STEP_MM = 5.0
MAX_POST_MM = 95.0
FAN_DEG = 70.0
FAN_STEP_DEG = 10.0
BASE_SIZE = (40.0, 30.0, 6.0)
TOL_MM3 = 0.5


def mf():
    return I._manifold()


def dist_to(point, tri):
    return mesh_geom.signed_distance(np.asarray(point, float), tri)


def clear_of(point, obstacles, clear):
    """Smallest margin: (distance - required clearance) over all obstacles, and which one."""
    worst, which = None, None
    for name, tri, need in obstacles:
        m = dist_to(point, tri) - (clear if need is None else need)
        if worst is None or m < worst:
            worst, which = m, name
    return worst, which


def push_clear(path, obstacles, iterations=60, step=0.8, lam=0.35):
    """Push interior points out to clearance, then smooth - the opposite of SP.relax."""
    pts = [np.asarray(p, float) for p in path]
    for _ in range(iterations):
        moved = list(pts)
        for k in range(1, len(pts) - 1):
            m, _ = clear_of(pts[k], obstacles, CLEAR_MM)
            if m < 0.0:
                best = None
                for name, tri, need in obstacles:
                    d, n = mesh_geom.signed_distance_and_normal(pts[k], tri)
                    need_ = CLEAR_MM if need is None else need
                    if d - need_ < 0.0 and (best is None or d - need_ < best[0]):
                        best = (d - need_, n)
                if best is not None:
                    moved[k] = pts[k] + best[1] * (-best[0]) * step
        for k in range(1, len(pts) - 1):
            moved[k] = moved[k] + lam * (moved[k - 1] + moved[k + 1] - 2.0 * moved[k])
        pts = moved
        if all(clear_of(p, obstacles, CLEAR_MM)[0] >= 0.0 for p in pts[1:-1]):
            break
    return pts


def plan_post(start, heading, obstacles, step=STEP_MM, max_len=MAX_POST_MM):
    """Greedy fan search: step while staying clear; return the path or (None, reason)."""
    path = [np.asarray(start, float)]
    h = frames.unit(np.asarray(heading, float))
    total = 0.0
    perp1 = frames.unit(np.cross(h, np.array([0.0, 0.0, 1.0])))
    perp2 = frames.unit(np.cross(h, perp1))
    while total < max_len:
        candidates = []
        for a in np.arange(-FAN_DEG, FAN_DEG + 1e-9, FAN_STEP_DEG):
            rad = np.radians(a)
            for b in np.arange(-FAN_DEG, FAN_DEG + 1e-9, FAN_STEP_DEG * 2):
                r2 = np.radians(b)
                cand = (np.cos(rad) * h
                        + np.sin(rad) * (np.cos(r2) * perp1 + np.sin(r2) * perp2))
                cand = frames.unit(cand)
                nxt = path[-1] + cand * step
                # Sample along the whole step, not just its ends: a 5 mm step whose endpoints
                # satisfy 4.5 mm can still graze a bone edge with the side of a radius-4 post,
                # which is what left 1.1 and 13.6 mm3 of contact in the first routed version.
                m = min(clear_of(path[-1] + cand * (f * step), obstacles, POST_CLEAR_MM)[0]
                        for f in (0.2, 0.4, 0.6, 0.8, 1.0))
                candidates.append((m, float(a) ** 2 + float(b) ** 2, cand, nxt))
        candidates.sort(key=lambda c: (-c[0], c[1]))
        best = candidates[0]
        if best[0] < 0.0:
            return path, 'blocked after %.0f mm' % total
        _, _, h, nxt = best
        path.append(nxt)
        total += step
        if all(clear_of(p, obstacles, POST_CLEAR_MM * 2.0)[0] >= 0.0 for p in [nxt]) \
                and total >= 45.0:
            return path, 'clear at %.0f mm' % total
    return path, 'reached max %.0f mm' % total


def tube(path, r):
    pieces = []
    for a, b in zip(path, path[1:]):
        s = FX.rod(a, b, r)
        if s is not None:
            pieces.append(s)
    for p in path:
        pieces.append(mf().Manifold.sphere(r, 24).translate(list(np.asarray(p, float))))
    return mf().Manifold.batch_boolean(pieces, mf().OpType.Add)


def main(argv=None):
    global CLEAR_MM
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', dest='json_out')
    ap.add_argument('--export', metavar='DIR',
                    help='write the beam, both posts and both bases as STL into this directory')
    ap.add_argument('--clear-mm', type=float, default=CLEAR_MM,
                    help='clearance the route must keep from every obstacle [mm]')
    ap.add_argument('--post-r-mm', type=float, default=FX.POST_R_MM,
                    help='post radius [mm]; the lateral zone may need a slimmer post')
    args = ap.parse_args(argv)
    CLEAR_MM = args.clear_mm

    gh = frames.joint_centres('source')['GH']
    axes = frames.gh_probe_axes()
    scap = SP.load_solid(FZ.BONE_DIR / SP.BONES['scapula'])
    clav = SP.load_solid(FZ.BONE_DIR / SP.BONES['clavicle'])
    humerus = SP.load_solid(FZ.BONE_DIR / SP.BONES['humerus'])
    SV, SF = mesh_geom.load_obj(str(FZ.BONE_DIR / SP.BONES['scapula']))
    SW, SG = mesh_geom.weld(SV, SF)
    scap_tri = SW[SG]
    clav_tri = (lambda a: a[0][a[1]])(mesh_geom.weld(*mesh_geom.load_obj(
        str(FZ.BONE_DIR / SP.BONES['clavicle']))))
    hum_tri = (lambda a: a[0][a[1]])(mesh_geom.weld(*mesh_geom.load_obj(
        str(FZ.BONE_DIR / SP.BONES['humerus']))))

    zones = FZ.zones(SW, gh)
    picks = {k: p for k, _, p, _ in zones
             if k in ('medial_mid', 'superior_angle', 'lateral_upper')}

    poses = []
    for name in ('flexion', 'abduction'):
        for sign in (+1, -1):
            for a in np.arange(0.0, FX.SWEEP_MAX + 1e-9, FX.SWEEP_STEP):
                T = I.rotation_about(gh, axes[name], sign * float(a))
                poses.append(humerus.transform(T))
    envelope = mf().Manifold.batch_boolean(poses, mf().OpType.Add)
    EV, EF = mesh_geom.manifold_mesh_to_arrays(envelope.to_mesh())
    EW, EG = mesh_geom.weld(EV, EF)
    env_tri = EW[EG]
    print('sweep envelope %.0f mm3, %d triangles' % (envelope.volume(), len(env_tri)))

    obstacles = [('scapula', scap_tri, None), ('clavicle', clav_tri, None),
                 ('humerus_sweep', env_tri, None)]

    pm = np.asarray(picks['medial_mid'], float)
    ps = np.asarray(picks['superior_angle'], float)
    pl = np.asarray(picks['lateral_upper'], float)

    # --- medial beam: follow the border's own vertices instead of chording them ---
    # Pushing a chord outwards (the first version) still cut 203.5 mm3 into the blade, because
    # the blade curves away faster than the push converges.  Taking the border's actual surface
    # points and offsetting each along its normal follows the curve by construction.
    x = SW[:, 0]
    band = SW[(x >= np.percentile(x, 97))
              & (SW[:, 2] >= min(pm[2], ps[2]) - 2.0)
              & (SW[:, 2] <= max(pm[2], ps[2]) + 2.0)]
    border = band[np.argsort(band[:, 2])] if len(band) > 4 else np.vstack([pm, ps])
    # thin out so the beam does not get a vertex every 0.2 mm
    keep = [border[0]]
    for p in border[1:]:
        if np.linalg.norm(p - keep[-1]) >= 3.0:
            keep.append(p)
    if np.linalg.norm(keep[-1] - ps) > 3.0:
        keep.append(ps)
    beam_path = []
    for p in keep:
        nrm = mesh_geom.signed_distance_and_normal(p, scap_tri)[1]
        beam_path.append(p + nrm * (FX.BEAM_R_MM + 0.6))
    beam_path = push_clear(beam_path, obstacles, iterations=40)
    beam = tube(beam_path, FX.BEAM_R_MM)
    # judge the solid, not the sampled points: the point check is conservative and flags a
    # clear beam because a vertex sits near the border it is following
    beam_bone = float(max(I.intersection_volume(beam, b).__float__() for b in (scap, clav)))
    beam_sweep = float(I.intersection_volume(beam, envelope))
    beam_ok = beam_bone <= TOL_MM3 and beam_sweep <= TOL_MM3
    print('\nbeam: %d points (follows the border), bone %.2f mm3, sweep %.2f mm3 -> %s'
          % (len(beam_path), beam_bone, beam_sweep, 'CLEAR' if beam_ok else 'HIT'))

    # --- pads first: the posts leave from the pads' outer faces ---
    pads_ok = []
    starts = {}
    for k, p in (('medial_mid', pm), ('superior_angle', ps), ('lateral_upper', pl)):
        n = mesh_geom.signed_distance_and_normal(p, scap_tri)[1]
        lift, _ = shoulder_pads.lift_to_clear(p, n, scap)
        lift = 0.0 if lift is None else lift
        pad = shoulder_pads.build_pad(p + n * lift, n)
        pads_ok.append({'zone': k, 'lift_mm': round(lift, 2),
                        'pad_bone_mm3': round(float(I.intersection_volume(pad, scap)), 3)})
        # the post leaves from the pad's OUTER face: starting on the bone surface leaves no
        # room for a 4 mm post, which is why the lateral route was blocked at step zero
        starts[k] = p + n * (lift + shoulder_pads.PAD_THICKNESS_MM)
    print('\npads: %s' % ', '.join('%s lift %.2f' % (p['zone'], p['lift_mm']) for p in pads_ok))

    # --- posts: greedy route, base at the end ---
    results = {}
    fixture_solids = {'medial_beam': beam}
    for tag, key, heading in (('medial_post', 'superior_angle', np.array([1.0, 0.0, -0.2])),
                              ('lateral_post', 'lateral_upper', np.array([-1.0, 0.0, -0.5]))):
        path, why = plan_post(starts[key], heading, obstacles)
        post = tube(path, args.post_r_mm) if len(path) > 1 else None
        base_c = path[-1]
        base = mf().Manifold.cube(list(BASE_SIZE)).translate(list(base_c))
        fixture_solids[tag] = post
        fixture_solids[tag.replace('post', 'base')] = base
        results[tag] = {
            'why': why, 'points': len(path),
            'length_mm': round(float(sum(np.linalg.norm(b - a)
                                         for a, b in zip(path, path[1:]))), 1),
            'base_centre': np.round(base_c, 1).tolist(),
            'post_bone_mm3': round(float(max(I.intersection_volume(post, b).__float__()
                                             for b in (scap, clav)) if post else 0.0), 3),
            'post_sweep_mm3': round(float(I.intersection_volume(post, envelope).__float__()
                                          if post else 0.0), 3),
            'base_sweep_mm3': round(float(I.intersection_volume(base, envelope)), 3),
            'base_bone_mm3': round(float(max(I.intersection_volume(base, b).__float__()
                                             for b in (scap, clav))), 3),
        }
        r = results[tag]
        print('%-14s %s  len %.0f mm  post/bone %.2f  post/sweep %.2f  base/sweep %.2f  '
              'base/bone %.2f'
              % (tag, why, r['length_mm'], r['post_bone_mm3'], r['post_sweep_mm3'],
                 r['base_sweep_mm3'], r['base_bone_mm3']))

    ok = (beam_ok and all(results[t]['points'] > 1
                          and v['post_bone_mm3'] <= TOL_MM3
                          and v['post_sweep_mm3'] <= TOL_MM3
                          and v['base_sweep_mm3'] <= TOL_MM3
                          and v['base_bone_mm3'] <= TOL_MM3
                          for t, v in results.items()))
    print('fixture exits all clear: %s' % ('YES' if ok else 'NO'))

    # Export every fixture solid as STL so the parts exist as files, not just as routes.
    exported = {}
    if args.export:
        dest = Path(args.export)
        dest.mkdir(parents=True, exist_ok=True)
        for tag, geom in fixture_solids.items():
            if geom is None:
                continue
            V, F = mesh_geom.manifold_mesh_to_arrays(geom.to_mesh())
            target = dest / ('%s.stl' % tag)
            mesh_geom.write_stl(str(target), V, F, header='S1 fixture %s' % tag)
            RV, RF = mesh_geom.read_stl(str(target))
            RW, RG = mesh_geom.weld(RV, RF)
            e, bnd, nm = mesh_geom.edge_census(RG)
            exported[tag] = {'file': target.name, 'volume_mm3': round(geom.volume(), 1),
                             'triangles': int(len(F)),
                             'boundary_edges': int(bnd), 'nonmanifold_edges': int(nm)}
            print('  exported %-14s %9.1f mm3  %6d tris  boundary %d  non-manifold %d'
                  % (tag, geom.volume(), len(F), bnd, nm))

    if args.json_out:
        I.write_report(args.json_out, {
            'scope': 'D17 v3 fixture exits routed with obstacle avoidance: beam follows the '
                     'medial border, posts planned by greedy fan search against scapula, clavicle '
                     'and the humerus sweep envelope, base placed at the route end',
            'beam_clearance_mm': CLEAR_MM, 'post_clearance_mm': POST_CLEAR_MM,
            'post_r_mm': args.post_r_mm, 'step_mm': STEP_MM, 'max_post_mm': MAX_POST_MM,
            'sweep_envelope_mm3': round(envelope.volume(), 1),
            'beam': {'clear': bool(beam_ok), 'points': len(beam_path),
                     'bone_mm3': round(beam_bone, 3), 'sweep_mm3': round(beam_sweep, 3)},
            'posts': results, 'pads': pads_ok, 'exits_all_clear': bool(ok),
            'exported': exported,
            'not_verified': ['no jaws, screws or rails', 'wedge gaps not filled',
                             'single-axis sweep union, not combined poses'],
        })
        print('\nreport written:', args.json_out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
