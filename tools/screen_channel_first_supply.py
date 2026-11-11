"""Screen supported-end supply channel candidates from the existing backboard."""
import argparse
from itertools import combinations

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save
from build_channel_first_inventory import provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expanded', action='store_true')
    args = parser.parse_args()
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    inventory_file = folder/'route_inventory.json'
    inventory = read(inventory_file)
    assembly_file = ROOT/inventory['original_elbow_assembly']['file']
    assert provenance(assembly_file) == inventory['original_elbow_assembly']
    sources = read(assembly_file)['sources']
    rigid = {key: load_component(source) for key, source in sources.items()
             if not key.startswith('tpu:')}
    board = rigid['bench:backboard']
    lo, hi = np.asarray(board.bounding_box()).reshape(2, 3)
    radius = 1.6
    parameters = np.linspace(0., 1., 81)
    rows = []
    for index, route in enumerate(inventory['routes']):
        local = np.asarray(route['local_route']['points_source_mm'])
        end = local[0]
        direction = A.D.unit(local[1]-local[0])
        start = np.array([lo[0]+25.+index*22., (lo[1]+hi[1])/2, hi[2]+12.])
        trials = []
        candidate = None
        lifts = (25., 50., 80., 120., 180.) if args.expanded else (25., 50., 80., 120.)
        approaches = (40., 60., 90., 130., 180.) if args.expanded else (8., 16., 30.)
        for lift in lifts:
            for approach in approaches:
                controls = np.array([start, start+[0., 0., lift], end-direction*approach, end])
                t = parameters[:, None]
                path = ((1-t)**3*controls[0]+3*(1-t)**2*t*controls[1]
                        +3*(1-t)*t**2*controls[2]+t**3*controls[3])
                derivative = (3*(1-t)**2*(controls[1]-controls[0])
                              +6*(1-t)*t*(controls[2]-controls[1])
                              +3*t**2*(controls[3]-controls[2]))
                second = (6*(1-t)*(controls[2]-2*controls[1]+controls[0])
                          +6*t*(controls[3]-2*controls[2]+controls[1]))
                curvature = np.linalg.norm(np.cross(derivative, second), axis=1)/np.linalg.norm(derivative, axis=1)**3
                bend_radius = float(1./max(curvature))
                if args.expanded and bend_radius < 8.:
                    trials.append(dict(lift_mm=lift, approach_handle_mm=approach,
                                       min_sampled_bend_radius_mm=bend_radius,
                                       rejection='bend_radius', rigid_hits_mm3=None))
                    continue
                channel = A.D.outer_envelope(path, radius)
                hits = {key: volume for key, solid in rigid.items()
                        if (volume := A.I.intersection_volume(channel, solid)) > .001}
                row = dict(lift_mm=lift, approach_handle_mm=approach, rigid_hits_mm3=hits,
                           obstacle_kinds={key: ('wire_reference_or_tendon' if 'wire' in key else 'rigid')
                                           for key in hits},
                           min_sampled_bend_radius_mm=bend_radius,
                           length_mm=float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum()))
                trials.append(row)
                if not hits and row['min_sampled_bend_radius_mm'] >= 8.:
                    candidate = dict(controls_source_mm=controls.tolist(), points_source_mm=path.tolist(), **row)
                    break
            if candidate:
                break
        rows.append(dict(id=route['id'], supply_start_source_mm=start.tolist(),
                         supply_end_source_mm=end.tolist(), trials=trials, selected_candidate=candidate,
                         mechanical_support=dict(
                             start='new rail/post attached to existing backboard top; mounting hardware pending',
                             end=route['proximal_conduit_reaction'],
                             intermediate='continuous TPU/Bowden conduit with its own bending stiffness; no fixed via sites'),
                         candidate_is_historical_local_inlet=route['local_route_status']
                         if 'local_route_status' in route else False))
        print(route['id'], 'clear candidate', candidate is not None, flush=True)
    refined_channels = {}
    if args.expanded:
        for route in rows:
            candidate = route['selected_candidate']
            if candidate is None:
                continue
            controls = np.asarray(candidate['controls_source_mm'])
            t = np.linspace(0., 1., 801)[:, None]
            path = ((1-t)**3*controls[0]+3*(1-t)**2*t*controls[1]
                    +3*(1-t)*t**2*controls[2]+t**3*controls[3])
            channel = A.D.outer_envelope(path, radius)
            derivative = (3*(1-t)**2*(controls[1]-controls[0])
                          +6*(1-t)*t*(controls[2]-controls[1])+3*t**2*(controls[3]-controls[2]))
            second = (6*(1-t)*(controls[2]-2*controls[1]+controls[0])
                      +6*t*(controls[3]-2*controls[2]+controls[1]))
            curvature = np.linalg.norm(np.cross(derivative, second), axis=1)/np.linalg.norm(derivative, axis=1)**3
            hits = {key: volume for key, solid in rigid.items()
                    if (volume := A.I.intersection_volume(channel, solid)) > .001}
            candidate['refined_check'] = dict(
                samples=801, min_sampled_bend_radius_mm=float(1./max(curvature)),
                obstacle_hits_mm3=hits,
                passed=not hits and float(1./max(curvature)) >= 8.)
            refined_channels[route['id']] = channel
    pair_checks = [dict(left=left, right=right,
                       overlap_mm3=A.I.intersection_volume(refined_channels[left], refined_channels[right]),
                       gap_mm_capped2=A.I.gap(refined_channels[left], refined_channels[right], 2.))
                   for left, right in combinations(refined_channels, 2)]
    output = 'supply_channel_expanded_screen.json' if args.expanded else 'supply_channel_screen.json'
    save(folder/output, dict(
        route_inventory=provenance(inventory_file), assembly=provenance(assembly_file),
        channel_outer_radius_mm=radius, min_candidate_bend_radius_mm=8.,
        cubic_samples=81, attachment_reference='bench:backboard top face', routes=rows,
        supply_pair_checks=pair_checks,
        scope='Neutral candidate rest shapes against non-TPU inventory, including recorded wire references. End supports require design. '
              'Pair checks cover only selected expanded candidates. No moving compliance, '
              'friction, motor dimensions or Y2 route adoption.'))
    print('refined checks', [(row['id'], row['selected_candidate']['refined_check']['passed'])
                             for row in rows if row['selected_candidate'] and args.expanded], pair_checks)


if __name__ == '__main__':
    main()
