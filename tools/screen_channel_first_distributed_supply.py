"""Compare backboard-edge supply starts while retaining accepted neutral routes."""
from itertools import combinations, product

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save
from build_channel_first_inventory import provenance


def sample_curve(controls, count):
    t = np.linspace(0., 1., count)[:, None]
    a, b, c, d = controls
    path = (1-t)**3*a+3*(1-t)**2*t*b+3*(1-t)*t**2*c+t**3*d
    velocity = 3*(1-t)**2*(b-a)+6*(1-t)*t*(c-b)+3*t**2*(d-c)
    acceleration = 6*(1-t)*(c-2*b+a)+6*t*(d-2*c+b)
    speed = np.linalg.norm(velocity, axis=1)
    if speed.min() < 1e-8:
        return path, 0.
    curvature = np.linalg.norm(np.cross(velocity, acceleration), axis=1)/speed**3
    return path, float(1./max(curvature.max(), 1e-12))


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    previous_file = folder/'supply_channel_expanded_screen.json'
    inventory_file = folder/'route_inventory.json'
    previous, inventory = read(previous_file), read(inventory_file)
    assert previous['route_inventory'] == provenance(inventory_file)
    assembly_file = ROOT/inventory['original_elbow_assembly']['file']
    assert previous['assembly'] == provenance(assembly_file)
    solids = {key: load_component(source) for key, source in read(assembly_file)['sources'].items()
              if not key.startswith('tpu:')}
    radius = previous['channel_outer_radius_mm']
    routes = {row['id']: row for row in inventory['routes']}
    accepted = {}
    selected = {}
    for route in previous['routes']:
        candidate = route['selected_candidate']
        if candidate and candidate['refined_check']['passed']:
            path, bend_radius = sample_curve(np.asarray(candidate['controls_source_mm']), 801)
            accepted[route['id']] = A.D.outer_envelope(path, radius)
            selected[route['id']] = dict(inherited=True, controls_source_mm=candidate['controls_source_mm'],
                                         min_sampled_bend_radius_mm=bend_radius)
    lo, hi = np.array(solids['bench:backboard'].bounding_box()).reshape(2, 3)
    anchors = [(np.array([lo[0]-20., (lo[1]+hi[1])/2, height]), np.array([-1., 0., 0.]), 'lateral_edge')
               for height in np.linspace(lo[2]+15., hi[2]-10., 6)]
    anchors += [(np.array([x, (lo[1]+hi[1])/2, hi[2]+12.]), np.array([0., 0., 1.]), 'top_edge')
                for x in np.linspace(lo[0]+15., hi[0]-15., 8)]
    trials = {}
    for head, route in routes.items():
        if head in accepted:
            continue
        local = np.asarray(route['local_route']['points_source_mm'])
        end, incoming = local[0], A.D.unit(local[1]-local[0])
        proposals = []
        for (start, normal, edge), reach, approach, launch in product(
                anchors, (20., 40., 70., 110., 160.), (20., 40., 70., 110., 160.), ('normal', 'toward_end')):
            initial = normal if launch == 'normal' else A.D.unit(end-start)
            controls = np.array([start, start+initial*reach, end-incoming*approach, end])
            path, bend_radius = sample_curve(controls, 81)
            if bend_radius >= 8.:
                proposals.append(dict(controls_source_mm=controls.tolist(), anchor_edge=edge,
                                      initial_handle_mm=reach, inlet_handle_mm=approach, launch=launch,
                                      coarse_length_mm=float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum())))
        proposals.sort(key=lambda row: row['coarse_length_mm'])
        trials[head] = []
        for proposal in proposals[:80]:
            controls = np.asarray(proposal['controls_source_mm'])
            path, bend_radius = sample_curve(controls, 81)
            row = dict(**proposal, min_sampled_bend_radius_mm=bend_radius)
            trials[head].append(row)
            channel = A.D.outer_envelope(path, radius)
            hits = {key: volume for key, solid in solids.items()
                    if (volume := A.I.intersection_volume(channel, solid)) > .001}
            row['obstacle_hits_mm3'] = hits
            if hits:
                continue
            gaps = {key: A.I.gap(channel, solid, 2.) for key, solid in accepted.items()}
            row['existing_supply_gap_mm_capped2'] = gaps
            if min(gaps.values(), default=2.) < .5:
                continue
            path, bend_radius = sample_curve(controls, 801)
            row['refined_min_bend_radius_mm'] = bend_radius
            if bend_radius < 8.:
                row['rejection'] = 'refined_bend_radius'
                continue
            channel = A.D.outer_envelope(path, radius)
            fine_hits = {key: volume for key, solid in solids.items()
                         if (volume := A.I.intersection_volume(channel, solid)) > .001}
            fine_gaps = {key: A.I.gap(channel, solid, 2.) for key, solid in accepted.items()}
            row['refined_obstacle_hits_mm3'] = fine_hits
            row['refined_supply_gap_mm_capped2'] = fine_gaps
            if fine_hits or min(fine_gaps.values(), default=2.) < .5:
                continue
            accepted[head] = channel
            selected[head] = row
            break
        print(head, 'numerical proposals', len(proposals), 'tested', len(trials[head]),
              'selected', head in selected, flush=True)
    pair_checks = [dict(left=left, right=right, gap_mm_capped2=A.I.gap(accepted[left], accepted[right], 2.),
                       overlap_mm3=A.I.intersection_volume(accepted[left], accepted[right]))
                   for left, right in combinations(accepted, 2)]
    save(folder/'supply_channel_distributed_screen.json', dict(
        inherited_screen=provenance(previous_file), route_inventory=provenance(inventory_file),
        assembly=provenance(assembly_file), channel_outer_radius_mm=radius,
        minimum_bend_radius_assumption_mm=8., minimum_supply_pair_gap_mm=.5,
        refined_samples=801, shortlist_limit_per_unresolved_route=80,
        selected=selected, trials=trials, supply_pair_checks=pair_checks,
        scope='Neutral rest-shape geometry only. New top/side support brackets, actuator envelopes, '
              'moving TPU compliance and reaction force capacity require design and verification.'))


if __name__ == '__main__':
    main()
