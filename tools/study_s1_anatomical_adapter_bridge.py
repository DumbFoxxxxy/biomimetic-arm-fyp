"""Reconnect relieved elbow adapter stock around protected muscle volumes."""
import argparse
from itertools import combinations
import numpy as np
from scipy.spatial import cKDTree

from build_s1_teres_origin_mount import A, ROOT, export, load_component, mesh_geom, read, rod, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flexor', action='store_true')
    flexor = parser.parse_args().flexor
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    cap_file = folder/'anatomical_elbow_cap_relief.json'
    cap_audit = read(cap_file)
    matched = next(row['matched_interface']['sources'] for row in cap_audit['adapter_ear_trials']
                   if 'sources' in row['matched_interface'])
    baseline = read(folder/'anatomical_elbow_baseline.json')
    sources = read(folder/'assembly_check.json')['sources']
    shoulder = read(folder/'anatomical_elbow_shoulder_clearance.json')
    deltoid = load_component(next(row['candidate_source'] for row in shoulder['trials'] if 'candidate_source' in row))
    head = 'FJ1486' if flexor else 'FJ1477'
    adapter = load_component(sources[f'elbow:{head}_bone_adapter'] if flexor else matched['adapter'])
    cap = load_component(sources[f'elbow:{head}_proximal_rigid_cap'] if flexor else matched['cap'])
    suffix = '_flexor' if flexor else ''
    obstacles = {'deltoid': deltoid, 'cap': cap}
    for label, item in baseline['muscles'].items():
        obstacles[label] = load_component(item['collision_proxy'])
        obstacles[f'{label}_wire'] = A.D.outer_envelope(np.array(item['wire_path_source_mm']), 2.)
    mf = A.I._manifold()
    clips = [adapter ^ deltoid.translate(direction*.12)
             for direction in [np.zeros(3), *np.eye(3), *(-np.eye(3))]]
    trimmed = (adapter-mf.Manifold.batch_boolean(clips, mf.OpType.Add)).as_original().simplify(1e-6)
    components = sorted((part for part in trimmed.decompose() if part.volume() > 1e-8),
                        key=lambda part: -part.volume())
    vertices = [mesh_geom.manifold_mesh_to_arrays(part.to_mesh64())[0] for part in components]
    rows, links = [], []
    for i, j in combinations(range(len(components)), 2):
        left, right = vertices[i], vertices[j]
        distance, nearest = cKDTree(right).query(left)
        passed = None
        for index in np.argsort(distance)[:8]:
            a, b = left[index], right[nearest[index]]
            axis = A.D.unit(b-a)
            u = A.D.unit(np.cross(axis, np.eye(3)[np.argmin(abs(axis))]))
            v = np.cross(axis, u)
            for reach in (4., 8., 12., 16.):
                for turn in np.arange(8)*np.pi/4:
                    mid = (a+b)/2+reach*(u*np.cos(turn)+v*np.sin(turn))
                    link = (rod(a, mid, 1.7)+rod(mid, b, 1.7)
                            +mf.Manifold.sphere(1.7, 16).translate(mid))
                    if any(A.I.intersection_volume(link, part) > .001 for part in obstacles.values()):
                        continue
                    if min(A.I.intersection_volume(link, components[k]) for k in (i, j)) <= .01:
                        continue
                    passed = dict(points_mm=[p.tolist() for p in (a, mid, b)], radius_mm=1.7)
                    links.append(link)
                    break
                if passed:
                    break
            if passed:
                break
        rows.append(dict(components=[i, j], nearest_gap_mm=float(min(distance)), link=passed))
    candidate = mf.Manifold.batch_boolean([trimmed, *links], mf.OpType.Add).as_original().simplify(1e-6)
    report = dict(cap_audit_sha256=sha(cap_file), trimmed_components_mm3=[part.volume() for part in components],
                  trials=rows, final_positive_components=sum(part.volume() > 1e-8 for part in candidate.decompose()),
                  remaining_deltoid_overlap_mm3=float(A.I.intersection_volume(candidate, deltoid)),
                  scope='Local neutral stock and muscle avoidance. Full fasteners, anatomy and motion pending.')
    # Separate the test of the main load-path candidate from the failed
    # attempt to retain every small fragment created by the Boolean cut.
    scraps = sum(float(part.volume()) for part in components[1:])
    if report['final_positive_components'] != 1:
        primary = components[0]
        interface_parts = {'proximal_cap': cap, 'humerus': load_component(sources['bone:humerus'])}
        for key, source in sources.items():
            if key.startswith((f'elbow:{head}_proximal_', f'elbow:{head}_bone_')) and any(
                    term in key for term in ('washer', 'bolt', 'nut')):
                interface_parts[key] = load_component(source)
        report['separated_stock_interface_checks'] = [
            dict(volume_mm3=float(part.volume()),
                 gaps_mm_capped2={key: float(A.I.gap(part, solid, 2.))
                                  for key, solid in interface_parts.items()},
                 overlaps_mm3={key: volume for key, solid in interface_parts.items()
                               if (volume := float(A.I.intersection_volume(part, solid))) > .001})
            for part in components[1:]]
        guard_loss = {}
        removed = adapter-primary
        for key, source in sources.items():
            if key.startswith(f'elbow:{head}_bone_') and any(term in key for term in ('washer', 'bolt', 'nut')):
                hardware = load_component(source)
                lo, hi = np.array(hardware.bounding_box()).reshape(2, 3)
                guard = mf.Manifold.cube(hi-lo+.5).translate(lo-.25)
                lost = float(A.I.intersection_volume(removed, guard))
                if lost > .001:
                    guard_loss[key] = lost
        report['main_stock_trial'] = dict(separated_stock_mm3=scraps, fragment_volume_limit_mm3=.3,
                                          hardware_guard_loss_mm3=guard_loss,
                                          residual_deltoid_overlap_mm3=float(
                                              A.I.intersection_volume(primary, deltoid)))
        isolated_from_interfaces = all(min(row['gaps_mm_capped2'].values()) >= 2.
                                       for row in report['separated_stock_interface_checks'])
        if ((scraps < .3 or isolated_from_interfaces) and not guard_loss
                and report['main_stock_trial']['residual_deltoid_overlap_mm3'] <= .001):
            origin = vertices[0][0]
            v, faces = mesh_geom.manifold_mesh_to_arrays(primary.translate(-origin).to_mesh64())
            welded, remap = np.unique(v, axis=0, return_inverse=True)
            indexed = remap[faces]
            edges = np.sort(np.vstack([indexed[:, [0, 1]], indexed[:, [1, 2]], indexed[:, [2, 0]]]), axis=1)
            unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
            bad = welded[unique_edges[counts > 2]]+origin
            trial = report['main_stock_trial']
            trial['removed_isolated_stock_mm3'] = scraps
            trial['stock_isolated_from_checked_interfaces'] = isolated_from_interfaces
            trial['nonmanifold_edges_after_exact_weld'] = bad.tolist()
            trial['repair_attempts'] = []
            for radius in (0., .005, .015, .03, -.005, -.015, -.03):
                pieces = [A.D.outer_envelope(edge, abs(radius)) for edge in bad] if radius else []
                operation = mf.OpType.Subtract if radius < 0 else mf.OpType.Add
                repaired = mf.Manifold.batch_boolean([primary, *pieces], operation)
                repaired = repaired.as_original().simplify(1e-6)
                added = repaired-primary
                item = dict(radius_mm=radius, added_mm3=max(0., float(added.volume())),
                            removed_mm3=max(0., float((primary-repaired).volume())),
                            added_obstacle_overlap_mm3={
                                key: float(A.I.intersection_volume(added, part)) for key, part in obstacles.items()})
                trial['repair_attempts'].append(item)
                if (item['added_mm3']+item['removed_mm3'] > .01
                        or max(item['added_obstacle_overlap_mm3'].values()) > .001):
                    continue
                try:
                    actual, trial['source'] = export(
                        repaired, folder/f'anatomical_elbow_trimmed_adapter{suffix}.obj', origin)
                except (AssertionError, ValueError) as error:
                    item['export_error'] = str(error) or type(error).__name__
                    continue
                trial['readback_deltoid_overlap_mm3'] = float(A.I.intersection_volume(actual, deltoid))
                break
    if report['final_positive_components'] == 1 and report['remaining_deltoid_overlap_mm3'] <= .001:
        try:
            _, report['source'] = export(
                candidate, folder/f'anatomical_elbow_rebridged_adapter{suffix}.obj', vertices[0][0])
        except (AssertionError, ValueError) as error:
            report['export_error'] = str(error)
    save(folder/f'anatomical_elbow_adapter_bridge{suffix}.json', report)
    print(report)


if __name__ == '__main__':
    main()
