"""Screen shoulder-side relief with the original elbow anatomy held fixed."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, save, sha


def main():
    base = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    baseline_file = base/'anatomical_elbow_baseline.json'
    baseline = read(baseline_file)
    assembly_file = base/'assembly_check.json'
    assert baseline['original_assembly_sha256'] == sha(assembly_file)
    sources = read(assembly_file)['sources']
    deltoid = load_component(sources['tpu:deltoid'])
    extensor = load_component(baseline['muscles']['extensor']['collision_proxy'])
    overlap = deltoid ^ extensor
    route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes = [row for row in read(route_file)['routes'] if row['head'] in ('FJ1468', 'FJ1513', 'FJ1467')]
    mf = A.I._manifold()
    bounds = np.array(overlap.bounding_box()).reshape(2, 3)
    lo, hi = bounds
    crop = mf.Manifold.cube(hi-lo+2.).translate(lo-1.)
    local_extensor = extensor ^ crop
    axis = int(np.argmax(hi-lo))
    guards = {f"{route['head']}_r{radius}": A.D.outer_envelope(
        np.array(route['full_wire_path_source_mm']), radius)
        for route in routes for radius in (.3, 2.)}
    candidates = [('exact_overlap', 0., local_extensor)]
    margin = .25
    for slices in (1, 3, 6):
        cutters = []
        edges = np.linspace(lo[axis], hi[axis], slices+1)
        for lower, upper in zip(edges, edges[1:]):
            start, stop = lo-1., hi+1.
            start[axis], stop[axis] = lower-.15, upper+.15
            band = overlap ^ mf.Manifold.cube(stop-start).translate(start)
            if not band.is_empty():
                cutters.append(band.hull().minkowski_sum(mf.Manifold.cube([2*margin]*3, True)))
        candidates.append((f'convex_bands_{slices}', margin,
                           mf.Manifold.batch_boolean(cutters, mf.OpType.Add)))
    rows = []
    for name, margin, cutter in candidates:
        revised = (deltoid-cutter).as_original().simplify(1e-6)
        removed = deltoid-revised
        remaining = float(A.I.intersection_volume(revised, extensor))
        guard_loss = {key: float(A.I.intersection_volume(removed, guard)) for key, guard in guards.items()}
        volumes = sorted((float(part.volume()) for part in revised.decompose()), reverse=True)
        row = dict(method=name, cutter_padding_mm=margin, removed_deltoid_mm3=float(removed.volume()),
                         clearance_capped_1mm=float(A.I.gap(revised, extensor, 1.)),
                         remaining_overlap_mm3=remaining, deltoid_component_volumes_mm3=volumes,
                         positive_material_components=sum(volume > 1e-8 for volume in volumes),
                         wire_region_removed_mm3=guard_loss)
        rows.append(row)
        if (row['positive_material_components'] == 1 and row['clearance_capped_1mm'] > .02
                and remaining <= .001 and max(guard_loss.values()) <= .001):
            try:
                actual, row['candidate_source'] = export(
                    revised, base/f'deltoid_original_elbow_{name}.obj', (lo+hi)/2)
            except (ValueError, AssertionError) as error:
                row['export_error'] = str(error)
                continue
            row['readback_clearance_mm_capped1'] = float(A.I.gap(actual, extensor, 1.))
            row['readback_overlap_mm3'] = float(A.I.intersection_volume(actual, extensor))
            comparison = []
            for key, source in sources.items():
                if key.startswith('tpu:elbow_'):
                    muscle = key.removeprefix('tpu:elbow_')
                    source = baseline['muscles'][muscle]['collision_proxy']
                elif not key.startswith('elbow:') and key != 'bone:forearm':
                    continue
                part = load_component(source)
                before = float(A.I.intersection_volume(deltoid, part))
                after = float(A.I.intersection_volume(actual, part))
                comparison.append(dict(part=key, before_mm3=before, after_mm3=after))
            row['elbow_inventory_comparison'] = comparison
    save(base/'anatomical_elbow_shoulder_clearance.json', dict(
        baseline_sha256=sha(baseline_file), route_sha256=sha(route_file),
        overlap_mm3=float(overlap.volume()), overlap_bounds_source_mm=bounds.tolist(),
        trials=rows, adopted=False,
        scope='Neutral shoulder-only relief trials. Elbow display geometry, anchors and wire remain the baseline.'))
    print([(row['method'], row['removed_deltoid_mm3'], row['clearance_capped_1mm'],
            len(row['deltoid_component_volumes_mm3']), max(row['wire_region_removed_mm3'].values())) for row in rows])


if __name__ == '__main__':
    main()
