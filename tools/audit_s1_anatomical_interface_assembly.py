"""Check the revised interfaces in the original anatomical elbow assembly."""
from itertools import combinations
import argparse

import numpy as np

from build_s1_teres_origin_mount import ROOT, A, export, load_component, read, save, sha


def relieve_deltoid_bolt(folder, sources, route_file):
    deltoid = load_component(sources['tpu:deltoid'])
    bolt = load_component(sources['elbow:FJ1477_bone_0_bolt'])
    mf = A.I._manifold()
    overlap = deltoid ^ bolt
    lo, hi = np.asarray(overlap.bounding_box()).reshape(2, 3)
    routes = [row for row in read(route_file)['routes'] if row['head'] in ('FJ1468', 'FJ1513', 'FJ1467')]
    assert len(routes) == 3
    guards = {row['head']: A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), 2.)
              for row in routes}
    guards.update({key: load_component(source) for key, source in sources.items()
                   if key.startswith('cap:') and key.split('/')[0][4:] in guards})
    trials = []
    for margin in (.20, .12, .08, .04):
        cutter = overlap.hull().minkowski_sum(mf.Manifold.cube([2*margin]*3, True))
        candidate = (deltoid-cutter).as_original().simplify(1e-6)
        removed = deltoid-candidate
        row = dict(cutter_padding_mm=margin, removed_mm3=float(removed.volume()),
                   protected_region_removed_mm3={
                       name: float(A.I.intersection_volume(removed, solid)) for name, solid in guards.items()})
        trials.append(row)
        if max(row['protected_region_removed_mm3'].values()) > .001:
            continue
        try:
            actual, source = export(candidate, folder/'deltoid_original_elbow_bolt_relief.obj', (lo+hi)/2)
        except (AssertionError, ValueError) as error:
            row['export_error'] = str(error) or type(error).__name__
            continue
        row['readback_bolt_gap_mm'] = float(A.I.gap(actual, bolt, 1.))
        row['readback_bolt_overlap_mm3'] = float(A.I.intersection_volume(actual, bolt))
        if row['readback_bolt_gap_mm'] > .02 and row['readback_bolt_overlap_mm3'] <= .001:
            return source, trials
    return None, trials


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shoulder-relief', action='store_true')
    parser.add_argument('--bolt-relief', action='store_true')
    args = parser.parse_args()
    if args.bolt_relief and not args.shoulder_relief:
        parser.error('--bolt-relief requires --shoulder-relief')
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    original_file = folder/'assembly_check.json'
    baseline_file = folder/'anatomical_elbow_baseline.json'
    relief_file = folder/'anatomical_elbow_shoulder_clearance.json'
    cap_file = folder/'anatomical_elbow_cap_relief.json'
    adapter_file = folder/'anatomical_elbow_adapter_bridge.json'
    flexor_file = folder/'anatomical_elbow_adapter_bridge_flexor.json'
    baseline, relief, cap, adapter = map(read, (baseline_file, relief_file, cap_file, adapter_file))
    assert baseline['original_assembly_sha256'] == sha(original_file)
    assert relief['baseline_sha256'] == cap['baseline_sha256'] == sha(baseline_file)
    assert cap['shoulder_clearance_sha256'] == sha(relief_file)
    assert adapter['cap_audit_sha256'] == sha(cap_file)
    flexor = read(flexor_file)
    assert flexor['cap_audit_sha256'] == sha(cap_file)
    original = read(original_file)['sources']
    matched = next(row['matched_interface']['sources'] for row in cap['adapter_ear_trials']
                   if 'sources' in row['matched_interface'])
    changed = {'tpu:deltoid': next(row['candidate_source'] for row in relief['trials']
                                  if 'candidate_source' in row),
               'elbow:FJ1477_proximal_rigid_cap': matched['cap'],
               'elbow:FJ1477_bone_adapter': adapter['main_stock_trial']['source'],
               'elbow:FJ1486_bone_adapter': flexor['main_stock_trial']['source']}
    output = folder/'anatomical_elbow_interface_assembly.json'
    reports = [relief_file, cap_file, adapter_file, flexor_file]
    if args.shoulder_relief:
        pocket_file = folder/'anatomical_elbow_adapter_shoulder_relief.json'
        pockets = read(pocket_file)
        assert pockets['interface_audit_sha256'] == sha(output)
        assert read(output)['changed_sources'] == changed
        for head, rows in pockets['trials'].items():
            changed[f'elbow:{head}_bone_adapter'] = next(
                row['source'] for row in rows if 'source' in row and row['positive_clearance'])
        reports.extend((output, pocket_file))
        output = folder/'anatomical_elbow_interface_with_shoulder_relief.json'
    bolt_trials = None
    if args.bolt_relief:
        route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
        assert sha(route_file) == relief['route_sha256']
        assert read(output)['changed_sources'] == changed
        reports.extend((output, route_file))
        output = folder/'anatomical_elbow_interface_with_bolt_relief.json'
        source, bolt_trials = relieve_deltoid_bolt(folder, {**original, **changed}, route_file)
        if source is None:
            save(output, dict(bolt_relief_trials=bolt_trials, selected=False))
            raise SystemExit('No bolt pocket passed the local checks')
        changed['tpu:deltoid'] = source
    proxies = {f'tpu:elbow_{label}': row['collision_proxy'] for label, row in baseline['muscles'].items()}
    # The old left fastener stack is removed for this provisional interface
    # check. Its angled replacement is still an unmodelled assembly item.
    removed = {key for key in original if key.startswith('elbow:FJ1477_proximal_')
               and key.endswith('_left')}
    inventory = {key: source for key, source in {**original, **proxies}.items() if key not in removed}
    old_solids = {key: load_component(source) for key, source in inventory.items()}
    new_solids = {**old_solids, **{key: load_component(row) for key, row in changed.items()}}
    rows = []
    tested = 0
    for left, right in combinations(new_solids, 2):
        if left not in changed and right not in changed:
            continue
        tested += 1
        before = float(A.I.intersection_volume(old_solids[left], old_solids[right]))
        after = float(A.I.intersection_volume(new_solids[left], new_solids[right]))
        if max(before, after) > .001:
            rows.append(dict(left=left, right=right, before_mm3=before, after_mm3=after))
    rows.sort(key=lambda row: -row['after_mm3'])
    save(output, dict(
        original_sha256=sha(original_file), baseline_sha256=sha(baseline_file),
        input_reports_sha256={p.name: sha(p) for p in reports},
        changed_sources=changed, omitted_original_left_fasteners=sorted(removed), tested_pairs=tested,
        threshold_mm3=.001, overlaps=rows, provisional=True,
        scope='Neutral pairs involving at least one changed solid; original elbow TPU collision proxies. Replacement left fastener and movement not covered.',
        **({'bolt_relief_trials': bolt_trials} if bolt_trials is not None else {})))
    print('checked pairs', tested, 'remaining', sum(row['after_mm3'] > .001 for row in rows),
          'new', sum(row['after_mm3'] > .001 >= row['before_mm3'] for row in rows))
    for row in rows[:8]:
        print(row)


if __name__ == '__main__':
    main()
