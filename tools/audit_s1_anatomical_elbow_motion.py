"""Screen the original elbow rigid chain in the assembled S1 neutral shoulder."""
from collections import defaultdict

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def main():
    design = ROOT/'02_人体骨骼晶格肌肉方案/输出'
    folder = design/'integration_review/s1_elbow_shared_humerus'
    assembly_file = folder/'assembly_anatomical_elbow_angled_candidate.json'
    fastener_file = folder/'anatomical_elbow_angled_fastener.json'
    owner_file = (design/'anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair'
                  /'paired_swivel_support/captive_keepers/continuous_channel_candidate/assembly.json')
    assembly = read(assembly_file)
    assert assembly['fastener_report_sha256'] == sha(fastener_file)
    assert read(fastener_file)['neutral_screen_passed']
    moving_keys = {'bone:forearm'}
    moving_keys.update(f"elbow:{row['key']}" for row in read(owner_file)['instances']
                       if row['owner'] == 'elbow' or row['owner'].endswith('_distal_pivot'))
    sources = assembly['sources']
    # The combined radius/ulna solid replaces the original FJ3391 forearm.
    moving_keys.remove('elbow:FJ3391')
    assert 'elbow:FJ3391' not in sources
    assert moving_keys <= sources.keys()
    excluded = {'tpu:elbow_flexor', 'tpu:elbow_extensor'}
    moving = {key: load_component(sources[key]) for key in sorted(moving_keys)}
    fixed = {key: load_component(source) for key, source in sources.items()
             if key not in moving_keys | excluded}
    fixed_keys = list(fixed)
    bounds = np.array([solid.bounding_box() for solid in fixed.values()])
    centre, axis = A.frames.elbow_axis()
    threshold = .001
    gap_cap = .5
    exact_checks = 0
    near_pairs = {}

    def pose(angle):
        nonlocal exact_checks
        transform = A.I.rotation_about(centre, axis, angle)
        hits = {}
        for key, solid in moving.items():
            posed = solid.transform(transform)
            lo, hi = np.asarray(posed.bounding_box()).reshape(2, 3)
            separation = np.maximum(lo-bounds[:, 3:], bounds[:, :3]-hi)
            candidates = np.flatnonzero(np.all(separation < gap_cap, axis=1))
            for index in candidates:
                target = fixed_keys[index]
                pair = key, target
                gap = A.I.gap(posed, fixed[target], gap_cap)
                if gap < near_pairs.get(pair, (gap_cap, 0))[0]:
                    near_pairs[pair] = (gap, angle)
                if np.any(separation[index] >= 0.):
                    continue
                exact_checks += 1
                volume = float(A.I.intersection_volume(posed, fixed[target]))
                if volume > threshold:
                    hits[pair] = volume
        return hits

    by_pair = defaultdict(dict)
    for angle in range(121):
        for pair, volume in pose(angle).items():
            by_pair[pair][angle] = volume
    refined = {}
    rows = []
    for (moving_key, fixed_key), hits in by_pair.items():
        pair = moving_key, fixed_key
        first = min(hits)
        refined_first = float(first)
        if first:
            for tenth in range(1, 10):
                angle = round(first-1+tenth/10, 1)
                if angle not in refined:
                    refined[angle] = pose(angle)
                if pair in refined[angle]:
                    refined_first = angle
                    break
        peak_angle = max(hits, key=hits.get)
        kind = ('tpu_geometry_contact' if fixed_key.startswith('tpu:') else
                'wire_geometry_contact' if 'wire' in fixed_key else 'rigid')
        rows.append(dict(moving=moving_key, fixed=fixed_key, kind=kind,
                         first_sample_deg=refined_first, neutral_mm3=hits.get(0, 0.),
                         peak_mm3=hits[peak_angle], peak_sample_deg=peak_angle,
                         integer_samples_with_overlap=len(hits)))
    rows.sort(key=lambda row: (row['kind'], row['first_sample_deg'], -row['peak_mm3']))
    save(folder/'anatomical_elbow_rigid_motion.json', dict(
        assembly_sha256=sha(assembly_file), owner_reference_sha256=sha(owner_file),
        axis_source=axis.tolist(), centre_source_mm=centre.tolist(), moving_keys=list(moving),
        fixed_count=len(fixed), excluded_unposed_elbow_tpu=sorted(excluded),
        path_deg=[0, 120, 0], symmetric_return_uses_same_geometry=True,
        sample_step_deg=1, onset_refinement_step_deg=.1, exact_checks=exact_checks,
        gap_cap_mm=gap_cap, close_pairs=[
            dict(moving=key, fixed=target, min_gap_mm=gap, angle_deg=angle)
            for (key, target), (gap, angle) in sorted(near_pairs.items(), key=lambda item: item[1])],
        threshold_mm3=threshold, overlaps=rows,
        scope='Elbow flexion at neutral shoulder with distal swivels locked to their neutral angle. '
              'Moving versus fixed pairs only. Neutral overlaps are listed without assuming legal contact. '
              'No soft deformation, elbow cable path motion, coupled shoulder path or load verification.'))
    print(f'{len(moving)} moving, {len(fixed)} fixed, {exact_checks} Boolean checks, '
          f'{len(rows)} overlap pairs, {len(near_pairs)} pairs below {gap_cap} mm')
    for row in rows:
        print(row)


if __name__ == '__main__':
    main()
