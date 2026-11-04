"""Screen an early upper jaw against the verified deltoid entry corridor."""
import numpy as np

from build_s1_teres_origin_mount import A, export, load_component, read, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    proof_file = folder / 'verification_checks.json'
    proof = read(proof_file)
    assert proof['assembly_sha256'] == sha(source)
    assert proof['conditional_bundle_placement_geometry_pass']
    sources = read(source)['sources']
    jaw = load_component(sources['collar:jaw_upper'])
    direction = np.array(proof['direction_source'])
    cube = A.I._manifold().Manifold.cube
    low, high = np.array(jaw.bounding_box()).reshape(2, 3)
    sweep_low = low + np.minimum(0, -direction * 80) - .1
    sweep_high = high + np.maximum(0, -direction * 80) + .1
    mask = cube(sweep_high - sweep_low).translate(sweep_low)
    offsets = np.unique(np.r_[np.arange(0, 80.001, 1), np.arange(0, 3.001, .1)])
    hits = []
    cutters = []
    overlaps = []
    cropped_parts = {}
    for key in proof['moving_keys']:
        moving = load_component(sources[key]) ^ mask
        cropped_parts[key] = moving
        for offset in offsets:
            overlap = jaw ^ moving.translate(direction * offset)
            volume = float(overlap.volume())
            if volume <= .001:
                continue
            hits.append(dict(moving=key, offset_mm=float(offset), volume_mm3=volume))
            overlaps.append(overlap)
            lower, upper = np.array(overlap.bounding_box()).reshape(2, 3)
            cutters.append(cube(upper - lower + .1).translate(lower - .05))
    trial = jaw
    for cutter in cutters:
        trial -= cutter
    removed = jaw - trial
    exact_cut = A.I._manifold().Manifold.batch_boolean(overlaps, A.I._manifold().OpType.Add)
    exact_trial = jaw - exact_cut
    manifold = A.I._manifold().Manifold
    swept_cutters = []
    for hit, overlap in zip(hits, overlaps):
        half_step = .06 if hit['offset_mm'] < 3 else .51
        shift = direction * half_step
        swept_cutters.append(manifold.batch_hull([
            overlap.translate(-shift), overlap.translate(shift),
        ]))
    swept_trial = jaw - manifold.batch_boolean(swept_cutters, A.I._manifold().OpType.Add)
    refinement_cutters = []
    swept_hits = []
    for offset in np.unique(np.r_[np.arange(0, 80.001, .25), np.arange(0, 3.001, .05)]):
        for key, moving in cropped_parts.items():
            volume = float(A.I.intersection_volume(swept_trial, moving.translate(direction * offset)))
            if volume > .001:
                swept_hits.append(dict(moving=key, offset_mm=float(offset), volume_mm3=volume))
                overlap = swept_trial ^ moving.translate(direction * offset)
                refinement_cutters.append(manifold.batch_hull([
                    overlap.translate(-direction * .13), overlap.translate(direction * .13),
                ]))
    islands = sorted(swept_trial.decompose(), key=lambda part: part.volume(), reverse=True)
    island_relief_boxes = []
    for island in islands[1:]:
        if island.volume() <= 1e-8:
            continue
        lower, upper = np.array(island.bounding_box()).reshape(2, 3)
        island_relief_boxes.append([lower.tolist(), upper.tolist()])
        refinement_cutters.append(cube(upper - lower + .02).translate(lower - .01))
    refined_trial = swept_trial - manifold.batch_boolean(refinement_cutters, A.I._manifold().OpType.Add)
    candidate_metadata = None
    if sum(p.volume() > 1e-8 for p in refined_trial.decompose()) == 1:
        refined_trial, candidate_metadata = export(
            refined_trial, folder / 'collar_first_relief_candidate.obj', (low + high) / 2,
        )
    refined_removed = jaw - refined_trial
    refined_hits = []
    for offset in np.unique(np.r_[np.arange(0, 80.001, .125), np.arange(0, 3.001, .025)]):
        for key, moving in cropped_parts.items():
            volume = float(A.I.intersection_volume(refined_trial, moving.translate(direction * offset)))
            if volume > .001:
                refined_hits.append(dict(moving=key, offset_mm=float(offset), volume_mm3=volume))
    midpoint_hits = []
    for offset in (offsets[:-1] + offsets[1:]) / 2:
        for key, moving in cropped_parts.items():
            volume = float(A.I.intersection_volume(exact_trial, moving.translate(direction * offset)))
            if volume > .001:
                midpoint_hits.append(dict(moving=key, offset_mm=float(offset), volume_mm3=volume))
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    opening = np.array(read(frame_file)['frame_source']['basis_columns'])[:, 1]
    washer_contacts = {}
    for index in range(2):
        key = f'collar:top_washer_{index}'
        probe = load_component(sources[key]).translate(-opening * .2)
        washer_contacts[key] = {
            'original_mm3': float(A.I.intersection_volume(jaw, probe)),
            'trial_mm3': float(A.I.intersection_volume(exact_trial, probe)),
            'swept_trial_mm3': float(A.I.intersection_volume(swept_trial, probe)),
            'refined_trial_mm3': float(A.I.intersection_volume(refined_trial, probe)),
        }
    protected = {}
    exact_protected = {}
    swept_protected = {}
    for key in ('collar:bolt_0', 'collar:bolt_1', 'link:FJ1500/bone_bolt', 'collar:liner_upper'):
        hardware = load_component(sources[key])
        lower, upper = np.array(hardware.bounding_box()).reshape(2, 3)
        guard = cube(upper - lower + 2).translate(lower - 1)
        protected[key] = float(A.I.intersection_volume(removed, guard))
        exact_protected[key] = float(A.I.intersection_volume(exact_cut, guard))
        swept_protected[key] = float(A.I.intersection_volume(refined_removed, guard))
    result = dict(assembly_sha256=sha(source), placement_proof_sha256=sha(proof_file),
                  direction=direction.tolist(), offsets_mm=offsets.tolist(), hits=hits,
                  trial_removed_mm3=float(removed.volume()),
                  trial_positive_components=sum(p.volume() > 1e-8 for p in trial.decompose()),
                  removed_within_expanded_hardware_boxes_mm3=protected,
                  exact_sample_union_mm3=float(exact_cut.volume()),
                  exact_trial_positive_components=sum(p.volume() > 1e-8 for p in exact_trial.decompose()),
                  exact_trial_component_volumes_mm3=sorted(
                      (float(p.volume()) for p in exact_trial.decompose() if abs(p.volume()) > 1e-8),
                      reverse=True),
                  midpoint_hits=midpoint_hits,
                  swept_removed_mm3=float((jaw - swept_trial).volume()),
                  refined_removed_mm3=float(refined_removed.volume()),
                  candidate_mesh=candidate_metadata,
                  swept_positive_component_volumes_mm3=sorted(
                      (float(p.volume()) for p in swept_trial.decompose() if p.volume() > 1e-8),
                      reverse=True),
                  swept_fine_grid_hits=swept_hits,
                  refined_hits=refined_hits,
                  refined_positive_components=sum(p.volume() > 1e-8 for p in refined_trial.decompose()),
                  island_relief_boxes_mm=island_relief_boxes,
                  refined_removed_within_expanded_hardware_boxes_mm3=swept_protected,
                  washer_contact_probe_shift_mm=.2,
                  washer_contact_probes=washer_contacts,
                  exact_union_within_expanded_hardware_boxes_mm3=exact_protected,
                  geometry_adopted=False, whole_shoulder_assembly_verified=False)
    save(folder / 'collar_first_screen.json', result)
    print('hits', len(hits), 'peak', max((h['volume_mm3'] for h in hits), default=0),
          'removed', result['trial_removed_mm3'], 'exact', exact_cut.volume(), exact_protected)
    print(washer_contacts, result['exact_trial_component_volumes_mm3'])
    print('midpoint hits', len(midpoint_hits), 'peak', max((h['volume_mm3'] for h in midpoint_hits), default=0))
    print('swept', result['swept_removed_mm3'], result['swept_positive_component_volumes_mm3'],
          'hits', len(swept_hits), 'peak', max((h['volume_mm3'] for h in swept_hits), default=0))
    print('refined', result['refined_positive_components'], 'hits', len(refined_hits),
          'peak', max((h['volume_mm3'] for h in refined_hits), default=0), swept_protected)


if __name__ == '__main__':
    main()
