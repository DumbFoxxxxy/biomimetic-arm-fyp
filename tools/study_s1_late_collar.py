"""Screen upper collar closure and fastener insertion after deltoid placement."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    proof_file = folder / 'verification_checks.json'
    model, proof = read(source), read(proof_file)
    assert proof['assembly_sha256'] == sha(source)
    assert proof['conditional_bundle_placement_geometry_pass']
    solids = {key: load_component(row) for key, row in model['sources'].items()}
    fixed = set(proof['fixed_keys'] + proof['moving_keys'])
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    basis = np.array(read(frame_file)['frame_source']['basis_columns'])
    opening = basis[:, 1]
    stages = []

    def screen(keys, direction, span):
        assert not fixed.intersection(keys)
        offsets = np.unique(np.r_[np.arange(0, span + .01, .5), np.arange(0, 3.01, .1)])
        for count, offset in enumerate(offsets, 1):
            hits = []
            for key in keys:
                posed = solids[key].translate(direction * offset)
                for target in sorted(fixed):
                    volume = float(A.I.intersection_volume(posed, solids[target]))
                    if volume > .001:
                        hits.append(dict(moving=key, target=target, volume_mm3=volume))
            if hits:
                return dict(keys=keys, direction=direction.tolist(), poses=count,
                            first_collision_mm=float(offset), hits=hits, clear=False)
        return dict(keys=keys, direction=direction.tolist(), poses=len(offsets),
                    first_collision_mm=None, hits=[], clear=True)

    upper = ['collar:jaw_upper', 'collar:liner_upper']
    original_fixed = sorted(fixed)
    nominal = screen(upper, opening, 30)
    anchor_parts = {key for key in fixed if key.startswith('link:FJ1500/bone_')}
    anchor_parts.add('link:FJ1500/tendon_wire')
    assert anchor_parts <= fixed
    fixed.difference_update(anchor_parts)
    without_anchor = screen(upper, opening, 30)
    original_liner = solids['collar:liner_upper']
    liner = original_liner
    for distance in np.arange(.1, 2.001, .1):
        liner = liner - solids['bone:humerus'].translate(-opening * distance)
    solids['collar:liner_upper'] = liner
    liner_removed = float((original_liner - liner).volume())
    liner_components = sum(part.volume() > 1e-8 for part in liner.decompose())
    assert liner_components == 1
    jaw = solids['collar:jaw_upper']
    deltoid = solids['tpu:deltoid']
    low, high = np.array(jaw.bounding_box()).reshape(2, 3)
    sweep_low = low + np.minimum(0, opening * 30) - .01
    sweep_high = high + np.maximum(0, opening * 30) + .01
    manifold = A.I._manifold().Manifold
    deltoid = deltoid ^ manifold.cube(sweep_high - sweep_low).translate(sweep_low)
    interference_samples = []
    relief_boxes = []
    exact_overlaps = []
    for offset in np.unique(np.r_[np.arange(0, 30.01, .5), np.arange(0, 3.01, .1)]):
        obstacle = deltoid.translate(-opening * offset)
        overlap = jaw ^ obstacle
        volume = float(overlap.volume())
        if volume > .001:
            interference_samples.append(dict(offset_mm=float(offset), volume_mm3=volume))
            exact_overlaps.append(overlap)
            box_low, box_high = np.array(overlap.bounding_box()).reshape(2, 3)
            relief_boxes.append(manifold.cube(box_high - box_low + .1).translate(box_low - .05))
    jaw_clearance_trial = jaw
    for box in relief_boxes:
        jaw_clearance_trial -= box
    jaw_relief = float((jaw - jaw_clearance_trial).volume())
    removed = jaw - jaw_clearance_trial
    protected_hardware = sorted(anchor_parts - {'link:FJ1500/tendon_wire'})
    protected_hardware += [f'collar:{part}_{index}' for index in range(2)
                           for part in ('bolt', 'top_washer')]
    nearby_material = {}
    exact_trial = jaw
    for overlap in exact_overlaps:
        exact_trial -= overlap
    exact_removed = jaw - exact_trial
    exact_nearby = {}
    for key in protected_hardware:
        box_low, box_high = np.array(solids[key].bounding_box()).reshape(2, 3)
        guard = manifold.cube(box_high - box_low + 2).translate(box_low - 1)
        nearby_material[key] = float(A.I.intersection_volume(removed, guard))
        exact_nearby[key] = float(A.I.intersection_volume(exact_removed, guard))
    save(folder / 'late_collar_exact_relief.json', dict(
        assembly_sha256=sha(source), placement_proof_sha256=sha(proof_file),
        script_sha256=sha(Path(__file__)), direction=opening.tolist(),
        method='Subtract intersections at sampled poses, with zero added clearance',
        removed_mm3=float(exact_removed.volume()),
        original_volume_mm3=float(jaw.volume()),
        positive_components=sum(part.volume() > 1e-8 for part in exact_trial.decompose()),
        removed_within_hardware_boxes_mm3=exact_nearby,
        samples=interference_samples, geometry_adopted=False))
    approaches = []
    for x, z in [(0, 0), (.25, 0), (-.25, 0), (0, .25), (0, -.25),
                 (.5, .25), (.5, -.25), (-.5, .25), (-.5, -.25)]:
        direction = A.D.unit(opening + x * basis[:, 0] + z * basis[:, 2])
        result = screen(upper, direction, 30)
        approaches.append(result)
        print(result, flush=True)
        if result['clear']:
            fixed.update(upper)
            for index in range(2):
                for keys, axis in [
                    ([f'collar:bolt_{index}', f'collar:top_washer_{index}'], opening),
                    ([f'collar:nut_{index}', f'collar:bottom_washer_{index}'], -opening),
                ]:
                    stage = screen(keys, axis, 20)
                    stages.append(stage)
                    if stage['clear']:
                        fixed.update(keys)
                    else:
                        break
                if not stages[-1]['clear']:
                    break
            break
    save(folder / 'late_collar_relief_screen.json', dict(
        scope=__doc__, assembly_sha256=sha(source), placement_proof_sha256=sha(proof_file),
        frame_sha256=sha(frame_file), script_sha256=sha(Path(__file__)),
        original_fixed=original_fixed, nominal_upper=nominal, upper_without_anchor=without_anchor,
        deferred_anchor_parts=sorted(anchor_parts), trial_liner_removed_mm3=liner_removed,
        trial_liner_exported=False,
        jaw_deltoid_interference_samples=interference_samples,
        unadopted_jaw_relief_volume_mm3=jaw_relief,
        removed_material_within_expanded_hardware_boxes_mm3=nearby_material,
        upper_approaches=approaches, fastener_stages=stages,
        installed_keys=sorted(fixed), remaining_keys=sorted(set(solids)-fixed),
        late_collar_geometry_pass=len(stages) == 4 and all(row['clear'] for row in stages),
        whole_shoulder_assembly_verified=False))


if __name__ == '__main__':
    main()
