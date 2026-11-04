"""Compare hole ligaments and screen fasteners after deltoid placement."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    proof_file = folder / 'verification_checks.json'
    study_file = folder / 'collar_first_screen.json'
    proof, study = read(proof_file), read(study_file)
    assert study['assembly_sha256'] == proof['assembly_sha256'] == sha(source)
    assert study['placement_proof_sha256'] == sha(proof_file) and not study['refined_hits']
    solids = {key: load_component(row) for key, row in read(source)['sources'].items()}
    original = solids['collar:jaw_upper']
    candidate = load_component(study['candidate_mesh'])
    solids['collar:jaw_upper'] = candidate
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    basis = np.array(read(frame_file)['frame_source']['basis_columns'])
    opening = basis[:, 1]
    thickness = []
    for index in range(2):
        bolt = solids[f'collar:bolt_{index}']
        local_bolt = bolt.transform(np.c_[basis.T, np.zeros(3)])
        low, high = np.array(local_bolt.bounding_box()).reshape(2, 3)
        centre = basis @ ((low + high) / 2)
        washer = solids[f'collar:top_washer_{index}']
        washer_low, _ = np.array(washer.transform(np.c_[basis.T, np.zeros(3)]).bounding_box()).reshape(2, 3)
        samples = []
        for depth in np.arange(.25, 5.001, .25):
            point = centre + opening * (washer_low[1] - centre @ opening - depth)
            for degrees in range(0, 360, 5):
                angle = np.deg2rad(degrees)
                radial = np.cos(angle) * basis[:, 0] + np.sin(angle) * basis[:, 2]
                baseline = original.ray_cast(point, point + radial * 30)
                if len(baseline) < 2 or not 1 < np.linalg.norm(np.array(baseline[0].position)-point) < 3:
                    continue
                start = np.linalg.norm(np.array(baseline[0].position)-point)
                before = np.linalg.norm(np.array(baseline[1].position)-point) - start
                after = candidate.ray_cast(point, point + radial * 30)
                remaining = 0.
                if len(after) >= 2 and abs(np.linalg.norm(np.array(after[0].position)-point)-start) < .001:
                    remaining = np.linalg.norm(np.array(after[1].position)-point) - start
                samples.append(dict(depth_mm=float(depth), angle_deg=degrees,
                                    original_mm=float(before), candidate_mm=float(remaining)))
        thickness.append(dict(bolt=index, samples=len(samples),
                              minimum_original_mm=min(s['original_mm'] for s in samples),
                              minimum_candidate_mm=min(s['candidate_mm'] for s in samples),
                              worst_loss=max(samples, key=lambda s: s['original_mm']-s['candidate_mm']),
                              zero_contiguous_samples=sum(s['candidate_mm'] <= .001 for s in samples)))
    fixed = set(proof['fixed_keys'] + proof['moving_keys'])
    fixed |= {'collar:jaw_upper', 'collar:liner_upper'}
    fixed -= {k for k in fixed if k.startswith('link:FJ1500/bone_') or k == 'link:FJ1500/tendon_wire'}
    stages = []
    for index in range(2):
        for keys, axis in [
            ([f'collar:bolt_{index}', f'collar:top_washer_{index}'], opening),
            ([f'collar:nut_{index}', f'collar:bottom_washer_{index}'], -opening),
        ]:
            assert not fixed.intersection(keys)
            hits = []
            for offset in np.arange(20., -.001, -.25):
                for key in keys:
                    posed = solids[key].translate(axis * offset)
                    for target in sorted(fixed):
                        volume = float(A.I.intersection_volume(posed, solids[target]))
                        if volume > .001:
                            hits.append(dict(moving=key, target=target, volume_mm3=volume))
                if hits:
                    break
            stages.append(dict(keys=keys, clear=not hits, hits=hits,
                               first_collision_from_entry_mm=float(offset) if hits else None))
            if not hits:
                fixed.update(keys)
    result = dict(assembly_sha256=sha(source), study_sha256=sha(study_file),
                  script_sha256=sha(Path(__file__)), thickness=thickness, fastener_stages=stages,
                  upper_jaw_and_liner_preinstallation_verified=False,
                  whole_shoulder_assembly_verified=False)
    save(folder / 'collar_first_fastener_review.json', result)
    print(thickness, stages)


if __name__ == '__main__':
    main()
