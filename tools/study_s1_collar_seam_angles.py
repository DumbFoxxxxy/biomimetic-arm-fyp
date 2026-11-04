"""Screen access for rotated collar fasteners before redrawing the collar."""
import argparse
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stagger', action='store_true')
    stagger = parser.parse_args().stagger
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    proof_file = folder / 'verification_checks.json'
    proof = read(proof_file)
    assert proof['assembly_sha256'] == sha(source)
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    frame = read(frame_file)['frame_source']
    basis = np.array(frame['basis_columns'])
    solids = {key: load_component(row) for key, row in read(source)['sources'].items()}
    fixed = set(proof['fixed_keys'] + proof['moving_keys'])
    fixed -= {key for key in fixed if key.startswith('collar:')}
    keys = [f'collar:{part}_{index}' for index in range(2)
            for part in ('bolt', 'top_washer', 'nut', 'bottom_washer')]
    obstacles = ['tpu:deltoid'] + sorted(fixed - {'tpu:deltoid'})
    rows = []
    for angle in range(0, 360, 20 if stagger else 10):
        transform = A.I.rotation_about(frame['origin_mm'], basis[:, 2], angle)
        axis = transform[:, :3] @ basis[:, 1]
        for shift in (range(-24, 25, 4) if stagger else [0]):
            shifted = transform.copy()
            shifted[:, 3] += basis[:, 2] * shift
            stages = screen_stages(keys, solids, obstacles, shifted, axis)
            rows.append(dict(angle_deg=angle, axial_shift_mm=shift, stages=stages,
                             clear_groups=sum(stage['collision'] is None for stage in stages)))
    passing = [row['angle_deg'] for row in rows if row['clear_groups'] == 4]
    paired = []
    if stagger:
        for angle in sorted({row['angle_deg'] for row in rows}):
            feasible = [[row['axial_shift_mm'] for row in rows if row['angle_deg'] == angle
                         and all(stage['collision'] is None for stage in row['stages'][i:i+2])]
                        for i in (0, 2)]
            if all(feasible):
                paired.append(dict(angle_deg=angle, group_0_shifts_mm=feasible[0], group_1_shifts_mm=feasible[1]))
    report = 'collar_staggered_fastener_screen.json' if stagger else 'collar_seam_angle_screen.json'
    save(folder / report, dict(
        assembly_sha256=sha(source), proof_sha256=sha(proof_file), frame_sha256=sha(frame_file),
        script_sha256=sha(Path(__file__)), rows=rows, passing_angles_deg=passing,
        staggered_pairs=paired,
        scope='Hardware-only 20 mm straight approach, sampled every 0.5 mm. Collar solids and tool envelopes pending.'))
    print('passing', passing, 'staggered pairs', paired, flush=True)


def screen_stages(keys, solids, obstacles, transform, axis):
        stages = []
        for start in range(0, len(keys), 2):
            group = keys[start:start+2]
            sign = 1 if start % 4 == 0 else -1
            bodies = [solids[key].transform(transform) for key in group]
            collision = None
            for offset in np.arange(20., -.001, -.5):
                for key, body in zip(group, bodies):
                    posed = body.translate(axis * sign * offset)
                    for target in obstacles:
                        volume = float(A.I.intersection_volume(posed, solids[target]))
                        if volume > .001:
                            collision = dict(key=key, target=target, offset_mm=float(offset), volume_mm3=volume)
                            break
                    if collision:
                        break
                if collision:
                    break
            stages.append(dict(keys=group, collision=collision))
        return stages


if __name__ == '__main__':
    main()
