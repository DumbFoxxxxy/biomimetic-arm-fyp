"""Stage-specific whole-scapula withdrawal and inverse insertion probes.

Uses actual exported files. GH/humerus/clavicle bone and muscles are not yet
installed in this stated stage; their absence is explicitly enumerated.
Co-moving loose liners are a temporary handling assumption, not a bond or
validated complete installation sequence. No force/strength/actuation result.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_clavicle_mount as K
import build_shoulder_front_liner_keepers as L

DEST = L.DEST/'insertion_review'
SOURCE = L.DEST/'assembly_check.json'


def path_check(moving, fixed, direction, maximum):
    def evaluate(distance):
        hits = []
        for name, body in moving.items():
            posed = body.translate(direction*distance)
            for other, obstacle in fixed.items():
                volume = A.I.intersection_volume(posed, obstacle)
                if volume > .001:
                    hits.append(dict(moving=name, fixed=other, intersection_mm3=volume))
        return sorted(hits, key=lambda p:-p['intersection_mm3'])
    initial = evaluate(0.)
    if initial:
        return dict(last_sample_below_threshold_mm=None, first_interference_mm=0., hits=initial,
                    sampled_path_clear=False, checked_samples_mm=[0.])
    previous = 0.
    checked = [0.]
    for distance in np.arange(1., maximum+.01, 1.):
        checked.append(float(distance))
        hits = evaluate(float(distance))
        if hits:
            fine_samples = []
            for fine in np.arange(previous+.1, distance+.001, .1):
                fine_samples.append(round(float(fine), 1))
                detailed = evaluate(float(fine))
                if detailed:
                    return dict(last_sample_below_threshold_mm=round(float(fine-.1), 1),
                                first_interference_mm=round(float(fine), 1), hits=detailed,
                                sampled_path_clear=False, checked_samples_mm=checked,
                                boundary_samples_mm=fine_samples)
            raise ValueError('Fine scan missed coarse collision')
        previous = float(distance)
    return dict(last_sample_below_threshold_mm=maximum, first_interference_mm=None,
                hits=[], sampled_path_clear=True, checked_samples_mm=checked)


def run(source, dest, maximum):
    dest.mkdir(parents=True, exist_ok=True)
    assembly = json.loads(source.read_text(encoding='utf8'))
    # Load/hash every source even if it is explicitly not yet installed.
    solids = {key:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for key, r in assembly['sources'].items()}
    front = {k for k in solids if k.startswith('clamp:') and k.endswith('/front_liner')}
    moving = {'bone:scapula'}
    keepers = {k for k in solids if k.startswith('liner_keeper:')}
    fixed = {k for k in solids if k.startswith(('fixture:', 'bench:')) or k in front
             or k == 'clavicle:fixed_support' or k.startswith('clavicle:M4_')} | keepers
    definitions = [
        dict(name='bone_against_preplaced_front_liners', moving_keys=sorted(moving), fixed_keys=sorted(fixed),
             handling='front liners at nominal positions with any declared keeper hardware included; actual retaining force/fit unverified'),
        dict(name='bone_and_temporarily_co_moving_front_liners', moving_keys=sorted(moving|front), fixed_keys=sorted(fixed-front-keepers),
             handling='keeper hardware not yet installed; operator temporarily coordinates 3 liners with bone; NOT bonded, handling unverified'),
        dict(name='bone_before_any_front_liners', moving_keys=sorted(moving), fixed_keys=sorted(fixed-front-keepers),
             handling='front liners not yet installed; later liner insertion must be verified separately')]
    rows = []
    for stage in definitions:
        present = set(stage['moving_keys'])|set(stage['fixed_keys'])
        stage['not_yet_installed_keys'] = sorted(set(solids)-present)
        bodies = {k:solids[k] for k in stage['moving_keys']}
        obstacles = {k:solids[k] for k in stage['fixed_keys']}
        trials = []
        for x, z in itertools.product([0., -.3, .3], repeat=2):
            direction = A.D.unit(np.array([x, -1., z]))
            result = path_check(bodies, obstacles, direction, maximum)
            trials.append(dict(direction_source=direction.tolist(), **result))
        row = dict(**stage, trials=trials,
                   clear_trial_indices=[i for i, r in enumerate(trials) if r['sampled_path_clear']])
        rows.append(row)
        print(stage['name'], 'clear', row['clear_trial_indices'], 'first limits',
              [(r['first_interference_mm'], r['hits'][0]['fixed'] if r['hits'] else None) for r in trials], flush=True)
    report = dict(scope=__doc__, assembly_source=source.relative_to(A.ROOT).as_posix(),
                  assembly_sha256=A.hashed(source), script_sha256=A.hashed(Path(__file__).resolve()),
                  source_count=len(solids), all_source_hashes_match=True, source_manifest=assembly['sources'],
                  direction_family='normalized [x,-1,z], x,z in {0,-0.3,+0.3}; no orientation changes',
                  maximum_translation_mm=maximum, coarse_step_mm=1., boundary_step_mm=.1,
                  interference_threshold_mm3=.001, stages=rows,
                  inverse_path_note='A sampled clear withdrawal has the same geometric samples in reverse; not yet an independent return or full assembly test.',
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['complete bone/liner/hardware handling, tooling and multi-part order',
                                'unexamined translations/rotations, nonlinear deformation and physical assembly',
                                'GH assembly and subsequent muscle installation',
                                'no shoulder movement, material load or actuation acceptance'])
    (dest/'scapula_insertion_audit.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser = argparse.ArgumentParser()
    parser.add_argument('--assembly', type=Path, default=SOURCE)
    parser.add_argument('--out', type=Path, default=DEST)
    parser.add_argument('--maximum', type=float, default=20.)
    args = parser.parse_args()
    run(args.assembly, args.out, args.maximum)
