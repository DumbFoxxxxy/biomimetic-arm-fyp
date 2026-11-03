"""Actual-neck-axis insertion and capture contrast for GH manufacturing meshes.

The bare staged fixture is explicit: scapula clamped before GH/humerus,
clavicle bone and muscles. Closed-cap straight escape is diagnostic only;
no arbitrary escape, strength, preload or actuation acceptance.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_front_liner_keepers as L

SOURCE = L.DEST/'assembly_check.json'
DEST = A.DESIGN/'s1_GH_assembly_review'


def scan(moving, direction, fixed, maximum=25., step=.1, continue_after_hit=False):
    rows = []
    for distance in np.linspace(0., maximum, int(round(maximum/step))+1):
        posed = moving.translate(direction*distance)
        volumes = {k:A.I.intersection_volume(posed, s) for k, s in fixed.items()}
        hits = [dict(part=k, volume_mm3=v) for k, v in volumes.items() if v > .001]
        rows.append(dict(distance_mm=round(float(distance), 6), maximum_overlap_mm3=max(volumes.values()),
                         pair_volumes_mm3=volumes, hits_over_threshold=hits))
        if hits and not continue_after_hit:
            break
    first = next((r for r in rows if r['hits_over_threshold']), None)
    return dict(direction_source=direction.tolist(), step_mm=step, requested_maximum_mm=maximum,
                checked_sample_count=len(rows), samples=rows, first_interference=first,
                entire_sampled_path_clear=len(rows) == int(round(maximum/step))+1 and first is None)


def run(source, dest):
    source, dest = source.resolve(), dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    data = json.loads(source.read_text(encoding='utf8'))
    joint_file = A.G1/'assembly.json'
    joint = json.loads(joint_file.read_text(encoding='utf8'))
    n = np.array(joint['neck_axis'])
    G = np.array(joint['centre_source_mm'])
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k, r in data['sources'].items()}
    bare_keys = {k for k in solids if k.startswith(('fixture:', 'bench:', 'clamp:', 'liner_keeper:'))
                 or k in ['bone:scapula', 'clavicle:fixed_support'] or k.startswith('clavicle:M4_')}
    gh = {k for k in solids if k.startswith('GH_')}
    open_context = {k:solids[k] for k in bare_keys}
    closed_context = {k:solids[k] for k in bare_keys|gh}
    ball = A.I._manifold().Manifold.sphere(joint['ball_radius_mm'], 40).translate(G)
    cases = []
    for label, body, fixed in [('actual_humerus_open_caps', solids['bone:humerus'], open_context),
                               ('actual_humerus_closed_caps', solids['bone:humerus'], closed_context),
                               ('diagnostic_ball_closed_caps', ball, closed_context)]:
        result = scan(body, n, fixed)
        cases.append(dict(name=label, fixed_keys=list(fixed),
                          not_yet_installed_keys=sorted(set(solids)-set(fixed)-{'bone:humerus'}),
                          probe_is_manufacturing_humerus=label.startswith('actual'), path=result))
        print(label, 'clear', result['entire_sampled_path_clear'],
              'first', None if result['first_interference'] is None else
              (result['first_interference']['distance_mm'], result['first_interference']['hits_over_threshold']), flush=True)
    depth = data.get('GH_retention_revision', {}).get('retaining_depth_mm', joint['cap_depth_mm'])
    result = dict(scope=__doc__, assembly_source=source.relative_to(A.ROOT).as_posix(), assembly_sha256=A.hashed(source),
                  script_sha256=A.hashed(Path(__file__).resolve()), all_source_hashes_match=True,
                  source_object_count=len(solids), joint_reference_sha256=A.hashed(joint_file),
                  actual_neck_axis_source=n.tolist(), ball_radius_mm=joint['ball_radius_mm'],
                  socket_inner_radius_mm=joint['socket_inner_radius_mm'], tested_cap_depth_mm=depth,
                  nominal_unslotted_aperture_radius_mm=float(np.sqrt(joint['socket_inner_radius_mm']**2-depth**2)),
                  staged_open_scene_objects=len(bare_keys)+1, staged_closed_scene_objects=len(bare_keys|gh)+1,
                  cases=cases, stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['cap removal and fastener/tool sequence must be checked separately',
                                'non-axial or curved escape, printed fits, strength and dynamic stability',
                                'rear clamp preload and full GH/muscle installation',
                                'no new arm drive or motion-range acceptance'])
    (dest/'axial_capture_audit.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8', newline='\n')
    return result


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser = argparse.ArgumentParser()
    parser.add_argument('--assembly', type=Path, default=SOURCE)
    parser.add_argument('--out', type=Path, default=DEST/'parent_187')
    args = parser.parse_args()
    run(args.assembly, args.out)
