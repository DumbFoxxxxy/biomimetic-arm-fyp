"""Trial an engineered humerus with integral, stationary tendon seats."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, export, load_component, read, save, sha


def main():
    base = A.DESIGN / 's1_deltoid_support_draft'
    source = base / 'assembly_check.json'
    proof_file = base / 'verification_checks.json'
    relief_file = base / 'collar_first_screen.json'
    model, proof, relief = read(source), read(proof_file), read(relief_file)
    assert proof['assembly_sha256'] == relief['assembly_sha256'] == sha(source)
    assert proof['conditional_bundle_placement_geometry_pass'] and not relief['refined_hits']
    solids = {key: load_component(row) for key, row in model['sources'].items()}
    bone_key = 'bone:humerus'
    original_bone = solids[bone_key]
    retired = {key for key in solids if key.startswith('collar:')}
    upper = load_component(relief['candidate_mesh'])
    lower = solids['collar:jaw_lower']
    liners = [solids['collar:liner_upper'], solids['collar:liner_lower']]
    mf = A.I._manifold()
    union = lambda pieces: mf.Manifold.batch_boolean(pieces, mf.OpType.Add)
    direction = np.array(proof['direction_source'])
    output = A.DESIGN / 's1_integrated_humerus_anchors'
    output.mkdir(exist_ok=True)
    trials = []
    selected = None
    # The liners become rigid transition stock. A small overlap fuses the
    # interfaces rather than relying on imported coincident triangle faces.
    for allowance in [.025, .05, .1]:
        transition = union([liner.translate(shift * allowance)
                            for liner in liners for shift in [np.zeros(3), *np.eye(3), *(-np.eye(3))]])
        trial = union([original_bone, upper, lower, transition])
        components = sorted((float(piece.volume()) for piece in trial.decompose()
                             if piece.volume() > 1e-8), reverse=True)
        fusion = {name: float(A.I.intersection_volume(transition, body))
                  for name, body in [('bone', original_bone), ('upper', upper), ('lower', lower)]}
        row = dict(transition_allowance_mm=allowance, positive_component_volumes_mm3=components,
                   transition_overlap_mm3=fusion)
        trials.append(row)
        if len(components) == 1 and min(fusion.values()) > .001:
            selected = trial.as_original().simplify(1e-6)
            break
    result = dict(scope=__doc__, assembly_sha256=sha(source), placement_sha256=sha(proof_file),
                  relief_sha256=sha(relief_file), script_sha256=sha(Path(__file__)),
                  trials=trials, retired_keys=sorted(retired), geometry_adopted=False,
                  whole_shoulder_assembly_verified=False)
    report_path = output / 'integration_screen.json'
    if selected is None:
        save(report_path, result)
        print('no connected fusion', trials, flush=True)
        return
    added = selected - original_bone
    offsets = np.unique(np.r_[np.arange(0, 80.001, .5), np.arange(0, 3.001, .05)])
    retained = set(solids) - retired - {bone_key}
    static_hits = []
    for key in sorted(retained):
        volume = float(A.I.intersection_volume(added, solids[key]))
        if volume > .001:
            static_hits.append(dict(target=key, overlap_mm3=volume))
    entry_hits = []
    for offset in offsets:
        for key in proof['moving_keys']:
            volume = float(A.I.intersection_volume(added, solids[key].translate(direction*offset)))
            if volume > .001:
                entry_hits.append(dict(target=key, offset_mm=float(offset), overlap_mm3=volume))
    result.update(added_mm3=float(added.volume()), bone_loss_mm3=float((original_bone-selected).volume()),
                  static_hits=static_hits, entry_hits=entry_hits, entry_samples=len(offsets))
    save(report_path, result)
    print('added', result['added_mm3'], 'static', static_hits, 'entry hits', len(entry_hits), flush=True)
    if static_hits or entry_hits:
        return
    actual, metadata = export(selected, output/'humerus_integral_anchors.obj',
                              np.array(model['sources'][bone_key]['translation_source_mm'] or [0, 0, 0]))
    result['candidate_mesh'] = metadata
    result['roundtrip_material_difference_mm3'] = float((actual-selected).volume() + (selected-actual).volume())
    save(report_path, result)


if __name__ == '__main__':
    main()
