"""Verify the exported integral-anchor humerus and write its assembly manifest."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def main():
    base = A.DESIGN / 's1_deltoid_support_draft'
    source = base / 'assembly_check.json'
    placement = base / 'verification_checks.json'
    folder = A.DESIGN / 's1_integrated_humerus_anchors'
    screen = folder / 'integration_screen.json'
    parent, proof, trial = read(source), read(placement), read(screen)
    assert proof['assembly_sha256'] == trial['assembly_sha256'] == sha(source)
    assert trial['placement_sha256'] == sha(placement)
    bone_key = 'bone:humerus'
    retired = set(trial['retired_keys'])
    sources = {key: row for key, row in parent['sources'].items() if key not in retired}
    original = load_component(sources[bone_key])
    sources[bone_key] = trial['candidate_mesh']
    solids = {key: load_component(row) for key, row in sources.items()}
    bone = solids[bone_key]
    added = bone-original
    loss = float((original-bone).volume())
    static_hits = []
    for key, solid in solids.items():
        if key == bone_key:
            continue
        volume = float(A.I.intersection_volume(bone, solid))
        if volume > .001:
            static_hits.append(dict(target=key, overlap_mm3=volume))
    direction = np.array(proof['direction_source'])
    offsets = np.unique(np.r_[np.arange(0, 80.001, .25), np.arange(0, 3.001, .05)])
    moving = proof['moving_keys']
    fixed = set(proof['fixed_keys']) - retired
    pending = set(proof['later_installation_keys']) - retired
    assert not (fixed & pending or fixed.intersection(moving) or pending.intersection(moving))
    assert fixed | pending | set(moving) == set(solids)
    entry_hits = []
    for offset in offsets:
        for key in moving:
            volume = float(A.I.intersection_volume(added, solids[key].translate(direction*offset)))
            if volume > .001:
                entry_hits.append(dict(target=key, offset_mm=float(offset), overlap_mm3=volume))
    route_file = A.DESIGN / 's1_port_path_revision/routes_compact_distal_candidate.json'
    assert parent['routes_sha256'] == sha(route_file)
    wires = {}
    for route in read(route_file)['routes']:
        wire = A.D.outer_envelope(np.array(route['full_wire_path_source_mm']), .3)
        wires[route['head']] = float(A.I.intersection_volume(added, wire))
    component_count = sum(piece.volume() > 1e-8 for piece in bone.decompose())
    passed = component_count == 1 and loss <= .001 and not static_hits and not entry_hits and max(wires.values()) <= .001
    report = dict(
        scope=__doc__, parent_report=source.relative_to(ROOT).as_posix(), parent_sha256=sha(source),
        placement_sha256=sha(placement), screen_sha256=sha(screen), script_sha256=sha(Path(__file__)),
        routes_sha256=sha(route_file), sources=sources, n_instances=len(sources),
        intended_body_owners={key: owner for key, owner in parent['intended_body_owners'].items() if key in sources},
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        replaced_keys=[bone_key], retired_keys=sorted(retired), bone_loss_mm3=loss,
        positive_components=component_count, added_mm3=float(added.volume()),
        final_pose_humerus_hits=static_hits, added_stock_wire_overlap_mm3=wires,
        changed_geometry_entry_hits=entry_hits, entry_samples=len(offsets),
        fixed_keys=sorted(fixed), moving_keys=moving, later_installation_keys=sorted(pending),
        integrated_geometry_screen_pass=passed, stage_pass=False, adopted_as_current=False,
        complete_assembly_verified=False, actual_drive_range_deg=None,
        pending_checks=['integral stock load path and printable material selection',
                        'later infraspinatus and deltoid fixation sequence',
                        'temporary supports, tendon and wire feeding, shoulder motion and actuation'])
    save(folder / 'assembly_check.json', report)
    print('sources', len(sources), 'components', component_count, 'bone loss', loss,
          'static hits', static_hits, 'entry hits', len(entry_hits), 'wires', wires, 'pass', passed, flush=True)


if __name__ == '__main__':
    main()
