"""Splice matching shaft stock while preserving proximal S1 and distal elbow features."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, save, sha


def main():
    shoulder_dir = A.DESIGN / 's1_integrated_humerus_anchors'
    shoulder_file = shoulder_dir / 'assembly_check.json'
    seam_file = shoulder_dir / 'shared_humerus_seam_screen.json'
    old_preview = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_preview'
    elbow_file = old_preview / 'integration_audit.json'
    overlay_file = old_preview / 'spatial_overlay.json'
    shoulder, seam, elbow, overlay = map(read, (shoulder_file, seam_file, elbow_file, overlay_file))
    assert seam['shoulder_sha256'] == sha(shoulder_file) and seam['elbow_audit_sha256'] == sha(elbow_file)
    assert elbow['overlay_sha256'] == sha(overlay_file)
    assert shoulder['integrated_geometry_screen_pass']
    bone_key = 'bone:humerus'
    shoulder_bone = load_component(shoulder['sources'][bone_key])
    elbow_bone = load_component(elbow['original_elbow_humerus'])
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    assert seam['frame_sha256'] == sha(frame_file)
    frame = read(frame_file)['frame_source']
    basis = np.array(frame['basis_columns'])
    origin = np.array(frame['origin_mm'])
    local = np.c_[basis.T, -basis.T @ origin]
    world = np.c_[basis, origin]
    shoulder_local, elbow_local = (body.transform(local) for body in (shoulder_bone, elbow_bone))
    bounds = np.array([shoulder_local.bounding_box(), elbow_local.bounding_box()]).reshape(2, 2, 3)
    low, high = bounds[:, 0].min(axis=0)-1, bounds[:, 1].max(axis=0)+1
    cube = A.I._manifold().Manifold.cube
    elbow_sources = {key: value for key, value in overlay['sources'].items()
                     if key.startswith('elbow:') or key in ('bone:forearm', 'tpu:elbow_flexor', 'tpu:elbow_extensor')}
    assert len(elbow_sources) == overlay['elbow_objects_added']
    assert not set(elbow_sources).intersection(shoulder['sources'])
    sources = {**shoulder['sources'], **elbow_sources}
    parts, readback_errors = {}, {}
    for key, row in sources.items():
        if key == bone_key:
            continue
        assert sha(ROOT/row['file']) == row['sha256']
        try:
            parts[key] = load_component(row)
        except (ValueError, AssertionError) as error:
            if key not in ('tpu:elbow_flexor', 'tpu:elbow_extensor'):
                raise
            readback_errors[key] = str(error) or type(error).__name__

    def fit(candidate):
        hits = []
        for key, part in parts.items():
            volume = float(A.I.intersection_volume(candidate, part))
            if volume > .001:
                hits.append(dict(target=key, overlap_mm3=volume))
        return hits

    trials = []
    selected = None
    closest = None
    chosen_station = None
    transferred = None
    for row in sorted(seam['rows'], key=lambda item: item['station_mm']):
        if row['symmetric_difference_mm3'] > .001:
            continue
        station = row['station_mm']
        upper_end = np.array([high[0], high[1], station+1])
        lower_start = np.array([low[0], low[1], station-1])
        upper = shoulder_local ^ cube(upper_end-low).translate(low)
        lower = elbow_local ^ cube(high-lower_start).translate(lower_start)
        candidate = (upper+lower).transform(world).as_original().simplify(1e-6)
        count = sum(part.volume() > 1e-8 for part in candidate.decompose())
        hits = fit(candidate)
        trials.append(dict(station_mm=station, components=count, bone_interface_hits=hits))
        if count == 1 and not hits:
            selected = candidate
            chosen_station = station
            break
        if count == 1 and (closest is None or len(hits) < len(closest[2])):
            closest = (station, candidate, hits)
    if selected is None and closest is not None:
        station, candidate, hits = closest
        elbow_hits = [hit['target'] for hit in hits if hit['target'] in elbow_sources]
        if len(elbow_hits) == len(hits):
            # Bring over the original elbow bone's local voids around its
            # interfaces. Keep the remaining proximal S1 material intact.
            mf = A.I._manifold()
            boxes = []
            for key in elbow_hits:
                minimum, maximum = np.array(parts[key].bounding_box()).reshape(2, 3)
                boxes.append(cube(maximum-minimum+1).translate(minimum-.5))
            region = mf.Manifold.batch_boolean(boxes, mf.OpType.Add)
            removed = (candidate-elbow_bone) ^ region
            local_transfer = (candidate-removed).as_original().simplify(1e-6)
            count = sum(part.volume() > 1e-8 for part in local_transfer.decompose())
            remaining = fit(local_transfer)
            transferred = dict(station_mm=station, targets=elbow_hits, margin_mm=.5,
                               removed_mm3=float(removed.volume()), components=count,
                               bone_interface_hits=remaining)
            if count == 1 and not remaining:
                selected = local_transfer
                chosen_station = station
    folder = old_preview.parent / 's1_elbow_shared_humerus'
    folder.mkdir(exist_ok=True)
    result = dict(scope=__doc__, shoulder_sha256=sha(shoulder_file), elbow_audit_sha256=sha(elbow_file),
                  elbow_overlay_sha256=sha(overlay_file), seam_sha256=sha(seam_file),
                  script_sha256=sha(Path(__file__)), trials=trials, readback_errors=readback_errors,
                  local_elbow_void_transfer=transferred,
                  common_humerus_geometry_unified=False, combined_assembly_verified=False,
                  public_release_ready=False, stage_pass=False)
    save(folder/'shared_humerus_audit.json', result)
    if selected is None:
        print('no splice candidate', len(trials), 'local transfer', transferred, flush=True)
        return
    export_attempts = []
    seating_relief = None
    for tolerance in [1e-6, 1e-5, 1e-4, .001]:
        prepared = selected.as_original().simplify(tolerance)
        try:
            actual, metadata = export(prepared, folder/'shared_humerus.obj', origin)
        except (AssertionError, ValueError) as error:
            export_attempts.append(dict(simplification_mm=tolerance, error=str(error) or type(error).__name__))
            continue
        hits = fit(actual)
        change = float((actual-selected).volume() + (selected-actual).volume())
        export_attempts.append(dict(simplification_mm=tolerance, change_mm3=change, interface_hits=hits))
        if not hits and abs(change) <= .01:
            sources[bone_key] = metadata
            break
    else:
        # Resolve only small residual seating interferences detected after
        # exact indexed import. Record this cut separately from void transfer.
        remaining = export_attempts[-1].get('interface_hits', [])
        if not remaining or any(hit['target'] not in elbow_sources or hit['overlap_mm3'] > .1 for hit in remaining):
            result['export_attempts'] = export_attempts
            save(folder/'shared_humerus_audit.json', result)
            print('splice readback failed', export_attempts, flush=True)
            return
        mf = A.I._manifold()
        clearance = .02
        guards = [parts[hit['target']].translate(delta * clearance)
                  for hit in remaining for delta in [np.zeros(3), *np.eye(3), *(-np.eye(3))]]
        guard = mf.Manifold.batch_boolean(guards, mf.OpType.Add)
        prepared = (selected - guard).as_original().simplify(1e-6)
        actual, metadata = export(prepared, folder/'shared_humerus.obj', origin)
        hits = fit(actual)
        seating_relief = dict(targets=[hit['target'] for hit in remaining], axis_shifts_mm=clearance,
                              removed_mm3=float((selected-actual).volume()), readback_hits=hits)
        if hits:
            result.update(export_attempts=export_attempts, seating_relief=seating_relief)
            save(folder/'shared_humerus_audit.json', result)
            print('seating relief failed', seating_relief, flush=True)
            return
        sources[bone_key] = metadata
    result.update(common_humerus_geometry_unified=True, selected_station_mm=chosen_station,
                  roundtrip_bone_interface_hits=hits, bone_interfaces_screened=len(parts),
                  export_attempts=export_attempts, seating_relief=seating_relief)
    save(folder/'shared_humerus_audit.json', result)
    save(folder/'assembly_check.json', dict(
        scope='Local merged geometry with one shared humerus. Whole assembly review pending.',
        sources=sources, n_instances=len(sources), audit_sha256=sha(folder/'shared_humerus_audit.json'),
        common_humerus_geometry_unified=True, combined_assembly_verified=False,
        public_release_ready=False, stage_pass=False))
    print('shared bone', chosen_station, 'objects', len(sources),
          'bone fits', len(parts), 'readback errors', readback_errors, 'hits', hits, flush=True)


if __name__ == '__main__':
    main()
