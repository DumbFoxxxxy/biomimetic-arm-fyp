"""Screen local elbow adapter pockets for the remaining shoulder hardware."""
import numpy as np
from tempfile import TemporaryDirectory
from pathlib import Path

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, save, sha


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    audit_file = folder/'anatomical_elbow_interface_assembly.json'
    audit = read(audit_file)
    assembly_file = folder/'assembly_check.json'
    assert sha(assembly_file) == audit['original_sha256']
    sources = {**read(assembly_file)['sources'], **audit['changed_sources']}
    mf = A.I._manifold()
    results = {}
    for head in ('FJ1486', 'FJ1477'):
        key = f'elbow:{head}_bone_adapter'
        adapter = load_component(sources[key])
        blockers = {}
        for row in audit['overlaps']:
            pair = (row['left'], row['right'])
            if key in pair and row['after_mm3'] > audit['threshold_mm3']:
                other = pair[1] if pair[0] == key else pair[0]
                blockers[other] = load_component(sources[other])
        interfaces = {'bone:humerus': load_component(sources['bone:humerus'])}
        for name, source in sources.items():
            if name == key or name in audit['omitted_original_left_fasteners']:
                continue
            if name.startswith((f'elbow:{head}_bone_', f'elbow:{head}_proximal_')):
                interfaces[name] = load_component(source)

        # The translated solid cuts create relief in six coordinate directions.
        # Direct measured clearance is recorded separately.
        offsets = [np.zeros(3), *np.eye(3), *(-np.eye(3))]
        rows = []
        best = None
        for clearance in (.15, .12, .10, .08, .06, .04, .02, 0.):
            cuts = [adapter ^ solid.translate(clearance*offset)
                    for solid in blockers.values() for offset in offsets]
            candidate = (adapter-mf.Manifold.batch_boolean(cuts, mf.OpType.Add))
            candidate = candidate.as_original().simplify(1e-6)
            removed = adapter-candidate
            row = dict(offset_mm=clearance, removed_mm3=float(removed.volume()),
                       positive_components=sum(part.volume() > 1e-8 for part in candidate.decompose()),
                       removal_to_interface_mm_capped1={
                           name: float(A.I.gap(removed, solid, 1.)) for name, solid in interfaces.items()})
            rows.append(row)
            if row['positive_components'] != 1 or min(row['removal_to_interface_mm_capped1'].values()) < .1:
                continue
            origin = np.asarray(sources[key]['translation_source_mm'])
            try:
                with TemporaryDirectory(dir=folder) as trial_folder:
                    actual, _ = export(candidate, Path(trial_folder)/'trial.obj', origin)
            except (AssertionError, ValueError) as error:
                row['export_error'] = str(error) or type(error).__name__
                continue
            row['readback_overlap_mm3'] = {
                name: float(A.I.intersection_volume(actual, solid)) for name, solid in blockers.items()}
            row['readback_gap_mm_capped1'] = {
                name: float(A.I.gap(actual, solid, 1.)) for name, solid in blockers.items()}
            if max(row['readback_overlap_mm3'].values()) <= audit['threshold_mm3']:
                row['positive_clearance'] = min(row['readback_gap_mm_capped1'].values()) > 0
                if best is None or row['positive_clearance']:
                    best = (candidate, origin, row)
                if row['positive_clearance']:
                    break
        if best:
            candidate, origin, row = best
            _, row['source'] = export(candidate, folder/f'anatomical_elbow_shoulder_relief_{head}.obj', origin)
        results[head] = rows
    save(folder/'anatomical_elbow_adapter_shoulder_relief.json', dict(
        interface_audit_sha256=sha(audit_file), trials=results,
        scope='Local neutral pocket screening with a 0.1 mm exclusion near checked elbow interfaces. '
              'Missing replacement fastener, remaining assembly and motion require separate checks.'))
    for head, rows in results.items():
        selected = next((row for row in rows if 'source' in row), None)
        print(head, selected)


if __name__ == '__main__':
    main()
