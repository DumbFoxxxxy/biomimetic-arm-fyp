"""Audit both repaired elbow lattice sources against the complete neutral inventory."""
from itertools import combinations
from pathlib import Path

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_check.json'
    repair_file = folder / 'elbow_tpu_readback_repair.json'
    neutral_file = folder / 'inter_module_neutral_audit.json'
    model, repair, neutral = map(read, (manifest, repair_file, neutral_file))
    assert repair['assembly_sha256'] == neutral['assembly_sha256'] == sha(manifest)
    replacements = {row['key']: row['candidate_source'] for row in repair['rows'] if 'candidate_source' in row}
    assert set(replacements) == set(neutral['unchecked_nonmanifold_sources'])
    sources = {**model['sources'], **replacements}
    solids = {key: load_component(source) for key, source in sources.items()}
    shoulder_file = A.DESIGN / 's1_integrated_humerus_anchors/assembly_check.json'
    shoulder = set(read(shoulder_file)['sources'])
    count = 0
    hits = []
    for left, right in combinations(sorted(solids), 2):
        if left not in replacements and right not in replacements:
            continue
        count += 1
        volume = float(A.I.intersection_volume(solids[left], solids[right]))
        if volume > .001:
            hits.append(dict(left=left, right=right, overlap_mm3=volume,
                             cross_module=(left in shoulder) != (right in shoulder),
                             includes_shared_humerus='bone:humerus' in (left, right)))
    hits.sort(key=lambda row: -row['overlap_mm3'])
    report = folder/'repaired_elbow_tpu_neutral_audit.json'
    save(report, dict(
        original_assembly_sha256=sha(manifest), repair_sha256=sha(repair_file),
        previous_neutral_sha256=sha(neutral_file), script_sha256=sha(Path(__file__)),
        replacement_keys=sorted(replacements), checked_pairs=count, threshold_mm3=.001,
        overlaps=hits, manufacturing_verified=False, combined_assembly_verified=False,
        public_release_ready=False))
    save(folder/'assembly_tpu_readback_candidate.json', dict(
        **{key: value for key, value in model.items() if key not in ('scope', 'sources')},
        scope='Shared humerus with locally repaired elbow TPU collision sources. Neutral conflicts pending.',
        sources=sources, tpu_readback_audit_sha256=sha(report)))
    print('pairs', count, 'overlaps', len(hits),
          'cross module', sum(row['cross_module'] for row in hits), 'largest', hits[:4], flush=True)


if __name__ == '__main__':
    main()
