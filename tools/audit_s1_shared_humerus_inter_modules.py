"""Check cross-module neutral collisions in the shared-humerus candidate."""
from pathlib import Path

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_check.json'
    audit = folder / 'shared_humerus_audit.json'
    model, bone_audit = read(manifest), read(audit)
    assert model['audit_sha256'] == sha(audit)
    assert model['common_humerus_geometry_unified'] and not bone_audit['roundtrip_bone_interface_hits']
    shoulder_file = A.DESIGN / 's1_integrated_humerus_anchors/assembly_check.json'
    assert bone_audit['shoulder_sha256'] == sha(shoulder_file)
    shoulder = set(read(shoulder_file)['sources'])
    elbow = set(model['sources']) - shoulder
    unchecked = set(bone_audit['readback_errors'])
    solids = {}
    for key, row in model['sources'].items():
        assert sha(ROOT / row['file']) == row['sha256']
        if key not in unchecked and key != 'bone:humerus':
            solids[key] = load_component(row)
    hits = []
    count = 0
    for upper in sorted(shoulder - unchecked - {'bone:humerus'}):
        for lower in sorted(elbow - unchecked):
            count += 1
            overlap = float(A.I.intersection_volume(solids[upper], solids[lower]))
            if overlap > .001:
                hits.append(dict(shoulder=upper, elbow=lower, overlap_mm3=overlap))
    hits.sort(key=lambda row: -row['overlap_mm3'])
    save(folder / 'inter_module_neutral_audit.json', dict(
        assembly_sha256=sha(manifest), bone_audit_sha256=sha(audit), script_sha256=sha(Path(__file__)),
        sampled_pose='neutral', checked_cross_module_pairs=count, overlaps=hits,
        unchecked_nonmanifold_sources=sorted(unchecked), threshold_mm3=.001,
        combined_assembly_verified=False, public_release_ready=False))
    print('pairs', count, 'hits', len(hits), 'largest', hits[:4], 'unchecked', sorted(unchecked), flush=True)


if __name__ == '__main__':
    main()
