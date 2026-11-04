"""Record the original elbow anatomy, attachment geometry and wire routes for S1 integration."""
from build_s1_teres_origin_mount import ROOT, read, save, sha


def main():
    base = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    original_file = base/'assembly_check.json'
    repaired_file = base/'assembly_tpu_readback_candidate.json'
    rigid_audit = base/'inter_module_neutral_audit.json'
    soft_audit = base/'repaired_elbow_tpu_neutral_audit.json'
    original, repaired = read(original_file)['sources'], read(repaired_file)['sources']
    rigid, soft = read(rigid_audit), read(soft_audit)
    assert rigid['assembly_sha256'] == soft['original_assembly_sha256'] == sha(original_file)
    assert soft['previous_neutral_sha256'] == sha(rigid_audit)
    assert read(repaired_file)['tpu_readback_audit_sha256'] == sha(soft_audit)
    muscles = {}
    for muscle, head in [('flexor', 'FJ1486'), ('extensor', 'FJ1477')]:
        key = f'tpu:elbow_{muscle}'
        source = original[key]
        mesh = ROOT/source['file']
        assert sha(mesh) == source['sha256']
        design = mesh.parent/'design.json'
        definition = read(design)
        shift = source['translation_source_mm']
        muscles[muscle] = dict(
            anatomical_id=head, display_source=source, collision_proxy=repaired[key],
            definition=design.relative_to(ROOT).as_posix(), definition_sha256=sha(design),
            pad_centres_source_mm=[[value+offset for value, offset in zip(point, shift, strict=True)]
                                   for point in definition['pad_centres_display_mm']],
            wire_path_source_mm=[[value+offset for value, offset in zip(point, shift, strict=True)]
                                 for point in definition['wire_polyline_display_mm']])
    hard = rigid['overlaps']
    two_muscle = [row for row in soft['overlaps'] if row['cross_module']
                  and row['left'].startswith('tpu:') and row['right'].startswith('tpu:')]
    save(base/'anatomical_elbow_baseline.json', dict(
        scope='Frozen integration inputs. Original elbow TPU display geometry, attachment centres and wire route. No combined motion approval.',
        original_assembly_sha256=sha(original_file), muscles=muscles,
        neutral_contacts=dict(rigid_inventory_audit_sha256=sha(rigid_audit),
                              tpu_proxy_audit_sha256=sha(soft_audit), previously_checked_pairs=rigid['checked_cross_module_pairs'],
                              original_inventory_cross_module_hits=hard, tpu_to_tpu=two_muscle),
        design_priority='Keep elbow lattice geometry, muscle anchors and cable routes. Resolve rigid interfaces first. Treat remaining shoulder/elbow TPU overlap as a separate geometry and pose problem.',
        combined_assembly_verified=False))
    print('Protected elbow muscles:', list(muscles), 'rigid inventory hits:', len(hard),
          'TPU pair hits:', two_muscle)


if __name__ == '__main__':
    main()
