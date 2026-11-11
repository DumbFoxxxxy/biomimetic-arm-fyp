"""Compare actual solids before transferring a route result across arm baselines."""
from build_s1_teres_origin_mount import A, ROOT, read, save, load_component
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    reference_file = folder/'thorax_arm_reference.json'
    selection_file = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    current_file = ROOT/selection['assembly_report']
    assert provenance(current_file)['sha256'] == selection['assembly_sha256']
    current, reference = read(current_file), read(reference_file)
    rows = []
    for key in ('bone:clavicle', 'bone:scapula', 'bone:humerus'):
        source = dict(current['sources'][key])
        source['translation_source_mm'] = source.get('translation_source_mm') or [0., 0., 0.]
        selected = load_component(source)
        baseline = load_component(reference['sources'][key])
        added, removed = selected-baseline, baseline-selected
        rows.append(dict(bone=key, reference_source=reference['sources'][key], selected_source=source,
            reference_volume_mm3=float(baseline.volume()), selected_volume_mm3=float(selected.volume()),
            selected_only_volume_mm3=float(added.volume()), reference_only_volume_mm3=float(removed.volume()),
            selected_only_bounds_mm=list(added.bounding_box()) if added.volume() > .001 else None,
            reference_only_bounds_mm=list(removed.bounding_box()) if removed.volume() > .001 else None,
            same_geometry_within_0_001_mm3=added.volume() <= .001 and removed.volume() <= .001))
        print(key, 'selected only', rows[-1]['selected_only_volume_mm3'],
              'reference only', rows[-1]['reference_only_volume_mm3'], flush=True)
    save(folder/'reference_basis_comparison.json', dict(reference=provenance(reference_file),
        current_selection=provenance(selection_file), selected_assembly=provenance(current_file), bones=rows,
        decision='Keep current channel results scoped to thorax_arm_reference. Do not transfer collision, fit, mass or drive claims to S1 selection without geometry reconciliation and revalidation. No source geometry or current pointer replaced.'))


if __name__ == '__main__':
    main()
