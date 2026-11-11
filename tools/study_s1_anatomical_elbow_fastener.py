"""Close the angled elbow cap joint with a nominal screw and bonded sleeve."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    assembly_file = folder/'assembly_check.json'
    cap_file = folder/'anatomical_elbow_cap_relief.json'
    audit_file = folder/'anatomical_elbow_interface_with_bolt_relief.json'
    baseline_file = folder/'anatomical_elbow_baseline.json'
    audit = read(audit_file)
    assert sha(assembly_file) == audit['original_sha256']
    assert sha(cap_file) == audit['input_reports_sha256'][cap_file.name]
    assert sha(baseline_file) == audit['baseline_sha256']
    cap_audit = read(cap_file)
    ear = next(row for row in cap_audit['adapter_ear_trials'] if 'sources' in row['matched_interface'])
    fit = next(row for row in cap_audit['angled_blind_fastener_trials']
               if all(row[k] == ear[k] for k in ('shift_mm', 'azimuth_deg', 'tilt_deg')))
    face, axis = np.array(fit['face_mm']), np.array(fit['axis'])

    def axial(lo, hi, radius):
        return rod(face+axis*lo, face+axis*hi, radius)

    # Negative coordinates run into the cap from its nominal entry plane.
    # The screw solid uses the thread core, with engagement evaluated separately
    # against the major-diameter envelope. No helical surfaces are assumed.
    sleeve_bottom, sleeve_top = -3.15, -1.25
    screw_tip = -3.10
    sleeve = axial(sleeve_bottom, sleeve_top, 1.64)-axial(-3.2, -1.2, 1.25)
    screw = (axial(screw_tip, -1.23, 1.20) + axial(sleeve_top, 2.41, 1.50)
             + axial(2.40, 3.90, 2.50))
    thread_envelope = axial(screw_tip, sleeve_top, 1.50)-axial(-3.2, -1.2, 1.25)
    sleeve_insertion = axial(sleeve_bottom, 20., 1.64)
    screw_insertion = axial(screw_tip, 20., 1.50)+axial(2.40, 20., 2.50)
    driver = axial(3.90, 20., 2.)

    assembly = read(assembly_file)
    sources = {**assembly['sources'], **audit['changed_sources']}
    baseline = read(baseline_file)
    for label, muscle in baseline['muscles'].items():
        sources[f'tpu:elbow_{label}'] = muscle['collision_proxy']
    inventory = {key: load_component(source) for key, source in sources.items()
                 if key not in audit['omitted_original_left_fasteners']}
    cap = inventory['elbow:FJ1477_proximal_rigid_cap']
    adapter_key = 'elbow:FJ1477_bone_adapter'
    old_adapter = inventory[adapter_key]
    bore = axial(-3.2, 2.41, 1.65)
    adapter, adapter_source = export(
        (old_adapter-bore).as_original().simplify(1e-6),
        folder/'anatomical_elbow_angled_bored_adapter.obj', face)
    inventory[adapter_key] = adapter
    for label, muscle in baseline['muscles'].items():
        inventory[f'wire_guard:{label}'] = A.D.outer_envelope(np.array(muscle['wire_path_source_mm']), 2.)

    hardware = {}
    for label, solid in (('sleeve', sleeve), ('screw', screw)):
        hardware[label], sources[label] = export(solid, folder/f'anatomical_elbow_angled_{label}.obj', face)
    checks = {**hardware, 'sleeve_insertion': sleeve_insertion,
              'screw_insertion': screw_insertion, 'driver': driver}
    hits = {label: {key: volume for key, part in inventory.items()
                    if (volume := float(A.I.intersection_volume(solid, part))) > .001}
            for label, solid in checks.items()}
    backing = axial(sleeve_bottom, sleeve_top, 1.85)-axial(-3.2, -1.2, 1.65)
    bearing = axial(2.20, 2.40, 2.50)-axial(2.19, 2.41, 1.65)
    engagement_fraction = A.I.intersection_volume(thread_envelope, hardware['sleeve'])/thread_envelope.volume()
    sleeve_backing_fraction = A.I.intersection_volume(backing, cap)/backing.volume()
    screw_bearing_fraction = A.I.intersection_volume(bearing, adapter)/bearing.volume()
    internal_overlap = A.I.intersection_volume(hardware['sleeve'], hardware['screw'])
    added = adapter-old_adapter
    adapter_added_mm3 = max(0., float(added.volume()))
    passed = (not any(hits.values()) and internal_overlap <= .001
              and adapter_added_mm3 <= .001
              and min(engagement_fraction, sleeve_backing_fraction, screw_bearing_fraction) >= .99)
    result = dict(interface_audit_sha256=sha(audit_file), source_axis=axis.tolist(),
                  source_face_mm=face.tolist(), solids={key: sources[key] for key in hardware},
                  bored_adapter_source=adapter_source,
                  adapter_stock_removed_mm3=float((old_adapter-adapter).volume()),
                  adapter_added_mm3=adapter_added_mm3,
                  sleeve_depth_range_mm=[-sleeve_top, -sleeve_bottom],
                  nominal_thread_engagement_mm=sleeve_top-screw_tip,
                  thread_envelope_engagement_fraction=float(engagement_fraction),
                  sleeve_backing_fraction=float(sleeve_backing_fraction),
                  head_bearing_fraction=float(screw_bearing_fraction),
                  hardware_internal_overlap_mm3=float(internal_overlap),
                  checked_obstacles=len(inventory), overlap_mm3=hits, neutral_screen_passed=passed,
                  scope='Nominal neutral geometry and straight insertion. Bonded sleeve, idealized thread core, '
                        'no helical thread, assembly tolerances, pull-out, or motion verification.')
    report = folder/'anatomical_elbow_angled_fastener.json'
    save(report, result)
    if not passed:
        raise SystemExit('Nominal fastener geometry failed the neutral screen')
    assembled_sources = {key: source for key, source in sources.items()
                         if key not in {*hardware, *audit['omitted_original_left_fasteners']}}
    assembled_sources[adapter_key] = adapter_source
    assembled_sources.update({f'elbow:FJ1477_proximal_angled_{key}': sources[key] for key in hardware})
    save(folder/'assembly_anatomical_elbow_angled_candidate.json', dict(
        scope='Shared humerus with original elbow TPU collision proxies and nominal angled fastener.',
        sources=assembled_sources, n_instances=len(assembled_sources),
        common_humerus_geometry_unified=assembly['common_humerus_geometry_unified'],
        fastener_report_sha256=sha(report), combined_assembly_verified=False, public_release_ready=False))
    print(result)


if __name__ == '__main__':
    main()
