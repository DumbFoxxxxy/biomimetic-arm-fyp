"""Place computed free-space beam states back into anatomy for contact screening."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, read, save, load_component
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    source_file = folder/'elastic_conduit_bench/FJ1468_E15_T5_ds5_dt2e-05.json'
    gravity_file = folder/'elastic_conduit_bench/FJ1468_E15_T0_ds5_dt2e-05.json'
    source, gravity = read(source_file), read(gravity_file)
    for report in (source, gravity):
        assert provenance(ROOT/report['generator']['file']) == report['generator']
        assert provenance(ROOT/report['source']['file']) == report['source']
    reference_file = folder/'thorax_arm_reference.json'
    cervical_file = folder/'endpoint_neutral_clearance.json'
    sources = dict(read(reference_file)['sources'])
    sources.update(read(cervical_file)['cervical_context'])
    solids = {key: load_component(value) for key, value in sources.items()}
    trace = source['result']['trace']
    cases = [('maximum_displacement', max(trace, key=lambda r: r['max_node_displacement_mm']), source),
             ('loaded_0_24s', min(trace, key=lambda r: abs(r['time_s']-.24)), source),
             ('recovery_last_sample', trace[-1], source),
             ('gravity_only_0_24s', min(gravity['result']['trace'], key=lambda r: abs(r['time_s']-.24)), gravity)]
    rows = []
    for label, state, report in cases:
        points = (np.array(state['nodes_local_m'])+np.array(report['specification']['source_origin_m']))*1000.
        wire, envelope = A.D.outer_envelope(points, .3), A.D.outer_envelope(points, 1.6)
        wire_hits, guide_hits = {}, {}
        for key, solid in solids.items():
            volume = A.I.intersection_volume(envelope, solid)
            if volume > .001:
                guide_hits[key] = volume
                inner = A.I.intersection_volume(wire, solid)
                if inner > .001:
                    wire_hits[key] = inner
        rows.append(dict(label=label, time_s=state['time_s'],
            tension_command_N=state['tension_command_N'], points_source_mm=points.tolist(),
            wire_bone_overlap_mm3=wire_hits, constant_radius_guide_envelope_overlap_mm3=guide_hits,
            anatomy_contact_needs_resolution=bool(wire_hits or guide_hits)))
        print(label, 'wire hits', list(wire_hits), 'guide envelope hits', list(guide_hits), flush=True)
    save(folder/'elastic_conduit_bench/anatomy_snapshot_screen.json', dict(
        loaded_report=provenance(source_file), gravity_report=provenance(gravity_file),
        anatomy_reference=provenance(reference_file), cervical_context=provenance(cervical_file),
        parts_checked=len(solids), wire_radius_mm=.3, assumed_guide_outer_radius_mm=1.6,
        threshold_mm3=.001, snapshots=rows,
        scope='Four computed free-space beam states only, not continuous collision verification. Fixed endpoint boundary bench, not whole-arm motion. Guide envelope is a constant-radius capsule chain, not solved TPU contact or Poisson deformation. Overlap is a contact-resolution flag, not automatic soft-material failure.'))


if __name__ == '__main__':
    main()
