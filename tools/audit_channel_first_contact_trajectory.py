"""Audit complete wire capsules at every saved contact-bench trajectory state."""
import argparse

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, read, save, load_component
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stem', default='FJ1468_T5_k200000_mu0_u5')
    args = parser.parse_args()
    report_file = folder/'elastic_conduit_bone_contact'/f'{args.stem}.json'
    report = read(report_file)
    for key in ('source', 'generator', 'surface_query', 'beam_generator', 'model', 'anatomy_reference'):
        assert provenance(ROOT/report[key]['file']) == report[key]
    reference_file = folder/'thorax_arm_reference.json'
    context_file = folder/'endpoint_neutral_clearance.json'
    sources = dict(read(reference_file)['sources'])
    sources.update(read(context_file)['cervical_context'])
    anatomy = {name: load_component(source) for name, source in sources.items()}
    origin = np.array(report['specification']['source_origin_m'])
    rows = []
    for i, state in enumerate(report['result']['trace']):
        points = (np.array(state['nodes_local_m'])+origin)*1000.
        wire = A.D.outer_envelope(points, .3)
        hits = {}
        gap, nearest = 1.5, None
        for name, solid in anatomy.items():
            volume = A.I.intersection_volume(wire, solid)
            if volume > .001:
                hits[name] = volume
            value = A.I.gap(wire, solid, gap)
            if value < gap:
                gap, nearest = value, name
        rows.append(dict(time_s=state['time_s'], wire_bone_overlap_mm3=hits,
                         wire_gap_mm_capped1_5=gap, nearest_part=nearest))
        if i % 50 == 0:
            print('Audited', i+1, 'of', len(report['result']['trace']), 'saved states', flush=True)
    output = dict(source=provenance(report_file), anatomy_reference=provenance(reference_file),
        cervical_context=provenance(context_file), auditor=provenance(ROOT/'tools/audit_channel_first_contact_trajectory.py'),
        checked_parts=len(anatomy), checked_states=len(rows), wire_radius_mm=.3, threshold_mm3=.001,
        maximum_saved_time_interval_s=float(np.diff([r['time_s'] for r in rows]).max()),
        states_with_wire_bone_overlap=sum(bool(r['wire_bone_overlap_mm3']) for r in rows),
        minimum_wire_gap_mm_capped1_5=min(r['wire_gap_mm_capped1_5'] for r in rows), states=rows,
        scope='Whole piecewise-linear wire capsule at every saved ~2 ms state; not a mathematical all-time proof or a deformable wall/contact validation. All 64 reference anatomy parts considered, only humerus/scapula exert contact forces in the simulated bench. No moving bones or mounted final interfaces.')
    save(report_file.with_name(f'{args.stem}_trajectory_audit.json'), output)
    print('Overlap states', output['states_with_wire_bone_overlap'], 'minimum wire gap mm',
          output['minimum_wire_gap_mm_capped1_5'], flush=True)


if __name__ == '__main__':
    main()
