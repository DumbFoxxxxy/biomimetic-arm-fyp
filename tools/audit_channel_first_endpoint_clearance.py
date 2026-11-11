"""Audit interior cable/channel envelopes independently of force feasibility."""
from itertools import combinations

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    endpoint_file = folder/'arm_girdle_equivalent_drive_geometry_mass_check.json'
    endpoint = read(endpoint_file)
    reference_file = folder/'thorax_arm_reference.json'
    reference = read(reference_file)
    solids, errors = {}, {}
    context_sources = dict(reference['sources'])
    support_folder = ROOT/'02_人体骨骼晶格肌肉方案/输入/骨骼/支撑参考'
    support_inventory_file = support_folder/'sources.json'
    cervical = {}
    for row in read(support_inventory_file)['models']:
        if 'cervical_spine' in row['categories']:
            file = support_folder/f"{row['id']}.obj"
            assert provenance(file)['sha256'] == row['sha256']
            cervical[f"cervical_reference:{row['id']}"] = dict(
                file=file.relative_to(ROOT).as_posix(), sha256=row['sha256'], translation_source_mm=[0., 0., 0.])
    context_sources.update(cervical)
    for key, source in context_sources.items():
        try:
            solids[key] = load_component(source)
        except (AssertionError, ValueError) as error:
            errors[key] = str(error) or type(error).__name__
    channels, rows = {}, []
    for route in endpoint['channels']:
        head = route['id']
        start, end = np.asarray(route['endpoints_source_mm'])
        direction = end-start
        length = np.linalg.norm(direction)
        direction /= length
        terminal_length = min(5., .1*length)
        interior = np.array([start+direction*terminal_length, end-direction*terminal_length])
        wire = A.D.outer_envelope(interior, .3)
        channel = A.D.outer_envelope(interior, 1.6)
        channels[head] = channel
        wire_hits, channel_hits = {}, {}
        minimum_gap = 2.
        nearest = None
        for key, part in solids.items():
            volume = A.I.intersection_volume(channel, part)
            if volume > .001:
                channel_hits[key] = volume
                wire_volume = A.I.intersection_volume(wire, part)
                if wire_volume > .001:
                    wire_hits[key] = wire_volume
            gap = A.I.gap(channel, part, minimum_gap)
            if gap < minimum_gap:
                minimum_gap, nearest = gap, key
        rows.append(dict(
            id=head, origin_body=route['origin_body'], distal_body=route['distal_body'],
            endpoints_source_mm=route['endpoints_source_mm'], straight_span_length_mm=float(length),
            excluded_terminal_zone_each_end_mm=float(terminal_length),
            interior_wire_hits_mm3=wire_hits, interior_channel_hits_mm3=channel_hits,
            channel_gap_mm_capped2=minimum_gap, nearest_part=nearest,
            neutral_anatomy_screen_clear=not channel_hits and not errors,
            terminal_attachment_review_pending=True))
    pairs = []
    for left, right in combinations(channels, 2):
        volume = A.I.intersection_volume(channels[left], channels[right])
        if volume > .001:
            pairs.append(dict(left=left, right=right, channel_envelope_overlap_mm3=volume))
    save(folder/'endpoint_neutral_clearance.json', dict(
        endpoint_report=provenance(endpoint_file), anatomy_reference=provenance(reference_file),
        wire_radius_mm=.3, channel_outer_radius_mm=1.6, threshold_mm3=.001,
        checked_anatomy_parts=len(solids), unreadable_solids=errors, channels=rows,
        cervical_context=cervical, support_inventory=provenance(support_inventory_file),
        interior_channel_pairs_with_overlap=pairs,
        scope='Neutral straight endpoint spans only. End zones excluded for separate attachment review. '
              'Anatomical context and engineered bones only; no source fasteners or moving contact. '
              'Channel envelope overlap may allow TPU compression but bare-wire/bone overlap needs routing.'))
    print('readable anatomy', len(solids), 'unreadable', len(errors), 'clear interior channels',
          sum(row['neutral_anatomy_screen_clear'] for row in rows), 'of', len(rows),
          'wire/bone flags', sum(bool(row['interior_wire_hits_mm3']) for row in rows),
          'interior channel pair overlaps', len(pairs))
    for row in rows:
        if row['interior_wire_hits_mm3']:
            print(row['id'], list(row['interior_wire_hits_mm3']))


if __name__ == '__main__':
    main()
