"""Screen inherited supported curves; do not transfer their dynamic approval."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save
from build_channel_first_inventory import provenance


def points(values):
    return np.array([np.fromstring(v, sep=' ') if isinstance(v, str) else v for v in values])


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    reference_file = folder/'thorax_arm_reference.json'
    previous_file = folder/'endpoint_neutral_clearance.json'
    route_file = ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_port_path_revision/routes_compact_distal_candidate.json'
    sources = dict(read(reference_file)['sources'])
    sources.update(read(previous_file)['cervical_context'])
    solids = {key: load_component(source) for key, source in sources.items()}
    rows = []
    for route in read(route_file)['routes']:
        segments = {
            'origin_transition': points(route['connectors']['origin']['path_source_mm']),
            'core': points(route['core_path_source_mm']),
            'distal_transition': points(route['connectors']['distal']['path_source_mm'])[::-1],
        }
        results = {}
        for name, path in segments.items():
            assert np.isfinite(path).all() and len(path) > 1
            wire, channel = A.D.outer_envelope(path, .3), A.D.outer_envelope(path, 1.6)
            hits, wire_hits = {}, {}
            gap, nearest = 2., None
            for key, solid in solids.items():
                volume = A.I.intersection_volume(channel, solid)
                if volume > .001:
                    hits[key] = volume
                    wire_volume = A.I.intersection_volume(wire, solid)
                    if wire_volume > .001:
                        wire_hits[key] = wire_volume
                distance = A.I.gap(channel, solid, gap)
                if distance < gap:
                    gap, nearest = distance, key
            results[name] = dict(points_source_mm=path.tolist(),
                length_mm=float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum()),
                wire_bone_overlap_mm3=wire_hits, channel_bone_overlap_mm3=hits,
                minimum_channel_gap_mm_capped2=gap, nearest_part=nearest)
        changed = route['head'] in ('FJ1500', 'FJ1504')
        junctions = [float(np.linalg.norm(segments['origin_transition'][-1]-segments['core'][0])),
                     float(np.linalg.norm(segments['core'][-1]-segments['distal_transition'][0]))]
        assert max(junctions) < 1e-6, 'Unscreened gap between route segments'
        row = dict(id=route['head'], segments=results,
            segment_junction_gaps_mm=junctions,
            current_topology_compatible=not changed,
            status='historical_only_Y2_terminal_topology_changed' if changed else 'inherited_curve_neutral_screen_only',
            physical_support='Existing TPU core and end-transition concept; fixed attachments, reaction capacity and moving channel law pending',
            core_wire_clear=not results['core']['wire_bone_overlap_mm3'],
            core_channel_clear=not results['core']['channel_bone_overlap_mm3'],
            interface_note='Transitions checked against solid bones without subtraction of unverified bores. Any overlap requires an explicit supported interface, not a bone-pair exclusion.')
        rows.append(row)
        print(row['id'], 'core wire hits', list(results['core']['wire_bone_overlap_mm3']),
              'core channel hits', list(results['core']['channel_bone_overlap_mm3']), flush=True)
    save(folder/'curved_shoulder_neutral_clearance.json', dict(
        source_routes=provenance(route_file), anatomy_reference=provenance(reference_file),
        cervical_context_source=provenance(previous_file), checked_parts=len(solids),
        wire_radius_mm=.3, channel_outer_radius_mm=1.6, threshold_mm3=.001,
        routes=rows, scope='Neutral capsule-chain Boolean screen of inherited curves only; no motor supply, bore approval, motion clearance, friction or force validation. Backboard excluded.'))


if __name__ == '__main__':
    main()
