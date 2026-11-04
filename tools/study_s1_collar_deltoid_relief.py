"""Locate material needed for group-1 insertion relative to TPU end zones and wires."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, rod, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    model = read(source)
    sources = model['sources']
    muscle = load_component(sources['tpu:deltoid'])
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    opening = np.array(read(frame_file)['frame_source']['basis_columns'])[:, 1]
    route_file = A.DESIGN / 's1_port_path_revision/routes_compact_distal_candidate.json'
    assert model['routes_sha256'] == sha(route_file)
    routes = {row['head']: row for row in read(route_file)['routes']}
    own = [routes[head] for head in ('FJ1468', 'FJ1513', 'FJ1467')]
    mf = A.I._manifold()
    zones = []
    for route in own:
        port = route['ports']['origin']
        centre = np.array(port['inner_face_source_mm'])
        normal = A.D.unit(port['normal_source'])
        zones.append(rod(centre-normal*18, centre+normal*.2, 12))
    port = own[0]['ports']['distal']
    centre = np.array(port['inner_face_source_mm']) - np.array(port['offset_from_seat_mm'])
    normal = A.D.unit(port['normal_source'])
    zones.append(rod(centre-normal*18, centre+normal*.2, 18))
    terminal_zone = mf.Manifold.batch_boolean(zones, mf.OpType.Add)
    guards = {row['head']: A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), 2)
              for row in own}
    rows = []
    for part, sign in [('bolt', 1), ('top_washer', 1), ('nut', -1), ('bottom_washer', -1)]:
        key = f'collar:{part}_1'
        hardware = load_component(sources[key])
        axis = opening * sign
        swept_envelope = mf.Manifold.batch_hull([hardware, hardware.translate(axis*20)])
        local_muscle = muscle ^ swept_envelope
        hits = []
        for offset in np.arange(0, 20.001, .5):
            conflict = local_muscle ^ hardware.translate(axis*offset)
            volume = float(conflict.volume())
            if volume <= .001:
                continue
            hits.append(dict(offset_mm=float(offset), overlap_mm3=volume,
                             outside_terminal_mm3=float((conflict-terminal_zone).volume()),
                             wire_guard_mm3={head: float(A.I.intersection_volume(conflict, guard))
                                             for head, guard in guards.items()}))
        rows.append(dict(key=key, hits=hits,
                         peak=max(hits, key=lambda hit: hit['overlap_mm3']) if hits else None,
                         outside_terminal_peak_mm3=max((hit['outside_terminal_mm3'] for hit in hits), default=0),
                         wire_guard_peak_mm3={head: max((hit['wire_guard_mm3'][head] for hit in hits), default=0)
                                              for head in guards}))
    save(folder / 'collar_deltoid_relief_screen.json', dict(
        assembly_sha256=sha(source), frame_sha256=sha(frame_file), routes_sha256=sha(route_file),
        script_sha256=sha(Path(__file__)), rows=rows, threshold_mm3=.001,
        approach_mm=20, step_mm=.5, material_removed=False,
        scope='Per-pose cut lower bounds for existing group-1 hardware at final deltoid placement.'))
    for row in rows:
        print(row['key'], len(row['hits']), row['peak'], row['outside_terminal_peak_mm3'],
              row['wire_guard_peak_mm3'], flush=True)


if __name__ == '__main__':
    main()
