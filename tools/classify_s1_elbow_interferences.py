"""Map shoulder-elbow intersections to deltoid terminals and wire protection."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, rod, save, sha


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_check.json'
    audit = folder / 'inter_module_neutral_audit.json'
    model, collisions = read(manifest), read(audit)
    assert collisions['assembly_sha256'] == sha(manifest)
    route_file = A.DESIGN / 's1_port_path_revision/routes_compact_distal_candidate.json'
    shoulder = read(A.DESIGN / 's1_integrated_humerus_anchors/assembly_check.json')
    assert shoulder['routes_sha256'] == sha(route_file)
    routes = [row for row in read(route_file)['routes'] if row['head'] in ('FJ1468', 'FJ1513', 'FJ1467')]
    guards = {row['head']: A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), 2.)
              for row in routes}
    distal = routes[0]['ports']['distal']
    centre = np.array(distal['inner_face_source_mm']) - np.array(distal['offset_from_seat_mm'])
    axis = A.D.unit(distal['normal_source'])
    terminal = rod(centre-axis*18, centre+axis*.2, 18)
    keys = {row[side] for row in collisions['overlaps'] for side in ('shoulder', 'elbow')}
    solids = {key: load_component(model['sources'][key]) for key in keys}
    rows = []
    for row in collisions['overlaps']:
        overlap = solids[row['shoulder']] ^ solids[row['elbow']]
        volume = float(overlap.volume())
        assert abs(volume-row['overlap_mm3']) < .001
        protected = {key: float(A.I.intersection_volume(overlap, guard)) for key, guard in guards.items()}
        rows.append(dict(**row, bounds_source_mm=list(overlap.bounding_box()),
                         outside_distal_terminal_mm3=max(0., float((overlap-terminal).volume())),
                         protected_wire_overlap_mm3=protected,
                         fraction_of_elbow_part=volume/solids[row['elbow']].volume()))
    save(folder/'interference_material_map.json', dict(
        assembly_sha256=sha(manifest), neutral_audit_sha256=sha(audit), routes_sha256=sha(route_file),
        script_sha256=sha(Path(__file__)), distal_terminal_radius_mm=18, wire_guard_radius_mm=2,
        rows=rows, unchecked_sources=collisions['unchecked_nonmanifold_sources']))
    for row in rows:
        print(row['shoulder'], row['elbow'], round(row['outside_distal_terminal_mm3'], 3),
              {key: round(value, 3) for key, value in row['protected_wire_overlap_mm3'].items()},
              round(row['fraction_of_elbow_part'], 3))


if __name__ == '__main__':
    main()
