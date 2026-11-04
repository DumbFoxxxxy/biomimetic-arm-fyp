"""Round the elbow bypass waypoint and check its neutral routing envelope."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def rounded_path(start, corner, end, radius):
    incoming = A.D.unit(corner-start)
    outgoing = A.D.unit(end-corner)
    normal = A.D.unit(np.cross(incoming, outgoing))
    turn = float(np.arccos(np.clip(incoming @ outgoing, -1, 1)))
    setback = radius*np.tan(turn/2)
    if setback >= min(np.linalg.norm(corner-start), np.linalg.norm(end-corner)):
        raise ValueError('Fillet exceeds available straight sections')
    tangent = corner-incoming*setback
    centre = tangent+np.cross(normal, incoming)*radius
    count = max(2, int(np.ceil(np.rad2deg(turn)/5)))
    angles = np.linspace(0, turn, count+1)
    radial = tangent-centre
    arc = centre + np.outer(np.cos(angles), radial) + np.outer(np.sin(angles), np.cross(normal, radial))
    sagitta = radius*(1-np.cos(turn/(2*count)))
    return np.vstack([start, arc, end]), float(sagitta), float(turn)


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder/'assembly_tpu_readback_candidate.json'
    chords = folder/'elbow_relocated_wire_chord_screen.json'
    stations = folder/'elbow_proximal_rotated_station_screen.json'
    model, routes, screen = map(read, (manifest, chords, stations))
    assert routes['assembly_sha256'] == screen['assembly_sha256'] == sha(manifest)
    assert routes['station_screen_sha256'] == sha(stations)
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    assert screen['frame_sha256'] == sha(frame_file)
    frame = read(frame_file)['frame_source']
    axis = np.array(frame['basis_columns'])[:, 2]
    solids = {key: load_component(row) for key, row in model['sources'].items()}
    mobile = {key for group in screen['moving_groups'].values() for key in group}
    fixed = set(solids)-mobile-{'tpu:elbow_flexor', 'tpu:elbow_extensor'}
    rows = []
    for pair in routes['waypoint_pairs']:
        if pair['relocated_hardware_hits']:
            continue
        obstacles = {key: solids[key] for key in sorted(fixed)}
        endpoints = {}
        for label, head in [('flexor', 'FJ1486'), ('extensor', 'FJ1477')]:
            angle, shift = pair[f'{label}_angle_deg'], pair[f'{label}_shift_mm']
            row = next(row for row in routes['rows'][label] if row['angle_deg'] == angle and row['shift_mm'] == shift)
            endpoints[label] = (np.array(row['proximal_source_mm']), np.array(row['distal_source_mm']))
            transform = A.I.rotation_about(frame['origin_mm'], axis, angle)
            transform[:, 3] += axis*shift
            for key in screen['moving_groups'][head]:
                obstacles[key] = solids[key].transform(transform)
        ext_start, ext_end = endpoints['extensor']
        ext_direction = A.D.unit(ext_end-ext_start)
        obstacles['probe:extensor'] = A.D.outer_envelope(np.stack([
            ext_start+ext_direction*8, ext_end-ext_direction*8]), 2.)
        proximal, distal = endpoints['flexor']
        waypoint = np.array(pair['flexor_waypoint']['waypoint_source_mm'])
        start = proximal+A.D.unit(waypoint-proximal)*8
        end = distal-A.D.unit(distal-waypoint)*8
        for radius in [4., 6., 8., 12., 16., 20.]:
            path, sagitta, turn = rounded_path(start, waypoint, end, radius)
            envelope = A.D.outer_envelope(path, 2.+sagitta)
            hits = []
            for key, solid in obstacles.items():
                volume = float(A.I.intersection_volume(envelope, solid))
                if volume > .001:
                    hits.append(dict(target=key, overlap_mm3=volume))
            rows.append(dict(pair=pair, radius_mm=radius, arc_turn_deg=float(np.rad2deg(turn)),
                             sampling_sagitta_mm=sagitta, trimmed_path_source_mm=path.tolist(),
                             flexor_endpoints_source_mm=[proximal.tolist(), distal.tolist()],
                             extensor_endpoints_source_mm=[ext_start.tolist(), ext_end.tolist()],
                             hits=hits))
    passing = [row for row in rows if not row['hits']]
    save(folder/'elbow_route_fillet_screen.json', dict(
        assembly_sha256=sha(manifest), chord_sha256=sha(chords), station_sha256=sha(stations),
        script_sha256=sha(Path(__file__)), rows=rows, passing=passing, threshold_mm3=.001,
        scope='Neutral route fillet with 2 mm envelope plus chord sagitta. Terminal 8 mm and guide structure pending.'))
    print('fillet trials', len(rows), 'clear', len(passing), 'largest clear radius',
          max((row['radius_mm'] for row in passing), default=None), flush=True)


if __name__ == '__main__':
    main()
