"""Measure mesh-edge distortion for an imposed elbow lattice displacement."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, mesh_geom, read, save, sha


def main():
    base = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    report = base/'extensor_late_reaction_audit.json'
    bridge = read(base/'elbow_relocated_bridge_screen.json')
    assembly = read(base/'assembly_tpu_readback_candidate.json')['sources']
    joint = next(row['combined_key_and_seat'] for row in read(report)['planar_bearing_trials']
                 if row.get('combined_key_and_seat', {}).get('sources'))
    channels = read(base/'elbow_relocated_channel_cores/channel_core_check.json')
    centre, axis = A.frames.elbow_axis()
    results = {}
    for label, head in [('flexor', 'FJ1486'), ('extensor', 'FJ1477')]:
        source = (joint['sources']['tpu'] if label == 'extensor' else next(
            terminal['source'] for row in bridge['flexor_lattice_midspan']
            for terminal in row.get('terminal_trials', []) if 'source' in terminal))
        mesh = load_component(source)
        vertices, faces = mesh_geom.manifold_mesh_to_arrays(mesh.to_mesh64())
        edges = np.unique(np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1), axis=0)
        lengths = np.linalg.norm(vertices[edges[:, 1]]-vertices[edges[:, 0]], axis=1)
        kept = lengths >= .2
        start, end = np.array(channels['choice'][f'{label}_endpoints_source_mm'])
        if label == 'extensor':
            start = np.array(read(base/'extensor_guide_body_audit.json')['guide_source_mm'])
            end -= A.D.unit(end-start)*7.
        else:
            neutral_path = np.array(next(row['path_source_mm'] for row in channels['channels']
                                         if row['label'] == 'flexor'))
            start += A.D.unit(neutral_path[0]-start)
            end += A.D.unit(neutral_path[-1]-end)
        chord = end-start
        fraction = np.clip((vertices-start)@chord/(chord@chord), 0., 1.)
        rotor = assembly[f'elbow:{head}_cap_rotor']
        screen = read((ROOT/rotor['file']).parent/'screen.json')
        pivot = np.array(screen['pivot_zero_display_mm'])+np.array(rotor['translation_source_mm'])
        hinge = screen['trunnion_axis_zero_display']
        rows = []
        for pose in joint['revised_elbow_motion']['tracked_swivels']['rows'][label]:
            angle = pose['elbow_deg']
            if angle % 5:
                continue
            turn = A.I.rotation_about(pivot, hinge, pose['local_cable_only_equilibrium']['swivel_deg'])
            bend = A.I.rotation_about(centre, axis, angle)
            tip = bend[:, :3]@(turn[:, :3]@end+turn[:, 3])+bend[:, 3]
            posed = vertices+np.outer(fraction, tip-end)
            edge_lengths = np.linalg.norm(posed[edges[kept, 1]]-posed[edges[kept, 0]], axis=1)
            ratios = edge_lengths/lengths[kept]
            sampled_edges = edges[kept]
            hotspots = {}
            for name, index in [('shortest_ratio', np.argmin(ratios)), ('longest_ratio', np.argmax(ratios))]:
                edge = sampled_edges[index]
                hotspots[name] = dict(ratio=float(ratios[index]),
                                      neutral_midpoint_mm=vertices[edge].mean(axis=0).tolist(),
                                      normalized_station=float(fraction[edge].mean()))
            jacobian = 1.+float(chord@(tip-end)/(chord@chord))
            rows.append(dict(elbow_deg=angle, axial_map_determinant=jacobian,
                             hotspots=hotspots,
                             edge_length_ratio_quantiles=np.quantile(ratios, [0., .05, .5, .95, 1.]).tolist()))
        results[label] = dict(source=source, mesh_edges_used=int(kept.sum()), poses=rows)
    save(base/'elbow_tpu_geometry_demand.json', dict(
        motion_report_sha256=sha(report), rows=results, step_deg=5., minimum_sampled_edge_mm=.2,
        scope='Imposed axial interpolation from fixed proximal reference to moving distal tip. Edge distortion is a geometric proxy, not TPU material strain or an equilibrium result.'))
    print({key: row['poses'][-1] for key, row in results.items()})


if __name__ == '__main__':
    main()
