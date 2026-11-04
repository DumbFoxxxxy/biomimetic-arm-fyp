"""Pair independent fastener directions and check mutual insertion clearance."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    screening = folder / 'collar_staggered_fastener_screen.json'
    report = read(screening)
    assert report['assembly_sha256'] == sha(source)
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    assert report['frame_sha256'] == sha(frame_file)
    frame = read(frame_file)['frame_source']
    basis = np.array(frame['basis_columns'])
    sources = read(source)['sources']
    manifold = A.I._manifold().Manifold
    candidates = [[], []]
    for index in range(2):
        keys = [f'collar:{name}_{index}' for name in ('bolt', 'top_washer', 'nut', 'bottom_washer')]
        parts = [load_component(sources[key]) for key in keys]
        for row in report['rows']:
            if any(stage['collision'] for stage in row['stages'][index*2:index*2+2]):
                continue
            angle, shift = row['angle_deg'], row['axial_shift_mm']
            transform = A.I.rotation_about(frame['origin_mm'], basis[:, 2], angle)
            transform[:, 3] += basis[:, 2] * shift
            axis = transform[:, :3] @ basis[:, 1]
            bodies = [part.transform(transform) for part in parts]
            sweep = [manifold.batch_hull([body, body.translate(axis * (20 if i < 2 else -20))])
                     for i, body in enumerate(bodies)]
            candidates[index].append((angle, shift, bodies, sweep))
    rows = []
    for angle0, shift0, bodies0, sweep0 in candidates[0]:
        for angle1, shift1, bodies1, sweep1 in candidates[1]:
            clearance = []
            for fixed, moving in [(bodies0, sweep1), (bodies1, sweep0)]:
                clearance.append(max(float(A.I.intersection_volume(a, b)) for a in fixed for b in moving))
            rows.append(dict(angles_deg=[angle0, angle1], axial_shifts_mm=[shift0, shift1],
                             angle_difference_deg=abs((angle1-angle0+180) % 360-180),
                             peak_conservative_overlap_mm3=clearance,
                             clear_orders=[order for order, volume in zip(([0, 1], [1, 0]), clearance)
                                           if volume <= .001]))
    passing = sorted((row for row in rows if row['clear_orders']),
                     key=lambda row: (row['angle_difference_deg'],
                                      sum(abs(value) for value in row['axial_shifts_mm'])))
    save(folder / 'collar_independent_ears_screen.json', dict(
        assembly_sha256=sha(source), screening_sha256=sha(screening), script_sha256=sha(Path(__file__)),
        rows=rows, passing_count=len(passing), best_candidates=passing[:20],
        scope='Convex swept hardware envelopes against the other installed group. Ear stock, joint contact and tool access pending.'))
    print('pairs', len(rows), 'passing', len(passing), 'best', passing[:3], flush=True)


if __name__ == '__main__':
    main()
