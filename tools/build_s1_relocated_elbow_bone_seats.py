"""Build and screen four integral screw seats on the shared humerus."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha


def fit_seat(choice, normal, bone, fixed):
    origin = np.array(choice['axis_point_source_mm'])
    front = origin + normal * choice['front_boss_face_relative_mm']
    back = origin + normal * choice['back_boss_face_relative_mm']
    collar = rod(back, front, 4.5)
    hole = rod(back-normal, front+normal, 1.7)
    added = (collar-hole)-bone
    hits = []
    for key, part in fixed.items():
        volume = float(A.I.intersection_volume(added, part))
        if volume > .001:
            hits.append(dict(key=key, overlap_mm3=volume))
    return collar, hole, dict(
        transverse_shift_mm=choice['transverse_shift_mm'], axial_shift_mm=choice['axial_shift_mm'],
        front_face_mm=front.tolist(), back_face_mm=back.tolist(),
        added_mm3=float(added.volume()), fixed_hits=hits)


def main():
    base = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = base / 'assembly_tpu_readback_candidate.json'
    seat_report = base / 'elbow_relocated_seat_audit.json'
    station_report = base / 'elbow_proximal_rotated_station_screen.json'
    channels = base / 'elbow_relocated_channel_cores/channel_core_check.json'
    seats, stations, channel_data = map(read, (seat_report, station_report, channels))
    assert all(data['assembly_sha256'] == sha(manifest) for data in (seats, stations, channel_data))
    assert seats['channel_report_sha256'] == sha(channels)
    sources = read(manifest)['sources']
    bone = load_component(sources['bone:humerus'])
    excluded = set(sum(stations['moving_groups'].values(), []))
    excluded.update(('bone:humerus', 'tpu:elbow_flexor', 'tpu:elbow_extensor'))
    fixed = {key: load_component(source) for key, source in sources.items() if key not in excluded}
    for row in channel_data['channels']:
        fixed[f"new_core:{row['label']}"] = load_component(row['source'])
    folder = base / 'elbow_relocated_bone_seats'
    folder.mkdir(exist_ok=True)
    stock, holes, rows = [], [], []
    for row in seats['anchors']:
        normal = np.array(row['axis'])
        rejected = []
        for choice in row['footprint_candidates']:
            collar, hole, result = fit_seat(choice, normal, bone, fixed)
            if not result['fixed_hits']:
                break
            rejected.append(result)
        else:
            raise ValueError(f"No clear seat for {row['muscle']} anchor {row['anchor']}: {rejected}")
        stock.append(collar)
        holes.append(hole)
        rows.append(dict(muscle=row['muscle'], anchor=row['anchor'],
                         **result, rejected_candidates=rejected))
    mf = A.I._manifold()
    combined = mf.Manifold.batch_boolean([bone, *stock], mf.OpType.Add)
    bores = mf.Manifold.batch_boolean(holes, mf.OpType.Add)
    revised = (combined-bores).as_original().simplify(1e-6)
    actual, source = export(revised, folder / 'humerus_with_elbow_seats.obj', np.zeros(3))
    save(folder / 'bone_seat_check.json', dict(
        assembly_sha256=sha(manifest), seat_report_sha256=sha(seat_report),
        script_sha256=sha(Path(__file__)), source=source, seats=rows,
        removed_from_original_bone_mm3=float((bone-actual).volume()),
        added_to_original_bone_mm3=float((actual-bone).volume()),
        neutral_added_stock_pass=all(not row['fixed_hits'] for row in rows),
        stage_pass=False,
        scope='Neutral added stock against fixed solids and channel cores. No bridge, fasteners, wire sweep or motion verification.'))
    print([(row['muscle'], row['anchor'], row['added_mm3'], row['fixed_hits']) for row in rows])


if __name__ == '__main__':
    main()
