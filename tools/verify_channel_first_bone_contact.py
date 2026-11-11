"""Verify contact distances/signs and distributed quadrature independently of runs."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, read, save, load_component
from build_channel_first_inventory import provenance
from channel_first_bone_contact import BoneContactSurface
from extract_attachment_regions import projector
from simulate_channel_first_conduit_bone_contact import contact_stations


def check(surface, points, unsigned, inside):
    maximum_error, near_count, rejected_deep = 0., 0, 0
    for point, exact, is_inside in zip(points, unsigned, inside):
        try:
            distance, normal = surface.nearest(point)
        except ValueError:
            assert is_inside and exact > surface.band-1e-7
            rejected_deep += 1
            continue
        if exact < surface.band-1e-7:
            expected = -exact if is_inside else exact
            error = abs(distance-expected)
            maximum_error = max(maximum_error, error)
            assert error < 1e-7
            assert abs(np.linalg.norm(normal)-1.) < 1e-7
            near_count += 1
        else:
            assert not is_inside and abs(distance-surface.band) < 1e-7
    return dict(points=len(points), near_points=near_count, rejected_deep_inside_points=rejected_deep,
                maximum_signed_distance_error_mm=maximum_error)


def main():
    rng = np.random.default_rng(67103)
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    surface = BoneContactSurface(A.I._manifold().Manifold.cube([10., 10., 10.], center=True))
    points = rng.uniform(-8., 8., (500, 3))
    q = abs(points)-5.
    sdf = np.linalg.norm(np.maximum(q, 0.), axis=1)+np.minimum(q.max(axis=1), 0.)
    box = check(surface, points, abs(sdf), sdf < 0.)
    source_file = folder/'elastic_conduit_bench/FJ1468_E15_T5_ds5_dt2e-05.json'
    source = read(source_file)
    trace = source['result']['trace']
    pool = np.concatenate([np.array(trace[i]['nodes_local_m']) for i in (0, 25, 75, 120, 200)])
    samples = pool[rng.integers(0, len(pool), size=300)]*1000.
    samples += np.array(source['specification']['source_origin_m'])*1000.+rng.normal(0., 1.5, (300, 3))
    reference_file = folder/'thorax_arm_reference.json'
    reference = read(reference_file)
    rows = []
    for name in ('bone:humerus', 'bone:scapula'):
        surface = BoneContactSurface(load_component(reference['sources'][name]))
        _, unsigned, _ = projector(surface.vertices, surface.faces)(samples)
        inside = np.array([surface.ray_inside(point) for point in samples])
        row = check(surface, samples, unsigned, inside)
        rows.append(dict(bone=name, **row))
    lengths = np.array(source['specification']['segment_rest_lengths_m'])
    stations = contact_stations(lengths)
    assert abs(sum(s[3] for s in stations)-lengths.sum()) < 1e-12
    assert all(0. <= s[2] <= 1. and s[3] > 0. for s in stations)
    dest = folder/'elastic_conduit_bone_contact'
    dest.mkdir(exist_ok=True)
    save(dest/'surface_query_verification.json', dict(seed=67103, analytic_box=box,
        anatomy_checks=rows, contact_weight_sum_m=sum(s[3] for s in stations),
        beam_rest_length_m=float(lengths.sum()), anatomy_reference=provenance(reference_file),
        sample_trajectory=provenance(source_file), surface_module=provenance(ROOT/'tools/channel_first_bone_contact.py'),
        verification_script=provenance(ROOT/'tools/verify_channel_first_bone_contact.py'),
        scope='Analytic box signed distances; actual-mesh candidate-cache coverage compared with untruncated global triangle search sharing the point/triangle kernel; signs independently compared with native solid ray parity. This verifies these sampled queries, not all possible geometry or continuum contact.'))
    print('Analytic box:', box, 'anatomy:', rows, flush=True)


if __name__ == '__main__':
    main()
