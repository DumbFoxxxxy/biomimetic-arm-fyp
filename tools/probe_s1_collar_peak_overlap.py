"""Measure the largest recorded collision against fastener neighborhoods."""
import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    screening_file = folder / 'late_collar_relief_screen.json'
    screening = read(screening_file)
    assert screening['assembly_sha256'] == sha(source)
    sources = read(source)['sources']
    peak = max(screening['jaw_deltoid_interference_samples'], key=lambda row: row['volume_mm3'])
    direction = np.array(screening['nominal_upper']['direction'])
    jaw = load_component(sources['collar:jaw_upper'])
    muscle = load_component(sources['tpu:deltoid'])
    overlap = jaw ^ muscle.translate(-direction * peak['offset_mm'])
    assert abs(overlap.volume() - peak['volume_mm3']) < .001
    cube = A.I._manifold().Manifold.cube
    neighborhoods = {}
    for key in ('collar:bolt_0', 'collar:bolt_1', 'link:FJ1500/bone_bolt'):
        low, high = np.array(load_component(sources[key]).bounding_box()).reshape(2, 3)
        guard = cube(high - low + 2).translate(low - 1)
        neighborhoods[key] = float(A.I.intersection_volume(overlap, guard))
    save(folder / 'late_collar_peak_overlap.json', dict(
        assembly_sha256=sha(source), screening_sha256=sha(screening_file),
        jaw_volume_mm3=float(jaw.volume()), peak=peak,
        peak_overlap_within_hardware_boxes_mm3=neighborhoods,
        interpretation='Single-pose lower bound for any subtraction-only correction on this path.'))
    print(peak, neighborhoods)


if __name__ == '__main__':
    main()
