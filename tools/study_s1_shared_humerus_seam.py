"""Locate matching shaft regions between integrated shoulder and elbow humeri."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def main():
    folder = A.DESIGN / 's1_integrated_humerus_anchors'
    shoulder_file = folder / 'assembly_check.json'
    elbow_file = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_preview/integration_audit.json'
    shoulder = read(shoulder_file)
    assert shoulder['integrated_geometry_screen_pass']
    shoulder_bone = load_component(shoulder['sources']['bone:humerus'])
    elbow_bone = load_component(read(elbow_file)['original_elbow_humerus'])
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    frame = read(frame_file)['frame_source']
    basis = np.array(frame['basis_columns'])
    origin = np.array(frame['origin_mm'])
    local = np.c_[basis.T, -basis.T @ origin]
    shoulder_local, elbow_local = (body.transform(local) for body in (shoulder_bone, elbow_bone))
    cube = A.I._manifold().Manifold.cube
    rows = []
    for station in range(100, 241, 5):
        slab = cube([120, 120, 2]).translate([-60, -60, station-1])
        upper = shoulder_local ^ slab
        lower = elbow_local ^ slab
        common = upper ^ lower
        upper_volume, lower_volume = float(upper.volume()), float(lower.volume())
        if min(upper_volume, lower_volume) <= .001:
            continue
        difference = float((upper-lower).volume() + (lower-upper).volume())
        rows.append(dict(station_mm=station, shoulder_mm3=upper_volume, elbow_mm3=lower_volume,
                         common_mm3=float(common.volume()), symmetric_difference_mm3=difference,
                         relative_difference=difference / max(upper_volume, lower_volume)))
    best = sorted(rows, key=lambda row: row['relative_difference'])[:8]
    save(folder / 'shared_humerus_seam_screen.json', dict(
        shoulder_sha256=sha(shoulder_file), elbow_audit_sha256=sha(elbow_file), frame_sha256=sha(frame_file),
        script_sha256=sha(Path(__file__)), slab_thickness_mm=2, rows=rows, best=best,
        common_humerus_geometry_unified=False, complete_assembly_verified=False))
    print('best', best[:3], flush=True)


if __name__ == '__main__':
    main()
