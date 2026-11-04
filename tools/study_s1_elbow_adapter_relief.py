"""Measure clearance cuts on elbow adapters against attachment material."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, read, save, sha


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_check.json'
    map_file = folder / 'interference_material_map.json'
    model, material = read(manifest), read(map_file)
    assert material['assembly_sha256'] == sha(manifest)
    mf = A.I._manifold()
    rows = []
    for head in ('FJ1477', 'FJ1486'):
        prefix = f'elbow:{head}'
        key = prefix+'_bone_adapter'
        adapter = load_component(model['sources'][key])
        conflicts = [row for row in material['rows'] if row['elbow'] == key]
        # Intersect with the adapter before union to keep the Boolean local.
        cuts = []
        for row in conflicts:
            obstacle = load_component(model['sources'][row['shoulder']])
            for shift in [np.zeros(3), *np.eye(3), *(-np.eye(3))]:
                cut = adapter ^ obstacle.translate(shift*.025)
                if not cut.is_empty():
                    cuts.append(cut)
        cut = mf.Manifold.batch_boolean(cuts, mf.OpType.Add)
        candidate = adapter-cut
        guards = {}
        for name in ('bone_0_bolt', 'bone_1_bolt', 'cap_rotor', 'pin_-1', 'pin_1'):
            part = load_component(model['sources'][f'{prefix}_{name}'])
            low, high = np.array(part.bounding_box()).reshape(2, 3)
            guard = mf.Manifold.cube(high-low+2).translate(low-1)
            guards[name] = dict(stock_mm3=float(A.I.intersection_volume(adapter, guard)),
                                removed_mm3=float(A.I.intersection_volume(cut, guard)))
        rows.append(dict(key=key, targets=[row['shoulder'] for row in conflicts],
                         removed_mm3=float(cut.volume()), adapter_mm3=float(adapter.volume()),
                         original_positive_components=sum(part.volume() > 1e-8 for part in adapter.decompose()),
                         candidate_positive_components=sum(part.volume() > 1e-8 for part in candidate.decompose()),
                         candidate_component_volumes_mm3=sorted(
                             (float(part.volume()) for part in candidate.decompose() if part.volume() > 1e-8),
                             reverse=True),
                         stock_in_expanded_interface_boxes=guards, candidate_exported=False))
    save(folder/'elbow_adapter_relief_screen.json', dict(
        assembly_sha256=sha(manifest), material_map_sha256=sha(map_file), script_sha256=sha(Path(__file__)),
        axis_shifts_mm=.025, box_margin_mm=1, rows=rows, geometry_adopted=False))
    print(rows, flush=True)


if __name__ == '__main__':
    main()
