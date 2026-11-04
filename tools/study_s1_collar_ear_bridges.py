"""Build and screen solid bridges for paired independent collar fasteners."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, rod, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    placement = folder / 'verification_checks.json'
    pairing = folder / 'collar_independent_ears_screen.json'
    proof, pairs = read(placement), read(pairing)
    assert proof['assembly_sha256'] == pairs['assembly_sha256'] == sha(source)
    sources = read(source)['sources']
    solids = {key: load_component(row) for key, row in sources.items()}
    solids['collar:jaw_upper'] = load_component(read(folder/'collar_first_screen.json')['candidate_mesh'])
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    frame = read(frame_file)['frame_source']
    basis = np.array(frame['basis_columns'])
    axis = basis[:, 1]
    local = np.c_[basis.T, np.zeros(3)]
    fixed = set(proof['fixed_keys'] + proof['moving_keys']) | {'collar:liner_upper'}
    fixed -= {'collar:jaw_upper', 'collar:jaw_lower'}
    manifold = A.I._manifold().Manifold
    rows = []
    for pair in pairs['best_candidates']:
        additions = {'collar:jaw_upper': [], 'collar:jaw_lower': []}
        for index, (angle, shift) in enumerate(zip(pair['angles_deg'], pair['axial_shifts_mm'])):
            transform = A.I.rotation_about(frame['origin_mm'], basis[:, 2], angle)
            transform[:, 3] += basis[:, 2] * shift
            for sign, side, washer in [(1, 'upper', 'top_washer'), (-1, 'lower', 'bottom_washer')]:
                low, high = np.array(solids[f'collar:{washer}_{index}'].transform(local).bounding_box()).reshape(2, 3)
                centre = (low + high) / 2
                centre[1] = (low if sign > 0 else high)[1] - sign * .15
                seat = basis @ centre
                stock = rod(seat - axis*sign*4, seat, 3.7)
                root = rod(seat - axis*sign*3.5, seat - axis*sign*.5, 2.5)
                bridge = manifold.batch_hull([root, stock.transform(transform)])
                bore = rod(seat-axis*10, seat+axis*10, 1.7).transform(transform)
                additions[f'collar:jaw_{side}'].append((bridge, bore))
        hits = []
        components = {}
        revised = {}
        for key, items in additions.items():
            jaw = solids[key]
            for bridge, _ in items:
                jaw += bridge
            for _, bore in items:
                jaw -= bore
            revised[key] = jaw
            components[key] = sum(part.volume() > 1e-8 for part in jaw.decompose())
            added = jaw - solids[key]
            for target in sorted(fixed):
                overlap = float(A.I.intersection_volume(added, solids[target]))
                if overlap > .001:
                    hits.append(dict(jaw=key, target=target, overlap_mm3=overlap))
        jaw_overlap = float(A.I.intersection_volume(revised['collar:jaw_upper'], revised['collar:jaw_lower']))
        rows.append(dict(pair=pair, jaw_components=components, fixed_hits=hits,
                         jaw_pair_overlap_mm3=jaw_overlap))
    save(folder / 'collar_ear_bridge_screen.json', dict(
        assembly_sha256=sha(source), placement_sha256=sha(placement), pairing_sha256=sha(pairing),
        script_sha256=sha(Path(__file__)), rows=rows, geometry_adopted=False,
        scope='Straight convex bridges from original ears to rotated seats. Final-pose added-material check.'))
    print('screened', len(rows), 'clear', sum(not r['fixed_hits'] and r['jaw_pair_overlap_mm3'] <= .001
          and all(n == 1 for n in r['jaw_components'].values()) for r in rows), flush=True)


if __name__ == '__main__':
    main()
