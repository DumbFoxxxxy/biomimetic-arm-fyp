"""Screen parallel relocation of collar fastener group 1."""
import argparse

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reverse', action='store_true')
    reverse = parser.parse_args().reverse
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    proof_file = folder / 'verification_checks.json'
    study_file = folder / 'collar_first_screen.json'
    proof, study = read(proof_file), read(study_file)
    assert proof['assembly_sha256'] == study['assembly_sha256'] == sha(source)
    solids = {key: load_component(row) for key, row in read(source)['sources'].items()}
    solids['collar:jaw_upper'] = load_component(study['candidate_mesh'])
    frame = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    basis = np.array(read(frame)['frame_source']['basis_columns'])
    opening = basis[:, 1]
    keys = ['collar:bolt_1', 'collar:top_washer_1', 'collar:nut_1', 'collar:bottom_washer_1']
    top_washer = solids[keys[1]]
    bottom_washer = solids[keys[3]]
    if reverse:
        local_frame = np.c_[basis.T, np.zeros(3)]
        top_low = top_washer.transform(local_frame).bounding_box()[1]
        bottom_high = bottom_washer.transform(local_frame).bounding_box()[4]
        reflection = np.eye(3) - 2 * np.outer(opening, opening)
        transform = np.c_[reflection, opening * (top_low + bottom_high)]
        for key in keys:
            solids[key] = solids[key].transform(transform)
    jaws = {'collar:jaw_upper', 'collar:jaw_lower'}
    fixed = (set(proof['fixed_keys'] + proof['moving_keys']) - jaws) | {'collar:liner_upper'}
    fixed |= {key.replace('_1', '_0') for key in keys}
    assert not fixed.intersection(keys)
    muscle = solids['tpu:deltoid']
    manifold = A.I._manifold().Manifold
    bolt_local = solids[keys[0]].transform(np.c_[basis.T, np.zeros(3)])
    low, high = np.array(bolt_local.bounding_box()).reshape(2, 3)
    centre = (low + high) / 2
    rows = []
    for x in range(-8, 9, 2):
        for z in range(-8, 9, 2):
            shift = basis @ [x, 0, z]
            obstruction = 0.
            for index, key in enumerate(keys):
                axis = opening if (index < 2) != reverse else -opening
                for offset in np.arange(0, 20.001, 1):
                    volume = A.I.intersection_volume(solids[key].translate(shift + axis * offset), muscle)
                    obstruction = max(obstruction, float(volume))
            row = dict(x_mm=x, z_mm=z, shift_source_mm=shift.tolist(), peak_muscle_overlap_mm3=obstruction)
            rows.append(row)
            if obstruction > .001:
                continue
            hits = []
            for index, key in enumerate(keys):
                axis = opening if (index < 2) != reverse else -opening
                for offset in np.arange(0, 20.001, .5):
                    posed = solids[key].translate(shift + axis * offset)
                    for target in sorted(fixed):
                        volume = float(A.I.intersection_volume(posed, solids[target]))
                        if volume > .001:
                            hits.append(dict(key=key, target=target, offset_mm=float(offset), volume_mm3=volume))
                    if hits:
                        break
                if hits:
                    break
            row['other_path_hits'] = hits
            # The 3.7 mm radius is a screening witness for a 1.7 mm bore plus 2 mm of stock.
            # The cylinder is evaluated before drilling, over 2 mm inward from each washer.
            stock = {}
            for key, washer_part, side in [('collar:jaw_upper', top_washer, 1), ('collar:jaw_lower', bottom_washer, -1)]:
                washer = washer_part.transform(np.c_[basis.T, np.zeros(3)])
                wlow, whigh = np.array(washer.bounding_box()).reshape(2, 3)
                seat = wlow[1] if side == 1 else whigh[1]
                bottom = seat - 2 if side == 1 else seat
                witness = manifold.cylinder(2, 3.7, circular_segments=64).transform(
                    np.c_[basis[:, [0, 2, 1]], basis @ [centre[0]+x, bottom, centre[2]+z]])
                stock[key] = float(A.I.intersection_volume(solids[key], witness) / witness.volume())
            row['stock_witness_coverage'] = stock
    clear = [row for row in rows if row.get('other_path_hits') == []]
    clear.sort(key=lambda row: min(row['stock_witness_coverage'].values()), reverse=True)
    report_name = 'collar_bolt_reverse_relocation_screen.json' if reverse else 'collar_bolt_relocation_screen.json'
    save(folder / report_name, dict(
        assembly_sha256=sha(source), study_sha256=sha(study_file), placement_proof_sha256=sha(proof_file),
        bolt_head_on_lower_side=reverse, grid_step_mm=2, stock_radius_mm=3.7, stock_depth_mm=2, rows=rows,
        clear_hardware_path_candidates=clear, geometry_adopted=False))
    print('clear paths', len(clear), 'best', clear[:3], flush=True)


if __name__ == '__main__':
    main()
