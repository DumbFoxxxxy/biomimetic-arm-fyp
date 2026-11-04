"""Screen an axial pin and interleaved knuckles for the group-1 collar seam."""
import argparse
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, rod, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offset-bridge', action='store_true')
    bridge = parser.parse_args().offset_bridge
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    proof_file = folder / 'verification_checks.json'
    study_file = folder / 'collar_first_screen.json'
    proof, study = read(proof_file), read(study_file)
    assert proof['assembly_sha256'] == study['assembly_sha256'] == sha(source)
    sources = read(source)['sources']
    solids = {key: load_component(row) for key, row in sources.items()}
    solids['collar:jaw_upper'] = load_component(study['candidate_mesh'])
    frame = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    basis = np.array(read(frame)['frame_source']['basis_columns'])
    pin_axis = basis[:, 2]
    direction = np.array(proof['direction_source'])
    local = np.c_[basis.T, np.zeros(3)]
    low, high = np.array(solids['collar:bolt_1'].transform(local).bounding_box()).reshape(2, 3)
    centre = (low + high) / 2
    centre[1] = (solids['collar:top_washer_1'].transform(local).bounding_box()[1]
                 + solids['collar:bottom_washer_1'].transform(local).bounding_box()[4]) / 2
    jaws = {'collar:jaw_upper', 'collar:jaw_lower'}
    fixed = set(proof['fixed_keys']) - jaws
    fixed.add('collar:liner_upper')
    rows = []
    manifold = A.I._manifold().Manifold
    for x in range(-6 if bridge else -2, 7, 2):
        for y in ([-4, -2, 0, 2, 4] if bridge else [-2, 0, 2]):
            for z in [-4, -2, 0, 2, 4]:
                point = basis @ (centre + [x, y, z])
                bore = rod(point - pin_axis*5, point + pin_axis*5, 1.6)
                knuckles = []
                for name, start, stop in [('upper', -1.4, 1.4), ('lower', -4.5, -1.6),
                                          ('lower', 1.6, 4.5)]:
                    key = f'collar:jaw_{name}'
                    sleeve = rod(point + pin_axis*start, point + pin_axis*stop, 3)
                    if bridge:
                        root = basis @ (centre + [0, 2 if name == 'upper' else -1, 0])
                        neck = rod(root + pin_axis*start, root + pin_axis*stop, 2.4)
                        sleeve = manifold.batch_hull([neck, sleeve])
                    sleeve -= bore
                    overlap = float(A.I.intersection_volume(sleeve, solids[key]))
                    knuckles.append((key, sleeve, overlap))
                row = dict(offset_local_mm=[x, y, z], pin_centre_source_mm=point.tolist(),
                           root_overlaps_mm3=[item[2] for item in knuckles])
                rows.append(row)
                if min(row['root_overlaps_mm3']) < 1:
                    row['rejected'] = 'knuckle attachment overlap below 1 mm3 screening threshold'
                    continue
                pin = rod(point-pin_axis*4.5, point+pin_axis*4.5, 1.5)
                # The pin bore would be cut in both jaws. Opposite-jaw knuckle
                # intrusion outside that bore remains a collision.
                stationary_hits = []
                for part in [pin] + [sleeve for _, sleeve, _ in knuckles]:
                    for target in sorted(fixed):
                        volume = float(A.I.intersection_volume(part, solids[target]))
                        if volume > .001:
                            stationary_hits.append(dict(target=target, volume_mm3=volume))
                row['stationary_hits'] = stationary_hits
                if stationary_hits:
                    continue
                revised = {}
                for key in sorted(jaws):
                    body = solids[key] - bore
                    for owner, sleeve, _ in knuckles:
                        if owner != key:
                            body -= sleeve
                    for owner, sleeve, _ in knuckles:
                        if owner == key:
                            body += sleeve
                    revised[key] = body
                row['revised_jaw_components'] = {
                    key: sum(piece.volume() > 1e-8 for piece in body.decompose())
                    for key, body in revised.items()
                }
                if any(count != 1 for count in row['revised_jaw_components'].values()):
                    continue
                additions = [body - solids[key] for key, body in revised.items()]
                upper, lower = (revised[key] for key in ('collar:jaw_upper', 'collar:jaw_lower'))
                row['revised_jaw_pair_overlap_mm3'] = float(A.I.intersection_volume(upper, lower))
                if row['revised_jaw_pair_overlap_mm3'] > .001:
                    continue
                moving_hit = None
                for offset in np.unique(np.r_[np.arange(0, 80.001, 1), np.arange(0, 3.001, .1)]):
                    for key in proof['moving_keys']:
                        posed = solids[key].translate(direction*offset)
                        volume = max(float(A.I.intersection_volume(part, posed)) for part in additions+[pin])
                        if volume > .001:
                            moving_hit = dict(target=key, offset_mm=float(offset), volume_mm3=volume)
                            break
                    if moving_hit:
                        break
                row['bundle_first_hit'] = moving_hit
    passing = [row for row in rows if row.get('stationary_hits') == [] and row.get('bundle_first_hit', True) is None]
    report = 'collar_offset_bridge_screen.json' if bridge else 'collar_side_pin_screen.json'
    save(folder / report, dict(
        assembly_sha256=sha(source), proof_sha256=sha(proof_file), study_sha256=sha(study_file),
        script_sha256=sha(Path(__file__)), frame_sha256=sha(frame),
        pin_diameter_mm=3, bore_diameter_mm=3.2, sleeve_diameter_mm=6,
        offset_bridge=bridge, bridge_root_radius_mm=2.4 if bridge else None, rows=rows,
        passing_envelope_sites=passing, geometry_adopted=False,
        scope='Preinstalled envelope screen. Pin approach, retention, load and full jaw installation pending.'))
    print('sites', len(rows), 'attached', sum('stationary_hits' in r for r in rows),
          'static clear', sum(r.get('stationary_hits') == [] for r in rows), 'passing', passing, flush=True)


if __name__ == '__main__':
    main()
