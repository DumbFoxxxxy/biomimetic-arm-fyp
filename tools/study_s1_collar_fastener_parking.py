"""Check whether group 1 can be parked during the deltoid entry sweep."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, load_component, read, save, sha


def main():
    folder = A.DESIGN / 's1_deltoid_support_draft'
    source = folder / 'assembly_check.json'
    proof_file = folder / 'verification_checks.json'
    study_file = folder / 'collar_first_screen.json'
    proof, study = read(proof_file), read(study_file)
    assert proof['assembly_sha256'] == study['assembly_sha256'] == sha(source)
    sources = read(source)['sources']
    frame_file = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    axis = np.array(read(frame_file)['frame_source']['basis_columns'])[:, 1]
    direction = np.array(proof['direction_source'])
    offsets = np.unique(np.r_[np.arange(0, 80.001, .5), np.arange(0, 3.001, .1)])
    cube = A.I._manifold().Manifold.cube
    parts = {key: load_component(sources[key]) for key in proof['moving_keys']}
    rows = []
    for part, sign in [('bolt', 1), ('top_washer', 1), ('nut', -1), ('bottom_washer', -1)]:
        key = f'collar:{part}_1'
        solid = load_component(sources[key])
        low, high = np.array(solid.bounding_box()).reshape(2, 3)
        mask_low = low + np.minimum(0, axis*sign*20) + np.minimum(0, -direction*80) - .1
        mask_high = high + np.maximum(0, axis*sign*20) + np.maximum(0, -direction*80) + .1
        mask = cube(mask_high-mask_low).translate(mask_low)
        cropped = {name: body ^ mask for name, body in parts.items()}
        trials = []
        for parking in np.arange(0, 20.001, .5):
            parked = solid.translate(axis*sign*parking)
            hits = []
            for name, body in cropped.items():
                for offset in offsets:
                    volume = float(A.I.intersection_volume(parked, body.translate(direction*offset)))
                    if volume > .001:
                        hits.append(dict(part=name, offset_mm=float(offset), volume_mm3=volume))
            trials.append(dict(parking_mm=float(parking), hits=len(hits),
                               peak=max(hits, key=lambda hit: hit['volume_mm3']) if hits else None))
        clear = [trial['parking_mm'] for trial in trials if not trial['hits']]
        rows.append(dict(key=key, direction=(axis*sign).tolist(), clear_parking_mm=clear, trials=trials))
        print(key, 'at final', trials[0], 'clear parking', clear, flush=True)
    save(folder / 'collar_fastener_parking_screen.json', dict(
        assembly_sha256=sha(source), proof_sha256=sha(proof_file), script_sha256=sha(Path(__file__)),
        frame_sha256=sha(frame_file), offsets_mm=offsets.tolist(), rows=rows,
        scope='Moving bundle versus each static parked part. Engagement, fixed obstacles and final tightening are pending.'))


if __name__ == '__main__':
    main()
