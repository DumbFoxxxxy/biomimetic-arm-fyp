"""Export open continuous-wall channel cores for one screened relocated pair."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, read, save, sha


def main():
    base = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = base/'assembly_tpu_readback_candidate.json'
    fillets = base/'elbow_route_fillet_screen.json'
    repairs = base/'elbow_tpu_readback_repair.json'
    screen, repair = read(fillets), read(repairs)
    assert screen['assembly_sha256'] == sha(manifest)
    choice = max(screen['passing'], key=lambda row: row['radius_mm'])
    assert choice['radius_mm'] == 20.
    ext_start, ext_end = np.array(choice['extensor_endpoints_source_mm'])
    ext_axis = A.D.unit(ext_end-ext_start)
    paths = {'flexor': np.array(choice['trimmed_path_source_mm']),
             'extensor': np.stack([ext_start+ext_axis*8, ext_end-ext_axis*8])}
    folder = base/'elbow_relocated_channel_cores'
    folder.mkdir(exist_ok=True)
    rows = []
    actual = {}
    for label, path in paths.items():
        key = f'tpu:elbow_{label}'
        source = next(row['original_source'] for row in repair['rows'] if row['key'] == key)
        design_file = (ROOT/source['file']).parent/'design.json'
        design = read(design_file)
        inner_radius, outer_radius = design['tube_inner_radius_mm'], design['tube_outer_radius_mm']
        assert 0 < inner_radius < outer_radius < 2.
        # Open both ends by extending the lumen beyond the round-ended stock.
        lumen_points = np.vstack([path[0]-A.D.unit(path[1]-path[0])*3, path,
                                   path[-1]+A.D.unit(path[-1]-path[-2])*3])
        outer = A.D.outer_envelope(path, outer_radius)
        lumen = A.D.outer_envelope(lumen_points, inner_radius)
        wire = A.D.outer_envelope(lumen_points, .3)
        wall = (outer-lumen).as_original().simplify(1e-6)
        actual[label], metadata = export(wall, folder/f'{label}_channel_core.obj', path[0])
        rows.append(dict(label=label, source=metadata, path_source_mm=path.tolist(),
                         original_design_sha256=sha(design_file), outer_radius_mm=outer_radius,
                         inner_radius_mm=inner_radius, nominal_wall_mm=outer_radius-inner_radius,
                         volume_mm3=float(actual[label].volume()),
                         wire_overlap_mm3=float(A.I.intersection_volume(actual[label], wire)),
                         outside_screen_envelope_mm3=max(0., float((actual[label]-A.D.outer_envelope(
                             path, 2.+choice['sampling_sagitta_mm'] if label == 'flexor' else 2.)).volume()))))
    pair_overlap = float(A.I.intersection_volume(actual['flexor'], actual['extensor']))
    passed = all(row['wire_overlap_mm3'] <= .001 and row['outside_screen_envelope_mm3'] <= .001
                 for row in rows) and pair_overlap <= .001
    save(folder/'channel_core_check.json', dict(
        assembly_sha256=sha(manifest), fillet_sha256=sha(fillets), script_sha256=sha(Path(__file__)),
        choice=choice, channels=rows, pair_overlap_mm3=pair_overlap, geometry_screen_pass=passed,
        full_lattice_designed=False, terminal_fittings_verified=False, stage_pass=False,
        scope='Nominal mid-span channel cores with open ends. Terminal 8 mm on each side excluded.'))
    print([(row['label'], row['volume_mm3'], row['wire_overlap_mm3'], row['outside_screen_envelope_mm3'])
           for row in rows], 'pair', pair_overlap, 'pass', passed, flush=True)


if __name__ == '__main__':
    main()
