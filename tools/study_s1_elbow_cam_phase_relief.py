"""Screen shoulder-side relief against nominal and offset cam phase paths."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, load_component, mesh_geom, read, save, sha
from study_s1_anatomical_elbow_tpu_demand import bypass_displacement, screen_deltoid_relief


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    combined_file = folder/'anatomical_elbow_combined_bypass.json'
    bypass_file = folder/'anatomical_elbow_extensor_bypass_screen.json'
    assembly_file = folder/'assembly_anatomical_elbow_angled_candidate.json'
    motion_file = folder/'anatomical_elbow_rigid_motion.json'
    baseline_file = folder/'anatomical_elbow_baseline.json'
    combined, bypass, assembly, motion, baseline = map(
        read, (combined_file, bypass_file, assembly_file, motion_file, baseline_file))
    assert combined['bypass_screen_sha256'] == sha(bypass_file)
    assert motion['assembly_sha256'] == sha(assembly_file)
    muscle = baseline['muscles']['extensor']
    rotor = assembly['sources']['elbow:FJ1477_cap_rotor']
    screen_file = (ROOT/rotor['file']).parent/'screen.json'
    screen = read(screen_file)
    pivot = np.asarray(screen['pivot_zero_display_mm'])+rotor['translation_source_mm']
    hinge = screen['trunnion_axis_zero_display']
    definition_file = ROOT/muscle['definition']
    assert sha(definition_file) == muscle['definition_sha256']
    longitudinal = np.asarray(read(definition_file)['basis_columns_display'])[:, 2]
    proximal, distal = np.asarray(muscle['pad_centres_source_mm'])
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(
        load_component(muscle['collision_proxy']).to_mesh64())
    faces = faces.astype(np.uint32)
    stations = (vertices-proximal-longitudinal*2.)@longitudinal/((distal-proximal)@longitudinal-4.)
    weights = np.clip(stations, 0., 1.)[:, None]
    mf = A.I._manifold()
    test_files = {
        -0.25: folder/'anatomical_elbow_cam_phase_-0.25.json',
        0.: folder/'anatomical_elbow_cam_pitch_geometry.json',
        .25: folder/'anatomical_elbow_cam_phase_+0.25.json'}
    cases = {}
    for phase, path in test_files.items():
        report = read(path)
        assert report['combined_bypass_sha256'] == sha(combined_file)
        assert report['anatomical_baseline_sha256'] == sha(baseline_file)
        assert report['swivel_screen_sha256']['extensor'] == sha(screen_file)
        cases[phase] = []
        for row in report['poses']:
            angle = row['elbow_deg']
            turn = A.I.rotation_about(pivot, hinge, row['muscles']['extensor']['swivel_deg'])
            bend = A.I.rotation_about(motion['centre_source_mm'], motion['axis_source'], angle)
            far = (vertices@turn[:, :3].T+turn[:, 3])@bend[:, :3].T+bend[:, 3]
            posed = vertices+weights*(far-vertices)
            posed[:, 1] += bypass_displacement(angle, stations, bypass['poses'])[0]
            solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), faces))
            assert solid.status() == mf.Error.NoError and not solid.is_empty()
            cases[phase].append(solid)
    deltoid = load_component(combined['deltoid_source'])
    caps = {key: load_component(source) for key, source in assembly['sources'].items()
            if key.startswith(('cap:FJ1468/', 'cap:FJ1513/', 'cap:FJ1467/'))}
    filename = 'deltoid_anatomical_elbow_cam_phase_relief.obj'
    result = screen_deltoid_relief(folder, deltoid, [solid for solids in cases.values() for solid in solids],
                                  caps, filename=filename,
                                  capture_threshold_mm3=1e-7, near_capture_mm=.15)
    result['phase_report_sha256'] = {str(phase): sha(path) for phase, path in test_files.items()}
    selected = next((row for row in result['trials'] if 'source' in row), None)
    if selected is None:
        (folder/filename).unlink(missing_ok=True)
    if selected:
        actual = load_component(selected['source'])
        result['phase_readback'] = {}
        for phase, poses in cases.items():
            gaps = [A.I.gap(actual, solid, 1.) for solid in poses]
            result['phase_readback'][str(phase)] = dict(
                minimum_gap_mm=min(gaps), minimum_gap_elbow_deg=int(np.argmin(gaps)),
                peak_overlap_mm3=max(A.I.intersection_volume(actual, solid) for solid in poses))
    save(folder/'anatomical_elbow_cam_phase_near_relief.json', result)
    print(result.get('phase_readback'), result['trials'])


if __name__ == '__main__':
    main()
