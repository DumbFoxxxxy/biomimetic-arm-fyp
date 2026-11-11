"""Trace the original elbow wires with a prescribed distal swivel motion."""
import argparse

import numpy as np
from scipy.optimize import brentq

from build_s1_teres_origin_mount import A, ROOT, load_component, mesh_geom, read, save, sha
from study_s1_anatomical_elbow_tpu_demand import bypass_displacement


def screen_swivel_tension(points, fraction, pivot, hinge, bend, swing, offset):
    """Evaluate the distal cap wrench from the incoming span at its bore."""
    pivot_world = bend[:, :3]@pivot+bend[:, 3]
    hinge_world = bend[:, :3]@hinge

    def torque_per_newton(degrees):
        turn = A.I.rotation_about(pivot, hinge, degrees)
        distal = points[-3:]@turn[:, :3].T+turn[:, 3]
        distal = distal@bend[:, :3].T+bend[:, 3]
        end_path = points[-3:]+fraction[-3:, None]*(distal-points[-3:])
        end_path[:, 1] += offset[-3:]
        incoming = A.D.unit(end_path[0]-end_path[1])
        return float(np.cross(end_path[1]-pivot_world, incoming)@hinge_world)

    candidates = np.arange(swing-30., swing+30.1, 2.)
    torques = [torque_per_newton(angle) for angle in candidates]
    roots = [brentq(torque_per_newton, lo, hi, xtol=1e-7)
             for lo, hi, left, right in zip(candidates[:-1], candidates[1:], torques[:-1], torques[1:])
             if left >= 0. and right <= 0. and left > right]
    root = min(roots, key=lambda angle: abs(angle-swing)) if roots else None
    return dict(prescribed_swivel_deg=swing, prescribed_torque_per_N_Nmm=torque_per_newton(swing),
                cable_only_equilibrium_deg=root,
                deviation_from_prescribed_deg=None if root is None else root-swing)


def cam_pitch_profile(angles_deg, target_deg):
    """Compute the working pitch locus for a fixed tangent and circular output drum."""
    from numpy.polynomial import Polynomial

    output_radius = 6.
    angles_rad = np.radians(angles_deg)
    fit = Polynomial.fit(angles_rad, -np.radians(target_deg), 5).convert()
    fit -= fit(0.)
    sample_deg = np.linspace(0., angles_deg[-1], 1201)
    sample_rad = np.radians(sample_deg)
    support = output_radius*fit.deriv()(sample_rad)
    curvature_radius = support+output_radius*fit.deriv(3)(sample_rad)
    # The cable tangent is fixed in the follower frame. In the cam frame
    # its outward normal rotates by minus the input shaft angle.
    normals = np.c_[np.cos(sample_rad), -np.sin(sample_rad)]
    tangents = np.c_[np.sin(sample_rad), np.cos(sample_rad)]
    pitch_points = support[:, None]*normals-output_radius*fit.deriv(2)(sample_rad)[:, None]*tangents
    assert min(support) > 0 and min(curvature_radius) > 0
    swivel_deg = -np.degrees(fit(angles_rad))
    return swivel_deg, dict(
        output_drum_pitch_radius_mm=output_radius, fit_coefficients_rad=fit.coef.tolist(),
        max_fit_error_deg=float(np.max(np.abs(swivel_deg-target_deg))),
        cable_take_up_mm=float(output_radius*fit(sample_rad[-1])),
        effective_radius_range_mm=[float(support.min()), float(support.max())],
        pitch_curve_min_curvature_radius_mm=float(curvature_radius.min()),
        pitch_curve_max_radial_extent_mm=float(np.linalg.norm(pitch_points, axis=1).max()),
        pitch_curve_step_deg=.1, working_pitch_points_xy_mm=pitch_points.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--combined-bypass', action='store_true')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--extensor-equilibrium', action='store_true')
    mode.add_argument('--fixed-ratio', action='store_true')
    mode.add_argument('--variable-ratio', action='store_true')
    mode.add_argument('--cam-pitch', action='store_true')
    parser.add_argument('--cam-phase-deg', type=float, default=0.)
    args = parser.parse_args()
    if args.cam_phase_deg and not args.cam_pitch:
        parser.error('--cam-phase-deg requires --cam-pitch')
    coupled = args.fixed_ratio or args.variable_ratio or args.cam_pitch
    if (args.extensor_equilibrium or coupled) and not args.combined_bypass:
        parser.error('the selected mode requires --combined-bypass')
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    assembly_file = folder/'assembly_anatomical_elbow_angled_candidate.json'
    motion_file = folder/'anatomical_elbow_rigid_motion.json'
    baseline_file = folder/'anatomical_elbow_baseline.json'
    motion = read(motion_file)
    assert motion['assembly_sha256'] == sha(assembly_file)
    sources = read(assembly_file)['sources']
    bypass = None
    if args.combined_bypass:
        combined_file = folder/'anatomical_elbow_combined_bypass.json'
        bypass_file = folder/'anatomical_elbow_extensor_bypass_screen.json'
        combined = read(combined_file)
        assert combined['bypass_screen_sha256'] == sha(bypass_file)
        bypass = read(bypass_file)['poses']
        sources['tpu:deltoid'] = combined['deltoid_source']
    equilibrium = None
    ratios = None
    cam_profiles = None
    if coupled:
        from scipy.interpolate import PchipInterpolator
        from scipy.optimize import linprog

        prescribed_file = folder/('anatomical_elbow_variable_ratio_geometry.json' if args.cam_pitch
                                  else 'anatomical_elbow_combined_bypass_wire.json')
        prescribed = read(prescribed_file)
        assert prescribed['combined_bypass_sha256'] == sha(combined_file)
        angles = np.array([row['elbow_deg'] for row in prescribed['poses']])
        ratios = {}
        cam_profiles = {}
        for label in ('flexor', 'extensor'):
            target = np.array([row['muscles'][label]['swivel_deg'] for row in prescribed['poses']])
            if args.cam_pitch:
                cam_profiles[label], ratios[label] = cam_pitch_profile(angles, target)
                if args.cam_phase_deg:
                    from numpy.polynomial import Polynomial

                    fit = Polynomial(ratios[label]['fit_coefficients_rad'])
                    cam_profiles[label] = -np.degrees(fit(np.radians(angles+args.cam_phase_deg)))
                continue
            if args.variable_ratio:
                knots = np.arange(0, 121, 30)
                profile = PchipInterpolator(knots, target[knots])
                cam_profiles[label] = profile
                ratios[label] = dict(elbow_knots_deg=knots.tolist(), swivel_knots_deg=target[knots].tolist(),
                                     max_fit_error_deg=float(np.max(np.abs(profile(angles)-target))),
                                     local_ratio_range=[float(profile.derivative()(angles).min()),
                                                        float(profile.derivative()(angles).max())])
                continue
            # Minimize the maximum angular error subject to a fixed zero offset.
            coefficients = np.c_[np.r_[angles, -angles], -np.ones(2*len(angles))]
            fit = linprog([0., 1.], A_ub=coefficients, b_ub=np.r_[target, -target],
                          bounds=[(None, None), (0., None)])
            assert fit.success
            ratios[label] = dict(swivel_per_elbow_degree=float(fit.x[0]),
                                 max_fit_error_deg=float(fit.x[1]))
    if args.extensor_equilibrium:
        tension_file = folder/'anatomical_elbow_swivel_tension_screen.json'
        tension = read(tension_file)
        assert tension['combined_wire_report_sha256'] == sha(folder/'anatomical_elbow_combined_bypass_wire.json')
        equilibrium = {row['elbow_deg']: row['muscles']['extensor']['cable_only_equilibrium_deg']
                       for row in tension['poses']}
        assert len(equilibrium) == 121 and all(value is not None for value in equilibrium.values())
    moving = set(motion['moving_keys'])
    solids = {key: load_component(source) for key, source in sources.items()
              if key not in motion['excluded_unposed_elbow_tpu']}
    original = read(baseline_file)['muscles']
    routes = {}
    for label, muscle in original.items():
        head = muscle['anatomical_id']
        source = sources[f'elbow:{head}_cap_rotor']
        screen_file = (ROOT/source['file']).parent/'screen.json'
        screen = read(screen_file)
        definition_file = ROOT/muscle['definition']
        assert sha(definition_file) == muscle['definition_sha256']
        definition = read(definition_file)
        longitudinal = np.array(definition['basis_columns_display'])[:, 2]
        pads = np.array(muscle['pad_centres_source_mm'])
        exit_point = pads[0]-longitudinal*2.15
        distal_entry = pads[1]+longitudinal*2.15
        guides = np.array(muscle['wire_path_source_mm'])
        guide_stations = (guides-exit_point)@longitudinal
        span = (distal_entry-exit_point)@longitudinal
        guides = guides[(guide_stations > 1e-6) & (guide_stations < span-1e-6)]
        points = np.vstack([pads[0]-longitudinal*5.15, exit_point, guides,
                            distal_entry, pads[1]+longitudinal*4.6])
        if bypass and label == 'extensor':
            mf = A.I._manifold()
            mesh = load_component(muscle['collision_proxy']).to_mesh64()
            ext_vertices, ext_faces = mesh_geom.manifold_mesh_to_arrays(mesh)
            ext_faces = ext_faces.astype(np.uint32)
            pad_origin = pads[0]+longitudinal*2.
            pad_span = (pads[1]-pads[0])@longitudinal-4.
            ext_stations = (ext_vertices-pad_origin)@longitudinal/pad_span
            wire_stations = (points-pad_origin)@longitudinal/pad_span
        fraction = np.clip((points-exit_point)@longitudinal/span, 0., 1.)
        pivot = np.array(screen['pivot_zero_display_mm'])+np.array(source['translation_source_mm'])
        hinge = np.array(screen['trunnion_axis_zero_display'])
        pack = {key for key in moving if key.startswith(f'elbow:{head}_distal_')
                or key in {f'elbow:{head}_cap_rotor', f'elbow:{head}_active_wire_crimp_barrel'}}
        routes[label] = (points, fraction, pivot, hinge, pack, sha(screen_file))
    assert not routes['flexor'][4] & routes['extensor'][4]
    swivel_keys = set().union(*(route[4] for route in routes.values()))
    fixed = {key: solid for key, solid in solids.items() if key not in moving}
    reference_lengths = {}
    rows = []
    tension_rows = []
    for angle in range(121):
        bend = A.I.rotation_about(motion['centre_source_mm'], motion['axis_source'], angle)
        posed = {key: solids[key].transform(bend) for key in moving-swivel_keys}
        wires, samples = {}, {}
        tension_sample = {}
        for label, (points, fraction, pivot, hinge, pack, _) in routes.items():
            end = bend[:, :3]@points[-1]+bend[:, 3]
            preceding = points[-3]+fraction[-3]*(end-points[-1])
            reference = points[-3]-points[-1]
            direction = bend[:, :3].T@(preceding-end)
            reference -= hinge*(reference@hinge)
            direction -= hinge*(direction@hinge)
            swing = float(np.degrees(np.arctan2(hinge@np.cross(reference, direction), reference@direction)))
            if equilibrium and label == 'extensor':
                swing = equilibrium[angle]
            if args.fixed_ratio:
                swing = angle*ratios[label]['swivel_per_elbow_degree']
            if args.variable_ratio:
                swing = float(cam_profiles[label](angle))
            if args.cam_pitch:
                swing = float(cam_profiles[label][angle])
            turn = A.I.rotation_about(pivot, hinge, swing)
            distal = points@turn[:, :3].T+turn[:, 3]
            distal = distal@bend[:, :3].T+bend[:, 3]
            path = points+fraction[:, None]*(distal-points)
            offset = np.zeros(len(points))
            if bypass and label == 'extensor':
                offset, amplitude = bypass_displacement(angle, wire_stations, bypass)
                assert amplitude == combined['poses'][angle]['imposed_bypass_mm']
                path[:, 1] += offset
            if bypass:
                tension_sample[label] = screen_swivel_tension(points, fraction, pivot, hinge,
                                                             bend, swing, offset)
            segments = np.diff(path, axis=0)
            lengths = np.linalg.norm(segments, axis=1)
            length = float(lengths.sum())
            reference_lengths.setdefault(label, length)
            tangent = segments/lengths[:, None]
            bends = np.degrees(np.arccos(np.clip(np.sum(tangent[:-1]*tangent[1:], axis=1), -1., 1.)))
            wires[f'wire:{label}'] = A.D.outer_envelope(path, .3)
            samples[label] = dict(swivel_deg=swing, length_mm=length,
                                  length_change_mm=length-reference_lengths[label],
                                  max_polyline_bend_deg=float(bends.max()))
            if bypass and label == 'extensor' and (equilibrium or ratios or angle == 0 or angle >= 105):
                distal_vertices = ext_vertices@turn[:, :3].T+turn[:, 3]
                distal_vertices = distal_vertices@bend[:, :3].T+bend[:, 3]
                vertices = ext_vertices+np.clip(ext_stations, 0., 1.)[:, None]*(distal_vertices-ext_vertices)
                vertices[:, 1] += bypass_displacement(angle, ext_stations, bypass)[0]
                tpu = mf.Manifold(mf.Mesh64(np.ascontiguousarray(vertices), ext_faces))
                assert tpu.status() == mf.Error.NoError and not tpu.is_empty()
                samples[label]['own_tpu_overlap_mm3'] = A.I.intersection_volume(wires[f'wire:{label}'], tpu)
                if equilibrium or ratios:
                    samples[label]['tpu_external_overlap_mm3'] = {
                        key: A.I.intersection_volume(tpu, fixed[key])
                        for key in ('bone:humerus', 'tpu:deltoid')}
                if cam_profiles:
                    samples[label]['tpu_external_clearance_mm_capped1'] = {
                        key: A.I.gap(tpu, fixed[key], 1.) for key in ('bone:humerus', 'tpu:deltoid')}
            posed.update({key: solids[key].transform(turn).transform(bend) for key in pack})
        obstacles = {**fixed, **posed}
        wire_hits = {label: {key: volume for key, part in obstacles.items()
                            if (volume := A.I.intersection_volume(wire, part)) > .001}
                     for label, wire in wires.items()}
        obstacles.update(wires)
        rotary_hits = {}
        for label, (_, _, _, _, pack, _) in routes.items():
            rotary_hits[label] = [
                dict(part=key, obstacle=target, overlap_mm3=volume)
                for key in sorted(pack) for target, part in obstacles.items()
                if target not in pack and not target.startswith('wire:')
                and (volume := A.I.intersection_volume(posed[key], part)) > .001]
        rows.append(dict(elbow_deg=angle, muscles=samples, wire_contacts=wire_hits,
                         rotary_contacts=rotary_hits,
                         wire_pair_overlap_mm3=A.I.intersection_volume(*wires.values())))
        if bypass:
            tension_rows.append(dict(elbow_deg=angle, muscles=tension_sample))
        if ratios:
            rows[-1]['torque_per_N_Nmm'] = {
                label: sample['prescribed_torque_per_N_Nmm'] for label, sample in tension_sample.items()}
    output = 'anatomical_elbow_constrained_wire_motion.json'
    if bypass:
        output = 'anatomical_elbow_combined_bypass_wire.json'
    if equilibrium:
        output = 'anatomical_elbow_extensor_equilibrium_geometry.json'
    if args.fixed_ratio:
        output = 'anatomical_elbow_fixed_ratio_geometry.json'
    if cam_profiles:
        output = 'anatomical_elbow_variable_ratio_geometry.json'
    if args.cam_pitch:
        output = 'anatomical_elbow_cam_pitch_geometry.json'
        if args.cam_phase_deg:
            output = f'anatomical_elbow_cam_phase_{args.cam_phase_deg:+g}.json'
    save(folder/output, dict(
        rigid_motion_sha256=sha(motion_file), anatomical_baseline_sha256=sha(baseline_file),
        swivel_screen_sha256={label: route[-1] for label, route in routes.items()},
        sample_step_deg=1, wire_radius_mm=.3, overlap_threshold_mm3=.001, poses=rows,
        **({'combined_bypass_sha256': sha(combined_file),
            'own_extensor_tpu_screen_angles_deg': (list(range(121)) if equilibrium or ratios else [0, *range(105, 121)])}
           if bypass else {}),
        **({'tension_screen_sha256': sha(tension_file), 'extensor_uses_cable_only_equilibrium': True}
           if equilibrium else {}),
        **({'prescribed_reference_sha256': sha(prescribed_file),
            ('cam_pitch_profiles' if args.cam_pitch else
             'variable_ratio_profiles' if cam_profiles else 'fixed_ratios'): ratios} if ratios else {}),
        **({'imposed_input_phase_error_deg': args.cam_phase_deg,
            'polynomial_extrapolation_beyond_nominal_travel': True} if args.cam_phase_deg else {}),
        scope='Prescribed neutral-shoulder flexion. Rigid cap bore endpoints are fixed to their respective caps; intermediate guide centres are blended. Distal swivels '
              'track the neutral-rotor terminal direction. No force equilibrium, TPU deformation, or friction. '
              'Contact flags include designed seats and do not establish failure.'))
    if bypass and not equilibrium and not ratios:
        save(folder/'anatomical_elbow_swivel_tension_screen.json', dict(
            combined_wire_report_sha256=sha(folder/output), sample_step_deg=1, search_half_span_deg=30,
            poses=tension_rows,
            scope='Ideal unit tension in the last free span at the distal cap entrance. Stable roots '
                  'use the prescribed guide interpolation and exclude TPU reaction, friction and stops. '
                  'Root poses have not been screened for contact.'))
    for label in routes:
        values = [row['muscles'][label] for row in rows]
        print(label, 'length change range', min(v['length_change_mm'] for v in values),
              max(v['length_change_mm'] for v in values),
              'wire contact poses', sum(bool(row['wire_contacts'][f'wire:{label}']) for row in rows),
              'rotary contact poses', sum(bool(row['rotary_contacts'][label]) for row in rows))


if __name__ == '__main__':
    main()
