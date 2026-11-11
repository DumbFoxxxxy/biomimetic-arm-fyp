"""Axially compliant articulated guide bench: fixed end seats, free internal nodes.

This is a boundary-condition study, not an installed whole-arm drive model.
Each beam segment has axial elasticity; successive segments have rotational
elasticity about their inherited rest pose. Cable tension acts through sites
on these moving segments. No internal path point is fixed to the world.
"""
import argparse
import math
import xml.etree.ElementTree as ET

import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector


def sample_path(points, spacing):
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    stations = np.linspace(0., arc[-1], int(np.ceil(arc[-1]/spacing))+1)
    return np.column_stack([np.interp(stations, arc, points[:, i]) for i in range(3)])


def build_model(route, modulus_MPa=15., segment_length_mm=5., timestep=.00002):
    points = sample_path(np.array(route['path_source_mm']), segment_length_mm)*.001
    origin = points[0].copy()
    points -= origin
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    area = np.pi*(.0016**2-.0008**2)
    second_moment = np.pi/4*(.0016**4-.0008**4)
    modulus = modulus_MPa*1e6
    root = ET.Element('mujoco', model='elastic_conduit_fixed_end_seat_bench')
    ET.SubElement(root, 'compiler', angle='radian')
    ET.SubElement(root, 'option', timestep=str(timestep), integrator='implicitfast',
                  gravity='0 0 -9.81', jacobian='sparse', iterations='50', tolerance='1e-10')
    world = ET.SubElement(root, 'worldbody')
    ET.SubElement(world, 'light', pos='0 0 .5')
    anchor = ET.SubElement(world, 'body', name='distal_seat', pos=vector(points[-1]))
    ET.SubElement(anchor, 'site', name='distal_anchor', size='.001')
    parent = ET.SubElement(world, 'body', name='proximal_clamp')
    ET.SubElement(parent, 'site', name='node_0', size='.0005')
    bending, axial = [], []
    deltas = np.diff(points, axis=0)
    for i, delta in enumerate(deltas):
        length = lengths[i]
        direction = delta/length
        # A slide BEFORE the ball at each node extends the previous segment
        # along the parent-frame segment direction. One body per node avoids
        # the unnecessary two-body nesting used in the first implementation.
        end = ET.SubElement(parent, 'body', name=f'end_{i}', pos=vector(delta))
        stiffness = modulus*area/length
        ET.SubElement(end, 'joint', name=f'axial_{i}', type='slide', axis=vector(direction),
                      stiffness=str(stiffness), damping='.5', limited='true',
                      range=vector([-.5*length, .5*length]))
        axial.append(stiffness)
        if i < len(lengths)-1:
            lever = .5*(lengths[i]+lengths[i+1])
            stiffness = modulus*second_moment/lever
            # Rotational damping is a numerical/model parameter, not measured TPU.
            ET.SubElement(end, 'joint', name=f'bend_{i+1}', type='ball',
                          stiffness=str(stiffness), damping='.00003')
            bending.append(stiffness)
        mass_length = .5*length+(.5*lengths[i+1] if i < len(lengths)-1 else 0.)
        mass = 1200.*area*mass_length
        scalar_inertia = mass*(mass_length**2/12.+(.0016**2+.0008**2)/4.)
        ET.SubElement(end, 'inertial', pos='0 0 0', mass=str(mass),
                      diaginertia=vector([scalar_inertia]*3))
        ET.SubElement(end, 'site', name=f'node_{i+1}', size='.0005', rgba='.9 .1 .1 1')
        parent = end
    equality = ET.SubElement(root, 'equality')
    ET.SubElement(equality, 'connect', name='distal_translation_support',
                  site1=f'node_{len(lengths)}', site2='distal_anchor',
                  solref='.001 1', solimp='.999 .9999 .0001')
    tendons = ET.SubElement(root, 'tendon')
    path = ET.SubElement(tendons, 'spatial', name='wire', width='.0003')
    for i in range(len(points)):
        ET.SubElement(path, 'site', site=f'node_{i}')
    actuators = ET.SubElement(root, 'actuator')
    ET.SubElement(actuators, 'motor', name='pull', tendon='wire', gear='-1',
                  ctrllimited='true', ctrlrange='0 60')
    return root, dict(source_origin_m=origin.tolist(), reference_nodes_local_m=points.tolist(),
        segment_rest_lengths_m=lengths.tolist(), segments=len(lengths),
        assumed_modulus_MPa=modulus_MPa, assumed_density_kg_m3=1200.,
        annular_area_m2=area, annular_second_moment_m4=second_moment,
        bending_stiffness_Nm_per_rad=bending, axial_stiffness_N_per_m=axial,
        timestep_s=timestep, rotational_damping_Nms_per_rad=.00003,
        axial_damping_Ns_per_m=.5,
        inertia_model='Lumped node masses from adjacent half-segments, isotropic nodal rotary inertia; clamped root half-segment mass omitted from moving system',
        axial_joint_limits_fraction=[-.5, .5],
        limit_note='Numerical protection only, not a proposed physical stroke stop',
        boundary_conditions='Proximal position and first-segment orientation fixed; distal position soft-constrained, orientation free. No intermediate world-fixed nodes.')


def tension_profile(time, peak):
    if time < .05:
        return peak*(.5-.5*math.cos(np.pi*time/.05))
    if time < .25:
        return peak
    if time < .30:
        return peak*(.5+.5*math.cos(np.pi*(time-.25)/.05))
    return 0.


def simulate(root, specification, peak):
    import mujoco
    model = mujoco.MjModel.from_xml_string(ET.tostring(root, encoding='unicode'))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    assert np.max(np.abs(data.efc_pos)) < 1e-8
    node_ids = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, f'node_{i}')
                for i in range(specification['segments']+1)]
    axial_ids = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, f'axial_{i}')
                 for i in range(specification['segments'])]
    axial_addresses = model.jnt_qposadr[axial_ids]
    rest = np.array(specification['segment_rest_lengths_m'])
    neutral = data.site_xpos[node_ids].copy()
    assert np.allclose(neutral, specification['reference_nodes_local_m'], atol=1e-10)
    initial_length = float(data.ten_length[0])
    records, maximum_support_error, limit_hit = [], 0., False
    max_node_displacement, min_axial_strain = 0., 0.
    for step in range(int(round(.5/model.opt.timestep))):
        data.ctrl[0] = tension_profile(data.time, peak)
        mujoco.mj_step(model, data)
        if step % max(1, round(.002/model.opt.timestep)):
            continue
        mujoco.mj_forward(model, data)
        assert np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
        current = data.site_xpos[node_ids].copy()
        error = float(np.linalg.norm(current[-1]-neutral[-1]))
        maximum_support_error = max(maximum_support_error, error)
        strain = data.qpos[axial_addresses]/rest
        min_axial_strain = min(min_axial_strain, float(strain.min()))
        limit_hit |= bool(np.max(np.abs(strain)) >= .49)
        displacement = float(np.linalg.norm(current-neutral, axis=1).max())
        max_node_displacement = max(max_node_displacement, displacement)
        vectors = np.diff(current, axis=0)
        lengths = np.linalg.norm(vectors, axis=1)
        turns = np.arccos(np.clip(np.sum(vectors[:-1]*vectors[1:], axis=1)/(lengths[:-1]*lengths[1:]), -1., 1.))
        discrete_radius = .5*(lengths[:-1]+lengths[1:])/np.maximum(turns, 1e-12)
        point_jacobian = np.zeros((3, model.nv))
        mujoco.mj_jacSite(model, data, point_jacobian, np.zeros_like(point_jacobian), node_ids[-1])
        support_force, _, _, _ = np.linalg.lstsq(point_jacobian.T, data.qfrc_constraint, rcond=None)
        reconstruction = float(np.max(np.abs(point_jacobian.T@support_force-data.qfrc_constraint)))
        records.append(dict(time_s=float(data.time), tension_command_N=float(data.ctrl[0]),
            guide_length_mm=float(data.ten_length[0])*1000.,
            length_change_mm=(float(data.ten_length[0])-initial_length)*1000.,
            max_node_displacement_mm=displacement*1000., distal_support_error_mm=error*1000.,
            minimum_segment_axial_strain=float(strain.min()), maximum_segment_axial_strain=float(strain.max()),
            minimum_discrete_bend_radius_mm=float(discrete_radius.min())*1000.,
            distal_support_force_on_guide_N=support_force.tolist(),
            constraint_force_reconstruction_max_abs_mixed_units=reconstruction,
            nodes_local_m=current.tolist()))
    warnings = {str(i): int(w.number) for i, w in enumerate(data.warning) if w.number}
    return dict(peak_tension_N=peak, initial_length_mm=initial_length*1000.,
        maximum_node_displacement_mm=max_node_displacement*1000.,
        maximum_distal_support_error_mm=maximum_support_error*1000.,
        minimum_segment_axial_strain=min_axial_strain, numerical_axial_limit_hit=limit_hit,
        final_length_change_mm=(float(data.ten_length[0])-initial_length)*1000.,
        final_max_node_displacement_mm=records[-1]['max_node_displacement_mm'],
        minimum_loaded_length_change_mm=min(r['length_change_mm'] for r in records),
        maximum_distal_support_force_N=max(float(np.linalg.norm(r['distal_support_force_on_guide_N'])) for r in records),
        warning_counts=warnings, trace=records,
        finite_and_no_limit_or_warning=not warnings and not limit_hit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--route', default='FJ1468')
    parser.add_argument('--modulus', type=float, default=15.)
    parser.add_argument('--tension', type=float, default=5.)
    parser.add_argument('--segment-length', type=float, default=5.)
    parser.add_argument('--timestep', type=float, default=.00002)
    args = parser.parse_args()
    if not (args.modulus > 0. and 0. <= args.tension <= 60. and args.segment_length > 0. and args.timestep > 0.):
        parser.error('Modulus, segment length and timestep must be positive; tension must be 0..60 N')
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    source_file = folder/'minimal_shoulder_conduits.json'
    route = next(r for r in read(source_file)['conduits'] if r['id'] == args.route)
    root, specification = build_model(route, args.modulus, args.segment_length, args.timestep)
    destination = folder/'elastic_conduit_bench'
    destination.mkdir(exist_ok=True)
    stem = f'{args.route}_E{args.modulus:g}_T{args.tension:g}_ds{args.segment_length:g}_dt{args.timestep:g}'
    xml_file = destination/f'{stem}.xml'
    ET.ElementTree(root).write(xml_file, encoding='utf-8', xml_declaration=True)
    result = simulate(root, specification, args.tension)
    save(destination/f'{stem}.json', dict(source=provenance(source_file),
        generator=provenance(ROOT/'tools/simulate_channel_first_elastic_conduits.py'),
        model=provenance(xml_file), route_id=args.route, specification=specification, result=result,
        scope='Isolated supported-guide sensitivity bench, NOT whole-arm motion or final hardware approval. No friction, bone contact, cross-section ovalisation, shear, self-contact, measured TPU material or fitted end connectors. Active input is prescribed cable tension, not prescribed pose.',
        cable_model='Uniform nonnegative pulling tension on a tendon through dynamically moving beam nodes. Beam axial and rotational springs carry guide reaction. Intermediate nodes are never fixed to world. Slipping/friction and separate cable slack are not modelled in this bench.'))
    print(stem, {k: v for k, v in result.items() if k != 'trace'})


if __name__ == '__main__':
    main()
