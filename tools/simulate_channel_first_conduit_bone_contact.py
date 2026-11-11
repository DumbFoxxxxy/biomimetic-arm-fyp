"""Guide beam with unilateral distributed contact on actual fixed bone surfaces."""
import argparse
import xml.etree.ElementTree as ET

import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save, load_component
from build_channel_first_inventory import provenance
from channel_first_bone_contact import BoneContactSurface
from simulate_channel_first_elastic_conduits import build_model, tension_profile


def contact_stations(lengths, spacing=.0025):
    stations = {}
    for i, length in enumerate(lengths):
        count = int(np.ceil(length/spacing))
        for j in range(count+1):
            fraction = j/count
            key = (i, i, 0.) if j == 0 else ((i+1, i+1, 0.) if j == count else (i, i+1, fraction))
            stations[key] = stations.get(key, 0.)+length/count*(.5 if j in (0, count) else 1.)
    return [(a, b, fraction, weight) for (a, b, fraction), weight in stations.items()]


def run(route, surfaces, peak, contact_density, friction, update_steps):
    import mujoco
    root, spec = build_model(route)
    model = mujoco.MjModel.from_xml_string(ET.tostring(root, encoding='unicode'))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    node_ids = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, f'node_{i}')
                for i in range(spec['segments']+1)]
    rest = np.array(spec['segment_rest_lengths_m'])
    neutral = data.site_xpos[node_ids].copy()
    axial_ids = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, f'axial_{i}') for i in range(spec['segments'])]
    axial_q = model.jnt_qposadr[axial_ids]
    stations = contact_stations(rest)
    aa = np.array([row[0] for row in stations])
    bb = np.array([row[1] for row in stations])
    fractions = np.array([row[2] for row in stations])
    weights = np.array([row[3] for row in stations])
    source_origin = np.array(spec['source_origin_m'])
    initial_length = float(data.ten_length[0])
    held_force = np.zeros(model.nv)
    trace, deepest = [], None
    max_penetration, min_distance, max_contact_force = 0., np.inf, 0.
    peak_displacement, minimum_length_change = 0., 0.
    limit_hit, maximum_support_error = False, 0.
    last_contacts = []
    for step in range(round(.5/model.opt.timestep)):
        mujoco.mj_step1(model, data)
        data.ctrl[0] = tension_profile(data.time, peak)
        nodes = data.site_xpos[node_ids]
        if step % update_steps == 0:
            held_force[:] = 0.
            jacobians = np.empty((len(node_ids), 3, model.nv))
            for i, site in enumerate(node_ids):
                mujoco.mj_jacSite(model, data, jacobians[i], np.zeros((3, model.nv)), site)
            jac = (1.-fractions[:, None, None])*jacobians[aa]+fractions[:, None, None]*jacobians[bb]
            points = (1.-fractions[:, None])*nodes[aa]+fractions[:, None]*nodes[bb]
            velocity = np.einsum('nij,j->ni', jac, data.qvel)
            world_mm = (points+source_origin)*1000.
            last_contacts = []
            total_normal = 0.
            for name, surface in surfaces.items():
                distance_mm, normals = surface.query(world_mm)
                min_distance = min(min_distance, float(distance_mm.min()))
                for i in np.where(distance_mm < 1.6)[0]:
                    penetration = (1.6-distance_mm[i])*.001
                    normal = normals[i]
                    normal_speed = float(velocity[i]@normal)
                    # Stiffness and damping per REST length, integrated by station weights.
                    force_n = max(0., weights[i]*(contact_density*penetration-50.*normal_speed))
                    tangent_velocity = velocity[i]-normal_speed*normal
                    force_t = -friction*force_n*tangent_velocity/np.sqrt(tangent_velocity@tangent_velocity+.001**2)
                    assert float(force_t@tangent_velocity) <= 1e-12
                    force = force_n*normal+force_t
                    held_force += jac[i].T@force
                    total_normal += force_n
                    last_contacts.append(dict(bone=name, station=int(i), penetration_mm=penetration*1000.,
                        normal_force_N=force_n, force_on_guide_N=force.tolist(),
                        source_position_mm=world_mm[i].tolist()))
                    if penetration*1000. > max_penetration:
                        max_penetration = penetration*1000.
                        deepest = dict(last_contacts[-1], time_s=float(data.time))
            max_contact_force = max(max_contact_force, total_normal)
        data.qfrc_applied[:] = held_force
        mujoco.mj_step2(model, data)
        if not np.isfinite(data.qpos).all() or not np.isfinite(data.qvel).all():
            raise ValueError('Nonfinite dynamics')
        if step % max(1, round(.002/model.opt.timestep)) == 0:
            mujoco.mj_forward(model, data)
            nodes = data.site_xpos[node_ids].copy()
            displacement = float(np.linalg.norm(nodes-neutral, axis=1).max())*1000.
            peak_displacement = max(peak_displacement, displacement)
            support_error = float(np.linalg.norm(nodes[-1]-neutral[-1]))*1000.
            maximum_support_error = max(maximum_support_error, support_error)
            strain = data.qpos[axial_q]/rest
            limit_hit |= bool(np.max(abs(strain)) >= .49)
            length_change = (float(data.ten_length[0])-initial_length)*1000.
            minimum_length_change = min(minimum_length_change, length_change)
            trace.append(dict(time_s=float(data.time), tension_command_N=float(data.ctrl[0]),
                max_node_displacement_mm=displacement, length_change_mm=length_change,
                distal_support_error_mm=support_error, minimum_axial_strain=float(strain.min()),
                maximum_axial_strain=float(strain.max()), nodes_local_m=nodes.tolist(),
                contact_samples=last_contacts))
    warnings = {str(i): int(w.number) for i, w in enumerate(data.warning) if w.number}
    result = dict(peak_tension_N=peak, peak_node_displacement_mm=peak_displacement,
        minimum_length_change_mm=minimum_length_change,
        maximum_guide_envelope_penetration_mm=max_penetration,
        minimum_sampled_centerline_signed_distance_mm=min_distance,
        minimum_sampled_wire_surface_gap_mm=min_distance-.3,
        maximum_sum_contact_normal_force_N=max_contact_force,
        maximum_distal_support_error_mm=maximum_support_error,
        numerical_axial_limit_hit=limit_hit, warning_counts=warnings,
        finite_no_warning_or_axial_limit=not warnings and not limit_hit,
        deepest_contact=deepest, trace=trace)
    spec.update(contact_radius_mm=1.6, station_max_rest_spacing_mm=2.5,
        contact_stations=len(stations), contact_update_steps=update_steps,
        contact_update_interval_s=update_steps*model.opt.timestep,
        contact_penalty_density_N_per_m2=contact_density,
        contact_damping_density_Ns_per_m2=50., bone_TPU_friction_coefficient=friction,
        friction_regularisation_speed_m_per_s=.001,
        contact_law='Unilateral penalty at dynamically moving beam stations on actual signed nearest-triangle bone surfaces. Point Jacobians distribute contact forces to both beam nodes by virtual work. No new fixed path points.',
        collision_scope='Humerus and scapula only in this fixed-end bench; full anatomy and continuum surface need separate verification')
    return root, spec, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tension', type=float, default=5.)
    parser.add_argument('--contact-density', type=float, default=200000.)
    parser.add_argument('--friction', type=float, default=0.)
    parser.add_argument('--update-steps', type=int, default=5)
    args = parser.parse_args()
    if args.contact_density <= 0 or args.friction < 0 or args.update_steps < 1 or not 0 <= args.tension <= 60:
        parser.error('Invalid contact parameters')
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    route_file = folder/'minimal_shoulder_conduits.json'
    reference_file = folder/'thorax_arm_reference.json'
    reference = read(reference_file)
    route = next(r for r in read(route_file)['conduits'] if r['id'] == 'FJ1468')
    sources = {name: reference['sources'][name] for name in ('bone:humerus', 'bone:scapula')}
    surfaces = {name: BoneContactSurface(load_component(source)) for name, source in sources.items()}
    root, spec, result = run(route, surfaces, args.tension, args.contact_density, args.friction, args.update_steps)
    dest = folder/'elastic_conduit_bone_contact'
    dest.mkdir(exist_ok=True)
    stem = f'FJ1468_T{args.tension:g}_k{args.contact_density:g}_mu{args.friction:g}_u{args.update_steps}'
    xml = dest/f'{stem}.xml'
    ET.ElementTree(root).write(xml, encoding='utf-8', xml_declaration=True)
    save(dest/f'{stem}.json', dict(source=provenance(route_file), anatomy_reference=provenance(reference_file),
        contact_sources=sources, model=provenance(xml), specification=spec, result=result,
        generator=provenance(ROOT/'tools/simulate_channel_first_conduit_bone_contact.py'),
        surface_query=provenance(ROOT/'tools/channel_first_bone_contact.py'),
        beam_generator=provenance(ROOT/'tools/simulate_channel_first_elastic_conduits.py'),
        scope='Fixed-end guide contact sensitivity only; no moving bone assembly, measured material, full lattice, steel/TPU internal friction or independent cable slack. XML alone omits Python contact forces. Bone friction here is TPU/bone, not cable/conduit friction.'))
    print(stem, {key: value for key, value in result.items() if key not in ('trace', 'deepest_contact')}, flush=True)


if __name__ == '__main__':
    main()
