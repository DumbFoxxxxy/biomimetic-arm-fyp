"""Attach tension-only actuators and passive elasticity to explicit endpoint spans."""
import argparse
import sys
import xml.etree.ElementTree as ET

import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector


def main():
    import mujoco

    sys.path.insert(0, str(ROOT))
    from project_runtime import load_mujoco_model
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--girdle-candidates', action='store_true')
    parser.add_argument('--upper-girdle', action='store_true')
    parser.add_argument('--geometry-mass', action='store_true')
    parser.add_argument('--arm-antagonists', action='store_true')
    parser.add_argument('--unified-sc-ac', action='store_true', help='Use current SC/AC candidate in a separate output folder')
    args = parser.parse_args()
    if args.arm_antagonists and not args.upper_girdle:
        parser.error('--arm-antagonists requires --upper-girdle')

    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    output_folder = folder/'unified_SC_AC' if args.unified_sc_ac else folder
    chain_stem = 'geometry_mass_body_chain' if args.geometry_mass else 'ideal_body_chain'
    chain_file = output_folder/f'{chain_stem}.xml'
    chain_check_file = output_folder/f'{chain_stem}_check.json'
    chain = read(chain_check_file)
    assert chain['generated_model'] == provenance(chain_file)
    inventory_file = folder/'route_inventory.json'
    inventory = read(inventory_file)
    root = ET.parse(chain_file).getroot()
    root.set('model', 'channel_first_endpoint_equivalent_drive')
    bodies = {body.get('name'): body for body in root.findall('.//body')}
    origins = {name: np.array(point) for name, point in chain['neutral_body_origins_m'].items()}
    origins['thorax'] = np.zeros(3)
    tendon = ET.SubElement(root, 'tendon')
    actuator = ET.SubElement(root, 'actuator')
    channels = []
    stiffness, damping = 200., 1.
    route_rows = list(inventory['routes'])
    with_girdle = args.girdle_candidates or args.upper_girdle
    if with_girdle:
        girdle_file = folder/('upper_girdle_endpoint_candidates.json' if args.upper_girdle
                              else 'girdle_endpoint_candidates.json')
        if args.arm_antagonists:
            girdle_file = folder/'arm_girdle_endpoint_candidates.json'
        route_rows.extend(read(girdle_file)['routes'])
    for route in route_rows:
        head = route['id']
        owner_a, owner_b = route['anatomical_origin_body'], route['distal_body']
        if 'candidate_endpoints_source_mm' in route:
            endpoints = route['candidate_endpoints_source_mm']
            reference_kind = 'unreviewed surface-patch sample'
        elif route['region'] == 'elbow':
            endpoints = route['pad_centres_source_mm']
            reference_kind = 'existing proximal and distal pad centres'
        else:
            ports = route['endpoint_ports_reference']
            endpoints = [ports['origin']['bone_seat_source_mm'], ports['distal']['bone_seat_source_mm']]
            reference_kind = 'historical origin/distal bone-seat reference; Y2 terminal details require remapping'
        for label, owner, point in zip(('origin', 'distal'), (owner_a, owner_b), endpoints):
            local = np.array(point)*.001-origins[owner]
            ET.SubElement(bodies[owner], 'site', name=f'{head}_{label}', pos=vector(local), size='.002',
                          rgba='.8 .1 .1 1')
        spatial = ET.SubElement(tendon, 'spatial', name=head, width='.0003', stiffness=str(stiffness),
                                damping=str(damping), springlength='-1')
        for label in ('origin', 'distal'):
            ET.SubElement(spatial, 'site', site=f'{head}_{label}')
        ET.SubElement(actuator, 'motor', name=f'pull_{head}', tendon=head, gear='-1',
                      ctrllimited='true', ctrlrange='0 60', forcelimited='true', forcerange='0 60')
        channels.append(dict(id=head, origin_body=owner_a, distal_body=owner_b,
                             endpoints_source_mm=endpoints, endpoint_reference=reference_kind,
                             straight_span=True, stiffness_N_per_m=stiffness, damping_Ns_per_m=damping,
                             max_commanded_tension_N=60., parameter_status='uncalibrated'))
    stem = 'girdle_candidate_equivalent_drive' if with_girdle else 'endpoint_equivalent_drive'
    if args.upper_girdle:
        stem = 'upper_girdle_equivalent_drive'
    if args.arm_antagonists:
        stem = 'arm_girdle_equivalent_drive'
    if args.geometry_mass:
        stem += '_geometry_mass'
    xml_file = output_folder/f'{stem}.xml'
    ET.ElementTree(root).write(xml_file, encoding='utf-8', xml_declaration=True)
    model = load_mujoco_model(xml_file)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    neutral = data.ten_length.copy()
    assert len(neutral) == len(route_rows) and np.all(neutral > 0)
    epsilon = 1e-6
    jacobian = np.zeros((model.ntendon, model.nv))
    for dof in range(model.nv):
        lengths = []
        for sign in (-1., 1.):
            probe = mujoco.MjData(model)
            velocity = np.zeros(model.nv)
            velocity[dof] = sign
            mujoco.mj_integratePos(model, probe.qpos, velocity, epsilon)
            mujoco.mj_forward(model, probe)
            lengths.append(probe.ten_length.copy())
        jacobian[:, dof] = (lengths[1]-lengths[0])/(2*epsilon)
    moment_arms = -jacobian
    blocks = {}
    for joint in ('SC', 'AC', 'GH', 'elbow'):
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, joint)
        count = 1 if joint == 'elbow' else 3
        start = model.jnt_dofadr[jid]
        block = moment_arms[:, start:start+count]
        blocks[joint] = dict(velocity_dof_indices=list(range(int(start), int(start)+count)),
                            rank=int(np.linalg.matrix_rank(block, tol=1e-7)),
                            moment_arms_mm=(block*1000.).tolist(),
                            note='neutral linear rank only; unilateral force closure is not established')
    pulses = {}
    for label, control in (('zero', 0.), ('flexor_10N', 10.)):
        state = mujoco.MjData(model)
        aid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, 'pull_FJ1486')
        state.ctrl[aid] = control
        for _ in range(200):
            mujoco.mj_step(model, state)
        assert np.isfinite(state.qpos).all() and np.isfinite(state.qvel).all()
        assert np.all(state.actuator_force >= -1e-8)
        elbow_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, 'elbow')
        pulses[label] = dict(elbow_angle_deg=float(np.degrees(state.qpos[model.jnt_qposadr[elbow_id]])),
                             qpos=state.qpos.tolist(), final_actuator_forces_N=state.actuator_force.tolist())
    from scipy.optimize import lsq_linear

    equilibrium_fit = lsq_linear(moment_arms.T, data.qfrc_bias-data.qfrc_passive,
                                 bounds=(0., 60.), tol=1e-10)
    sensitivity = []
    if args.upper_girdle:
        for bound in (60., 120., 240., np.inf):
            fit = lsq_linear(moment_arms.T, data.qfrc_bias-data.qfrc_passive,
                             bounds=(0., bound), tol=1e-10, max_iter=1000)
            sensitivity.append(dict(tension_upper_bound_N=None if np.isinf(bound) else bound,
                                    residual_norm_Nm=float(np.linalg.norm(fit.fun)),
                                    peak_tension_N=float(fit.x.max())))
    save(output_folder/f'{stem}_check.json', dict(
        chain_check=provenance(chain_check_file), route_inventory=provenance(inventory_file),
        generated_model=provenance(xml_file), mujoco_version=mujoco.__version__, channels=channels,
        neutral_span_lengths_mm=(neutral*1000.).tolist(), neutral_joint_blocks=blocks,
        full_neutral_moment_arm_rank=int(np.linalg.matrix_rank(moment_arms, tol=1e-7)),
        **({'girdle_candidates': provenance(girdle_file)} if with_girdle else {}),
        neutral_gravity_balance=dict(tension_bounds_N=[0., 60.], tensions_N=equilibrium_fit.x.tolist(),
                                     torque_residual_Nm=(moment_arms.T@equilibrium_fit.x
                                                        -(data.qfrc_bias-data.qfrc_passive)).tolist(),
                                     residual_norm_Nm=float(np.linalg.norm(equilibrium_fit.fun))),
        neutral_balance_residual_by_joint_Nm={
            name: equilibrium_fit.fun[block['velocity_dof_indices']].tolist()
            for name, block in blocks.items()},
        tension_bound_sensitivity=sensitivity,
        pulse_duration_seconds=.2, gravity_pulse_comparison=pulses,
        limitations=[
            'Straight spans are endpoint force surrogates, not routed steel or a collision-cleared layout',
            'Passive elasticity represents total endpoint compliance, with provisional parameters',
            'Only proximal and distal reaction forces are represented; thorax supply conduit reactions are absent',
            'Tension is force-commanded with a zero-release state; motor spool displacement and length-triggered slack are pending',
            'SC/AC ideal mobility, mass, inertia, damping and collision limitations remain from chain model',
            ('Thorax-to-girdle endpoints are unreviewed proximity candidates; upper support paths remain incomplete'
             if with_girdle else 'Eight local paths do not include thorax-to-shoulder-girdle actuators'),
            'No range acceptance, nonlinear lattice response, or friction prediction']))
    print('neutral rank', int(np.linalg.matrix_rank(moment_arms, tol=1e-7)),
          'joint ranks', {name: row['rank'] for name, row in blocks.items()})
    print('gravity versus flexor pulse', {name: row['elbow_angle_deg'] for name, row in pulses.items()})


if __name__ == '__main__':
    main()
