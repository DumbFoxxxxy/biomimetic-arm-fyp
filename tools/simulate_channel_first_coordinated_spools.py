"""Plan cable lengths, then simulate a free force-driven elbow return cycle."""
import sys

import numpy as np
from scipy.optimize import lsq_linear

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance
from simulate_channel_first_scapular_pads import pad_forces
from simulate_channel_first_spool_drive import cable_tensions


def endpoint_moment_matrix(model, data, channels):
    import mujoco

    matrix = []
    for channel in channels:
        jacobians, points = [], []
        for label in ('origin', 'distal'):
            site = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, f"{channel['id']}_{label}")
            jacobian = np.zeros((3, model.nv))
            rotation = np.zeros_like(jacobian)
            mujoco.mj_jacSite(model, data, jacobian, rotation, site)
            jacobians.append(jacobian)
            points.append(data.site_xpos[site].copy())
        direction = points[1]-points[0]
        direction /= np.linalg.norm(direction)
        matrix.append(-direction@(jacobians[1]-jacobians[0]))
    return np.array(matrix).T


def requested_angle(time):
    if .5 < time < 2.:
        fraction = (time-.5)/1.5
        return 30.*(.5-.5*np.cos(np.pi*fraction))
    if 2. <= time <= 3.:
        return 30.
    if 3. < time < 4.5:
        fraction = (time-3.)/1.5
        return 30.*(.5+.5*np.cos(np.pi*fraction))
    return 0.


def main():
    import mujoco

    sys.path.insert(0, str(ROOT))
    from project_runtime import load_mujoco_model

    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    pad_file = folder/'scapular_pad_dynamic_screen.json'
    pad_report = read(pad_file)
    endpoint_file = ROOT/pad_report['endpoint_report']['file']
    assert provenance(endpoint_file) == pad_report['endpoint_report']
    endpoint = read(endpoint_file)
    support_file = ROOT/pad_report['support_direction_screen']['file']
    assert provenance(support_file) == pad_report['support_direction_screen']
    support = read(support_file)
    model_file = ROOT/pad_report['model']['file']
    assert provenance(model_file) == pad_report['model']
    model = load_mujoco_model(model_file)
    elbow_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, 'elbow')
    elbow_q = model.jnt_qposadr[elbow_id]
    stiffness, damping = 10000., .5
    preferred = np.asarray(support['compact_force_solution']['cable_tensions_N'])
    plan = []
    regularizer = 1e-4
    for angle in range(31):
        state = mujoco.MjData(model)
        state.qpos[elbow_q] = np.radians(angle)
        mujoco.mj_forward(model, state)
        state.qfrc_applied[:] = 0.
        pad_forces(model, state, pad_report['pads'], 1., 3.)
        target = state.qfrc_bias-state.qfrc_passive-state.qfrc_applied
        moment = endpoint_moment_matrix(model, state, endpoint['channels'])
        if angle == 0:
            recorded = np.zeros_like(moment)
            for block in endpoint['neutral_joint_blocks'].values():
                recorded[block['velocity_dof_indices']] = np.array(block['moment_arms_mm']).T*.001
            assert np.allclose(moment, recorded, atol=1e-7)
        fit = lsq_linear(np.vstack([moment, regularizer*np.eye(model.nu)]),
                         np.r_[target, regularizer*preferred], bounds=(0., 60.), tol=1e-10)
        preferred = fit.x
        available = state.ten_length-preferred/stiffness
        plan.append(dict(reference_elbow_deg=angle, planned_tension_N=preferred.tolist(),
                         available_length_m=available.tolist(),
                         static_residual_Nm=float(np.linalg.norm(moment@preferred-target))))
    available_table = np.array([row['available_length_m'] for row in plan])
    data = mujoco.MjData(model)
    initial = data.qpos.copy()
    trace = []
    saturated = np.zeros(model.nu, dtype=int)
    bottomed, left_edge = False, False
    for step in range(round(5./model.opt.timestep)):
        mujoco.mj_step1(model, data)
        angle = requested_angle(float(data.time))
        low = min(29, int(np.floor(angle)))
        fraction = angle-low
        available = (1.-fraction)*available_table[low]+fraction*available_table[low+1]
        tension, slack = cable_tensions(data.ten_length, data.ten_velocity, available, stiffness, damping)
        saturated += tension > 60.
        data.ctrl[:] = np.minimum(tension, 60.)
        data.qfrc_applied[:] = 0.
        loads = pad_forces(model, data, pad_report['pads'], 1., 3.)
        bottomed |= any(row['pad_bottomed'] for row in loads)
        left_edge |= any(row['over_pad_edge'] for row in loads)
        if step % 50 == 0:
            delta = np.zeros(model.nv)
            mujoco.mj_differentiatePos(model, delta, 1., initial, data.qpos)
            trace.append(dict(time_s=float(data.time), planned_elbow_deg=angle,
                              actual_elbow_deg=float(np.degrees(data.qpos[elbow_q])),
                              joint_displacement_deg=np.degrees(delta).tolist(),
                              available_length_mm=(available*1000.).tolist(),
                              cable_tension_N=data.ctrl.tolist(), slack=slack.tolist(), pad_loads=loads))
        mujoco.mj_step2(model, data)
        assert np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
    actual = np.array([row['actual_elbow_deg'] for row in trace])
    planned = np.array([row['planned_elbow_deg'] for row in trace])
    displacement = np.array([row['joint_displacement_deg'] for row in trace])
    shoulder_drift = {name: float(np.linalg.norm(displacement[:, block['velocity_dof_indices']], axis=1).max())
                      for name, block in endpoint['neutral_joint_blocks'].items() if name != 'elbow'}
    save(folder/'coordinated_spool_cycle_screen.json', dict(
        pad_model=provenance(pad_file), endpoint_report=provenance(endpoint_file), model=provenance(model_file),
        plan=plan, trace=trace, duration_s=5., drive_stiffness_N_per_m=stiffness,
        drive_damping_Ns_per_m=damping, allocation_regularizer_m=regularizer,
        channel_ids=[row['id'] for row in endpoint['channels']],
        spool_travel_ranges_mm=((available_table.max(axis=0)-available_table.min(axis=0))*1000.).tolist(),
        summary=dict(actual_elbow_min_deg=float(actual.min()), actual_elbow_max_deg=float(actual.max()),
                     actual_elbow_final_deg=float(np.degrees(data.qpos[elbow_q])),
                     max_tracking_error_deg=float(np.max(abs(actual-planned))),
                     peak_shoulder_joint_rotation_vector_norm_deg=shoulder_drift,
                     saturated_steps=saturated.tolist(), pad_bottomed=bottomed, left_pad_footprint=left_edge),
        scope='Reference poses only generate offline cable length commands. Simulated joint positions '
              'are integrated from forces. Twenty-seven endpoint spans, uncalibrated compliance and pads. '
              'No routed supply, friction, bone collision, installed anchors or manufacturing acceptance.'))
    print('actual elbow min/max/final', actual.min(), actual.max(), np.degrees(data.qpos[elbow_q]),
          'tracking peak', np.max(abs(actual-planned)), 'saturation', saturated.sum(),
          'bottomed', bottomed, 'left footprint', left_edge)


if __name__ == '__main__':
    main()
