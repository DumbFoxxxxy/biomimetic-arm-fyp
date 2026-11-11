"""Run displacement-commanded unilateral cable trials with passive TPU and pads."""
import sys
import argparse

import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance
from simulate_channel_first_scapular_pads import pad_forces


def cable_tensions(length, velocity, available_length, stiffness, damping):
    extension = length-available_length
    force = np.where(extension > 0., np.maximum(0., stiffness*extension+damping*velocity), 0.)
    assert np.all(force >= 0.) and np.all(force[extension <= 0.] == 0.)
    return force, extension <= 0.


def spool_command(time, trial, neutral, flexor):
    result = neutral.copy()
    if trial == 'flexor_cycle':
        times = [.5, 1.5, 2., 3.]
        pull = np.interp(time, times, [0., .006, .006, 0.])
        result[flexor] -= pull
    elif trial == 'flexor_release' and time >= .5:
        result[flexor] += .010*min(1., (time-.5)/.25)
    return result


def main():
    import mujoco

    sys.path.insert(0, str(ROOT))
    from project_runtime import load_mujoco_model

    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--unified-sc-ac',action='store_true')
    args=parser.parse_args()
    if args.unified_sc_ac:folder=folder/'unified_SC_AC'
    pad_file = folder/'scapular_pad_dynamic_screen.json'
    pad_report = read(pad_file)
    model_file = ROOT/pad_report['model']['file']
    assert provenance(model_file) == pad_report['model']
    support_file = ROOT/pad_report['support_direction_screen']['file']
    assert provenance(support_file) == pad_report['support_direction_screen']
    support = read(support_file)
    endpoint_file = ROOT/pad_report['endpoint_report']['file']
    assert provenance(endpoint_file) == pad_report['endpoint_report']
    endpoint = read(endpoint_file)
    model = load_mujoco_model(model_file)
    reference = mujoco.MjData(model)
    mujoco.mj_forward(model, reference)
    flexor = next(index for index, row in enumerate(endpoint['channels']) if row['id'] == 'FJ1486')
    elbow = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, 'elbow')
    elbow_q = model.jnt_qposadr[elbow]
    stiffness, damping = 10000., .5
    nominal_tension = np.array(support['compact_force_solution']['cable_tensions_N'])
    neutral_available_length = reference.ten_length-nominal_tension/stiffness
    initial_force, _ = cable_tensions(reference.ten_length, reference.ten_velocity,
                                       neutral_available_length, stiffness, damping)
    assert np.allclose(initial_force, nominal_tension, atol=1e-8)
    # Explicit constitutive checks independent of a simulation trajectory.
    zero_force, slack = cable_tensions(np.array([.2, .2]), np.array([1., -1.]),
                                       np.array([.21, .21]), stiffness, damping)
    assert slack.all() and np.all(zero_force == 0.)
    trials = []
    for trial in ('hold', 'flexor_cycle', 'flexor_release'):
        data = mujoco.MjData(model)
        initial_qpos = data.qpos.copy()
        trace = []
        saturation_steps = np.zeros(model.nu, dtype=int)
        bottomed = False
        left_edge = False
        peak_tension = np.zeros(model.nu)
        for step in range(round(4./model.opt.timestep)):
            mujoco.mj_step1(model, data)
            available = spool_command(float(data.time), trial, neutral_available_length, flexor)
            requested, slack = cable_tensions(data.ten_length, data.ten_velocity, available, stiffness, damping)
            saturation_steps += requested > 60.
            data.ctrl[:] = np.minimum(requested, 60.)
            data.qfrc_applied[:] = 0.
            loads = pad_forces(model, data, pad_report['pads'], 1., 3.)
            peak_tension = np.maximum(peak_tension, data.ctrl)
            bottomed |= any(row['pad_bottomed'] for row in loads)
            left_edge |= any(row['over_pad_edge'] for row in loads)
            if step % 50 == 0:
                joint_delta=np.zeros(model.nv)
                mujoco.mj_differentiatePos(model,joint_delta,1.,initial_qpos,data.qpos)
                trace.append(dict(time_s=float(data.time), elbow_deg=float(np.degrees(data.qpos[elbow_q])),
                                  joint_displacement_rad=joint_delta.tolist(),
                                  cable_length_mm=(data.ten_length*1000.).tolist(),
                                  available_length_mm=(available*1000.).tolist(),
                                  cable_tension_N=data.ctrl.tolist(), slack=slack.tolist(),
                                  pad_loads=loads))
            mujoco.mj_step2(model, data)
            assert np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
        elbow_angles = [row['elbow_deg'] for row in trace]
        summary = dict(
            trial=trial, trace=trace, elbow_min_deg=min(elbow_angles), elbow_max_deg=max(elbow_angles),
            elbow_final_deg=float(np.degrees(data.qpos[elbow_q])),
            final_qvel_rad_per_s=data.qvel.tolist(), peak_tension_N=peak_tension.tolist(),
            saturated_step_counts=saturation_steps.tolist(),
            flexor_slack_sample_count=sum(row['slack'][flexor] for row in trace),
            pad_bottomed=bottomed, left_pad_footprint=left_edge)
        trials.append(summary)
        print(trial, 'elbow min/max/final', summary['elbow_min_deg'], summary['elbow_max_deg'],
              summary['elbow_final_deg'], 'saturated', int(saturation_steps.sum()),
              'bottomed', bottomed, 'left edge', left_edge, flush=True)
    save(folder/'spool_drive_dynamic_screen.json', dict(
        pad_dynamic_report=provenance(pad_file), endpoint_report=provenance(endpoint_file),
        model=provenance(model_file), channel_ids=[row['id'] for row in endpoint['channels']],
        nominal_available_lengths_mm=(neutral_available_length*1000.).tolist(),
        effective_drive_stiffness_N_per_m=stiffness, effective_drive_damping_Ns_per_m=damping,
        duration_s=4., constitutive_slack_check=True, trials=trials,
        scope='Available-length commands with effective cable/drive compliance, not measured steel modulus. '
              'Zero force when slack. Force limited to 60 N; motor dynamics under saturation unmodelled. '
              'Passive endpoint TPU and finite normal pads retained. No routed supply, friction or collision.'))


if __name__ == '__main__':
    main()
