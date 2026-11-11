"""Run finite-footprint compression-only scapulothoracic support trials."""
import sys
import argparse
import xml.etree.ElementTree as ET

import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector


def pad_forces(model, data, pads, stiffness_scale, damping):
    import mujoco

    body = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, 'scapula')
    loads = []
    for pad in pads:
        point = data.xpos[body]+data.xmat[body].reshape(3, 3)@np.array(pad['scapula_local_point_m'])
        normal = np.asarray(pad['normal_source'])
        offset = point-np.array(pad['rib_point_m'])
        separation = float(offset@normal)
        tangent = float(np.linalg.norm(offset-separation*normal))
        jacobian = np.zeros((3, model.nv))
        rotational = np.zeros_like(jacobian)
        mujoco.mj_jac(model, data, jacobian, rotational, point, body)
        velocity = float(normal@(jacobian@data.qvel))
        compression = pad['rest_thickness_m']-separation
        force = 0.
        if compression > 0. and tangent <= pad['contact_radius_m']:
            force = max(0., stiffness_scale*pad['stiffness_N_per_m']*compression-damping*velocity)
            data.qfrc_applied[:] += jacobian.T@(normal*force)
        loads.append(dict(normal_force_N=force, separation_mm=separation*1000.,
                          tangential_offset_mm=tangent*1000.,
                          over_pad_edge=tangent > pad['contact_radius_m'],
                          pad_bottomed=separation <= 0.))
    return loads


def main():
    import mujoco

    sys.path.insert(0, str(ROOT))
    from project_runtime import load_mujoco_model

    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--unified-sc-ac',action='store_true')
    args=parser.parse_args()
    if args.unified_sc_ac:folder=folder/'unified_SC_AC'
    support_file = folder/'scapulothoracic_support_direction_screen.json'
    support = read(support_file)
    endpoint_report_file = ROOT/support['endpoint_model']['file']
    assert provenance(endpoint_report_file) == support['endpoint_model']
    endpoint = read(endpoint_report_file)
    model_file = ROOT/endpoint['generated_model']['file']
    assert provenance(model_file) == endpoint['generated_model']
    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    body = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, 'scapula')
    compact = support['compact_force_solution']
    pads = []
    stiffness = 10000.
    for index in compact['candidate_indices']:
        item = support['candidates'][index]
        surface = np.array(item['scapular_surface_source_mm'])*.001
        rest = item['gap_mm']*.001+compact['support_forces_N'][index]/stiffness
        pads.append(dict(candidate_index=index, rib_point_m=(np.array(item['rib_surface_source_mm'])*.001).tolist(),
                         scapula_local_point_m=(data.xmat[body].reshape(3, 3).T@(surface-data.xpos[body])).tolist(),
                         normal_source=item['normal_force_on_scapula_source'],
                         rest_thickness_m=rest, contact_radius_m=.012, stiffness_N_per_m=stiffness,
                         nominal_preload_N=compact['support_forces_N'][index]))
    root = ET.parse(model_file).getroot()
    root.find('option').set('timestep', '.0002')
    thorax = root.find(".//body[@name='thorax']")
    for pad in pads:
        normal = np.asarray(pad['normal_source'])
        centre = np.array(pad['rib_point_m'])+normal*pad['rest_thickness_m']/2.
        z = np.array([0., 0., 1.])
        quaternion = np.r_[1.+z@normal, np.cross(z, normal)]
        quaternion /= np.linalg.norm(quaternion)
        ET.SubElement(thorax, 'geom', name=f"pad_envelope_{pad['candidate_index']}", type='cylinder',
                      pos=vector(centre), quat=vector(quaternion),
                      size=vector([pad['contact_radius_m'], pad['rest_thickness_m']/2.]),
                      rgba='.2 .7 .3 .35', contype='0', conaffinity='0', group='2')
    output = folder/'scapular_pad_equivalent_drive.xml'
    ET.ElementTree(root).write(output, encoding='utf-8', xml_declaration=True)
    model = load_mujoco_model(output)
    dt = model.opt.timestep
    trials = []
    for scale, damping in ((1., 3.), (.5, 3.), (2., 3.), (1., 1.), (1., 10.)):
        state = mujoco.MjData(model)
        initial = state.qpos.copy()
        state.ctrl[:] = compact['cable_tensions_N']
        trace = []
        peak_force = np.zeros(len(pads))
        bottomed = False
        left_edge = False
        for step in range(round(2./dt)):
            mujoco.mj_step1(model, state)
            state.qfrc_applied[:] = 0.
            loads = pad_forces(model, state, pads, scale, damping)
            peak_force = np.maximum(peak_force, [row['normal_force_N'] for row in loads])
            bottomed |= any(row['pad_bottomed'] for row in loads)
            left_edge |= any(row['over_pad_edge'] for row in loads)
            mujoco.mj_step2(model, state)
            assert np.isfinite(state.qpos).all() and np.isfinite(state.qvel).all()
            if step % 50 == 0:
                delta = np.zeros(model.nv)
                mujoco.mj_differentiatePos(model, delta, 1., initial, state.qpos)
                trace.append(dict(time_s=float(state.time), joint_displacement_rad=delta.tolist(),
                                  pad_loads=loads))
        delta = np.zeros(model.nv)
        mujoco.mj_differentiatePos(model, delta, 1., initial, state.qpos)
        trials.append(dict(stiffness_scale=scale, pad_damping_Ns_per_m=damping, trace=trace,
                           final_joint_displacement_deg=np.degrees(delta).tolist(),
                           final_qvel_rad_per_s=state.qvel.tolist(), peak_pad_force_N=peak_force.tolist(),
                           pad_bottomed=bottomed, left_pad_footprint=left_edge))
        print(scale, damping, 'max angle change', float(np.degrees(abs(delta)).max()),
              'bottomed', bottomed, 'left footprint', left_edge, flush=True)
    save(folder/'scapular_pad_dynamic_screen.json', dict(
        support_direction_screen=provenance(support_file), endpoint_report=provenance(endpoint_report_file),
        model=provenance(output), pads=pads, duration_s=2., timestep_s=dt, trials=trials,
        scope='Compression-only flat pad law applied at moving scapular material points. Fixed candidate '
              'cable tensions and ideal joints. Pads are uncalibrated, with no tangential friction, '
              'rolling contact, mechanical mounts, routed cable reaction or bone collision checks. '
              'XML alone does not apply pad forces; run this script to reproduce the force model.'))


if __name__ == '__main__':
    main()
