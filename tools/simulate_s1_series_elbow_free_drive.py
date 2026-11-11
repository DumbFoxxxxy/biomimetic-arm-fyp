"""Drive the fixed-girdle S1 and series elbow together with local tensile inputs."""

from itertools import combinations
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np
from scipy.optimize import lsq_linear
from scipy.spatial.transform import Rotation

from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from simulate_elbow_series_slider import SLIDE, SLIDE_DAMPING, SLIDE_STIFFNESS
from study_s1_local_force_geometry import load_reference


def vector(values):
    return " ".join(str(float(value)) for value in values)


def make_model(elbow_xml, shoulder):
    root = ET.parse(elbow_xml).getroot()
    root.set("model", "S1 and series elbow local equivalent drive")
    humerus = root.find(".//body[@name='humerus']")
    origin = np.fromstring(humerus.get("pos"), sep=" ")
    centre = np.array(shoulder["gh_centre_source_mm"]) / 1000.
    if not np.allclose(origin, centre, atol=1e-9, rtol=0.):
        raise ValueError("Shoulder and elbow use different GH frames")
    ET.SubElement(humerus, "joint", name="GH", type="ball", damping=".15")
    fixed = ET.SubElement(root.find("worldbody"), "body", name="fixed_shoulder_girdle")
    for head, row in shoulder["neutral"].items():
        proximal, distal = f"{head}_proximal", f"{head}_distal"
        ET.SubElement(fixed, "site", name=proximal, size=".002", pos=vector(np.array(row["proximal_port_source_mm"]) / 1000.))
        ET.SubElement(humerus, "site", name=distal, size=".002",
                      pos=vector(np.array(row["distal_port_source_mm"]) / 1000. - centre))
        line = ET.SubElement(root.find("tendon"), "spatial", name=head, width=".0006",
                             springlength=str(row["chord_mm"] / 1000.), stiffness="200", damping="1")
        for site in (proximal, distal):
            ET.SubElement(line, "site", site=site)
        ET.SubElement(root.find("actuator"), "motor", name=head, tendon=head, gear="-1",
                      ctrllimited="true", ctrlrange="0 25", forcelimited="true", forcerange="0 25")
    return ET.tostring(root, encoding="unicode")


def simulate(model, controls, axis_frame):
    data = mujoco.MjData(model)
    dt = model.opt.timestep
    steps = round(2. / dt)
    slide = model.joint(SLIDE)
    slide_q, slide_v = slide.qposadr.item(), slide.dofadr.item()
    elbow_q = model.joint("elbow").qposadr.item()
    gh_q = model.joint("GH").qposadr.item()
    model.dof_frictionloss[slide_v] = .15
    mujoco.mj_setConst(model, data)
    neutral_lengths = model.tendon_lengthspring[:, 0].copy()
    ids = [model.actuator(name).id for name in controls]
    magnitudes = np.array(list(controls.values()))
    rows, states = [], []
    sample = round(.01 / dt)
    for step in range(steps + 1):
        time = step * dt
        mujoco.mj_step1(model, data)
        command = float(np.interp(time, [0., .2, .65, 1.15, 1.55, 2.], [0., 0., 1., 1., 0., 0.]))
        data.ctrl[ids] = command * magnitudes
        data.qfrc_applied[slide_v] = -max(0., .25 + SLIDE_STIFFNESS * data.qpos[slide_q]
                                       + SLIDE_DAMPING * data.qvel[slide_v])
        if step % sample == 0 or step == steps:
            quaternion = data.qpos[gh_q:gh_q + 4]
            turn = Rotation.from_quat(np.roll(quaternion, -1))
            components = axis_frame @ np.rad2deg(turn.as_rotvec())
            direction = turn.apply(axis_frame[2])
            elevation = np.rad2deg(np.arccos(np.clip(direction @ axis_frame[2], -1., 1.)))
            rows.append([time, elevation, *components, np.rad2deg(data.qpos[elbow_q]), data.qpos[slide_q] * 1000.,
                         *((data.ten_length - neutral_lengths) * 1000.), *data.ctrl, *data.qpos])
            states.append(data.qpos.copy())
        if step < steps:
            mujoco.mj_step2(model, data)
    values = np.array(rows)
    if not np.isfinite(values).all() or any(item.number for item in data.warning):
        raise RuntimeError("Invalid shoulder–elbow simulation")
    return values, states


def neutral_shoulder_balance(model, channel_names):
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    start = model.joint("GH").dofadr.item()
    block = slice(start, start + 3)
    target = (data.qfrc_bias - data.qfrc_passive)[block].copy()
    columns = []
    for name in channel_names:
        data.ctrl[:] = 0.
        data.ctrl[model.actuator(name).id] = 1.
        mujoco.mj_forward(model, data)
        columns.append(data.qfrc_actuator[block].copy())
    influence = np.array(columns).T
    fit = lsq_linear(influence, target, bounds=(0., 25.), tol=1e-12)
    if not fit.success:
        raise RuntimeError(fit.message)
    return dict(channels=list(channel_names), influence_Nm_per_N=influence.tolist(),
                required_moment_Nm=target.tolist(), residual_moment_Nm=fit.fun.tolist(),
                required_norm_Nm=float(np.linalg.norm(target)), residual_norm_Nm=float(np.linalg.norm(fit.fun)),
                fitted_input_N=fit.x.tolist(), input_bound_N=[0., 25.])


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    shoulder_path = DEST.parent / "s1_local_muscle_load/gh_endpoint_geometry.json"
    shoulder = read(shoulder_path)
    load_reference(shoulder["shoulder_assembly"])
    series_path = DEST.parent / "elbow_local_muscle_load/elbow_series_free_slider.json"
    series = read(series_path)
    xml_input = ROOT / series["local_model"]["file"]
    if provenance(xml_input) != series["local_model"]:
        raise ValueError("Elbow model hash differs from the screened slider")
    assembly_path = folder / "assembly_FJ1504_local_interface_bidirectional.json"
    assembly = read(assembly_path)
    xml_path = folder / "shoulder_elbow_free_drive.xml"
    xml = make_model(xml_input, shoulder)
    xml_path.write_text(xml, encoding="utf-8")
    cases = dict(gravity={}, shoulder={"FJ1468": 10., "FJ1467": 10.},
                 elbow={"FJ1486": 25.}, combined={"FJ1468": 10., "FJ1467": 10., "FJ1486": 25.})
    balance = neutral_shoulder_balance(mujoco.MjModel.from_xml_string(xml), shoulder["neutral"])
    frame = np.array(shoulder["axes_source"])
    keys = list(assembly["sources"])
    owners = assembly["owners"]
    solids = [load_component(assembly["sources"][key]) for key in keys]
    pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                      if owners[keys[i]] != owners[keys[j]]
                      and (owners[keys[i]] in MOBILE or owners[keys[j]] in MOBILE)])
    outcomes = {}
    for label, controls in cases.items():
        model = mujoco.MjModel.from_xml_string(xml)
        data = mujoco.MjData(model)
        mujoco.mj_kinematics(model, data)
        origins = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
        columns = ["time_s", "humerus_elevation_deg", "GH_flexion_projection_deg", "GH_abduction_projection_deg",
                   "GH_axial_projection_deg", "elbow_deg", "carriage_mm",
                   *(f"{model.tendon(i).name}_span_change_mm" for i in range(model.ntendon)),
                   *(f"{model.actuator(i).name}_input_N" for i in range(model.nu)),
                   *(f"qpos_{i}" for i in range(model.nq))]
        rows, states = simulate(model, controls, frame)
        output = folder / f"free_drive_{label}.csv"
        write_csv(output, columns, rows)
        first = None
        for index, state in enumerate(states):
            data.qpos[:] = state
            mujoco.mj_kinematics(model, data)
            poses = transforms(model, data, origins)
            positioned = [part.transform(poses[owners[key]]) if owners[key] in poses else part
                          for key, part in zip(keys, solids, strict=True)]
            hits = find_hits(positioned, keys, pairs)
            if hits:
                first = dict(time_s=float(rows[index, 0]), elevation_deg=float(rows[index, 1]),
                             elbow_deg=float(rows[index, 5]), hits=hits)
                break
        outcomes[label] = dict(commands_N=controls, samples=len(states), first_rigid_interference=first,
                               model_peak_elevation_deg=float(rows[:, 1].max()),
                               model_peak_elbow_deg=float(rows[:, 5].max()), final_elevation_deg=float(rows[-1, 1]),
                               final_elbow_deg=float(rows[-1, 5]), model_peak_carriage_mm=float(rows[:, 6].max()),
                               trajectory=provenance(output))
        print(label, outcomes[label], flush=True)
    save(folder / "shoulder_elbow_free_drive.json", dict(
        generator=provenance(ROOT / "tools/simulate_s1_series_elbow_free_drive.py"),
        model=provenance(xml_path), shoulder_endpoints=provenance(shoulder_path), elbow=provenance(series_path),
        geometry_assembly=provenance(assembly_path), dof=model.nv, qpos_count=model.nq,
        audited_cross_body_pairs=len(pairs), sample_step_s=.01, neutral_shoulder_balance=balance, results=outcomes,
        scope="Fixed shoulder girdle with free GH ball joint, elbow hinge, two cap swivels and one series slider."
              " Four S1 endpoint force surrogates and two elbow channels."
              " The two Y muscles and external supply are excluded. The elbow extension stop uses original contact."
              " Other rigid intersections are identified by an independent sampling audit, not solved by dynamics."
              " Peaks after a first intersection are diagnostic only. TPU and wire geometry are not checked."))


if __name__ == "__main__":
    main()
