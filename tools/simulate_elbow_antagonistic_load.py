"""Evaluate local elbow muscle geometry and force-driven antagonistic response."""

import copy
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read


HEADS = ("FJ1486", "FJ1477")
JOINTS = ("elbow", *(f"{head}_swivel" for head in HEADS))


def build_local_model(source, check):
    humerus = copy.deepcopy(ET.parse(source).find(".//body[@name='humerus']"))
    humerus.set("pos", " ".join(map(str, check["neutral_origins_m"]["humerus"])))
    kept_sites = {f"{head}_{end}" for head in HEADS for end in ("origin", "distal")}
    for body in humerus.iter("body"):
        for node in list(body):
            keep = (node.tag in ("body", "inertial")
                    or node.tag == "joint" and node.get("name") in JOINTS
                    or node.tag == "site" and node.get("name") in kept_sites)
            if not keep:
                body.remove(node)
    root = ET.Element("mujoco", model="elbow local equivalent loads")
    ET.SubElement(root, "compiler", angle="radian")
    ET.SubElement(root, "option", timestep="0.0005", gravity="0 0 -9.81",
                  integrator="implicitfast")
    ET.SubElement(root, "worldbody").append(humerus)
    tendons = ET.SubElement(root, "tendon")
    for head in HEADS:
        tendon = ET.SubElement(tendons, "spatial", name=head, width=".001")
        for end in ("origin", "distal"):
            ET.SubElement(tendon, "site", site=f"{head}_{end}")
    return ET.tostring(root, encoding="unicode")


def measure(model, data, site_ids):
    lengths = []
    gradients = []
    for origin, distal in site_ids:
        separation = data.site_xpos[distal] - data.site_xpos[origin]
        length = np.linalg.norm(separation)
        jacobian = np.zeros((3, model.nv))
        mujoco.mj_jacSite(model, data, jacobian, None, distal)
        # Both origins are on the fixed humerus.
        lengths.append(length)
        gradients.append(separation @ jacobian / length)
    return np.asarray(lengths), np.asarray(gradients)


def screen_geometry(model, site_ids, neutral, elbow_qpos, elbow_dof):
    data = mujoco.MjData(model)
    rows = []
    max_error = 0.
    for angle in range(151):
        position = np.deg2rad(angle)
        data.qpos[elbow_qpos] = position
        mujoco.mj_forward(model, data)
        lengths, gradients = measure(model, data, site_ids)
        arms = -gradients[:, elbow_dof]
        gravity = -data.qfrc_bias[elbow_dof]
        # Compare with the middle stiffness scenario and static, neutral cap angles.
        spring_torques = 200. * (lengths - neutral) * arms
        required_tension = -(gravity + spring_torques.sum()) / arms[0] if arms[0] > 0. else None
        offset = 1e-6
        sampled = []
        for value in (position - offset, position + offset):
            data.qpos[elbow_qpos] = value
            mujoco.mj_kinematics(model, data)
            sampled.append([np.linalg.norm(data.site_xpos[b] - data.site_xpos[a])
                            for a, b in site_ids])
        finite_diff = np.diff(sampled, axis=0)[0] / (2. * offset)
        max_error = max(max_error, float(np.max(abs(finite_diff - gradients[:, elbow_dof]))))
        rows.append([angle, *((lengths - neutral) * 1000.), *(arms * 1000.),
                     gravity, *spring_torques, required_tension])
    if max_error > 1e-8:
        raise RuntimeError("Moment arm does not match the finite-difference length derivative")
    return rows, max_error


def simulate(model, site_ids, neutral, addresses, stiffness, peak, active, trace_interval=.01):
    data = mujoco.MjData(model)
    trace = []
    stops = set()
    contact_steps = 0
    minimum_contact_distance = None
    timestep = model.opt.timestep
    for step in range(round(4. / timestep) + 1):
        mujoco.mj_step1(model, data)
        lengths, gradients = measure(model, data, site_ids)
        length_rates = gradients @ data.qvel
        time = step * timestep
        tension = peak * np.interp(time, [0., .5, 1.5, 2.5, 3., 4.], [0., 0., 1., 1., 0., 0.])
        drive = np.zeros(len(HEADS))
        if active is not None:
            drive[active] = tension
        # Positive drive pulls. TPU restores both shortening and extension.
        passive = stiffness * (lengths - neutral) + length_rates
        data.qfrc_applied[:] = -(drive + passive) @ gradients
        for index in range(data.nefc):
            if data.efc_type[index] == mujoco.mjtConstraint.mjCNSTR_LIMIT_JOINT:
                stops.add(int(data.efc_id[index]))
        if data.ncon:
            contact_steps += 1
            distance = float(min(data.contact.dist[:data.ncon]))
            minimum_contact_distance = (distance if minimum_contact_distance is None
                                        else min(distance, minimum_contact_distance))
        if step % round(trace_interval / timestep) == 0:
            trace.append([time, *np.rad2deg(data.qpos[addresses]),
                          *((lengths - neutral) * 1000.), *drive, *passive])
        if time < 4. - timestep / 2.:
            mujoco.mj_step2(model, data)
    values = np.asarray(trace)
    if not np.isfinite(values).all() or any(w.number for w in data.warning):
        raise RuntimeError("Invalid elbow simulation")
    return trace, dict(
        active_muscle=HEADS[active] if active is not None else None, peak_tension_N=peak,
        stiffness_N_per_m=stiffness, damping_Ns_per_m=1.,
        elbow_min_deg=float(values[:, 1].min()), elbow_max_deg=float(values[:, 1].max()),
        elbow_final_deg=float(values[-1, 1]), cap_peak_abs_deg=abs(values[:, 2:4]).max(axis=0).tolist(),
        muscle_min_length_change_mm=values[:, 4:6].min(axis=0).tolist(),
        muscle_max_length_change_mm=values[:, 4:6].max(axis=0).tolist(),
        joint_limits_engaged=[model.joint(idx).name for idx in sorted(stops)],
        contact_steps=contact_steps,
        minimum_contact_distance_mm=(None if minimum_contact_distance is None
                                     else minimum_contact_distance * 1000.),
    )


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    check_path = DEST / "elbow_passive_model/model_check.json"
    check = read(check_path)
    source = ROOT / check["model"]["file"]
    if provenance(source) != check["model"]:
        raise RuntimeError("Passive elbow model fingerprint changed")
    output = DEST.parent / "elbow_local_muscle_load"
    output.mkdir(exist_ok=True)
    model_file = output / "elbow_local_force.xml"
    model_file.write_text(build_local_model(source, check), encoding="utf-8")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    site_ids = [[model.site(f"{head}_{end}").id for end in ("origin", "distal")] for head in HEADS]
    if model.nv != len(JOINTS):
        raise RuntimeError("Expected one elbow hinge and two passive cap hinges")
    for port in check["ports"]:
        site = model.site(f"{port['head']}_{port['role']}").id
        np.testing.assert_allclose(data.site_xpos[site] * 1000., port["source_position_mm"], atol=1e-6)
    neutral, neutral_gradient = measure(model, data, site_ids)
    np.testing.assert_allclose(neutral, data.ten_length, atol=1e-10)
    elbow_dof = model.joint("elbow").dofadr.item()
    geometry, error = screen_geometry(model, site_ids, neutral, addresses[0], elbow_dof)
    geometry_path = output / "elbow_geometry.csv"
    write_csv(geometry_path, ["elbow_deg", *(f"{h}_length_change_mm" for h in HEADS),
                              *(f"{h}_signed_moment_arm_mm" for h in HEADS), "gravity_torque_Nm",
                              *(f"{h}_spring_torque_at_k200_Nm" for h in HEADS),
                              "flexor_force_to_balance_at_k200_N"], geometry)
    cases = []
    for stiffness in (100., 200., 400.):
        for name, active, peak in (("baseline", None, 0.), ("flexor", 0, 5.), ("extensor", 1, 5.)):
            rows, case = simulate(model, site_ids, neutral, addresses, stiffness, peak, active)
            path = output / f"elbow_{name}_k{stiffness:g}.csv"
            write_csv(path, ["time_s", *(f"{j}_deg" for j in JOINTS),
                             *(f"{h}_length_change_mm" for h in HEADS),
                             *(f"{h}_drive_N" for h in HEADS),
                             *(f"{h}_TPU_restoring_N" for h in HEADS)], rows)
            cases.append(dict(**case, trace=provenance(path)))
            print(name, stiffness, "elbow min/max/final",
                  *(round(case[f"elbow_{field}_deg"], 3) for field in ("min", "max", "final")))
    report = dict(
        source=check["model"], source_interfaces=provenance(check_path), solver=mujoco.__version__,
        local_model=provenance(model_file), generator=provenance(Path(__file__).resolve()),
        order=list(HEADS), neutral_lengths_mm=(neutral * 1000.).tolist(),
        neutral_signed_moment_arms_mm=(-neutral_gradient[:, elbow_dof] * 1000.).tolist(),
        geometry=provenance(geometry_path), finite_difference_max_error_m=error,
        geometry_scope="Elbow angles sampled at 1 degree with both caps held neutral; no collision tests.",
        static_budget_scope="Fixed neutral cap angles, both muscle springs at 200 N/m, zero velocity. "
        "Flexor tension balances signed elbow gravity and both passive forces when flexor arm is positive.",
        dynamics_scope="Fixed humerus at neutral pose. Existing forearm/cap masses, gravity, hinge axes, "
        "joint damping, cap spring and software elbow limits retained. Assumed axial spring/damping "
        "and applied tension, no contact, cable, friction or lattice deformation. Elbow and caps free.",
        cases=cases,
    )
    (output / "elbow_antagonistic_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
