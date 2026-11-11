"""Compare the same free elbow drive with gravity at fixed S1 orientations."""

from pathlib import Path

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, save
from simulate_elbow_antagonistic_load import HEADS, JOINTS, write_csv
from simulate_elbow_series_slider import SLIDE, run


def gravity_moment(model, data, gh, row, rotation):
    for joint, angle in zip(JOINTS, row[1:4], strict=True):
        data.joint(joint).qpos[:] = np.deg2rad(angle)
    data.joint(SLIDE).qpos[:] = row[4] * .001
    mujoco.mj_kinematics(model, data)
    offsets = data.xipos - gh
    forces = model.body_mass[:, None] * model.opt.gravity
    return rotation @ np.cross(offsets, forces).sum(axis=0)


def main():
    elbow = DEST.parent / "elbow_local_muscle_load"
    folder = DEST.parent / "integration_review/s1_series_elbow"
    folder.mkdir(exist_ok=True)
    source_file = elbow / "elbow_series_free_slider.json"
    source = read(source_file)
    model_path = ROOT / source["local_model"]["file"]
    if provenance(model_path) != source["local_model"]:
        raise RuntimeError("Local slider model changed")
    selection = read(ROOT / "design/S1_CURRENT_CANDIDATE.json")
    bare_path = ROOT / selection["bare_motion_report"]
    if provenance(bare_path)["sha256"] != selection["bare_motion_report_sha256"]:
        raise RuntimeError("S1 coordinate source changed")
    bare = read(bare_path)
    xml = model_path.read_text(encoding="utf-8")
    gh = np.asarray(bare["joint_centre_source_mm"]) * .001
    g = np.array([0., 0., -9.81])
    cases = []
    for axis, angle in (("flexion", 0), ("flexion", 30), ("flexion", 60), ("flexion", 90),
                        ("abduction", 30), ("abduction", 60), ("abduction", 90)):
        model = mujoco.MjModel.from_xml_string(xml)
        rotation = A.I.rotation_about(gh * 1000., bare["axes_source"][axis], angle)[:, :3]
        model.opt.gravity[:] = rotation.T @ g
        rows, result = run(model, .25, .15, 25., 0.)
        if angle == 0:
            reference = source["cases"][2]
            for key in ("elbow_peak_deg", "elbow_final_deg", "carriage_peak_mm"):
                np.testing.assert_allclose(result[key], reference[key], atol=1e-8)
        check_data = mujoco.MjData(model)
        check_data.ctrl[:] = [25., 0.]
        mujoco.mj_forward(model, check_data)
        elbow_dof = model.joint("elbow").dofadr.item()
        drive_torque = float(check_data.qfrc_actuator[elbow_dof])
        gravity_torque = -float(check_data.qfrc_bias[elbow_dof])
        peak = max(rows, key=lambda row: row[1])
        m0 = gravity_moment(model, check_data, gh, rows[0], rotation)
        m1 = gravity_moment(model, check_data, gh, peak, rotation)
        path = folder / f"held_{axis}_{angle:03d}.csv"
        write_csv(path, ["time_s", *(f"{joint}_deg" for joint in JOINTS),
                         "carriage_mm", *(f"{head}_span_change_mm" for head in HEADS),
                         "flexor_force_N", "extensor_force_N", "series_return_force_N"], rows)
        cases.append(dict(**result, axis=axis, held_GH_deg=angle, gravity_local_m_per_s2=model.opt.gravity.tolist(),
                          neutral_elbow_25N_actuator_torque_Nm=drive_torque,
                          neutral_elbow_gravity_torque_Nm=gravity_torque,
                          GH_gravity_moment_neutral_elbow_Nm=m0.tolist(),
                          GH_gravity_moment_peak_elbow_Nm=m1.tolist(), trace=provenance(path)))
        print(axis, angle, "peak elbow", round(result["elbow_peak_deg"], 3),
              "peak carriage", round(result["carriage_peak_mm"], 3),
              "GH gravity moment", round(float(np.linalg.norm(m1)), 3), flush=True)
    save(folder / "held_shoulder_elbow_load.json", dict(
        generator=provenance(Path(__file__).resolve()), elbow_response=provenance(source_file),
        shoulder_frame=provenance(bare_path), cases=cases,
        scope="Gravity re-expressed in the fixed humerus frame at each held GH orientation. "
        "Four elbow/slider degrees of freedom solve dynamics, GH/AC/SC do not move. GH moment is "
        "the weight contribution only, including provisional humerus/forearm/cap masses. No "
        "shoulder muscle force, dynamics, combined pose contact or physical motion acceptance."))


if __name__ == "__main__":
    main()
