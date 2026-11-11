"""Extract clear-trajectory endpoint travel and local tension from the saved cycles."""

import csv

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    case_file = folder / "gravity_balanced_active_cycles.json"
    report = load_reference(provenance(case_file))
    baseline = load_reference(report["baseline"])
    model_file = ROOT / baseline["model"]["file"]
    if provenance(model_file) != baseline["model"]:
        raise ValueError("Source MuJoCo model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    joint = model.joint("GH").dofadr.item()
    names = [model.tendon(i).name for i in range(model.ntendon)]
    result = {}
    for case, row in report["cases"].items():
        csv_path = ROOT / row["trace"]["file"]
        if provenance(csv_path) != row["trace"]:
            raise ValueError(f"Trace changed: {case}")
        with csv_path.open(newline="", encoding="utf-8") as stream:
            frames = list(csv.DictReader(stream))[:row["clear_frames"]]
        if not frames:
            continue
        lengths, tensions, torque = [], [], []
        for frame in frames:
            data.qpos[:] = [float(frame[f"qpos_{i}"]) for i in range(model.nq)]
            data.ctrl[:] = [float(frame[f"{model.actuator(i).name}_input_N"]) for i in range(model.nu)]
            mujoco.mj_forward(model, data)
            lengths.append(data.ten_length.copy() * 1000.)
            tensions.append([float(data.actuator_force[model.actuator(name).id]) for name in names])
            torque.append(np.linalg.norm(data.qfrc_actuator[joint:joint + 3]))
        travel = np.array(lengths) - lengths[0]
        tensions = np.array(tensions)
        if np.min(tensions) < -1e-8:
            raise ValueError(f"Unexpected pushing actuator: {case}")
        result[case] = dict(
            checked_clear_frames=len(frames), max_active_GH_moment_Nm=max(torque),
            endpoint_channels={name: dict(shortening_peak_mm=max(0., -travel[:, i].min()),
                                         extension_peak_mm=max(0., travel[:, i].max()),
                                         active_tension_peak_N=tensions[:, i].max())
                               for i, name in enumerate(names)})
    save(folder / "gravity_balanced_active_stroke.json", dict(
        generator=provenance(ROOT / "tools/summarize_s1_series_elbow_active_stroke.py"),
        active_cycles=provenance(case_file), cases=result,
        scope="Endpoint travel and ideal force-controlled tension through the independently screened"
              " pre-interference samples. Endpoint distance is not wire routing or motor stroke."))


if __name__ == "__main__":
    main()
