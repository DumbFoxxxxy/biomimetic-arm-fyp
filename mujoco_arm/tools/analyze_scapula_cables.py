"""Report the task-space effect of each scapular cable."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import json
from project_runtime import WorkspacePath as Path

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "arm_tendon_demo.xml"
CONFIG = ROOT / "config" / "tendon_demo_16.json"

AC_MM = np.array([156.86, 233.24, -11.97])
P_MM = np.array([62.18, 210.10, 91.39])
W_MM = np.array([0.6556, 0.0, 0.7552])
U_MM = np.array([72.76 - 44.55, 268.45 - 112.87, 82.21 - 106.69])
U_MM /= np.linalg.norm(U_MM)
E_SLIDE_MM = np.array([0.7552, 0.0, -0.6556])
E_SLIDE_MM /= np.linalg.norm(E_SLIDE_MM)


def mj(point_mm):
    return np.array([point_mm[0], -point_mm[2], point_mm[1]]) * 0.001


def rotvec(matrix):
    angle = np.arccos(np.clip((np.trace(matrix) - 1.0) / 2.0, -1.0, 1.0))
    if angle < 1e-8:
        return np.zeros(3), 0.0
    axis = np.array(
        [
            matrix[2, 1] - matrix[1, 2],
            matrix[0, 2] - matrix[2, 0],
            matrix[1, 0] - matrix[0, 1],
        ]
    ) / (2.0 * np.sin(angle))
    return axis, angle


def simulate(model, neutral_lengths, pulls, seconds=4.0):
    out = mujoco.MjData(model)
    out.ctrl[:] = 0.0
    for cable_id, value in neutral_lengths.items():
        actuator = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{cable_id}"
        )
        out.ctrl[actuator] = value
    for cable_id, pull_mm in pulls.items():
        actuator = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{cable_id}"
        )
        out.ctrl[actuator] = neutral_lengths[cable_id] - pull_mm * 0.001
    out.eq_active[:] = 0
    for _ in range(int(seconds / model.opt.timestep)):
        mujoco.mj_step(model, out)
    return out


def main():
    model = _load_project_model(str(MODEL))
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    neutral = {
        route["id"]: float(route["neutral_length_m"])
        for route in config["routes"]
    }
    base = mujoco.MjData(model)
    base.eq_active[:] = 0
    base.ctrl[:] = 0.0
    for cable_id, value in neutral.items():
        actuator = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{cable_id}"
        )
        base.ctrl[actuator] = value
    for _ in range(int(4.0 / model.opt.timestep)):
        mujoco.mj_step(model, base)

    scapula = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "scapula")
    r0 = base.xmat[scapula].reshape(3, 3).copy()
    p0 = base.xpos[scapula] + r0 @ (mj(P_MM) - mj(AC_MM))
    w = mj(W_MM) / np.linalg.norm(mj(W_MM))
    u = mj(U_MM) / np.linalg.norm(mj(U_MM))
    slide = mj(E_SLIDE_MM) / np.linalg.norm(mj(E_SLIDE_MM))

    print("cable   slide_mm  up_deg  ext_deg  total_deg")
    for index in range(1, 7):
        cable_id = f"S{index:02d}"
        result = simulate(model, neutral, {cable_id: 30.0})
        r = result.xmat[scapula].reshape(3, 3)
        p = result.xpos[scapula] + r @ (mj(P_MM) - mj(AC_MM))
        axis, angle = rotvec(r @ r0.T)
        print(
            f"{cable_id}  {np.dot(p - p0, slide) * 1000:8.2f}"
            f"  {np.degrees(angle * np.dot(axis, w)):7.2f}"
            f"  {np.degrees(angle * np.dot(axis, u)):7.2f}"
            f"  {np.degrees(angle):9.2f}"
        )


if __name__ == "__main__":
    main()
