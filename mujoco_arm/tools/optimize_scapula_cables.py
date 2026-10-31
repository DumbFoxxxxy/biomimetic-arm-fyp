"""Find coordinated scapular cable pulls for the actual SC/AC chain."""

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
from scipy.optimize import differential_evolution

from analyze_scapula_cables import (
    AC_MM,
    CONFIG,
    E_SLIDE_MM,
    MODEL,
    P_MM,
    U_MM,
    W_MM,
    mj,
    rotvec,
)


def prepare(model):
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
    return neutral, base


def run(model, neutral, pulls_mm, seconds=3.0):
    data = mujoco.MjData(model)
    data.eq_active[:] = 0
    data.ctrl[:] = 0.0
    for cable_id, value in neutral.items():
        actuator = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{cable_id}"
        )
        data.ctrl[actuator] = value
    for index, amount_mm in enumerate(pulls_mm, start=1):
        actuator = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_S{index:02d}"
        )
        data.ctrl[actuator] = neutral[f"S{index:02d}"] - amount_mm * 0.001
    for _ in range(int(seconds / model.opt.timestep)):
        mujoco.mj_step(model, data)
    return data


def task_pose(model, baseline, data):
    scapula = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "scapula")
    r0 = baseline.xmat[scapula].reshape(3, 3)
    p0 = baseline.xpos[scapula] + r0 @ (mj(P_MM) - mj(AC_MM))
    r = data.xmat[scapula].reshape(3, 3)
    p = data.xpos[scapula] + r @ (mj(P_MM) - mj(AC_MM))
    axis, angle = rotvec(r @ r0.T)
    w = mj(W_MM) / np.linalg.norm(mj(W_MM))
    u = mj(U_MM) / np.linalg.norm(mj(U_MM))
    s = mj(E_SLIDE_MM) / np.linalg.norm(mj(E_SLIDE_MM))
    penetration = min((contact.dist for contact in data.contact), default=0.0)
    return {
        "slide_mm": float(np.dot(p - p0, s) * 1000.0),
        "up_deg": float(np.degrees(angle * np.dot(axis, w))),
        "ext_deg": float(np.degrees(angle * np.dot(axis, u))),
        "penetration_mm": float(penetration * 1000.0),
    }


TARGETS = {
    "up": ("up_deg", 1.0),
    "down": ("up_deg", -1.0),
    "slide_forward": ("slide_mm", 1.0),
    "slide_back": ("slide_mm", -1.0),
    "external": ("ext_deg", 1.0),
    "internal": ("ext_deg", -1.0),
}


def main():
    model = _load_project_model(str(MODEL))
    neutral, baseline = prepare(model)
    results = {}
    for name, (key, sign) in TARGETS.items():
        history = []

        def objective(values):
            data = run(model, neutral, values)
            if not np.isfinite(data.qpos).all():
                return 1e6
            pose = task_pose(model, baseline, data)
            primary = sign * pose[key]
            other_a, other_b = [k for k in ("slide_mm", "up_deg", "ext_deg") if k != key]
            coupling = 0.08 * abs(pose[other_a]) + 0.08 * abs(pose[other_b])
            penalty = 20.0 * max(0.0, -pose["penetration_mm"] - 0.5)
            penalty += 0.02 * float(np.sum(values))
            value = -primary + coupling + penalty
            history.append((value, pose, values))
            return value

        bounds = [(0.0, 40.0)] * 6
        result = differential_evolution(
            objective,
            bounds,
            maxiter=4,
            popsize=4,
            tol=0.5,
            seed=7,
            polish=False,
        )
        best_value = result.fun
        best_entry = min(history, key=lambda item: item[0])
        results[name] = {
            "score": float(best_value),
            "pull_mm": {
                f"S{index:02d}": float(value)
                for index, value in enumerate(result.x, start=1)
            },
            "pose": best_entry[1],
        }
        print(
            name,
            results[name]["pull_mm"],
            results[name]["pose"],
        )

    output = MODEL.parent / "config" / "scapula_cable_presets.json"
    output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
