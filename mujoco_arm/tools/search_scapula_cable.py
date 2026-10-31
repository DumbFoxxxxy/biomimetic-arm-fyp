"""Search alternate anchors for the scapular downward cable S02."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import copy
import json
from project_runtime import WorkspacePath as Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "arm_tendon_demo.xml"
CONFIG = ROOT / "config" / "tendon_demo_16.json"

AC_MM = np.array([156.86, 233.24, -11.97])
P_MM = np.array([62.18, 210.10, 91.39])
SUP_MM = np.array([72.76, 268.45, 82.21])
INF_MM = np.array([44.55, 112.87, 106.69])


def mj(point_mm):
    return np.array([point_mm[0], -point_mm[2], point_mm[1]]) * 0.001


def lp(u, v, w=10.0):
    U = SUP_MM - INF_MM
    U /= np.linalg.norm(U)
    V0 = np.array([159.06, 187.54, 7.30]) - P_MM
    V = V0 - np.dot(V0, U) * U
    V /= np.linalg.norm(V)
    W = np.array([0.6556, 0.0, 0.7552])
    return P_MM + u * U + v * V + w * W


def th(x, y):
    z = -210.69 + np.sqrt(400.0**2 - (x + 200.01) ** 2)
    return np.array([x, y, z])


def set_site(tree, name, position):
    for site in tree.iter("site"):
        if site.get("name") == name:
            site.set("pos", " ".join(f"{value:.9f}" for value in position))
            return
    raise KeyError(name)


def run_model(model, neutral, pull_mm):
    data = mujoco.MjData(model)
    data.eq_active[:] = 0
    data.ctrl[:] = 0.0
    for cable_id, value in neutral.items():
        actuator = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{cable_id}"
        )
        data.ctrl[actuator] = value
    cable_id = "S02"
    actuator = mujoco.mj_name2id(
        model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{cable_id}"
    )
    data.ctrl[actuator] = neutral[cable_id] - pull_mm * 0.001
    for _ in range(3000):
        mujoco.mj_step(model, data)
    values = {}
    for name in ("sc_elevation", "sc_protraction", "ac_ext"):
        joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
        values[name] = float(np.degrees(data.qpos[model.jnt_qposadr[joint]]))
    return values


def main():
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    neutral = {
        route["id"]: float(route["neutral_length_m"])
        for route in config["routes"]
    }
    base_tree = ET.parse(MODEL)
    base_model = _load_project_model(str(MODEL))
    base_data = mujoco.MjData(base_model)
    base_data.eq_active[:] = 0
    base_data.ctrl[:] = 0.0
    for cable_id, value in neutral.items():
        actuator = mujoco.mj_name2id(
            base_model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{cable_id}"
        )
        base_data.ctrl[actuator] = value
    for _ in range(3000):
        mujoco.mj_step(base_model, base_data)
    base = {}
    for name in ("sc_elevation", "sc_protraction", "ac_ext"):
        joint = mujoco.mj_name2id(base_model, mujoco.mjtObj.mjOBJ_JOINT, name)
        base[name] = float(
            np.degrees(base_data.qpos[base_model.jnt_qposadr[joint]])
        )

    a_candidates = {
        "T55_160": th(55, 160),
        "T95_160": th(95, 160),
        "T120_160": th(120, 160),
        "T55_210": th(55, 210),
        "T95_210": th(95, 210),
        "T120_210": th(120, 210),
        "T55_255": th(55, 255),
        "T95_255": th(95, 255),
        "T120_255": th(120, 255),
    }
    b_candidates = {
        "SUP": SUP_MM,
        "INF": INF_MM,
        "LAT_UP": lp(40, 80),
        "LAT_MID": lp(0, 60),
        "LAT_LOW": lp(-40, 80),
        "MED_MID": lp(0, 0),
    }
    rows = []
    for a_name, a_point in a_candidates.items():
        for b_name, b_point in b_candidates.items():
            tree = copy.deepcopy(base_tree)
            set_site(tree, "A_S02", mj(a_point))
            set_site(tree, "B_S02", mj(b_point) - mj(AC_MM))
            xml = ET.tostring(tree.getroot(), encoding="unicode")
            try:
                model = _load_project_xml(xml)
            except Exception:
                continue
            values = run_model(model, neutral, 40.0)
            delta = {key: values[key] - base[key] for key in base}
            score = -delta["sc_elevation"] - 0.2 * abs(delta["sc_protraction"])
            rows.append((score, a_name, b_name, delta))

    for score, a_name, b_name, delta in sorted(rows, reverse=True)[:20]:
        print(
            f"{a_name:9} -> {b_name:8} "
            f"delev={delta['sc_elevation']:7.2f} "
            f"dprot={delta['sc_protraction']:7.2f} "
            f"dext={delta['ac_ext']:7.2f}"
        )


if __name__ == "__main__":
    main()
