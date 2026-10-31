"""Open the V2 cable-design model with actuator sliders."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
from project_runtime import WorkspacePath as Path
import time

import mujoco
import mujoco.viewer

from v2_common import CONFIG_PATH, MODEL_PATH, load_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    args = parser.parse_args()
    config = load_json(args.config)
    model = _load_project_model(str(args.model))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    neutral = {}
    for route in config["routes"]:
        tendon = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_TENDON, f"cable_{route['id']}"
        )
        actuator = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{route['id']}"
        )
        length = float(data.ten_length[tendon])
        neutral[route["id"]] = length
        model.actuator_ctrlrange[actuator] = [
            max(0.0, length - 0.06),
            length + 0.01,
        ]
        data.ctrl[actuator] = length
    with mujoco.viewer.launch_passive(model, data) as viewer:
        while viewer.is_running():
            mujoco.mj_step(model, data)
            viewer.sync()
            time.sleep(model.opt.timestep)


if __name__ == "__main__":
    main()
