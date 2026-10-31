"""Open the K6/S01 engineering baseline with a pull slider."""

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
import numpy as np

from k6_common import DEFAULT_CONFIG, DEFAULT_XML, load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_XML)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    config = load_config(args.config)
    model = _load_project_model(str(args.model))
    data = mujoco.MjData(model)
    tendon = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_TENDON, "cable_S01")
    actuator = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, "pull_S01")
    mujoco.mj_forward(model, data)
    neutral = float(data.ten_length[tendon])
    model.actuator_ctrlrange[actuator] = [
        neutral - float(config["target"]["max_pull_mm"]) * 0.001,
        neutral + 0.01,
    ]
    data.ctrl[actuator] = neutral

    with mujoco.viewer.launch_passive(model, data) as viewer:
        while viewer.is_running():
            mujoco.mj_step(model, data)
            viewer.sync()
            time.sleep(model.opt.timestep)


if __name__ == "__main__":
    main()
