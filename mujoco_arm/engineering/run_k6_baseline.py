"""Open the K6/S01 engineering baseline with a pull slider."""

from __future__ import annotations

import argparse
from pathlib import Path
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
    model = mujoco.MjModel.from_xml_path(str(args.model))
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
