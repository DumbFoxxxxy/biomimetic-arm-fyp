"""Open the V2 cable-design model with actuator sliders."""

from __future__ import annotations

import argparse
from pathlib import Path
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
    model = mujoco.MjModel.from_xml_path(str(args.model))
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
