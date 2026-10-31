"""Interactive 16-tendon demo.

Each actuator slider is a target cable length in metres.  Dragging a slider to
a smaller value pulls that cable.  Number keys load basic presets.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import mujoco
import mujoco.viewer
import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_MODEL = ROOT / "arm_tendon_demo.xml"
DEFAULT_CONFIG = ROOT / "config" / "tendon_demo_16.json"


PRESETS_MM = {
    "1": {"E01": 40.0},
    "2": {"G01": 30.0},
    "3": {"G03": 40.0},
    "4": {"G05": 30.0},
    "5": {"G06": 30.0},
    "6": {"S01": 40.0},
    "7": {"G01": 30.0, "S01": 40.0},
    "8": {"S03": 30.0},
    "9": {"S04": 30.0},
    "a": {"S06": 40.0},
    "b": {"S05": 30.0},
    "c": {"S02": 40.0},
    "d": {"P01": 30.0, "P02": 40.0},
}


def actuator_names(model: mujoco.MjModel) -> list[str]:
    return [
        mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, index) or ""
        for index in range(model.nu)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--no-presets", action="store_true")
    args = parser.parse_args()

    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    neutral = {
        route["id"]: float(route["neutral_length_m"])
        for route in config["routes"]
    }
    names = actuator_names(model)
    actuator_by_id = {
        name.removeprefix("pull_"): index
        for index, name in enumerate(names)
        if name.startswith("pull_")
    }
    data.eq_active[:] = 0
    data.ctrl[:] = [
        neutral[name.removeprefix("pull_")]
        if name.startswith("pull_")
        else 0.0
        for name in names
    ]
    mujoco.mj_forward(model, data)

    def reset() -> None:
        data.qpos[:] = 0.0
        data.qvel[:] = 0.0
        data.ctrl[:] = [
            neutral[name.removeprefix("pull_")]
            if name.startswith("pull_")
            else 0.0
            for name in names
        ]

    def apply_preset(pull_mm: dict[str, float]) -> None:
        reset()
        for cable_id, amount_mm in pull_mm.items():
            actuator = actuator_by_id[cable_id]
            data.ctrl[actuator] = neutral[cable_id] - amount_mm * 0.001

    def key_callback(keycode: int) -> None:
        try:
            key = chr(keycode)
        except ValueError:
            return
        if key == "0":
            reset()
            print("preset: zero")
        elif not args.no_presets and key in PRESETS_MM:
            apply_preset(PRESETS_MM[key])
            print(f"preset: {key} -> {PRESETS_MM[key]}")

    print("Tendon demo controls")
    print("  0 : zero pose")
    for key, pulls in PRESETS_MM.items():
        print(f"  {key} : {pulls}")
    print("Right panel: drag pull_<cable> sliders to change cable length.")

    with mujoco.viewer.launch_passive(
        model,
        data,
        key_callback=key_callback,
    ) as viewer:
        while viewer.is_running():
            mujoco.mj_step(model, data)
            viewer.sync()
            time.sleep(model.opt.timestep)


if __name__ == "__main__":
    main()
