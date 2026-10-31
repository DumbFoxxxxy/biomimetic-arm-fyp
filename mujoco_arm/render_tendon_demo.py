"""Render representative poses of the tendon demo to PNG files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import struct

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_MODEL = ROOT / "arm_tendon_demo.xml"
DEFAULT_CONFIG = ROOT / "config" / "tendon_demo_16.json"
DEFAULT_OUT = ROOT / "exports" / "tendon_demo"

POSES = {
    "zero": {},
    "K1_E01": {"E01": 40.0},
    "K3_G03": {"G03": 40.0},
    "K6_S01": {"S01": 40.0},
    "P_both": {"P01": 30.0, "P02": 40.0},
}


def actuator_map(model: mujoco.MjModel) -> dict[str, int]:
    result = {}
    for index in range(model.nu):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, index) or ""
        if name.startswith("pull_"):
            result[name[5:]] = index
    return result


def write_bmp(path: Path, pixels: np.ndarray) -> None:
    """Write a bottom-up 24-bit BMP without an image-library dependency."""

    height, width, channels = pixels.shape
    if channels != 3:
        raise ValueError("Expected an RGB image")
    row_bytes = ((width * 3 + 3) // 4) * 4
    image_size = row_bytes * height
    file_size = 54 + image_size
    with path.open("wb") as handle:
        handle.write(b"BM")
        handle.write(struct.pack("<IHHI", file_size, 0, 0, 54))
        handle.write(
            struct.pack(
                "<IIIHHIIIIII",
                40,
                width,
                height,
                1,
                24,
                0,
                image_size,
                2835,
                2835,
                0,
                0,
            )
        )
        for row in pixels[::-1]:
            bgr = row[:, ::-1].astype(np.uint8).tobytes()
            handle.write(bgr)
            handle.write(b"\0" * (row_bytes - width * 3))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    model = mujoco.MjModel.from_xml_path(str(args.model))
    config = json.loads(args.config.read_text(encoding="utf-8"))
    neutral = {
        route["id"]: float(route["neutral_length_m"])
        for route in config["routes"]
    }
    actuators = actuator_map(model)
    renderer = mujoco.Renderer(model, height=480, width=640)
    option = mujoco.MjvOption()
    option.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = True
    option.flags[mujoco.mjtVisFlag.mjVIS_JOINT] = True
    camera = mujoco.MjvCamera()
    camera.type = mujoco.mjtCamera.mjCAMERA_FREE
    camera.lookat[:] = (0.10, 0.0, -0.02)
    camera.distance = 0.9
    camera.azimuth = 45.0
    camera.elevation = -12.0

    for name, pulls in POSES.items():
        data = mujoco.MjData(model)
        data.eq_active[:] = 0
        data.ctrl[:] = 0.0
        for cable_id, actuator in actuators.items():
            data.ctrl[actuator] = neutral[cable_id]
        for cable_id, pull_mm in pulls.items():
            data.ctrl[actuators[cable_id]] = neutral[cable_id] - pull_mm * 0.001
        for _ in range(int(5.0 / model.opt.timestep)):
            mujoco.mj_step(model, data)
        renderer.update_scene(data, camera=camera, scene_option=option)
        pixels = renderer.render()
        path = args.out / f"{name}.bmp"
        write_bmp(path, pixels)
        print(f"wrote {path}")

    renderer.close()


if __name__ == "__main__":
    main()
