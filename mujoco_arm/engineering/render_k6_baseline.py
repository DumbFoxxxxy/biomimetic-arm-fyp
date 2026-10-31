"""Render representative K6/S01 engineering-baseline poses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import mujoco
import numpy as np

from k6_audit import K6Evaluator
from k6_common import DEFAULT_CONFIG, DEFAULT_XML, MUJOCO_ROOT, load_config


sys.path.insert(0, str(MUJOCO_ROOT))
from render_tendon_demo import write_bmp  # noqa: E402


def find_pull_for_angle(
    evaluator: K6Evaluator,
    angle_deg: float,
    max_pull_mm: float,
) -> tuple[float, float]:
    previous_pull = 0.0
    previous_angle = evaluator.measure(0.0)["up_deg"]
    best_pull = 0.0
    best_angle = previous_angle
    for pull in np.arange(1.0, max_pull_mm + 1e-9, 1.0):
        sample = evaluator.measure(float(pull))
        if abs(sample["up_deg"] - angle_deg) < abs(best_angle - angle_deg):
            best_pull = float(pull)
            best_angle = float(sample["up_deg"])
        if previous_angle <= angle_deg <= sample["up_deg"]:
            fraction = (angle_deg - previous_angle) / (
                sample["up_deg"] - previous_angle
            )
            best_pull = previous_pull + fraction
            best_angle = angle_deg
            break
        previous_pull = float(pull)
        previous_angle = float(sample["up_deg"])
    return best_pull, best_angle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_XML)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--out",
        type=Path,
        default=MUJOCO_ROOT / "exports" / "k6_baseline" / "poses",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    config = load_config(args.config)
    evaluator = K6Evaluator(args.model, config)

    requested_angles = (0.0, 12.0, 24.0, 36.0)
    poses = []
    for target_angle in requested_angles:
        pull, angle = find_pull_for_angle(
            evaluator,
            target_angle,
            float(config["target"]["max_pull_mm"]),
        )
        poses.append(
            {
                "requested_angle_deg": target_angle,
                "rendered_angle_deg": angle,
                "pull_mm": pull,
            }
        )

    renderer = mujoco.Renderer(evaluator.model, height=480, width=640)
    option = mujoco.MjvOption()
    option.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = True
    option.flags[mujoco.mjtVisFlag.mjVIS_JOINT] = True
    camera = mujoco.MjvCamera()
    camera.type = mujoco.mjtCamera.mjCAMERA_FREE
    camera.lookat[:] = (0.11, -0.06, 0.20)
    camera.distance = 0.82
    camera.azimuth = 138.0
    camera.elevation = -9.0

    for index, pose in enumerate(poses):
        data = evaluator.simulate(pose["pull_mm"])
        renderer.update_scene(data, camera=camera, scene_option=option)
        pixels = renderer.render()
        path = args.out / f"{index:02d}_target_{pose['requested_angle_deg']:.0f}deg.bmp"
        write_bmp(path, pixels)
        print(
            f"wrote {path}  pull={pose['pull_mm']:.2f} mm "
            f"angle={pose['rendered_angle_deg']:.2f} deg"
        )
    renderer.close()
    (args.out / "poses.json").write_text(
        json.dumps(poses, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
