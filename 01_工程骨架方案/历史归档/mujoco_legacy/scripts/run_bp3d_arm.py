"""Open the BodyParts3D arm model with live cable sliders.

Every slider is one wire: ``pull_<wire>`` in metres.  Positive is a pull,
negative releases the wire so it cannot act as a spring against the others.
The tendons are one-sided, so a released wire carries nothing at all.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import time

import mujoco
import mujoco.viewer

import bp3d_common as B


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=B.MODEL_PATH)
    parser.add_argument("--release", action="store_true",
                        help="start with every wire released instead of taut")
    args = parser.parse_args()

    if not args.model.exists():
        raise FileNotFoundError(
            f"Missing {args.model}. Run build_bp3d_arm_model.py first."
        )

    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    data.qpos[:] = 0.0
    data.ctrl[:] = 0.0
    if args.release:
        for index in range(model.nu):
            name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, index)
            if name and name.startswith("pull_"):
                data.ctrl[index] = float(model.actuator_ctrlrange[index][0])

    with mujoco.viewer.launch_passive(model, data) as viewer:
        while viewer.is_running():
            B.apply_pull_only_tension(model, data)
            mujoco.mj_step(model, data)
            viewer.sync()
            time.sleep(model.opt.timestep)


if __name__ == "__main__":
    main()
