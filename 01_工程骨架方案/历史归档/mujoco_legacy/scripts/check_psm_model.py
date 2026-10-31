"""Compile and exercise the generated flexible-muscle model."""

import argparse
from pathlib import Path

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "arm_psm_generated.xml"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    args = parser.parse_args()
    model_path = args.model

    if not model_path.exists():
        raise FileNotFoundError(
            "Run build_arm_model.py before check_psm_model.py"
        )

    model = mujoco.MjModel.from_xml_path(str(model_path))
    data = mujoco.MjData(model)

    required = (
        "cable_supraspinatus_0",
        "cable_infraspinatus_0",
        "cable_deltoid_anterior_0",
        "cable_deltoid_middle_0",
        "cable_deltoid_posterior_0",
        "cable_elbow_flexors_0",
        "cable_triceps_0",
    )
    for tendon_name in required:
        if mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_TENDON, tendon_name
        ) < 0:
            raise RuntimeError(f"Required tendon missing: {tendon_name}")

    for _ in range(500):
        mujoco.mj_step(model, data)

    if not np.isfinite(data.qpos).all() or not np.isfinite(data.qvel).all():
        raise RuntimeError("The generated model produced non-finite state.")

    print(f"model: {model_path}")
    print(f"bodies: {model.nbody}")
    print(f"joints: {model.njnt}")
    print(f"flex bodies: {model.nflex}")
    print(f"tendons: {model.ntendon}")
    print(f"actuators: {model.nu}")
    print(f"DOF: {model.nv}")
    print("zero-load stability: OK")


if __name__ == "__main__":
    main()
