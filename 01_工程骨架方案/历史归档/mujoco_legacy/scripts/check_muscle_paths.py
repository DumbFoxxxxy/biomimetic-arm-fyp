"""Check that each enabled cable has a usable joint moment arm."""

import argparse
from pathlib import Path
import json

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "arm_psm_generated.xml"
MUSCLE_CONFIG = ROOT / "config" / "muscle_manifest.json"


def tendon_length(model, data, tendon_id: int) -> float:
    mujoco.mj_forward(model, data)
    return float(data.ten_length[tendon_id])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument(
        "--muscle-config", type=Path, default=MUSCLE_CONFIG
    )
    args = parser.parse_args()

    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    manifest = json.loads(
        args.muscle_config.read_text(encoding="utf-8")
    )

    enabled = [
        muscle
        for muscle in manifest["muscles"]
        if muscle.get("enabled", False)
    ]
    joint_addresses = [
        (joint_id, int(model.jnt_qposadr[joint_id]))
        for joint_id in range(model.njnt)
        if int(model.jnt_type[joint_id])
        == int(mujoco.mjtJoint.mjJNT_HINGE)
    ]

    failures = []
    epsilon = 1e-4
    for muscle in enabled:
        for cable_index in range(int(muscle["cables"])):
            tendon_name = f"cable_{muscle['id']}_{cable_index}"
            tendon_id = mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_TENDON, tendon_name
            )
            if tendon_id < 0:
                failures.append(f"{tendon_name}: missing")
                continue

            data.qpos[:] = 0.0
            base_length = tendon_length(model, data, tendon_id)
            sensitivities = []
            for joint_id, qpos_adr in joint_addresses:
                data.qpos[:] = 0.0
                data.qpos[qpos_adr] = epsilon
                shifted_length = tendon_length(model, data, tendon_id)
                sensitivities.append(
                    (
                        abs((shifted_length - base_length) / epsilon),
                        joint_id,
                        (shifted_length - base_length) / epsilon,
                    )
                )
            sensitivity, joint_id, signed_sensitivity = max(sensitivities)
            joint_name = mujoco.mj_id2name(
                model, mujoco.mjtObj.mjOBJ_JOINT, joint_id
            )
            print(
                f"{tendon_name:36s} "
                f"dL/dq={signed_sensitivity: .6f} m/rad via {joint_name}"
            )
            if not np.isfinite(sensitivity) or sensitivity < 1e-5:
                failures.append(f"{tendon_name}: no usable moment arm")

    data.qpos[:] = 0.0
    mujoco.mj_forward(model, data)
    if failures:
        raise RuntimeError("; ".join(failures))

    print("muscle path moment-arm check: OK")


if __name__ == "__main__":
    main()
