"""Check that shoulder elevation produces coordinated scapular motion."""

from pathlib import Path

import mujoco
import numpy as np

from shoulder_coupling import actuator_id, apply_shoulder_rhythm


MODEL_PATH = Path(__file__).with_name("arm_reference.xml")


def body_id(model: mujoco.MjModel, name: str) -> int:
    body = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, name)
    if body < 0:
        raise ValueError(f"Body not found: {name}")
    return int(body)


def joint_value_deg(model: mujoco.MjModel, data: mujoco.MjData, name: str) -> float:
    joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
    if joint < 0:
        raise ValueError(f"Joint not found: {name}")
    return float(np.rad2deg(data.qpos[int(model.jnt_qposadr[joint])]))


def quat_angle(first: np.ndarray, second: np.ndarray) -> float:
    dot = np.clip(abs(float(np.dot(first, second))), -1.0, 1.0)
    return 2.0 * np.arccos(dot)


def main() -> None:
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
    data = mujoco.MjData(model)

    scapula = body_id(model, "scapula_external_rotation")
    mujoco.mj_forward(model, data)
    rest_quat = np.array(data.xquat[scapula])

    data.ctrl[actuator_id(model, "shoulder_flexion_act")] = np.deg2rad(120.0)
    data.ctrl[actuator_id(model, "shoulder_abduction_act")] = 0.0

    for _ in range(1000):
        apply_shoulder_rhythm(model, data)
        mujoco.mj_step(model, data)

    elevated_quat = np.array(data.xquat[scapula])
    angle = quat_angle(rest_quat, elevated_quat)

    print(f"shoulder flexion: 120.0 deg")
    print(f"scapula rotation: {np.rad2deg(angle):.2f} deg")
    print(f"scapula position: {np.round(data.xpos[scapula], 5)}")
    print(
        "coupled joints: "
        f"clavicle_elevation={joint_value_deg(model, data, 'clavicle_elevation'):.1f} deg, "
        f"scapula_upward_rotation={joint_value_deg(model, data, 'scapula_upward_rotation'):.1f} deg"
    )

    if angle < np.deg2rad(20.0):
        raise RuntimeError("The scapula did not rotate enough during elevation.")

    print("scapulohumeral rhythm: OK")


if __name__ == "__main__":
    main()
