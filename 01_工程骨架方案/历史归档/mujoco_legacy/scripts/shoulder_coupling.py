"""Apply a simplified scapulohumeral rhythm to the arm model."""

import math

import mujoco
import numpy as np


def actuator_id(model: mujoco.MjModel, name: str) -> int:
    actuator = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, name)
    if actuator < 0:
        raise ValueError(f"Actuator not found: {name}")
    return int(actuator)


def apply_shoulder_rhythm(model: mujoco.MjModel, data: mujoco.MjData) -> None:
    """Couple clavicle and scapula motion to shoulder elevation.

    This is a first-pass rhythm model, not a clinical shoulder model:

    - positive shoulder flexion/abduction drives upward scapular rotation;
    - clavicle elevation and protraction follow with smaller gains;
    - depression/adduction produces the opposite residual motion.
    """

    flexion = float(data.ctrl[actuator_id(model, "shoulder_flexion_act")])
    abduction = float(data.ctrl[actuator_id(model, "shoulder_abduction_act")])
    _apply_rhythm_targets(model, data, flexion, abduction)


def joint_qpos(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    joint_name: str,
) -> float:
    joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, joint_name)
    if joint < 0:
        raise ValueError(f"Joint not found: {joint_name}")
    return float(data.qpos[int(model.jnt_qposadr[joint])])


def apply_shoulder_rhythm_from_state(
    model: mujoco.MjModel,
    data: mujoco.MjData,
) -> None:
    """Drive the temporary scapular rhythm from actual humeral joint angles."""

    flexion = joint_qpos(model, data, "shoulder_flexion")
    abduction = joint_qpos(model, data, "shoulder_abduction")
    _apply_rhythm_targets(model, data, flexion, abduction)


def _apply_rhythm_targets(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    flexion: float,
    abduction: float,
) -> None:
    """Set the temporary clavicle/scapula position actuator targets."""

    elevation = float(np.hypot(max(flexion, 0.0), max(abduction, 0.0)))
    depression = min(flexion, abduction, 0.0)

    targets = {
        "clavicle_protraction_act": 0.10 * elevation + 0.03 * depression,
        "clavicle_elevation_act": 0.20 * elevation + 0.05 * depression,
        "clavicle_axial_rotation_act": 0.10 * elevation + 0.03 * depression,
        "scapula_upward_rotation_act": 0.30 * elevation + 0.08 * depression,
        "scapula_posterior_tilt_act": 0.15 * elevation + 0.04 * depression,
        "scapula_external_rotation_act": 0.10 * elevation + 0.03 * depression,
    }

    for actuator_name, target in targets.items():
        actuator = actuator_id(model, actuator_name)
        lower, upper = model.actuator_ctrlrange[actuator]
        data.ctrl[actuator] = np.clip(target, lower, upper)


def shoulder_rhythm_example_degrees() -> tuple[float, float, float]:
    """Return a representative 120 degree humeral elevation example."""

    return (math.radians(120.0), math.radians(60.0), math.radians(30.0))
