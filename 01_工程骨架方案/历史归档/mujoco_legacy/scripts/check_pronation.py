"""Verify radial-head contact and distal-radius motion during pronation."""

import math
from pathlib import Path

import mujoco
import numpy as np


MODEL_PATH = Path(__file__).with_name("arm_reference.xml")

# BodyParts3D coordinates in metres. The proximal contact point lies on
# the pronation axis, while the distal-radius probe shows the motion.
RADIOCAPITELLAR_CONTACT = np.array([-0.233684, -0.071294, 1.036219])
RADIAL_HEAD_CENTRE = np.array([-0.23225, -0.07689, 1.03364])
RADIUS_DISTAL_PROBE = np.array([-0.27433, -0.11687, 0.80670])


def joint_address(model: mujoco.MjModel, joint_name: str) -> int:
    joint_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, joint_name)
    if joint_id < 0:
        raise ValueError(f"Joint not found: {joint_name}")
    return int(model.jnt_qposadr[joint_id])


def body_id(model: mujoco.MjModel, body_name: str) -> int:
    body_id_value = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, body_name)
    if body_id_value < 0:
        raise ValueError(f"Body not found: {body_name}")
    return int(body_id_value)


def main() -> None:
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
    data = mujoco.MjData(model)

    pronation_qpos = joint_address(model, "forearm_pronation")
    ulna_body = body_id(model, "forearm_ulna")
    radius_body = body_id(model, "radius")

    samples = (-90.0, -45.0, 0.0, 45.0, 90.0)
    ulna_quats = []
    radius_points = []

    print(
        "pronation_deg  radiocapitellar_contact_m"
        "                 radius_distal_position_m"
    )
    radial_heads = []
    for angle_deg in samples:
        data.qpos[:] = 0.0
        data.qpos[pronation_qpos] = math.radians(angle_deg)
        mujoco.mj_forward(model, data)

        ulna_quats.append(np.array(data.xquat[ulna_body]))
        radial_head = data.xpos[radius_body].copy()
        radius_rotation = data.xmat[radius_body].reshape(3, 3)
        radial_head_centre = (
            data.xpos[radius_body]
            + radius_rotation @ (RADIAL_HEAD_CENTRE - RADIOCAPITELLAR_CONTACT)
        )
        radius_distal = (
            data.xpos[radius_body]
            + radius_rotation @ (RADIUS_DISTAL_PROBE - RADIOCAPITELLAR_CONTACT)
        )
        radial_heads.append(radial_head)
        radius_points.append(radius_distal)
        print(
            f"{angle_deg:11.1f}  "
            f"[{radial_head[0]: .4f}, {radial_head[1]: .4f}, "
            f"{radial_head[2]: .4f}]  "
            f"[{radius_distal[0]: .4f}, {radius_distal[1]: .4f}, "
            f"{radius_distal[2]: .4f}]"
        )
        print(
            f"             radial_head_centre="
            f"[{radial_head_centre[0]: .4f}, {radial_head_centre[1]: .4f}, "
            f"{radial_head_centre[2]: .4f}]"
        )

    ulna_fixed = all(
        np.allclose(ulna_quats[0], quat, atol=1e-10)
        for quat in ulna_quats[1:]
    )
    radius_moved = not np.allclose(radius_points[0], radius_points[-1], atol=1e-6)
    radial_head_fixed = all(
        np.allclose(radial_heads[0], point, atol=1e-10)
        for point in radial_heads[1:]
    )

    if not ulna_fixed:
        raise RuntimeError("The ulna orientation changed during pronation.")
    if not radius_moved:
        raise RuntimeError("The radius did not move during pronation.")
    if not radial_head_fixed:
        raise RuntimeError("The radial head left the capitulum reference point.")

    print()
    print("ulna fixed during pronation: OK")
    print("radial head fixed at the capitulum reference point: OK")
    print("radius swings around the pronation axis: OK")


if __name__ == "__main__":
    main()
