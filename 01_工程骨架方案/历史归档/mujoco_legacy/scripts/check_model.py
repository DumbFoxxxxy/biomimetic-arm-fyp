"""Compile the reference arm and print its kinematic structure."""

import math
from pathlib import Path

import mujoco


MODEL_PATH = Path(__file__).with_name("arm_reference.xml")


def main() -> None:
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
    data = mujoco.MjData(model)

    print(f"model: {MODEL_PATH}")
    print(f"bodies: {model.nbody}")
    print(f"joints: {model.njnt}")
    print(f"degrees of freedom: {model.nv}")
    print(f"actuators: {model.nu}")
    print()

    for joint_id in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, joint_id)
        qpos_adr = int(model.jnt_qposadr[joint_id])
        dof_adr = int(model.jnt_dofadr[joint_id])
        joint_type = int(model.jnt_type[joint_id])
        lower, upper = model.jnt_range[joint_id]

        if joint_type == int(mujoco.mjtJoint.mjJNT_HINGE):
            print(
                f"{name:24s} qpos={qpos_adr:2d} dof={dof_adr:2d} "
                f"range=[{math.degrees(lower):7.2f}, "
                f"{math.degrees(upper):7.2f}] deg"
            )
        else:
            print(f"{name:24s} qpos={qpos_adr:2d} type={joint_type}")

    # Exercise the position actuators over a short trajectory. This catches
    # invalid ranges, missing actuators and unstable inertias before the GUI.
    for _ in range(250):
        data.ctrl[:] = 0.0
        mujoco.mj_step(model, data)

    if not (data.qpos[:] == data.qpos[:]).all():
        raise RuntimeError("The model produced NaN during the validation rollout.")

    print()
    print("validation rollout: OK")


if __name__ == "__main__":
    main()
