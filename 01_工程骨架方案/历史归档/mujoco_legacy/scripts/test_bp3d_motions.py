"""Scripted cable-pull motions for the BodyParts3D arm model.

Every motion pulls the wires of one muscle group along their travel and
releases every other channel, which is how the real device is commanded: the
antagonist has to be let out or it acts as a spring holding the limb back.

Ball-joint angles are reported about one named anatomical axis, taken from the
joint quaternion, because the magnitude of a ball rotation on its own says
nothing about which shoulder motion happened.

Reports land in ``exports/bp3d/``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mujoco
import numpy as np

import bp3d_common as B

# name, agonist groups, joint, anatomical axis, target degrees, co-contraction.
#
# The co-contraction fraction is the part of travel applied to the antagonist
# group at the same time.  An arm is never driven by one muscle at a time: the
# rotator cuff has to hold the humeral head in the glenoid plane while the
# deltoid elevates it, otherwise the ball joint takes the shortest path and
# the motion turns into a twist.  Leaving it out is what a single-cable test
# measures, and it is not how the device is commanded.
MOTIONS = (
    ("elbow_flexion", ("elbow_flexors",), "elbow", "flexion", 120.0, {}),
    ("shoulder_abduction", ("deltoid_middle", "supraspinatus"), "glenohumeral", "abduction", 60.0,
     {"infraspinatus_teres_minor": 0.35, "subscapularis": 0.35}),
    ("shoulder_flexion", ("deltoid_anterior", "coracobrachialis"), "glenohumeral", "flexion", 90.0,
     {"infraspinatus_teres_minor": 0.35, "subscapularis": 0.35}),
    ("shoulder_external_rotation", ("infraspinatus_teres_minor",), "glenohumeral", "internal_rotation", -30.0,
     {"deltoid_middle": 0.2}),
    ("shoulder_internal_rotation", ("subscapularis",), "glenohumeral", "internal_rotation", 30.0,
     {"deltoid_middle": 0.2}),
    ("scapular_upward_rotation", ("trapezius_levator", "serratus_anterior"), "acromioclavicular", "upward_rotation", 25.0, {}),
    ("clavicle_elevation", ("trapezius_levator",), "sternoclavicular", "elevation", 15.0, {}),
)


class ArmModel:
    def __init__(self, path: Path, contract: Path) -> None:
        self.model = mujoco.MjModel.from_xml_path(str(path))
        self.contract = B.load_json(contract)
        self.joint_address: dict[str, int] = {}
        self.joint_type: dict[str, int] = {}
        for index in range(self.model.njnt):
            name = mujoco.mj_id2name(self.model, mujoco.mjtObj.mjOBJ_JOINT, index)
            if name and not name.startswith("slider_"):
                self.joint_address[name] = int(self.model.jnt_qposadr[index])
                self.joint_type[name] = int(self.model.jnt_type[index])
        self.actuators: dict[str, int] = {}
        for index in range(self.model.nu):
            name = mujoco.mj_id2name(self.model, mujoco.mjtObj.mjOBJ_ACTUATOR, index)
            if name and name.startswith("pull_"):
                self.actuators[name[len("pull_") :]] = index
        first = next(iter(self.actuators.values()))
        self.release = abs(float(self.model.actuator_ctrlrange[first][0]))
        self.travel = {
            name: float(self.model.actuator_ctrlrange[index][1])
            for name, index in self.actuators.items()
        }

    # -- kinematics --------------------------------------------------------

    def rotation_vector(self, data: mujoco.MjData, joint: str) -> np.ndarray:
        address = self.joint_address[joint]
        if self.joint_type[joint] == mujoco.mjtJoint.mjJNT_BALL:
            quaternion = np.array(data.qpos[address : address + 4], dtype=float)
            if quaternion[0] < 0.0:
                quaternion = -quaternion
            vector = quaternion[1:]
            norm = float(np.linalg.norm(vector))
            if norm < 1e-12:
                return np.zeros(3)
            return 2.0 * np.arctan2(norm, quaternion[0]) * vector / norm
        entry = self.contract["joints"][joint]
        axis = np.asarray(entry["axes"][entry["primary_axis"]], dtype=float)
        return axis * float(data.qpos[address])

    def angle(self, data: mujoco.MjData, joint: str, axis_name: str) -> float:
        """Rotation about one named anatomical axis, in degrees."""

        axis = np.asarray(
            self.contract["joints"][joint]["axes"][axis_name], dtype=float
        )
        return float(np.degrees(self.rotation_vector(data, joint) @ axis))

    def magnitude(self, data: mujoco.MjData, joint: str) -> float:
        return float(np.degrees(np.linalg.norm(self.rotation_vector(data, joint))))

    # -- simulation --------------------------------------------------------

    def settle(self, data: mujoco.MjData, steps: int = 3000) -> None:
        for _ in range(steps):
            B.apply_pull_only_tension(self.model, data)
            mujoco.mj_step(self.model, data)

    def fresh(self) -> mujoco.MjData:
        data = mujoco.MjData(self.model)
        data.qpos[:] = 0.0
        data.qvel[:] = 0.0
        data.ctrl[:] = 0.0
        self.settle(data)
        return data

    def drive(
        self, data: mujoco.MjData, pulls: dict[str, float], seconds: float
    ) -> None:
        steps = max(1, int(seconds / self.model.opt.timestep))
        for step in range(steps):
            fraction = (step + 1) / steps
            data.ctrl[:] = 0.0
            for name, index in self.actuators.items():
                data.ctrl[index] = (
                    pulls[name] * fraction if name in pulls else -self.release
                )
            B.apply_pull_only_tension(self.model, data)
            mujoco.mj_step(self.model, data)

    def pulls_for(self, groups) -> dict[str, float]:
        return {
            channel: self.travel[channel]
            for channel in self.actuators
            if any(channel.startswith(group) for group in groups)
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=B.MODEL_PATH)
    parser.add_argument("--contract", type=Path, default=B.CONTRACT_PATH)
    parser.add_argument("--seconds", type=float, default=4.0)
    parser.add_argument("--gravity", choices=("on", "off"), default="on")
    parser.add_argument(
        "--report", type=Path, default=B.REPORT_DIR / "motion_test.json"
    )
    args = parser.parse_args()

    arm = ArmModel(args.model, args.contract)
    if args.gravity == "off":
        arm.model.opt.gravity[:] = 0.0
    print(f"model:    {args.model}")
    print(f"gravity:  {args.gravity}")
    print(f"channels: {len(arm.actuators)}")

    rest = arm.fresh()
    rest_angles = {
        joint: round(arm.magnitude(rest, joint), 1)
        for joint in ("sternoclavicular", "acromioclavicular", "glenohumeral", "elbow")
    }
    print(f"rest pose: {rest_angles}")

    report: dict[str, object] = {
        "seconds": args.seconds,
        "gravity": args.gravity,
        "rest_deg": rest_angles,
        "motions": {},
    }

    print(f"\n{'motion':28s} {'measured':>9s} {'target':>8s}  verdict   detail")
    for name, groups, joint, axis_name, target, brace in MOTIONS:
        pulls = arm.pulls_for(groups)
        for brace_group, fraction in brace.items():
            for channel, travel in arm.pulls_for((brace_group,)).items():
                pulls[channel] = travel * fraction
        if not pulls:
            print(f"{name:28s} {'-':>9s} {target:8.1f}  NO CHANNEL")
            continue
        data = arm.fresh()
        arm.drive(data, pulls, args.seconds)
        measured = arm.angle(data, joint, axis_name)
        detail = {
            key: round(arm.magnitude(data, key), 1)
            for key in ("sternoclavicular", "acromioclavicular", "glenohumeral", "elbow")
        }
        passed = abs(measured) >= 0.8 * abs(target) and np.sign(measured) == np.sign(target)
        report["motions"][name] = {
            "pulled": list(pulls),
            "co_contraction": brace,
            "measured_deg": round(measured, 2),
            "target_deg": target,
            "joint_magnitudes_deg": detail,
            "pass": bool(passed),
        }
        print(
            f"{name:28s} {measured:9.1f} {target:8.1f}  "
            f"{'PASS' if passed else 'below'}   {detail}"
        )

    # Antagonist test: flex the elbow, then straighten it with the triceps.
    flex = arm.fresh()
    arm.drive(flex, arm.pulls_for(("elbow_flexors",)), args.seconds)
    flexed = arm.angle(flex, "elbow", "flexion")
    arm.drive(flex, arm.pulls_for(("triceps",)), args.seconds)
    extended = arm.angle(flex, "elbow", "flexion")
    extension_pass = extended <= 20.0
    report["motions"]["elbow_extension_return"] = {
        "flexed_deg": round(flexed, 2),
        "returned_deg": round(extended, 2),
        "pass": bool(extension_pass),
    }
    print(
        f"{'elbow_extension_return':28s} {extended:9.1f} {20.0:8.1f}  "
        f"{'PASS' if extension_pass else 'below'}   (flexed to {flexed:.1f})"
    )

    # Shoulder rhythm: how much did the scapula follow the humerus?
    elevate = arm.fresh()
    arm.drive(
        elevate,
        arm.pulls_for(
            (
                "deltoid_middle",
                "deltoid_anterior",
                "trapezius_levator",
                "serratus_anterior",
            )
        ),
        args.seconds,
    )
    humerus = arm.angle(elevate, "glenohumeral", "abduction")
    scapula = arm.angle(elevate, "acromioclavicular", "upward_rotation")
    ratio = scapula / humerus if abs(humerus) > 1.0 else float("nan")
    report["rhythm"] = {
        "glenohumeral_abduction_deg": round(humerus, 2),
        "scapular_upward_rotation_deg": round(scapula, 2),
        "ratio": None if np.isnan(ratio) else round(float(ratio), 3),
    }
    print(
        f"{'shoulder_rhythm':28s} humerus {humerus:6.1f}  scapula {scapula:6.1f}  "
        f"ratio {ratio:5.2f} (reference 0.5)"
    )

    B.save_json(args.report, report)
    print(f"\nreport: {args.report}")


if __name__ == "__main__":
    main()
