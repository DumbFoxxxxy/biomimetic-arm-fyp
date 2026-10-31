"""Self-test for the paired-frame pose chain.

The test invents a target pose for every bone, derives the two coordinate
systems of every joint from that pose, writes them as a ``part_frames.txt``
table, and then checks that the exporter rebuilds the target pose exactly.
It covers the two naming styles the SolidWorks export can produce:

* the same frame name in both parts of a joint (``SC_J`` in thorax and
  clavicle)
* a suffixed parent frame (``SC_J`` in the clavicle, ``SC_J_thorax`` in the
  thorax)

Run it with the MuJoCo environment python:

    D:\\MujocoEnv\\Scripts\\python.exe D:\\XJTLU\\PSP\\mujoco_arm\\test_paired_frames.py
"""

from __future__ import annotations

from pathlib import Path
import tempfile

import numpy as np

from cable_drive_common import BONE_NAMES, JOINT_CHAIN
from prepare_cable_drive_exports import (
    build_bone_transforms,
    chain_joint_matrices,
    load_part_frames,
    resolve_frame,
)


def random_rigid(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(3, 3))
    rotation, upper = np.linalg.qr(matrix)
    rotation = rotation @ np.diag(np.sign(np.diag(upper)))
    if np.linalg.det(rotation) < 0.0:
        rotation[:, 0] *= -1.0
    transform = np.eye(4)
    transform[:3, :3] = rotation
    transform[:3, 3] = rng.uniform(-0.5, 0.5, 3)
    return transform


def write_part_frames(
    path: Path,
    frames: dict[str, dict[str, np.ndarray]],
) -> None:
    lines = [
        "name|component|v0|v1|v2|v3|v4|v5|v6|v7|v8|v9|v10|v11|v12|v13|v14|v15"
    ]
    for name, owners in frames.items():
        for bone, transform in owners.items():
            values = list(transform[:3, :3].reshape(-1)) + list(transform[:3, 3])
            payload = values + [1.0, 0.0, 0.0, 0.0]
            lines.append(
                name
                + "|"
                + f"{bone}-1"
                + "|"
                + "|".join(repr(float(value)) for value in payload)
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-16")


def build_case(
    rng: np.random.Generator,
    suffix_parent_frames: bool,
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, np.ndarray]]:
    targets = {bone: random_rigid(rng) for bone in BONE_NAMES}
    targets["thorax"] = np.eye(4)

    frames: dict[str, dict[str, np.ndarray]] = {}
    for step in JOINT_CHAIN:
        joint, parent, child = step["joint"], step["parent"], step["child"]
        joint_pose = random_rigid(rng)
        frames.setdefault(joint, {})[child] = (
            np.linalg.inv(targets[child]) @ joint_pose
        )
        parent_name = joint + f"_{parent}" if suffix_parent_frames else joint
        frames.setdefault(parent_name, {})[parent] = (
            np.linalg.inv(targets[parent]) @ joint_pose
        )
    return frames, targets


def check_case(label: str, suffix_parent_frames: bool) -> None:
    rng = np.random.default_rng(7)
    frames, targets = build_case(rng, suffix_parent_frames)

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "part_frames.txt"
        write_part_frames(path, frames)
        table = load_part_frames(path)

    usage_frames = {}
    for step in JOINT_CHAIN:
        joint, parent, child = step["joint"], step["parent"], step["child"]
        child_name, child_frame = resolve_frame(table, child, joint)
        parent_name, parent_frame = resolve_frame(table, parent, joint)
        assert child_frame is not None, f"{label}: missing child frame {joint}"
        assert parent_frame is not None, f"{label}: missing parent frame {joint}"
        usage_frames[joint] = (parent_name, child_name)

    world, usage = build_bone_transforms({}, table)
    assert all(item["mode"] == "paired_frames" for item in usage), (
        f"{label}: not every joint used the paired frames"
    )

    worst = 0.0
    for bone in BONE_NAMES:
        error = float(np.max(np.abs(world[bone] - targets[bone])))
        worst = max(worst, error)
    assert worst < 1e-9, f"{label}: bone pose error {worst}"

    matrices = chain_joint_matrices(table, world)
    joint_error = 0.0
    rng2 = np.random.default_rng(7)
    for step in JOINT_CHAIN:
        joint = step["joint"]
        rebuilt = matrices[joint]
        # The joint pose is frame-to-frame, so compare the frames themselves.
        _name, child_frame = resolve_frame(
            table, step["child"], joint
        )
        joint_error = max(
            joint_error,
            float(
                np.max(
                    np.abs(world[step["child"]] @ child_frame
                           - rebuilt)
                )
            ),
        )
    assert joint_error < 1e-9, f"{label}: joint frame error {joint_error}"

    print(f"PASS {label}: 5 bones and 4 joints rebuilt within {worst:.2e}")
    for joint, (parent_name, child_name) in usage_frames.items():
        print(f"     {joint:8s} parent frame {parent_name:16s} "
              f"child frame {child_name}")


def check_fallback() -> None:
    rng = np.random.default_rng(11)
    component_transforms = {bone: random_rigid(rng) for bone in BONE_NAMES}
    # The thorax is the fixed root component and always sits at the origin.
    component_transforms["thorax"] = np.eye(4)
    table = {
        joint: {step["child"]: random_rigid(rng) for step in JOINT_CHAIN
                if step["joint"] == joint}
        for joint in (step["joint"] for step in JOINT_CHAIN)
    }
    world, usage = build_bone_transforms(component_transforms, table)
    for item in usage:
        assert item["mode"] == "assembly_pose", item
    for bone in BONE_NAMES:
        error = float(
            np.max(np.abs(world[bone] - component_transforms[bone]))
        )
        assert error < 1e-9, f"fallback pose error on {bone}: {error}"
    print("PASS fallback: missing parent frames keep the assembly pose")


def main() -> None:
    check_case("same frame name in both parts", suffix_parent_frames=False)
    check_case("suffixed parent frame name", suffix_parent_frames=True)
    check_fallback()


if __name__ == "__main__":
    main()
