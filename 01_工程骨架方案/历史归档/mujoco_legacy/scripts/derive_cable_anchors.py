"""Re-derive the wire anchor points from the current bone geometry.

The anchors in ``cable_drive_manifest.json`` were tuned against an older
bone layout.  This tool puts every insertion point back on the surface of its
own bone, inside an anatomically sensible distance band from the joint it
drives, and picks the side that produces the required moment-arm sign.

The result is written as ``insertion_part_m`` -- a point in the bone's part
coordinates -- so the anchor survives any later change of the assembly pose.
The guide point stays procedural: the model builder derives it from the
insertion point and the joint axis.

    D:\\MujocoEnv\\Scripts\\python.exe D:\\XJTLU\\PSP\\mujoco_arm\\derive_cable_anchors.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from cable_drive_common import (
    BONE_NAMES,
    DEFAULT_EXPORT_DIR,
    DEFAULT_MANIFEST,
    JOINT_SPECS,
    load_bone_transforms,
    load_cable_manifest,
    load_joint_frames,
    read_binary_stl,
)


# Distance from the joint centre in which an attachment pad may sit.
ANCHOR_BAND_M = {
    "SC_J": (0.015, 0.045),
    "ACST_J": (0.015, 0.045),
    "GH_J": (0.025, 0.070),
    "ELBOW_J": (0.025, 0.070),
}

BODY_BONE = {"torso": "thorax"}
PARENT_BODY = {}
for _spec in JOINT_SPECS:
    PARENT_BODY[_spec["body"]] = _spec["parent"]
    if _spec["mesh"]:
        BODY_BONE[_spec["body"]] = _spec["mesh"]

JOINT_OF_BODY = {spec["body"]: spec for spec in JOINT_SPECS}


def upstream_joints(body: str) -> list[dict]:
    """Every joint that moves ``body``, nearest joint first."""

    result = []
    current = body
    while current in PARENT_BODY:
        spec = JOINT_OF_BODY.get(current)
        if spec is not None:
            result.append(spec)
        current = PARENT_BODY[current]
    return result


def unit(vectors: np.ndarray) -> np.ndarray:
    lengths = np.linalg.norm(vectors, axis=-1, keepdims=True)
    lengths[lengths < 1e-12] = 1.0
    return vectors / lengths


def load_bone_vertices(bone_dir: Path) -> dict[str, np.ndarray]:
    vertices = {}
    for bone in BONE_NAMES:
        path = bone_dir / f"{bone}.stl"
        points = read_binary_stl(path).reshape(-1, 3)
        if float(np.max(np.ptp(points, axis=0))) > 5.0:
            points = points * 0.001
        vertices[bone] = points
    return vertices


def joint_axes(frames: dict) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    axes = {}
    for spec in JOINT_SPECS:
        frame = frames[spec["frame"]]
        origin = np.asarray(frame["origin_m"], dtype=float)
        axis = np.asarray(
            frame["axes"][int(spec["axis_index"])], dtype=float
        )
        axes[spec["joint"]] = (origin, axis / np.linalg.norm(axis))
    return axes


def guide_point(
    origin: np.ndarray,
    axis: np.ndarray,
    point: np.ndarray,
    offset: float,
    sign: float,
) -> np.ndarray:
    """The builder's constant-moment-arm guide for one insertion point."""

    vector = point - origin
    along = float(np.dot(vector, axis))
    perpendicular = vector - axis * along
    length = float(np.linalg.norm(perpendicular))
    if length < 1e-9:
        return origin + axis * along
    perpendicular = perpendicular / length
    tangent = np.cross(axis, perpendicular)
    if sign < 0.0:
        tangent = -tangent
    return origin + axis * along + perpendicular * length + tangent * offset


def choose_anchor(
    points: np.ndarray,
    origin: np.ndarray,
    axis: np.ndarray,
    joints: list[dict],
    axes: dict[str, tuple[np.ndarray, np.ndarray]],
    target_joint: str,
    offset: float,
    sign: float,
    band: tuple[float, float],
) -> tuple[np.ndarray, float, float]:
    """Pick the surface point with the best moment arm for the target joint."""

    vector = points - origin
    distance = np.linalg.norm(vector, axis=1)
    inside = (distance >= band[0]) & (distance <= band[1])
    candidates = points[inside]
    if not len(candidates):
        raise ValueError("No bone vertices inside the anchor band.")

    along = (candidates - origin) @ axis
    perpendicular = (candidates - origin) - np.outer(along, axis)
    length = np.linalg.norm(perpendicular, axis=1)
    length[length < 1e-9] = 1.0
    perpendicular = perpendicular / length[:, None]
    tangent = np.cross(np.broadcast_to(axis, perpendicular.shape), perpendicular)
    if sign < 0.0:
        tangent = -tangent
    direction = unit(tangent)

    target_origin, target_axis = axes[target_joint]
    target_moment = -np.einsum(
        "ij,ij->i",
        direction,
        np.cross(np.broadcast_to(target_axis, candidates.shape), candidates - target_origin),
    )

    other = np.zeros(len(candidates))
    other_joint = ""
    for spec in joints:
        name = spec["joint"]
        if name == target_joint:
            continue
        joint_origin, joint_axis = axes[name]
        moment = -np.einsum(
            "ij,ij->i",
            direction,
            np.cross(
                np.broadcast_to(joint_axis, candidates.shape),
                candidates - joint_origin,
            ),
        )
        peak = float(np.max(np.abs(moment)))
        if peak > other.max(initial=0.0):
            other = np.abs(moment)
            other_joint = name

    want_negative = sign > 0.0
    usable = target_moment < 0.0 if want_negative else target_moment > 0.0
    if not np.any(usable):
        raise ValueError(
            f"No candidate on the bone drives {target_joint} in the required "
            "direction. Check the joint axis direction."
        )
    # Only the moment arm about the joint the cable is meant to drive is
    # scored.  A cable that runs along the limb also torques the joints above
    # it, which is physical, so those are reported but not penalised.
    score = np.where(usable, np.abs(target_moment), -1e9)
    best = int(np.argmax(score))
    return (
        candidates[best].copy(),
        float(target_moment[best]),
        float(other[best]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT_DIR)
    parser.add_argument(
        "--joint-frames",
        type=Path,
        default=DEFAULT_EXPORT_DIR / "joint_frames.json",
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    manifest = load_cable_manifest(args.manifest)
    frames = load_joint_frames(args.joint_frames)
    axes = joint_axes(frames)
    transforms = load_bone_transforms(
        args.export_dir / "bone_transforms.json"
    )
    vertices = load_bone_vertices(args.export_dir / "bones")

    output = args.output or args.manifest
    updated = json.loads(json.dumps(manifest))

    for cable in updated["cables"]:
        if not cable.get("enabled", True):
            continue
        insertion_body = cable["insertion_body"]
        bone = BODY_BONE[insertion_body]
        target_frame = next(
            spec["frame"]
            for spec in JOINT_SPECS
            if spec["joint"] == cable["target_joint"]
        )
        if bone not in transforms:
            raise ValueError(f"No bone transform for {bone}.")
        origin, axis = axes[cable["target_joint"]]
        sign = float(cable.get("guide_moment_sign", 1.0))
        offset = float(cable.get("guide_offset_m", 0.03))
        band = ANCHOR_BAND_M[target_frame]

        anchor_world, moment, other = choose_anchor(
            vertices[bone],
            origin,
            axis,
            upstream_joints(insertion_body),
            axes,
            cable["target_joint"],
            offset,
            sign,
            band,
        )
        bone_transform = transforms[bone]
        anchor_part = np.linalg.inv(bone_transform) @ np.append(
            anchor_world, 1.0
        )
        anchor_part = anchor_part[:3]
        local = anchor_world - np.asarray(frames[target_frame]["origin_m"])

        cable["insertion_part_m"] = anchor_part.tolist()
        cable["insertion_local_m"] = local.tolist()
        cable["guide_mode"] = "constant_moment_arm"
        cable.pop("guides", None)
        cable["via_body"] = insertion_body
        cable["via_local_m"] = local.tolist()
        print(
            f"{cable['id']:18s} -> {bone:9s} dL/dq={moment: .4f} "
            f"(next {other: .4f})"
        )

    output.write_text(
        json.dumps(updated, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(f"manifest: {output}")


if __name__ == "__main__":
    main()
