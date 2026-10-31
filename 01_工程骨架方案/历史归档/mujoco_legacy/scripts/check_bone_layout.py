"""Check that the exported bones and joint frames form a usable arm skeleton.

The MuJoCo model places every bone at the pose recorded in
``component_transforms.txt`` and every joint at the origin of the matching
named coordinate system.  If the SolidWorks assembly still has the bone parts
scattered, the simulation reproduces that scatter and the viewer looks broken.

The script reports the two quantities that decide whether the export is
usable:

* how far each named joint frame sits from the bone that owns it
* how far apart the two bones that share a joint are

It needs no CAD rebuild and no MuJoCo model, so it can run right after
``prepare_cable_drive_exports.py``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from cable_drive_common import (
    BONE_NAMES,
    DEFAULT_EXPORT_DIR,
    DEFAULT_REPORT_DIR,
    JOINT_CHAIN,
    load_joint_frames,
    read_binary_stl,
)
from prepare_cable_drive_exports import read_pipe_table


# A joint centre sits inside the bone end, so its distance to the nearest
# surface vertex is roughly the local bone radius.  Anything much larger means
# the coordinate system was not placed on the joint.
FRAME_ON_BONE_LIMIT_M = 0.030

# Bones that share a joint do not have to touch exactly, but a gap of more
# than a centimetre means the assembly pose is not usable for a simulation.
JOINT_GAP_LIMIT_M = 0.015

FRAME_OWNERS = {
    "SC_J": "clavicle",
    "ACST_J": "scapula",
    "GH_J": "humerus",
    "ELBOW_J": "ulna",
}

JOINTS = tuple(
    (step["joint"], step["parent"], step["child"]) for step in JOINT_CHAIN
)


def load_vertices(bone_dir: Path) -> dict[str, np.ndarray]:
    vertices: dict[str, np.ndarray] = {}
    for bone in BONE_NAMES:
        path = bone_dir / f"{bone}.stl"
        if not path.exists():
            raise FileNotFoundError(
                f"Missing bone mesh {path}. Run "
                "prepare_cable_drive_exports.py first."
            )
        points = read_binary_stl(path).reshape(-1, 3)
        # The exporter writes metres; anything larger is still millimetres.
        if float(np.max(np.ptp(points, axis=0))) > 5.0:
            points = points * 0.001
        vertices[bone] = points
    return vertices


def nearest_distance(
    points: np.ndarray,
    target: np.ndarray,
    coarse_step: int = 20,
) -> float:
    """Distance from ``points`` to the closest vertex of ``target``."""

    coarse = target[::coarse_step]
    best = np.inf
    best_point = np.zeros(3)
    for start in range(0, len(points), 2048):
        block = points[start : start + 2048]
        squared = ((block[:, None, :] - coarse[None, :, :]) ** 2).sum(axis=2)
        flat = int(np.argmin(squared))
        row, column = divmod(flat, len(coarse))
        if squared[row, column] < best:
            best = float(squared[row, column])
            best_point = block[row].copy()

    if not np.isfinite(best):
        raise ValueError("Cannot measure a distance on an empty mesh.")

    # Refine against the full target mesh around the coarse match.
    radius = 5.0 * float(np.sqrt(best)) + 0.02
    local = target[np.linalg.norm(target - best_point, axis=1) < radius]
    if len(local):
        squared = ((best_point[None, :] - local) ** 2).sum(axis=1)
        return float(np.sqrt(np.min(squared)))
    return float(np.sqrt(best))


def subsample(vertices: dict[str, np.ndarray], bone: str, step: int = 3) -> np.ndarray:
    return vertices[bone][::step]


def load_part_boxes(path: Path) -> dict[str, np.ndarray]:
    """bone name -> part bounding box in part coordinates (metres)."""

    boxes: dict[str, np.ndarray] = {}
    if not path.exists():
        return boxes
    for row in read_pipe_table(path):
        bone = row.get("bone", "").strip().lower()
        if bone not in BONE_NAMES:
            continue
        try:
            boxes[bone] = np.array(
                [float(row[key]) for key in
                 ("minx", "miny", "minz", "maxx", "maxy", "maxz")],
                dtype=float,
            )
        except (KeyError, TypeError, ValueError):
            continue
    return boxes


def stl_part_mismatch(
    stl_points: np.ndarray, part_box: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """Compare the STL size with the size of the part in part coordinates.

    The comparison is scale free, so it works whether the STL was written in
    millimetres or metres.  A rotation or a full solid-body offset survives the
    comparison and shows up as a size mismatch.
    """

    stl_extent = np.ptp(stl_points, axis=0)
    part_extent = part_box[3:] - part_box[:3]
    reference = float(np.max(part_extent))
    if reference <= 0.0 or float(np.max(stl_extent)) <= 0.0:
        return stl_extent, part_extent, 0.0
    scaled = stl_extent * (reference / float(np.max(stl_extent)))
    error = float(np.max(np.abs(scaled - part_extent))) / reference
    return stl_extent, part_extent, error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT_DIR)
    parser.add_argument(
        "--joint-frames",
        type=Path,
        default=DEFAULT_EXPORT_DIR / "joint_frames.json",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_REPORT_DIR / "bone_layout_check.json",
    )
    args = parser.parse_args()

    vertices = load_vertices(args.export_dir / "bones")
    frames = load_joint_frames(args.joint_frames)
    layout = json.loads(args.joint_frames.read_text(encoding="utf-8")).get(
        "layout", {}
    )

    print(f"Pose source: {layout.get('mode', 'assembly_pose')}")
    for item in layout.get("joints", []):
        parent = item.get("parent_frame") or "-"
        print(
            f"  {item['joint']:8s} {item['mode']:14s} "
            f"parent frame = {parent:16s} child frame = {item['child_frame']}"
        )

    print("Bone bounding boxes in assembly coordinates (mm)")
    for bone in BONE_NAMES:
        points = vertices[bone]
        low = points.min(axis=0) * 1000.0
        high = points.max(axis=0) * 1000.0
        print(f"  {bone:9s} min {np.round(low, 1)}  max {np.round(high, 1)}")

    part_boxes = load_part_boxes(args.export_dir / "part_boxes.txt")
    stl_errors: dict[str, float] = {}
    if part_boxes:
        print()
        print("STL geometry against the part's own bounding box")
        for bone in BONE_NAMES:
            if bone not in part_boxes:
                continue
            stl_extent, part_extent, error = stl_part_mismatch(
                vertices[bone], part_boxes[bone]
            )
            stl_errors[bone] = error
            flag = "OK  " if error <= 0.03 else "FAIL"
            print(
                f"  {flag} {bone:9s} stl size "
                f"{np.round(stl_extent * 1000.0, 1)} vs part size "
                f"{np.round(part_extent * 1000.0, 1)} mm"
            )
        if any(value > 0.03 for value in stl_errors.values()):
            print(
                "  The STL was not written in the part's own coordinate "
                "system. Re-export it with the output coordinate system set "
                "to the part origin."
            )

    results: dict[str, object] = {"frames": {}, "joints": {}, "passed": True}
    results["stl_vs_part_error"] = stl_errors

    print()
    print("Joint frame placement")
    for frame_name, owner in FRAME_OWNERS.items():
        origin = np.asarray(frames[frame_name]["origin_m"], dtype=float)
        distance = nearest_distance(origin[None, :], subsample(vertices, owner))
        ok = distance <= FRAME_ON_BONE_LIMIT_M
        results["frames"][frame_name] = {
            "owner": owner,
            "origin_m": origin.tolist(),
            "distance_to_bone_m": distance,
            "ok": ok,
        }
        flag = "OK  " if ok else "FAIL"
        print(
            f"  {flag} {frame_name:8s} on {owner:9s} "
            f"distance to bone = {distance * 1000.0:7.1f} mm"
        )

    print()
    print("Bone-to-bone gap at each joint")
    for frame_name, parent, child in JOINTS:
        gap = nearest_distance(
            subsample(vertices, child), subsample(vertices, parent)
        )
        ok = gap <= JOINT_GAP_LIMIT_M
        results["joints"][frame_name] = {
            "parent": parent,
            "child": child,
            "gap_m": gap,
            "ok": ok,
        }
        flag = "OK  " if ok else "FAIL"
        print(f"  {flag} {parent:9s} -> {child:9s} gap = {gap * 1000.0:7.1f} mm")

    bad_frames = [
        name for name, value in results["frames"].items() if not value["ok"]
    ]
    bad_joints = [
        name for name, value in results["joints"].items() if not value["ok"]
    ]
    results["passed"] = not (bad_frames or bad_joints)

    print()
    if not results["passed"]:
        if bad_frames:
            print("Frames not on their bone: " + ", ".join(bad_frames))
        if bad_joints:
            print("Bones not in joint contact: " + ", ".join(bad_joints))
        print(
            "The exported bone pose is not a usable arm. Place the bone parts "
            "in a proper joint pose inside biomimetic_structure.SLDASM, keep "
            "the joint frames on the bones, then export again."
        )
    else:
        print("Layout check passed.")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(results, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(f"report: {args.report}")
    raise SystemExit(0 if results["passed"] else 1)


if __name__ == "__main__":
    main()
