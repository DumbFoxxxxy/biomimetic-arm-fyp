"""Convert SolidWorks part STLs and joint frames into MuJoCo-ready exports.

Two pose sources are supported:

``paired frames``
    For every joint there is a coordinate system on the child bone and a
    matching one on the parent bone.  Each child bone is placed so that its
    frame coincides with the parent frame, which makes the pose inside the
    SolidWorks assembly irrelevant.

``assembly pose``
    Only the child-side frame exists, so the bone is placed at the pose
    recorded in ``component_transforms.txt``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from cable_drive_common import (
    BONE_NAMES,
    DEFAULT_EXPORT_DIR,
    FRAME_AXIS_NAMES,
    JOINT_CHAIN,
    component_base_name,
    load_bone_scales,
    read_binary_stl,
    transform_from_array,
    write_binary_stl,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_COMPONENTS = DEFAULT_EXPORT_DIR / "component_transforms.txt"
DEFAULT_LOCAL_STL_DIR = DEFAULT_EXPORT_DIR / "bones_local"
DEFAULT_JOINT_FRAMES_TXT = DEFAULT_EXPORT_DIR / "joint_frames.txt"
DEFAULT_JOINT_FRAMES_JSON = DEFAULT_EXPORT_DIR / "joint_frames.json"
DEFAULT_PART_FRAMES_TXT = DEFAULT_EXPORT_DIR / "part_frames.txt"

FRAME_RANGES = {
    "SC_J": [[-30, 40], [-20, 40], [-30, 30]],
    "ACST_J": [[-20, 50], [-15, 35], [-20, 30]],
    "GH_J": [[-90, 90], [-45, 180], [-45, 180]],
    "ELBOW_J": [[0, 145]],
}


def read_pipe_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-16").splitlines()
    if not lines:
        return []
    headers = lines[0].split("|")
    return [
        dict(zip(headers, line.split("|")))
        for line in lines[1:]
        if line.strip()
    ]


def find_component(
    rows: list[dict[str, str]], bone_name: str
) -> dict[str, str]:
    for row in rows:
        if component_base_name(row.get("name", "")) == bone_name:
            if row.get("suppressed", "").lower() != "true":
                return row
    raise ValueError(f"No active component found for bone: {bone_name}")


def local_stl_path(local_dir: Path, bone_name: str) -> Path:
    candidates = (
        local_dir / f"{bone_name}.stl",
        local_dir / f"{bone_name}.STL",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Missing local STL for {bone_name}: expected {candidates[0]}"
    )


def convert_stl_to_global(
    source: Path,
    destination: Path,
    transform: np.ndarray,
    offset: np.ndarray | None = None,
    bone_scale: float = 1.0,
) -> None:
    triangles = read_binary_stl(source)
    vertices = triangles.reshape(-1, 3)
    extent = np.ptp(vertices, axis=0)
    scale = 0.001 if float(np.max(extent)) > 5.0 else 1.0
    vertices_m = vertices * scale
    if offset is not None:
        vertices_m = vertices_m + offset
    if bone_scale != 1.0:
        # The offset puts the mesh into part coordinates first, so this scale
        # is taken about the part origin -- the proximal joint centre.
        vertices_m = vertices_m * bone_scale
    homogeneous = np.column_stack(
        (vertices_m, np.ones(len(vertices_m), dtype=float))
    )
    transformed = (transform @ homogeneous.T).T[:, :3]
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_binary_stl(
        destination,
        transformed.reshape(-1, 3, 3),
        f"global assembly STL for {destination.stem}",
    )


def orthonormalize(columns: np.ndarray) -> np.ndarray:
    basis: list[np.ndarray] = []
    for index in range(columns.shape[1]):
        vector = columns[:, index].astype(float)
        for existing in basis:
            vector -= existing * float(np.dot(existing, vector))
        length = float(np.linalg.norm(vector))
        if length > 1e-8:
            basis.append(vector / length)
    for candidate in np.eye(3):
        if len(basis) >= 3:
            break
        vector = candidate.copy()
        for existing in basis:
            vector -= existing * float(np.dot(existing, vector))
        length = float(np.linalg.norm(vector))
        if length > 1e-8:
            basis.append(vector / length)
    result = np.column_stack(basis[:3])
    if np.linalg.det(result) < 0.0:
        result[:, 2] *= -1.0
    return result


def frames_from_matrices(matrices: dict[str, np.ndarray]) -> dict:
    frames = {}
    for frame_name, axis_names in FRAME_AXIS_NAMES.items():
        transform = matrices.get(frame_name)
        if transform is None:
            raise ValueError(f"joint frame is missing {frame_name}.")
        rotation = orthonormalize(transform[:3, :3])
        origin = transform[:3, 3]
        axes = []
        for index, axis_name in enumerate(axis_names):
            axes.append(
                {
                    "name": axis_name,
                    "direction": rotation[:, index].tolist(),
                }
            )
        frames[frame_name] = {
            "origin_m": origin.tolist(),
            "axes": axes,
            "range_deg": FRAME_RANGES[frame_name],
        }
    return frames


def write_joint_frames(
    destination: Path,
    frames: dict,
    source: str,
    layout: dict | None = None,
) -> None:
    payload = {"units": "meters", "source": source, "frames": frames}
    if layout is not None:
        payload["layout"] = layout
    destination.write_text(
        json.dumps(payload, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def load_component_transforms(
    rows: list[dict[str, str]],
) -> dict[str, np.ndarray]:
    transforms = {}
    for bone_name in BONE_NAMES:
        row = find_component(rows, bone_name)
        values = [float(row[f"xform{index}"]) for index in range(16)]
        transforms[bone_name] = transform_from_array(values)
    return transforms


def load_part_probe(path: Path) -> dict[str, dict[str, object]]:
    """Centre of mass and volume per bone, in the part's own coordinates.

    The SolidWorks STL exporter writes the mesh in its own output coordinate
    system, which is not always the part origin.  The centre of mass is always
    reported in part coordinates, so the difference between it and the mesh
    centroid recovers the offset of the exported STL.
    """

    probe: dict[str, dict[str, object]] = {}
    if not path.exists():
        return probe
    for row in read_pipe_table(path):
        bone = component_base_name(row.get("component", ""))
        if bone not in BONE_NAMES:
            continue
        try:
            probe[bone] = {
                "centre_of_mass_m": np.array(
                    [float(row["comx"]), float(row["comy"]), float(row["comz"])],
                    dtype=float,
                ),
                "volume_m3": float(row["volume"]),
            }
        except (KeyError, TypeError, ValueError):
            continue
    return probe


def mesh_statistics(
    source: Path,
) -> tuple[np.ndarray, float, float]:
    """Return (centroid, volume, scale) of a binary STL in metres."""

    triangles = read_binary_stl(source)
    vertices = triangles.reshape(-1, 3)
    extent = np.ptp(vertices, axis=0)
    scale = 0.001 if float(np.max(extent)) > 5.0 else 1.0
    vertices = vertices * scale
    corner0 = triangles[:, 0] * scale
    corner1 = triangles[:, 1] * scale
    corner2 = triangles[:, 2] * scale
    signed = np.einsum(
        "ij,ij->i", corner0, np.cross(corner1, corner2)
    )
    total = float(signed.sum())
    if abs(total) < 1e-18:
        return vertices.mean(axis=0), 0.0, scale
    volume = total / 6.0
    centroid = (
        ((corner0 + corner1 + corner2) / 4.0) * signed[:, None]
    ).sum(axis=0) / total
    return centroid, volume, scale


def stl_part_offset(
    source: Path,
    probe: dict[str, object],
) -> np.ndarray:
    """Offset to add to the exported STL so its origin matches the part.

    Returns zero when no probe data is available, which keeps the previous
    behaviour for assemblies exported with the part origin as the STL output
    coordinate system.
    """

    if not probe:
        return np.zeros(3)
    centroid, volume, _scale = mesh_statistics(source)
    reference = float(probe["volume_m3"])
    if volume > 1e-15 and reference > 1e-15:
        ratio = reference / volume
        if not 0.9 < ratio < 1.1:
            raise ValueError(
                f"{source.name}: exported mesh volume {volume:.3e} m3 does "
                f"not match the SolidWorks volume {reference:.3e} m3."
            )
    return np.asarray(probe["centre_of_mass_m"], dtype=float) - centroid


def load_part_frames(
    part_frames_path: Path,
    bone_scales: dict[str, float] | None = None,
) -> dict[str, dict[str, np.ndarray]]:
    """frame name -> {bone name: transform of the frame in that part}."""

    table: dict[str, dict[str, np.ndarray]] = {}
    for row in read_pipe_table(part_frames_path):
        name = row.get("name", "").strip()
        if not name or name == "ERROR":
            continue
        bone = component_base_name(row.get("component", ""))
        if bone not in BONE_NAMES:
            continue
        values = [float(row[f"v{index}"]) for index in range(16)]
        matrix = transform_from_array(values)
        factor = 1.0 if bone_scales is None else float(bone_scales.get(bone, 1.0))
        if factor != 1.0:
            # A frame inside a scaled part scales with it: the frame origin is
            # measured from the part origin, so it simply moves out by the
            # same factor.  Its orientation is unchanged.
            matrix = matrix.copy()
            matrix[:3, 3] = matrix[:3, 3] * factor
        table.setdefault(name, {})[bone] = matrix
    if not table:
        raise ValueError(
            f"{part_frames_path} contains no usable part coordinate systems."
        )
    return table


def resolve_frame(
    table: dict[str, dict[str, np.ndarray]],
    bone: str,
    joint: str,
) -> tuple[str | None, np.ndarray | None]:
    """Pick the coordinate system that describes ``joint`` on ``bone``.

    An exact name match wins.  Otherwise a unique name that starts with the
    joint name is accepted, so ``SC_J``, ``SC_J_thorax`` and ``SC_J_P`` all
    work as long as they are unambiguous.
    """

    for name, owners in table.items():
        if bone in owners and name.upper() == joint.upper():
            return name, owners[bone]
    matches = [
        (name, owners[bone])
        for name, owners in table.items()
        if bone in owners and name.upper().startswith(joint.upper())
    ]
    if len(matches) == 1:
        return matches[0]
    return None, None


def frame_layout(
    table: dict[str, dict[str, np.ndarray]],
) -> tuple[dict[str, str], list[str]]:
    """Resolve the frame pair used by every joint."""

    resolved: dict[str, str] = {}
    problems: list[str] = []
    for step in JOINT_CHAIN:
        joint, parent, child = step["joint"], step["parent"], step["child"]
        parent_name, _parent_frame = resolve_frame(table, parent, joint)
        child_name, child_frame = resolve_frame(table, child, joint)
        if child_frame is None:
            problems.append(
                f"{joint}: no usable coordinate system on {child}"
            )
            continue
        resolved[joint] = child_name
        resolved[f"{joint}:parent"] = parent_name or ""
    return resolved, problems


def build_bone_transforms(
    component_transforms: dict[str, np.ndarray],
    table: dict[str, dict[str, np.ndarray]],
) -> tuple[dict[str, np.ndarray], list[dict]]:
    """Return part-to-world transforms for every bone plus a usage report."""

    world: dict[str, np.ndarray] = {"thorax": np.eye(4)}
    usage: list[dict] = []
    for step in JOINT_CHAIN:
        joint, parent, child = step["joint"], step["parent"], step["child"]
        parent_name, parent_frame = resolve_frame(table, parent, joint)
        child_name, child_frame = resolve_frame(table, child, joint)
        if child_frame is None:
            raise ValueError(
                f"{joint} needs a coordinate system inside {child}."
            )
        if parent_frame is not None:
            world[child] = (
                world[parent] @ parent_frame @ np.linalg.inv(child_frame)
            )
            mode = "paired_frames"
        else:
            world[child] = (
                world[parent]
                @ np.linalg.inv(component_transforms[parent])
                @ component_transforms[child]
            )
            mode = "assembly_pose"
        usage.append(
            {
                "joint": joint,
                "parent": parent,
                "child": child,
                "mode": mode,
                "parent_frame": parent_name,
                "child_frame": child_name,
            }
        )
    return world, usage


def chain_joint_matrices(
    table: dict[str, dict[str, np.ndarray]],
    bone_transforms: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    matrices: dict[str, np.ndarray] = {}
    for step in JOINT_CHAIN:
        joint, child = step["joint"], step["child"]
        _name, child_frame = resolve_frame(table, child, joint)
        if child_frame is None:
            raise ValueError(
                f"{joint} needs a coordinate system inside {child}."
            )
        matrices[joint] = bone_transforms[child] @ child_frame
    return matrices


def assembly_level_matrices(path: Path) -> dict[str, np.ndarray] | None:
    if not path.exists():
        return None
    rows = read_pipe_table(path)
    by_name = {row.get("name", ""): row for row in rows}
    matrices: dict[str, np.ndarray] = {}
    for frame_name in FRAME_AXIS_NAMES:
        row = by_name.get(frame_name)
        if row is None or row.get("v0") == "MISSING" or not row.get("v0"):
            return None
        matrices[frame_name] = transform_from_array(
            [float(row[f"v{i}"]) for i in range(16)]
        )
    return matrices


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--components", type=Path, default=DEFAULT_COMPONENTS)
    parser.add_argument(
        "--local-stl-dir", type=Path, default=DEFAULT_LOCAL_STL_DIR
    )
    parser.add_argument(
        "--joint-frames-txt", type=Path, default=DEFAULT_JOINT_FRAMES_TXT
    )
    parser.add_argument(
        "--part-frames-txt", type=Path, default=DEFAULT_PART_FRAMES_TXT
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_EXPORT_DIR
    )
    parser.add_argument(
        "--joint-frames-json",
        type=Path,
        default=DEFAULT_JOINT_FRAMES_JSON,
    )
    args = parser.parse_args()
    scales = load_bone_scales()
    for bone_name, factor in scales.items():
        if factor != 1.0:
            print(
                f"bone scale: {bone_name} x{factor:g} about the part origin "
                "(set in exports/bone_scale.json)"
            )

    if not args.components.exists():
        raise FileNotFoundError(f"Missing components table: {args.components}")

    component_transforms = load_component_transforms(
        read_pipe_table(args.components)
    )

    assembly_matrices = assembly_level_matrices(args.joint_frames_txt)
    if assembly_matrices is not None:
        bone_transforms = component_transforms
        matrices = assembly_matrices
        layout = {
            "mode": "assembly_level_frames",
            "joints": [],
        }
        source = str(args.joint_frames_txt)
    else:
        if not args.part_frames_txt.exists():
            raise FileNotFoundError(
                "Missing joint_frames.txt and part_frames.txt. Create the "
                "joint coordinate systems in SolidWorks and run the export "
                "macro."
            )
        table = load_part_frames(args.part_frames_txt, scales)
        bone_transforms, usage = build_bone_transforms(
            component_transforms, table
        )
        matrices = chain_joint_matrices(table, bone_transforms)
        modes = {item["mode"] for item in usage}
        if modes == {"paired_frames"}:
            mode_label = "paired_frames"
        elif modes == {"assembly_pose"}:
            mode_label = "assembly_pose"
        else:
            mode_label = "mixed"
        layout = {"mode": mode_label, "joints": usage}
        source = str(args.part_frames_txt)

    bone_dir = args.output_dir / "bones"
    probe = load_part_probe(args.output_dir / "part_probe.txt")
    for bone_name in BONE_NAMES:
        stl_source = local_stl_path(args.local_stl_dir, bone_name)
        probe_entry = probe.get(bone_name)
        offset = (
            stl_part_offset(stl_source, probe_entry) if probe_entry else None
        )
        bone_scale = scales[bone_name]
        convert_stl_to_global(
            stl_source,
            bone_dir / f"{bone_name}.stl",
            bone_transforms[bone_name],
            offset,
            bone_scale,
        )
        notes = []
        if offset is not None:
            notes.append(
                "stl offset corrected by "
                f"{np.round(offset * 1000.0, 1).tolist()} mm"
            )
        if bone_scale != 1.0:
            notes.append(f"scaled x{bone_scale:g} about the part origin")
        suffix = ("  " + "; ".join(notes)) if notes else ""
        print(f"mesh: {bone_name} -> {bone_dir / (bone_name + '.stl')}{suffix}")

    write_joint_frames(
        args.joint_frames_json,
        frames_from_matrices(matrices),
        source,
        layout,
    )
    print(f"joint frames: {args.joint_frames_json}")

    bone_transforms_path = args.output_dir / "bone_transforms.json"
    bone_transforms_path.write_text(
        json.dumps(
            {
                "units": "meters",
                "source": source,
                "comment": (
                    "part-to-world transform for every bone. Attach points "
                    "stored in part coordinates are mapped with these."
                ),
                "transforms": {
                    bone: bone_transforms[bone].tolist()
                    for bone in BONE_NAMES
                },
            },
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"bone transforms: {bone_transforms_path}")
    print(f"pose source: {layout['mode']}")
    for item in layout.get("joints", []):
        parent = item.get("parent_frame") or "-"
        print(
            f"  {item['joint']:8s} {item['mode']:14s} "
            f"parent frame = {parent:16s} child frame = {item['child_frame']}"
        )
    if layout["mode"] != "paired_frames":
        missing = [
            item["joint"]
            for item in layout.get("joints", [])
            if item["mode"] != "paired_frames"
        ]
        if missing:
            print(
                "  note: these joints still use the assembly pose: "
                + ", ".join(missing)
            )


if __name__ == "__main__":
    main()
