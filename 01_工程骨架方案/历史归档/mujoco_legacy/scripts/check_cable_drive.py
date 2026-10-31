"""Validate the generated shoulder/elbow cable-drive model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mujoco
import numpy as np

from cable_drive_common import (
    DEFAULT_EXPORT_DIR,
    DEFAULT_MANIFEST,
    DEFAULT_OUTPUT_XML,
    DEFAULT_REPORT_DIR,
    JOINT_SPECS,
    SKELETON_JOINT_NAMES,
    cable_guides,
    load_cable_manifest,
    read_binary_stl,
    task_envelope_states,
)


def cable_channels(manifest: dict) -> list[dict]:
    result = []
    default_travel = float(
        manifest.get("defaults", {}).get("travel_m", 0.0)
    )
    for muscle in manifest["cables"]:
        if not muscle.get("enabled", True):
            continue
        guides = cable_guides(muscle)
        for channel in range(int(muscle.get("channels", 1))):
            cable_id = f"{muscle['id']}_{channel}"
            result.append(
                {
                    "id": cable_id,
                    "muscle_id": muscle["id"],
                    "group": muscle.get("group", muscle["id"]),
                    "cable": f"cable_{cable_id}",
                    "slider_joint": f"slider_{cable_id}",
                    "actuator": f"pull_{cable_id}",
                    "insertion_body": muscle["insertion_body"],
                    "guides": [
                        {
                            "body": guide["body"],
                            "site": f"cable_{cable_id}_via_{index}",
                        }
                        for index, guide in enumerate(guides)
                    ],
                    "max_force_n": float(muscle["max_force_n"]),
                    "travel_m": float(
                        muscle.get("travel_m", default_travel)
                    ),
                }
            )
    return result


def body_to_bone() -> dict[str, str]:
    result = {"torso": "thorax"}
    for spec in JOINT_SPECS:
        if spec["mesh"]:
            result[spec["body"]] = spec["mesh"]
    return result


def tendon_id(model: mujoco.MjModel, name: str) -> int:
    value = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_TENDON, name)
    if value < 0:
        raise RuntimeError(f"Tendon missing: {name}")
    return int(value)


def qpos_address(model: mujoco.MjModel, joint_name: str) -> int:
    joint_id = mujoco.mj_name2id(
        model, mujoco.mjtObj.mjOBJ_JOINT, joint_name
    )
    if joint_id < 0:
        raise RuntimeError(f"Joint missing: {joint_name}")
    return int(model.jnt_qposadr[joint_id])


def sensitivity_matrix(
    model: mujoco.MjModel,
    channels: list[dict],
) -> dict[str, dict[str, float]]:
    data = mujoco.MjData(model)
    result = {}
    epsilon = 1e-5
    for joint_name in SKELETON_JOINT_NAMES:
        address = qpos_address(model, joint_name)
        sensitivities = {}
        for channel in channels:
            tid = tendon_id(model, channel["cable"])
            data.qpos[:] = 0.0
            data.qpos[address] = epsilon
            mujoco.mj_forward(model, data)
            plus = float(data.ten_length[tid])
            data.qpos[:] = 0.0
            data.qpos[address] = -epsilon
            mujoco.mj_forward(model, data)
            minus = float(data.ten_length[tid])
            sensitivities[channel["id"]] = (plus - minus) / (2.0 * epsilon)
        result[joint_name] = sensitivities
    return result


def sample_states(model: mujoco.MjModel, samples: int = 300) -> list[np.ndarray]:
    return task_envelope_states(model, random_samples=samples, seed=23)


def clearance_report(
    model: mujoco.MjModel,
    manifest: dict,
    export_dir: Path,
) -> list[dict]:
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return [
            {
                "severity": "warning",
                "message": "scipy is unavailable; STL clearance check skipped.",
            }
        ]

    bone_points = {}
    bone_meshes = {}
    for bone in ("thorax", "clavicle", "scapula", "humerus", "ulna"):
        path = export_dir / "bones" / f"{bone}.stl"
        if path.exists():
            triangles = read_binary_stl(path)
            points = triangles.reshape(-1, 3)
            if np.max(np.ptp(points, axis=0)) > 5.0:
                triangles = triangles * 0.001
                points = points * 0.001
            bone_points[bone] = points
            bone_meshes[bone] = triangles
    trees = {
        bone: cKDTree(points) for bone, points in bone_points.items()
    }
    body_bone = body_to_bone()
    findings = []
    for channel in cable_channels(manifest):
        insertion_name = f"cable_{channel['id']}_insertion"
        slider_name = f"slider_{channel['id']}_site"
        guide_names = [guide["site"] for guide in channel["guides"]]
        site_names = [insertion_name] + guide_names + [slider_name]
        site_ids = []
        for name in site_names:
            site_id = mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_SITE, name
            )
            if site_id < 0:
                findings.append(
                    {
                        "severity": "error",
                        "message": f"Missing site: {name}",
                    }
                )
                site_ids.append(None)
            else:
                site_ids.append(int(site_id))
        if any(site_id is None for site_id in site_ids):
            continue

        data = mujoco.MjData(model)
        data.qpos[:] = 0.0
        mujoco.mj_forward(model, data)
        points = np.array([data.site_xpos[site_id] for site_id in site_ids])
        path_bodies = [channel["insertion_body"]] + [
            guide["body"] for guide in channel["guides"]
        ] + [None]
        for segment_index, (start, end) in enumerate(
            zip(points[:-1], points[1:])
        ):
            ignore = {
                body_bone.get(path_bodies[segment_index]),
                body_bone.get(path_bodies[segment_index + 1]),
            }
            for ratio in np.linspace(0.0, 1.0, 80)[1:-1]:
                point = start + ratio * (end - start)
                for bone, tree in trees.items():
                    if bone in ignore:
                        continue
                    distance = float(tree.query(point)[0])
                    if distance < 0.006 and bool(
                        points_inside_mesh(
                            point[None, :], bone_meshes[bone]
                        )[0]
                    ):
                        findings.append(
                            {
                                "severity": "error",
                                "cable": channel["id"],
                                "segment": segment_index,
                                "bone": bone,
                                "message": (
                                    f"{channel['id']} segment {segment_index} "
                                    f"passes through {bone}."
                                ),
                            }
                        )
                        break
                    if distance < 0.003:
                        findings.append(
                            {
                                "severity": "warning",
                                "cable": channel["id"],
                                "segment": segment_index,
                                "bone": bone,
                                "distance_m": distance,
                                "message": (
                                    f"{channel['id']} segment {segment_index} "
                                    f"passes within {distance:.4f} m of {bone}."
                                ),
                            }
                        )
                        break
    return findings


ATTACH_SURFACE_LIMIT_M = 0.040


def points_inside_mesh(points: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Ray-parity inside test for a closed triangle mesh."""

    direction = np.array([1.0, 0.0, 0.0])
    corner0 = triangles[:, 0]
    edge1 = triangles[:, 1] - corner0
    edge2 = triangles[:, 2] - corner0
    normal = np.cross(direction, edge2)
    determinant = np.einsum("ij,ij->i", edge1, normal)
    usable = np.abs(determinant) > 1e-12
    inverse = np.zeros_like(determinant)
    inverse[usable] = 1.0 / determinant[usable]

    inside = np.zeros(len(points), dtype=bool)
    for index, point in enumerate(points):
        offset = point - corner0
        u = inverse * np.einsum("ij,ij->i", offset, normal)
        candidate = usable & (u >= 0.0) & (u <= 1.0)
        cross = np.cross(offset, edge1)
        v = inverse * (cross @ direction)
        candidate &= (v >= 0.0) & (u + v <= 1.0)
        distance = inverse * np.einsum("ij,ij->i", edge2, cross)
        candidate &= distance > 1e-9
        inside[index] = int(np.count_nonzero(candidate)) % 2 == 1
    return inside


def attach_point_report(
    model: mujoco.MjModel,
    manifest: dict,
    export_dir: Path,
) -> list[dict]:
    """Distance from every wire anchor to the surface of its own bone.

    The anchors in ``cable_drive_manifest.json`` are offsets from a joint
    frame, so they only stay on the bone while the bone keeps the pose they
    were tuned in.  After the bones are re-posed the anchors must be checked
    again.
    """

    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return []

    trees = {}
    for bone in ("thorax", "clavicle", "scapula", "humerus", "ulna"):
        path = export_dir / "bones" / f"{bone}.stl"
        if not path.exists():
            continue
        triangles = read_binary_stl(path)
        if float(np.max(np.ptp(triangles.reshape(-1, 3), axis=0))) > 5.0:
            triangles = triangles * 0.001
        # Triangle corners alone overstate the distance on a coarse mesh: the
        # middle of a 100 mm thorax triangle is far from every corner.  Sample
        # the surface instead.
        corner0, corner1, corner2 = (
            triangles[:, 0],
            triangles[:, 1],
            triangles[:, 2],
        )
        points = np.vstack(
            [
                (corner0 + corner1 + corner2) / 3.0,
                (corner0 + corner1) / 2.0,
                (corner1 + corner2) / 2.0,
                (corner2 + corner0) / 2.0,
            ]
        )
        trees[bone] = cKDTree(points)

    data = mujoco.MjData(model)
    data.qpos[:] = 0.0
    mujoco.mj_forward(model, data)
    body_bone = body_to_bone()

    rows = []
    for channel in cable_channels(manifest):
        entries = [
            (
                f"cable_{channel['id']}_insertion",
                channel["insertion_body"],
                "insertion",
            )
        ]
        for index, guide in enumerate(channel["guides"]):
            entries.append((guide["site"], guide["body"], f"guide{index}"))
        for site_name, body, kind in entries:
            site_id = mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_SITE, site_name
            )
            if site_id < 0:
                continue
            bone = body_bone.get(body)
            if bone not in trees:
                continue
            point = data.site_xpos[site_id]
            distance = float(trees[bone].query(point)[0])
            rows.append(
                {
                    "cable": channel["id"],
                    "kind": kind,
                    "bone": bone,
                    "distance_to_bone_m": distance,
                    "on_surface": distance <= ATTACH_SURFACE_LIMIT_M,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_OUTPUT_XML)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT_DIR)
    parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_REPORT_DIR / "cable_drive_check.json",
    )
    args = parser.parse_args()

    manifest = load_cable_manifest(args.manifest)
    channels = cable_channels(manifest)
    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)

    report = {
        "model": str(args.model),
        "counts": {
            "bodies": int(model.nbody),
            "joints": int(model.njnt),
            "tendons": int(model.ntendon),
            "actuators": int(model.nu),
            "dof": int(model.nv),
            "cable_channels": len(channels),
        },
        "findings": [],
        "cables": {},
    }

    expected_joints = len(SKELETON_JOINT_NAMES) + len(channels)
    if int(model.njnt) != expected_joints:
        report["findings"].append(
            {
                "severity": "error",
                "message": (
                    f"Expected {expected_joints} joints "
                    f"({len(SKELETON_JOINT_NAMES)} skeleton + {len(channels)} "
                    f"sliders), found {int(model.njnt)}."
                ),
            }
        )

    sensitivities = sensitivity_matrix(model, channels)
    for channel in channels:
        values = {
            joint: sensitivities[joint][channel["id"]]
            for joint in SKELETON_JOINT_NAMES
        }
        dominant_joint = max(values, key=lambda joint: abs(values[joint]))
        dominant = values[dominant_joint]
        report["cables"][channel["id"]] = {
            "dominant_joint": dominant_joint,
            "dominant_dldq": dominant,
            "sensitivities": values,
        }
        if abs(dominant) < 1e-5:
            report["findings"].append(
                {
                    "severity": "error",
                    "message": f"{channel['id']} has no usable moment arm.",
                }
            )

    texts = []
    one_way_failures = 0
    max_force_failures = 0
    states = sample_states(model)
    for channel in channels:
        tid = tendon_id(model, channel["cable"])
        stiffness = float(model.tendon_stiffness[tid])
        springlength = float(model.tendon_lengthspring[tid, 0])
        damping = float(model.tendon_damping[tid])
        minimum_margin = float("inf")
        maximum_force = 0.0
        for state in states:
            data.qpos[:] = state
            mujoco.mj_forward(model, data)
            length = float(data.ten_length[tid])
            velocity = float(data.ten_velocity[tid])
            force = stiffness * (length - springlength) + damping * velocity
            minimum_margin = min(minimum_margin, length - springlength)
            maximum_force = max(maximum_force, force)
        if minimum_margin <= 0.0:
            one_way_failures += 1
        if maximum_force > channel["max_force_n"] * 1.25:
            max_force_failures += 1
        texts.append(
            {
                "id": channel["id"],
                "minimum_length_margin_m": minimum_margin,
                "maximum_sampled_force_n": maximum_force,
                "max_force_n": channel["max_force_n"],
                "springlength_m": springlength,
            }
        )
    report["tension"] = texts
    if one_way_failures:
        report["findings"].append(
            {
                # With spring_rest_mode="reference" the tendon is taut at the
                # neutral pose, so a pose that shortens a cable makes the
                # linear spring push.  Physically that cable would simply go
                # slack; the force should just be ignored there.
                "severity": "warning",
                "message": (
                    f"{one_way_failures} cable(s) go slack over the sampled "
                    "state range (the linear spring pushes instead of "
                    "releasing)."
                ),
            }
        )
    if max_force_failures:
        report["findings"].append(
            {
                "severity": "warning",
                "message": (
                    f"{max_force_failures} cable(s) exceed 125% of the nominal "
                    "maximum force over the sampled state range."
                ),
            }
        )

    report["findings"].extend(clearance_report(model, manifest, args.export_dir))

    attach_rows = attach_point_report(model, manifest, args.export_dir)
    report["attach_points"] = attach_rows
    floating = [
        row
        for row in attach_rows
        if not row["on_surface"] and row["kind"] == "insertion"
    ]
    if floating:
        worst = max(row["distance_to_bone_m"] for row in floating)
        report["findings"].append(
            {
                "severity": "warning",
                "message": (
                    f"{len(floating)} wire anchor(s) are more than "
                    f"{ATTACH_SURFACE_LIMIT_M * 1000:.0f} mm from their own "
                    f"bone (worst {worst * 1000:.0f} mm). Re-derive the "
                    "anchors for the current bone pose."
                ),
            }
        )

    data.qpos[:] = 0.0
    data.qvel[:] = 0.0
    data.ctrl[:] = 0.0
    for _ in range(1000):
        mujoco.mj_step(model, data)
    report["finite_state"] = bool(
        np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
    )
    if not report["finite_state"]:
        report["findings"].append(
            {"severity": "error", "message": "Simulation produced non-finite state."}
        )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"model: {args.model}")
    print(f"bodies: {report['counts']['bodies']}")
    print(f"joints: {report['counts']['joints']}")
    print(f"tendons: {report['counts']['tendons']}")
    print(f"actuators: {report['counts']['actuators']}")
    print(f"cable channels: {len(channels)}")
    for channel in channels:
        item = report["cables"][channel["id"]]
        print(
            f"  {channel['id']:28s} "
            f"dL/dq={item['dominant_dldq']: .6f} via {item['dominant_joint']}"
        )
    if report["findings"]:
        for finding in report["findings"]:
            print(f"[{finding['severity'].upper()}] {finding['message']}")
    else:
        print("findings: none")
    print(f"report: {args.report}")

    fatal = any(
        finding["severity"] == "error" for finding in report["findings"]
    )
    raise SystemExit(1 if fatal else 0)


if __name__ == "__main__":
    main()
