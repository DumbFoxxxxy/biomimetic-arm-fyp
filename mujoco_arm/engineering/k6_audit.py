"""Audit the isolated K6/S01 engineering baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mujoco
import numpy as np

from k6_common import (
    DEFAULT_ASSETS,
    DEFAULT_CONFIG,
    DEFAULT_REPORT,
    DEFAULT_XML,
    MUJOCO_ROOT,
    load_config,
    r400_center,
    to_mujoco_design,
)


JOINT_NAMES = ("sc_elevation", "sc_protraction", "ac_ext", "gh_flex")


def design_from_mujoco(point_m: np.ndarray) -> np.ndarray:
    point_m = np.asarray(point_m, dtype=float)
    return np.array([point_m[0], point_m[2], -point_m[1]], dtype=float) * 1000.0


def rotation_axis_angle(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    angle = float(
        np.arccos(np.clip((np.trace(matrix) - 1.0) / 2.0, -1.0, 1.0))
    )
    if angle < 1e-9:
        return np.zeros(3), 0.0
    axis = np.array(
        [
            matrix[2, 1] - matrix[1, 2],
            matrix[0, 2] - matrix[2, 0],
            matrix[1, 0] - matrix[0, 1],
        ],
        dtype=float,
    )
    return axis / (2.0 * np.sin(angle)), angle


def geom_id(model: mujoco.MjModel, name: str) -> int:
    result = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, name)
    if result < 0:
        raise KeyError(name)
    return result


def body_id(model: mujoco.MjModel, name: str) -> int:
    result = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, name)
    if result < 0:
        raise KeyError(name)
    return result


def geom_mesh_vertices(model: mujoco.MjModel, geom_name: str) -> np.ndarray:
    geom = geom_id(model, geom_name)
    mesh = int(model.geom_dataid[geom])
    if mesh < 0:
        raise ValueError(f"{geom_name} is not a mesh geom")
    start = int(model.mesh_vertadr[mesh])
    count = int(model.mesh_vertnum[mesh])
    return model.mesh_vert[start : start + count].copy()


class K6Evaluator:
    def __init__(
        self,
        model_path: Path | None,
        config: dict,
        model: mujoco.MjModel | None = None,
    ):
        self.model_path = model_path
        self.config = config
        if model is None:
            if model_path is None:
                raise ValueError("model_path or model is required")
            self.model = mujoco.MjModel.from_xml_path(str(model_path))
        else:
            self.model = model
        self.actuator_id = mujoco.mj_name2id(
            self.model, mujoco.mjtObj.mjOBJ_ACTUATOR, "pull_S01"
        )
        self.tendon_id = mujoco.mj_name2id(
            self.model, mujoco.mjtObj.mjOBJ_TENDON, "cable_S01"
        )
        if self.actuator_id < 0 or self.tendon_id < 0:
            raise RuntimeError("K6 model is missing pull_S01 or cable_S01")
        self.scapula_id = body_id(self.model, "scapula")
        self.gh_id = body_id(self.model, "gh_flex")
        self.pad_geom_id = geom_id(self.model, "scapula_contact_pad")
        self.band_geom_id = geom_id(self.model, "thorax_contact_band")
        self.visual_geom_id = geom_id(self.model, "scapula_visual")
        self.pad_vertices = geom_mesh_vertices(
            self.model, "scapula_contact_pad"
        )
        self.visual_vertices = geom_mesh_vertices(
            self.model, "scapula_visual"
        )
        self.band_centre = r400_center(config)
        self.band_radius = float(config["contact"]["thorax_band"]["radius_mm"])
        self.band_x = (
            float(config["contact"]["thorax_band"]["x_min_mm"]),
            float(config["contact"]["thorax_band"]["x_max_mm"]),
        )
        self.band_y = (
            float(config["contact"]["thorax_band"]["y_min_mm"]),
            float(config["contact"]["thorax_band"]["y_max_mm"]),
        )
        self.zero_data = mujoco.MjData(self.model)
        mujoco.mj_forward(self.model, self.zero_data)
        self.zero_scapula_rotation = self.zero_data.xmat[
            self.scapula_id
        ].reshape(3, 3).copy()
        self.neutral_length_m = float(self.zero_data.ten_length[self.tendon_id])
        self.w_axis_mj = np.array([0.6556, -0.7552, 0.0])
        self.w_axis_mj /= np.linalg.norm(self.w_axis_mj)

    def source_contact_gap(self) -> float:
        """Zero-pose radial gap of the pad source mesh to the R400 band."""

        return self.contact_gap_mm(self.zero_data)

    def simulate(self, pull_mm: float, steps: int | None = None) -> mujoco.MjData:
        data = mujoco.MjData(self.model)
        data.qpos[:] = 0.0
        data.qvel[:] = 0.0
        target = self.neutral_length_m - pull_mm * 0.001
        data.ctrl[self.actuator_id] = target
        mujoco.mj_forward(self.model, data)
        count = steps if steps is not None else int(self.config["simulation"]["pose_steps"])
        for _ in range(count):
            mujoco.mj_step(self.model, data)
        mujoco.mj_forward(self.model, data)
        return data

    def joint_values_deg(self, data: mujoco.MjData) -> dict[str, float]:
        result = {}
        for name in JOINT_NAMES:
            joint = mujoco.mj_name2id(
                self.model, mujoco.mjtObj.mjOBJ_JOINT, name
            )
            result[name] = float(
                np.degrees(data.qpos[int(self.model.jnt_qposadr[joint])])
            )
        return result

    def task_angle_deg(self, data: mujoco.MjData) -> float:
        rotation = (
            data.xmat[self.scapula_id].reshape(3, 3)
            @ self.zero_scapula_rotation.T
        )
        axis, angle = rotation_axis_angle(rotation)
        return float(np.degrees(angle * np.dot(axis, self.w_axis_mj)))

    def gh_design_mm(self, data: mujoco.MjData) -> np.ndarray:
        return design_from_mujoco(data.xpos[self.gh_id])

    def contact_gap_mm(self, data: mujoco.MjData) -> float:
        rotation = data.geom_xmat[self.pad_geom_id].reshape(3, 3)
        position = data.geom_xpos[self.pad_geom_id]
        world = (rotation @ self.pad_vertices.T).T + position
        design = np.asarray([design_from_mujoco(point) for point in world])
        valid = (
            (design[:, 0] >= self.band_x[0])
            & (design[:, 0] <= self.band_x[1])
            & (design[:, 1] >= self.band_y[0])
            & (design[:, 1] <= self.band_y[1])
        )
        if not np.any(valid):
            return float("inf")
        points = design[valid]
        radial = (
            np.linalg.norm(points[:, [0, 2]] - self.band_centre, axis=1)
            - self.band_radius
        )
        return float(np.min(radial))

    def thorax_signed_clearance_mm(self, points_design: np.ndarray) -> np.ndarray:
        points = np.asarray(points_design, dtype=float)
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        inside_height = (y >= 0.0) & (y <= float(self.config["thorax"]["height_mm"]))
        radial_pad = (
            np.linalg.norm(points[:, [0, 2]] - self.band_centre, axis=1)
            - self.band_radius
        )
        in_pad_sector = (
            (x >= self.band_x[0])
            & (x <= self.band_x[1])
            & (z >= 0.0)
        )
        semi_x = float(self.config["thorax"]["semi_x_mm"])
        semi_z = float(self.config["thorax"]["semi_z_mm"])
        ellipse_radius = np.sqrt((x / semi_x) ** 2 + (z / semi_z) ** 2)
        ellipse_clearance = (ellipse_radius - 1.0) * min(semi_x, semi_z)
        clearance = np.where(in_pad_sector, radial_pad, ellipse_clearance)
        clearance = np.where(inside_height, clearance, 10.0)
        return clearance

    def visual_clearance_mm(self, data: mujoco.MjData) -> float:
        rotation = data.geom_xmat[self.visual_geom_id].reshape(3, 3)
        position = data.geom_xpos[self.visual_geom_id]
        world = (rotation @ self.visual_vertices.T).T + position
        design = np.asarray([design_from_mujoco(point) for point in world])
        return float(np.min(self.thorax_signed_clearance_mm(design)))

    def cable_clearance_mm(self, data: mujoco.MjData) -> float:
        address = int(data.ten_wrapadr[self.tendon_id])
        count = int(data.ten_wrapnum[self.tendon_id])
        flat_wrap = data.wrap_xpos.reshape(-1, 3)
        path_points = flat_wrap[address : address + count]
        design = np.asarray([design_from_mujoco(point) for point in path_points])
        sampled = []
        for start, end in zip(design[:-1], design[1:]):
            length = float(np.linalg.norm(end - start))
            count_segment = max(2, int(np.ceil(length / 1.5)) + 1)
            for fraction in np.linspace(0.0, 1.0, count_segment):
                sampled.append(start + fraction * (end - start))
        if not sampled:
            return float("inf")
        return float(np.min(self.thorax_signed_clearance_mm(np.asarray(sampled))))

    def measure(self, pull_mm: float) -> dict:
        data = self.simulate(pull_mm)
        gap = self.contact_gap_mm(data)
        contacts = [
            contact
            for contact in data.contact
            if {
                contact.geom1,
                contact.geom2,
            }
            == {self.pad_geom_id, self.band_geom_id}
        ]
        return {
            "pull_mm": float(pull_mm),
            "joints_deg": self.joint_values_deg(data),
            "up_deg": self.task_angle_deg(data),
            "gh_mm": self.gh_design_mm(data).tolist(),
            "cable_length_mm": float(data.ten_length[self.tendon_id] * 1000.0),
            "contact_gap_mm": gap,
            "mujoco_contact_count": len(contacts),
            "mujoco_min_contact_mm": float(
                min((contact.dist for contact in contacts), default=0.0) * 1000.0
            ),
            "visual_clearance_mm": self.visual_clearance_mm(data),
            "cable_clearance_mm": self.cable_clearance_mm(data),
        }

    def sweep(self, max_pull_mm: float, increment_mm: float) -> list[dict]:
        values = []
        pull = 0.0
        while pull <= max_pull_mm + 1e-9:
            values.append(self.measure(pull))
            pull += increment_mm
        return values


def build_report(
    model_path: Path,
    config: dict,
    route_path: Path,
    max_pull_mm: float,
    increment_mm: float,
) -> dict:
    evaluator = K6Evaluator(model_path, config)
    samples = evaluator.sweep(max_pull_mm, increment_mm)
    target = config["target"]
    acceptance = config["acceptance"]
    zero_gap = evaluator.source_contact_gap()
    target_pose = next(
        (
            sample
            for sample in samples
            if sample["up_deg"]
            >= float(target["up_deg"]) - float(target["up_tolerance_deg"])
        ),
        None,
    )

    moment_arm_mm = None
    if target_pose is not None:
        target_pull = target_pose["pull_mm"]
        if target_pull > increment_mm:
            before = evaluator.measure(target_pull - increment_mm)
            after = evaluator.measure(target_pull + increment_mm)
            delta_angle = np.radians(after["up_deg"] - before["up_deg"])
            delta_length = after["cable_length_mm"] - before["cable_length_mm"]
            if abs(delta_angle) > 1e-9:
                moment_arm_mm = float(-delta_length / delta_angle)

    contraction_percent = None
    if target_pose is not None:
        contraction_percent = float(
            (evaluator.neutral_length_m * 1000.0 - target_pose["cable_length_mm"])
            / (evaluator.neutral_length_m * 1000.0)
            * 100.0
        )

    max_pose_jump = max(
        (
            abs(samples[index + 1]["up_deg"] - samples[index]["up_deg"])
            for index in range(len(samples) - 1)
        ),
        default=0.0,
    )
    max_length_jump = max(
        (
            abs(
                samples[index + 1]["cable_length_mm"]
                - samples[index]["cable_length_mm"]
            )
            for index in range(len(samples) - 1)
        ),
        default=0.0,
    )
    min_contact_gap = min(sample["contact_gap_mm"] for sample in samples)
    max_separation = max(
        max(0.0, sample["contact_gap_mm"]) for sample in samples
    )
    min_visual_clearance = min(
        sample["visual_clearance_mm"] for sample in samples
    )
    min_cable_clearance = min(
        sample["cable_clearance_mm"] for sample in samples
    )
    max_pad_penetration = max(0.0, -min_contact_gap)
    max_visual_penetration = max(0.0, -min_visual_clearance)
    max_cable_penetration = max(0.0, -min_cable_clearance)

    failures = []
    zero_gap_range = acceptance["zero_gap_mm"]
    if not zero_gap_range[0] <= zero_gap <= zero_gap_range[1]:
        failures.append(
            f"zero contact gap {zero_gap:.3f} mm outside {zero_gap_range}"
        )
    if target_pose is None:
        failures.append(
            f"target {target['up_deg']} deg not reached within {max_pull_mm} mm"
        )
    else:
        if target_pose["pull_mm"] > float(target["max_pull_mm"]) + 1e-9:
            failures.append("target reached beyond maximum pull")
        if contraction_percent is None or contraction_percent > float(
            target["max_contraction_percent"]
        ):
            failures.append(f"contraction {contraction_percent} percent exceeds limit")
        arm_range = target["moment_arm_mm"]
        if moment_arm_mm is None or not arm_range[0] <= moment_arm_mm <= arm_range[1]:
            failures.append(
                f"effective moment arm {moment_arm_mm} mm outside {arm_range}"
            )
    if max_pad_penetration > float(acceptance["max_pad_penetration_mm"]):
        failures.append(
            f"pad penetration {max_pad_penetration:.3f} mm exceeds limit"
        )
    if max_visual_penetration > float(acceptance["max_visual_penetration_mm"]):
        failures.append(
            f"visual penetration {max_visual_penetration:.3f} mm exceeds limit"
        )
    if max_cable_penetration > float(acceptance["max_cable_penetration_mm"]):
        failures.append(
            f"cable penetration {max_cable_penetration:.3f} mm exceeds limit"
        )
    if max_pose_jump > float(acceptance["max_pose_jump_deg"]):
        failures.append(f"pose jump {max_pose_jump:.3f} deg exceeds limit")
    if max_length_jump > float(acceptance["max_length_jump_mm"]):
        failures.append(f"cable length jump {max_length_jump:.3f} mm exceeds limit")

    route = json.loads(route_path.read_text(encoding="utf-8"))
    return {
        "model": str(model_path),
        "route": route,
        "neutral_length_mm": evaluator.neutral_length_m * 1000.0,
        "zero_contact_gap_mm": zero_gap,
        "target_up_deg": float(target["up_deg"]),
        "target_pose": target_pose,
        "contraction_percent": contraction_percent,
        "effective_moment_arm_mm": moment_arm_mm,
        "min_contact_gap_mm": min_contact_gap,
        "max_pad_penetration_mm": max_pad_penetration,
        "max_free_separation_mm": max_separation,
        "min_visual_clearance_mm": min_visual_clearance,
        "max_visual_penetration_mm": max_visual_penetration,
        "min_cable_clearance_mm": min_cable_clearance,
        "max_cable_penetration_mm": max_cable_penetration,
        "max_pose_jump_deg": max_pose_jump,
        "max_length_jump_mm": max_length_jump,
        "samples": samples,
        "failures": failures,
        "pass": not failures,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_XML)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--route",
        type=Path,
        default=DEFAULT_ASSETS / "route_metadata.json",
    )
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--max-pull-mm", type=float, default=None)
    parser.add_argument("--increment-mm", type=float, default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    max_pull = (
        float(config["target"]["max_pull_mm"])
        if args.max_pull_mm is None
        else args.max_pull_mm
    )
    increment = (
        float(config["simulation"]["pull_increment_mm"])
        if args.increment_mm is None
        else args.increment_mm
    )
    report = build_report(
        args.model,
        config,
        args.route,
        max_pull,
        increment,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        "zero_gap={zero_contact_gap_mm:.3f} mm "
        "target_up={target_up_deg:.1f} deg "
        "reached={reached} "
        "moment_arm={effective_moment_arm_mm} "
        "pad_pen={max_pad_penetration_mm:.3f} mm "
        "visual_pen={max_visual_penetration_mm:.3f} mm".format(
            **report,
            reached=report["target_pose"] is not None,
        )
    )
    print(f"report: {args.report}")
    if report["failures"]:
        print("FAIL")
        for failure in report["failures"]:
            print(" ", failure)
        raise SystemExit(1)
    print("PASS: K6/S01 engineering baseline")


if __name__ == "__main__":
    main()
