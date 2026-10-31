"""Trajectory interpolation and cable-length playback for V6."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

from dataclasses import dataclass
import json
from project_runtime import WorkspacePath as Path
import sys

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
V5_ROOT = ROOT.parent / "cable_design_v5"
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V5_ROOT), str(V2_ROOT)]

from evaluate_v2 import JOINT_NAMES, dense_tendon_jacobian  # noqa: E402
from v2_common import design_to_mj, mj_to_design  # noqa: E402
from v6_evaluator import V6Evaluator  # noqa: E402


@dataclass
class TrajectorySample:
    time_s: float
    q_target_deg: dict[str, float]
    wrist_target_mm: np.ndarray
    wrist_actual_mm: np.ndarray
    path_error_mm: float
    visual_clearance_mm: float
    contact_gap_mm: float
    tension_ok: bool
    relative_residual: float
    selected_cables: list[str]


def smoothstep(value: float) -> float:
    value = float(np.clip(value, 0.0, 1.0))
    return value * value * (3.0 - 2.0 * value)


def interpolate_waypoint_q(
    trajectory: dict,
    time_s: float,
) -> dict[str, float]:
    waypoints = sorted(
        trajectory["waypoints"],
        key=lambda item: float(item["time_s"]),
    )
    if time_s <= float(waypoints[0]["time_s"]):
        return {
            name: float(value)
            for name, value in waypoints[0].get("q_target_deg", {}).items()
        }
    if time_s >= float(waypoints[-1]["time_s"]):
        return {
            name: float(value)
            for name, value in waypoints[-1].get("q_target_deg", {}).items()
        }
    for start, end in zip(waypoints[:-1], waypoints[1:]):
        start_t = float(start["time_s"])
        end_t = float(end["time_s"])
        if start_t <= time_s <= end_t:
            fraction = (time_s - start_t) / max(end_t - start_t, 1e-9)
            blend = smoothstep(fraction)
            names = set(start.get("q_target_deg", {})) | set(
                end.get("q_target_deg", {})
            )
            return {
                name: (
                    float(start.get("q_target_deg", {}).get(name, 0.0))
                    * (1.0 - blend)
                    + float(end.get("q_target_deg", {}).get(name, 0.0))
                    * blend
                )
                for name in names
            }
    return {}


class TrajectoryRunner:
    def __init__(
        self,
        model_path: Path,
        config: dict,
        matrix: dict,
        library: dict,
        model: mujoco.MjModel | None = None,
    ):
        self.model_path = model_path
        self.config = config
        self.matrix = matrix
        self.library = library
        self.evaluator = V6Evaluator(
            model_path,
            config,
            matrix,
            model=model,
        )
        self.model = self.evaluator.model
        self.data = mujoco.MjData(self.model)
        self.reference_data = mujoco.MjData(self.model)
        self.actuator_by_cable = {}
        for index in range(self.model.nu):
            name = mujoco.mj_id2name(
                self.model,
                mujoco.mjtObj.mjOBJ_ACTUATOR,
                index,
            )
            if name and name.startswith("pull_"):
                self.actuator_by_cable[name[5:]] = index
        self.scapula_body = self.evaluator.scapula_body
        self.roller_site = self.evaluator.roller_site
        self.wrist_site = self.evaluator.wrist_site
        self.zero_pose_offsets = {
            name: value
            for name, value in self.evaluator.zero_pose_offsets.items()
        }
        self.control_period = float(
            config.get("simulation", {}).get("control_period_s", 0.02)
        )
        self.reset()

    def reset(self) -> None:
        self.data.qpos[:] = 0.0
        self.data.qvel[:] = 0.0
        for name, value in self.zero_pose_offsets.items():
            self.data.qpos[self.evaluator.qpos_indices[name]] = value
        self.data.ctrl[:] = 0.0
        if self.config["simulation"].get("actuator_type") == "position":
            for cable_id, actuator in self.actuator_by_cable.items():
                self.data.ctrl[actuator] = self.evaluator.neutral_lengths[cable_id]
        mujoco.mj_forward(self.model, self.data)

    def q_for_trajectory(
        self,
        trajectory: dict,
        time_s: float,
    ) -> tuple[np.ndarray, dict[str, float]]:
        relative_deg = interpolate_waypoint_q(trajectory, time_s)
        qpos = np.zeros(self.model.nq)
        full_deg = {name: 0.0 for name in JOINT_NAMES}
        for name, value in self.zero_pose_offsets.items():
            qpos[self.evaluator.qpos_indices[name]] = value
            full_deg[name] = float(np.degrees(value))
        for name, value in relative_deg.items():
            if name not in self.evaluator.qpos_indices:
                continue
            qpos[self.evaluator.qpos_indices[name]] += np.deg2rad(value)
            full_deg[name] += float(value)
        return qpos, full_deg

    def reference_state(self, qpos: np.ndarray) -> mujoco.MjData:
        self.reference_data.qpos[:] = qpos
        self.reference_data.qvel[:] = 0.0
        mujoco.mj_forward(self.model, self.reference_data)
        return self.reference_data

    def wrist_design_mm(self, data: mujoco.MjData) -> np.ndarray:
        return mj_to_design(data.site_xpos[self.wrist_site])

    def apply_shoulder_preload(self) -> None:
        if self.roller_site < 0:
            return
        force = self.evaluator.preload_direction_mj * self.evaluator.preload_n
        self.data.xfrc_applied[self.scapula_body, :3] = force
        offset = self.data.site_xpos[self.roller_site] - self.data.xpos[
            self.scapula_body
        ]
        self.data.xfrc_applied[self.scapula_body, 3:] = np.cross(offset, force)

    def command_cable_tensions(
        self,
        q_target: np.ndarray,
        qd_target: np.ndarray,
        active_cables: list[str],
    ) -> dict:
        active = set(active_cables)
        if self.config["simulation"].get("actuator_type") == "position":
            q_target = np.asarray(q_target, dtype=float)
            pretension_m = float(
                self.config["simulation"].get(
                    "position_pretension_mm", 0.1
                )
            ) * 0.001
            mode = self.config["simulation"].get(
                "position_mode", "jacobian"
            )
            if mode == "feedforward":
                reference = self.reference_state(q_target)
                for cable_id, actuator in self.actuator_by_cable.items():
                    if cable_id not in active:
                        self.data.ctrl[actuator] = 1.0
                        continue
                    self.data.ctrl[actuator] = float(
                        reference.ten_length[
                            self.evaluator.tendon_ids[cable_id]
                        ]
                    ) - pretension_m
            else:
                q_error = q_target - self.data.qpos
                jacobian = dense_tendon_jacobian(self.model, self.data)
                dof_indices = [
                    int(self.evaluator.dof_indices[name])
                    for name in JOINT_NAMES
                ]
                projection = jacobian[:, dof_indices] @ q_error[dof_indices]
                gain = float(
                    self.config["simulation"].get(
                        "position_jacobian_gain", 5.0
                    )
                )
                for cable_id, actuator in self.actuator_by_cable.items():
                    current_length = float(
                        self.data.ten_length[
                            self.evaluator.tendon_ids[cable_id]
                        ]
                    )
                    index = self.evaluator.tendon_indices[cable_id]
                    if (
                        cable_id not in active
                        or projection[index] >= -1e-6
                    ):
                        self.data.ctrl[actuator] = 1.0
                        continue
                    self.data.ctrl[actuator] = (
                        current_length
                        + gain * float(projection[index])
                        - pretension_m
                    )
            return {
                "ok": True,
                "solver": f"position_actuator_{mode}_length",
                "selected": list(active_cables),
                "relative_residual": 0.0,
                "tensions_N": {},
            }
        dof_indices = [
            int(self.evaluator.dof_indices[name]) for name in JOINT_NAMES
        ]
        q = np.asarray(
            [self.data.qpos[self.evaluator.qpos_indices[name]] for name in JOINT_NAMES]
        )
        qd = np.asarray(
            [self.data.qvel[self.evaluator.dof_indices[name]] for name in JOINT_NAMES]
        )
        q_target = np.asarray(
            [q_target[self.evaluator.qpos_indices[name]] for name in JOINT_NAMES]
        )
        qd_target = np.asarray(
            [qd_target[self.evaluator.dof_indices[name]] for name in JOINT_NAMES]
        )
        kp = float(
            self.config.get("simulation", {}).get(
                "position_kp_n_per_rad", 3.0
            )
        )
        kd = float(
            self.config.get("simulation", {}).get(
                "position_kd_nms_per_rad", 0.18
            )
        )
        tracking = kp * (q_target - q) + kd * (qd_target - qd)
        action = {
            "id": "TRAJECTORY_CONTROL",
            "allowed_cables": list(active_cables),
            "active_cables": list(active_cables),
            "mandatory_cables": [],
            "solve_all_dofs": True,
            "torque_tracking": tracking.tolist(),
        }
        solution = self.evaluator._solve_tensions(action, self.data)
        for cable_id, actuator in self.actuator_by_cable.items():
            if cable_id not in active:
                self.data.ctrl[actuator] = 0.0
                continue
            self.data.ctrl[actuator] = float(
                solution.get("tensions_N", {}).get(cable_id, 0.0)
            )
        return solution

    def run(
        self,
        trajectory_id: str,
        speed: float = 1.0,
        include_samples: bool = True,
    ) -> dict:
        trajectory = next(
            item
            for item in self.library["trajectories"]
            if item["id"] == trajectory_id
        )
        speed = float(np.clip(speed, 0.25, 1.0))
        self.reset()
        duration = float(trajectory["duration_s"])
        sample_period = float(
            self.config["acceptance"].get(
                "trajectory_sample_period_s",
                0.02,
            )
        )
        sample_times = np.arange(0.0, duration + 1e-9, sample_period)
        if sample_times[-1] < duration:
            sample_times = np.r_[sample_times, duration]
        active_cables = list(trajectory.get("active_cables", []))
        samples = []
        tension_action = {
            "id": f"TRAJ_{trajectory_id}",
            "allowed_cables": active_cables,
            "active_cables": active_cables,
            "mandatory_cables": [],
            "solve_all_dofs": True,
        }
        previous_q = self.data.qpos.copy()
        max_joint_step_deg = 0.0
        max_wrist_step_mm = 0.0
        previous_wrist = None
        for index, time_s in enumerate(sample_times):
            qpos, q_deg = self.q_for_trajectory(trajectory, float(time_s))
            reference = self.reference_state(qpos)
            wrist_target = self.wrist_design_mm(reference)
            q_next, _ = self.q_for_trajectory(
                trajectory,
                min(duration, float(time_s) + sample_period),
            )
            q_previous, _ = self.q_for_trajectory(
                trajectory,
                max(0.0, float(time_s) - sample_period),
            )
            qd_target = (q_next - q_previous) / (
                2.0 * sample_period
            )
            self.command_cable_tensions(
                qpos,
                qd_target,
                active_cables,
            )
            substeps = max(
                1,
                int(round((sample_period / speed) / self.model.opt.timestep)),
            )
            for _ in range(substeps):
                self.apply_shoulder_preload()
                mujoco.mj_step(self.model, self.data)
            mujoco.mj_forward(self.model, self.data)
            wrist_actual = self.wrist_design_mm(self.data)
            connection = self.evaluator._solve_tensions(
                tension_action,
                self.data,
            )
            joint_step = max(
                abs(float(value))
                for value in np.degrees(self.data.qpos - previous_q)
            )
            max_joint_step_deg = max(max_joint_step_deg, joint_step)
            if previous_wrist is not None:
                max_wrist_step_mm = max(
                    max_wrist_step_mm,
                    float(np.linalg.norm(wrist_actual - previous_wrist)),
                )
            previous_wrist = wrist_actual.copy()
            previous_q = self.data.qpos.copy()
            sample = TrajectorySample(
                time_s=float(time_s),
                q_target_deg=q_deg,
                wrist_target_mm=wrist_target,
                wrist_actual_mm=wrist_actual,
                path_error_mm=float(
                    np.linalg.norm(wrist_actual - wrist_target)
                ),
                visual_clearance_mm=self.evaluator._visual_clearance(
                    self.data
                ),
                contact_gap_mm=self.evaluator._contact_gap(self.data),
                tension_ok=bool(connection.get("ok", False)),
                relative_residual=float(
                    connection.get("relative_residual", float("inf"))
                ),
                selected_cables=list(connection.get("selected", [])),
            )
            samples.append(sample)
        errors = np.asarray(
            [sample.path_error_mm for sample in samples],
            dtype=float,
        )
        min_visual = min(
            sample.visual_clearance_mm for sample in samples
        )
        return {
            "id": trajectory_id,
            "label": trajectory["label"],
            "duration_s": duration,
            "playback_speed": speed,
            "rms_error_mm": float(np.sqrt(np.mean(errors**2))),
            "p95_error_mm": float(np.percentile(errors, 95)),
            "max_error_mm": float(np.max(errors)),
            "max_joint_step_deg": float(max_joint_step_deg),
            "max_wrist_step_mm": float(max_wrist_step_mm),
            "min_visual_clearance_mm": float(min_visual),
            "all_tensions_ok": all(sample.tension_ok for sample in samples),
            "max_relative_residual": max(
                sample.relative_residual for sample in samples
            ),
            "pass_path": bool(
                float(np.sqrt(np.mean(errors**2)))
                <= float(
                    self.config["acceptance"]["wrist_path_rms_mm_max"]
                )
                and float(np.percentile(errors, 95))
                <= float(
                    self.config["acceptance"]["wrist_path_p95_mm_max"]
                )
            ),
            "pass_clearance": bool(
                min_visual
                >= float(self.config["acceptance"]["full_mesh_clearance_mm_min"])
            ),
            "samples": [
                {
                    "time_s": sample.time_s,
                    "q_target_deg": sample.q_target_deg,
                    "wrist_target_mm": sample.wrist_target_mm.tolist(),
                    "wrist_actual_mm": sample.wrist_actual_mm.tolist(),
                    "path_error_mm": sample.path_error_mm,
                    "visual_clearance_mm": sample.visual_clearance_mm,
                    "contact_gap_mm": sample.contact_gap_mm,
                    "tension_ok": sample.tension_ok,
                    "relative_residual": sample.relative_residual,
                    "selected_cables": sample.selected_cables,
                }
                for sample in samples
            ]
            if include_samples
            else [],
        }

    def step_toward_time(
        self,
        trajectory: dict,
        time_s: float,
    ) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
        qpos, q_deg = self.q_for_trajectory(trajectory, time_s)
        reference = self.reference_state(qpos)
        wrist_target = self.wrist_design_mm(reference)
        active = list(trajectory.get("active_cables", []))
        q_next, _ = self.q_for_trajectory(
            trajectory,
            min(float(trajectory["duration_s"]), time_s + self.control_period),
        )
        q_previous, _ = self.q_for_trajectory(
            trajectory,
            max(0.0, time_s - self.control_period),
        )
        qd_target = (q_next - q_previous) / (2.0 * self.control_period)
        self.command_cable_tensions(qpos, qd_target, active)
        substeps = max(
            1,
            int(round(self.control_period / self.model.opt.timestep)),
        )
        for _ in range(substeps):
            self.apply_shoulder_preload()
            mujoco.mj_step(self.model, self.data)
        mujoco.mj_forward(self.model, self.data)
        return q_deg, wrist_target, self.wrist_design_mm(self.data)


def load_inputs() -> tuple[dict, dict, dict]:
    config = json.loads(
        (ROOT / "cable_design_v6.json").read_text(encoding="utf-8")
    )
    matrix = json.loads(
        (ROOT / "action_cable_matrix_v6.json").read_text(encoding="utf-8")
    )
    library = json.loads(
        (ROOT / "trajectory_library_v6.json").read_text(encoding="utf-8")
    )
    return config, matrix, library
