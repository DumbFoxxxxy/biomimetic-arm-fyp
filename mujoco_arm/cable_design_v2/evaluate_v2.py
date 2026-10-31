"""Evaluate actions, cable tensions, wrap geometry and route deliverables."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import csv
import itertools
import json
import math
from project_runtime import WorkspacePath as Path

import mujoco
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import linprog, lsq_linear

from v2_common import (
    ACTION_MATRIX_PATH,
    CONFIG_PATH,
    EXPORT_ROOT,
    MODEL_PATH,
    POINTS_MM,
    PROJECT_ROOT,
    design_to_mj,
    load_json,
    mj_to_design,
    normalize,
)


JOINT_NAMES = (
    "sc_elevation",
    "sc_protraction",
    "ac_ext",
    "gh_flex",
    "gh_abd",
    "gh_rot",
    "elbow",
)


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
    return axis / (2.0 * math.sin(angle)), angle


def scalar_joint(model: mujoco.MjModel, name: str) -> int:
    joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
    if joint < 0:
        raise KeyError(name)
    return int(model.jnt_qposadr[joint])


def scalar_object(model: mujoco.MjModel, object_type: mujoco.mjtObj, name: str) -> int:
    value = mujoco.mj_name2id(model, object_type, name)
    if value < 0:
        raise KeyError(name)
    return int(value)


def dense_tendon_jacobian(
    model: mujoco.MjModel, data: mujoco.MjData
) -> np.ndarray:
    jacobian = np.zeros((model.ntendon, model.nv), dtype=float)
    for tendon in range(model.ntendon):
        start = int(model.ten_J_rowadr[tendon])
        count = int(model.ten_J_rownnz[tendon])
        columns = model.ten_J_colind[start : start + count]
        jacobian[tendon, columns] = data.ten_J[start : start + count]
    return jacobian


def design_vector_to_mj(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return normalize(np.array([values[0], -values[2], values[1]]))


def r400_center(config: dict) -> np.ndarray:
    band = config["contact"]["thorax_band"]
    radius = float(band["radius_mm"])
    x0 = float(band["x_min_mm"])
    x1 = float(band["x_max_mm"])
    semi_x = 140.0
    semi_z = 106.0

    def ellipse_z(x: float) -> float:
        return semi_z * math.sqrt(max(0.0, 1.0 - (x / semi_x) ** 2))

    p0 = np.array([x0, ellipse_z(x0)])
    p1 = np.array([x1, ellipse_z(x1)])
    midpoint = 0.5 * (p0 + p1)
    chord = p1 - p0
    chord_length = float(np.linalg.norm(chord))
    perpendicular = np.array([-chord[1], chord[0]]) / chord_length
    height = math.sqrt(radius**2 - (chord_length / 2.0) ** 2)
    candidate_a = midpoint + height * perpendicular
    candidate_b = midpoint - height * perpendicular
    return candidate_a if np.linalg.norm(candidate_a) < np.linalg.norm(candidate_b) else candidate_b


def geom_mesh_vertices(model: mujoco.MjModel, geom_name: str) -> np.ndarray:
    geom = scalar_object(model, mujoco.mjtObj.mjOBJ_GEOM, geom_name)
    mesh = int(model.geom_dataid[geom])
    if mesh < 0:
        raise ValueError(f"{geom_name} is not a mesh")
    start = int(model.mesh_vertadr[mesh])
    count = int(model.mesh_vertnum[mesh])
    return model.mesh_vert[start : start + count].copy()


class V2Evaluator:
    def __init__(
        self,
        model_path: Path,
        config: dict,
        actions: dict,
        tpu_scale: float = 1.0,
        guide_friction: float | None = None,
        model: mujoco.MjModel | None = None,
    ):
        self.model_path = model_path
        self.config = config
        self.actions = actions
        self.tpu_scale = float(tpu_scale)
        self.guide_friction = (
            float(config["guides"]["friction_coefficient"])
            if guide_friction is None
            else float(guide_friction)
        )
        if model is None:
            self.model = _load_project_model(str(model_path))
        else:
            self.model = model
        self.qpos_indices = {
            name: scalar_joint(self.model, name) for name in JOINT_NAMES
        }
        self.dof_indices = {}
        for name in JOINT_NAMES:
            joint = scalar_object(
                self.model, mujoco.mjtObj.mjOBJ_JOINT, name
            )
            self.dof_indices[name] = int(self.model.jnt_dofadr[joint])

        self.route_by_id = {
            route["id"]: route for route in config["routes"]
        }
        self.tendon_ids = {
            route_id: scalar_object(
                self.model, mujoco.mjtObj.mjOBJ_TENDON, f"cable_{route_id}"
            )
            for route_id in self.route_by_id
        }
        self.tendon_indices = {
            route_id: index
            for index, route_id in enumerate(self.route_by_id)
        }
        self.scapula_body = scalar_object(
            self.model, mujoco.mjtObj.mjOBJ_BODY, "scapula"
        )
        self.gh_body = scalar_object(
            self.model, mujoco.mjtObj.mjOBJ_BODY, "gh_flex"
        )
        self.scapula_visual_geom = scalar_object(
            self.model, mujoco.mjtObj.mjOBJ_GEOM, "scapula_mesh"
        )
        self.contact_pad_geom = scalar_object(
            self.model, mujoco.mjtObj.mjOBJ_GEOM, "v2_scapula_contact_pad"
        )
        self.contact_band_geom = scalar_object(
            self.model, mujoco.mjtObj.mjOBJ_GEOM, "v2_thorax_contact_band"
        )
        self.scapula_visual_vertices = geom_mesh_vertices(
            self.model, "scapula_mesh"
        )
        self.contact_pad_vertices = geom_mesh_vertices(
            self.model, "v2_scapula_contact_pad"
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

        self.u_axis = normalize(POINTS_MM["SUP"] - POINTS_MM["INF"])
        self.w_axis = normalize(np.array([0.6556, 0.0, 0.7552]))
        self.slide_axis = normalize(np.array([0.7552, 0.0, -0.6556]))
        self.u_axis_mj = design_vector_to_mj(self.u_axis)
        self.w_axis_mj = design_vector_to_mj(self.w_axis)
        self.slide_axis_mj = design_vector_to_mj(self.slide_axis)

        zero = mujoco.MjData(self.model)
        zero.qpos[:] = 0.0
        mujoco.mj_forward(self.model, zero)
        self.zero_scapula_rotation = zero.xmat[
            self.scapula_body
        ].reshape(3, 3).copy()
        self.zero_contact_point = self._point_on_body(
            zero, self.scapula_body, POINTS_MM["P"]
        )
        self.neutral_lengths = {
            route_id: float(zero.ten_length[tendon_id])
            for route_id, tendon_id in self.tendon_ids.items()
        }

    def _guide_efficiency(self, route_id: str) -> float:
        route = self.route_by_id[route_id]
        total_wrap = sum(
            math.radians(float(element.get("wrap_angle_deg", 0.0)))
            for element in route["path"]
            if element["kind"] == "guide"
        )
        return math.exp(-self.guide_friction * total_wrap)

    def _point_on_body(
        self, data: mujoco.MjData, body_id: int, design_world_mm: np.ndarray
    ) -> np.ndarray:
        # At zero all body frames are aligned with the world frame.
        del data
        return design_to_mj(design_world_mm)

    def _action_pose(self, action: dict, fraction: float) -> mujoco.MjData:
        data = mujoco.MjData(self.model)
        data.qpos[:] = 0.0
        for joint_name, degrees in action.get("q_target_deg", {}).items():
            data.qpos[self.qpos_indices[joint_name]] = math.radians(
                float(degrees) * fraction
            )
        mujoco.mj_forward(self.model, data)
        return data

    def _scapula_rotation(self, data: mujoco.MjData) -> tuple[np.ndarray, float, np.ndarray]:
        rotation = (
            data.xmat[self.scapula_body].reshape(3, 3)
            @ self.zero_scapula_rotation.T
        )
        axis, angle = rotation_axis_angle(rotation)
        return rotation, angle, axis

    def _task_metric(self, action: dict, data: mujoco.MjData) -> float:
        metric = action["metric"]
        if metric == "joint_abs":
            return abs(
                float(
                    np.degrees(
                        data.qpos[self.qpos_indices[action["metric_joint"]]]
                    )
                )
            )
        rotation, angle, axis = self._scapula_rotation(data)
        up_deg = float(np.degrees(angle * np.dot(axis, self.w_axis_mj)))
        if metric == "scapula_up":
            return up_deg
        if metric == "scapula_slide":
            body_rotation = data.xmat[self.scapula_body].reshape(3, 3)
            point = body_rotation @ design_to_mj(POINTS_MM["P"] - POINTS_MM["AC"])
            point += data.xpos[self.scapula_body]
            return float(
                np.dot(point - self.zero_contact_point, self.slide_axis_mj)
                * 1000.0
            )
        if metric == "scapula_external":
            return float(np.degrees(angle * np.dot(axis, self.u_axis_mj)))
        if metric == "combined_lift":
            gh_flex = abs(
                float(np.degrees(data.qpos[self.qpos_indices["gh_flex"]]))
            )
            return gh_flex + up_deg
        if metric == "pectoral":
            values = [
                abs(float(np.degrees(data.qpos[self.qpos_indices[name]])))
                for name in ("gh_flex", "gh_abd", "gh_rot")
            ]
            return float(np.linalg.norm(values))
        raise KeyError(metric)

    def _action_achievement(self, action: dict, data: mujoco.MjData) -> dict:
        value = self._task_metric(action, data)
        if "target_mm" in action:
            target = abs(float(action["target_mm"]))
            achieved = abs(value)
            visible = float(action["minimum_visible_mm"])
            unit = "mm"
        else:
            target = abs(float(action["target_deg"]))
            achieved = abs(value)
            visible = float(action["minimum_visible_deg"])
            unit = "deg"
        acceptance = abs(
            float(action.get("acceptance_deg", action.get("acceptance_mm", target * 0.8)))
        )
        return {
            "value": float(value),
            "achieved": achieved,
            "target": target,
            "acceptance_80": acceptance,
            "minimum_visible": visible,
            "unit": unit,
            "fraction_of_target": achieved / target if target > 0 else 0.0,
            "meets_80_percent": achieved >= acceptance - 1e-6,
            "meets_minimum_visible": achieved >= visible - 1e-6,
        }

    def _tpu_forces(self, action: dict, data: mujoco.MjData) -> dict[str, float]:
        forces = {}
        for route_id in action["allowed_cables"]:
            route = self.route_by_id[route_id]
            tpu = route["tpu"]
            length = float(data.ten_length[self.tendon_ids[route_id]])
            neutral = self.neutral_lengths[route_id]
            strain = (neutral - length) / neutral if neutral > 0 else 0.0
            if not tpu.get("enabled", False) or strain <= 0.0:
                forces[route_id] = 0.0
                continue
            reference_strain = float(tpu["reference_strain"])
            reference_force = float(tpu["reference_force_N"])
            strain_limit = (
                float(self.config["acceptance"]["cable_contraction_percent_max"])
                / 100.0
            )
            effective_strain = min(strain, strain_limit)
            forces[route_id] = self.tpu_scale * reference_force * (
                effective_strain / reference_strain
            ) ** 2
        return forces

    def _dofs_for_action(self, action: dict) -> list[int]:
        if action["id"] == "K1":
            names = ["elbow"]
        elif action["id"] in ("K2", "K3"):
            names = [action["metric_joint"]]
        elif action["id"] in ("K4", "K5"):
            names = ["gh_rot"]
        elif action["id"] in ("K6", "SLIDE_FWD", "SLIDE_BACK", "SCAP_EXT", "SCAP_INT"):
            names = ["sc_elevation", "sc_protraction", "ac_ext"]
        elif action["id"] == "K7":
            names = [
                "sc_elevation",
                "sc_protraction",
                "ac_ext",
                "gh_flex",
            ]
        elif action["id"] == "PECTORAL":
            names = ["gh_flex", "gh_abd", "gh_rot"]
        else:
            names = ["gh_flex", "gh_abd", "gh_rot"]
        return [self.dof_indices[name] for name in names]

    def _solve_tensions(self, action: dict, data: mujoco.MjData) -> dict:
        allowed = list(action["allowed_cables"])
        tendon_columns = [self.tendon_indices[route_id] for route_id in allowed]
        dofs = self._dofs_for_action(action)
        jacobian = dense_tendon_jacobian(self.model, data)
        full_moment = -jacobian[:, dofs][tendon_columns].T
        efficiencies = np.asarray(
            [self._guide_efficiency(route_id) for route_id in allowed]
        )
        full_moment = full_moment * efficiencies[None, :]
        tpu_map = self._tpu_forces(action, data)
        bias = np.asarray(data.qfrc_bias)[dofs]
        if len(allowed) < 2:
            return {
                "ok": False,
                "reason": "fewer than two allowed cables",
                "allowed": allowed,
                "tensions_N": {},
                "tpu_N": tpu_map,
                "residual": None,
            }
        maximum = float(self.config["acceptance"]["tension_max_N"])
        bounds = (0.0, maximum)
        mandatory = list(action.get("mandatory_cables", []))
        positive_spanning = self._positive_spanning(full_moment, allowed)
        if mandatory:
            for route_id in mandatory:
                if route_id not in allowed:
                    return {
                        "ok": False,
                        "reason": f"mandatory cable {route_id} is not allowed",
                        "allowed": allowed,
                        "tensions_N": {},
                        "tpu_N": tpu_map,
                        "residual": None,
                        "positive_spanning": positive_spanning,
                    }
            mandatory_indices = [allowed.index(route_id) for route_id in mandatory]
            optional_indices = [
                index for index in range(len(allowed))
                if index not in mandatory_indices
            ]
            best = None
            for optional_size in range(0, min(3, len(optional_indices)) + 1):
                for optional in itertools.combinations(
                    optional_indices, optional_size
                ):
                    subset = tuple(mandatory_indices + list(optional))
                    moment = full_moment[:, subset]
                    selected_ids = [allowed[index] for index in subset]
                    tpu = np.asarray(
                        [tpu_map[route_id] for route_id in selected_ids]
                    )
                    lower = np.zeros(len(subset))
                    for route_id in mandatory:
                        lower[selected_ids.index(route_id)] = float(
                            action.get("minimum_mandatory_tension_N", 2.0)
                        )
                    target = bias - moment @ tpu
                    solution = lsq_linear(
                        moment,
                        target,
                        bounds=(lower, np.full(len(subset), maximum)),
                        lsmr_tol="auto",
                        max_iter=3000,
                    )
                    tensions = np.asarray(solution.x)
                    residual = moment @ tensions - target
                    residual_norm = float(np.linalg.norm(residual))
                    relative = residual_norm / (
                        float(np.linalg.norm(target)) + 1e-9
                    )
                    selected_routes = [
                        self._route_geometry(self.route_by_id[route_id], data)
                        for route_id in selected_ids
                    ]
                    max_contraction = max(
                        route["contraction_percent"]
                        for route in selected_routes
                    )
                    min_clearance = min(
                        route["thorax_clearance_mm"]
                        for route in selected_routes
                    )
                    spanning = self._positive_spanning(moment, selected_ids)
                    score = (
                        relative
                        + 0.05 * float(np.max(tensions, initial=0.0)) / maximum
                        + (0.0 if spanning["pass"] else 0.5)
                        + 0.03
                        * max(
                            0.0,
                            max_contraction
                            - float(
                                self.config["acceptance"][
                                    "cable_contraction_percent_max"
                                ]
                            ),
                        )
                        + 0.05
                        * max(
                            0.0,
                            float(
                                self.config["acceptance"][
                                    "cable_clearance_mm_min"
                                ]
                            )
                            - min_clearance,
                        )
                    )
                    if max_contraction > float(
                        self.config["acceptance"][
                            "cable_contraction_percent_max"
                        ]
                    ):
                        score += 100.0
                    if best is None or score < best["score"]:
                        best = {
                            "score": score,
                            "selected_ids": selected_ids,
                            "tensions": tensions,
                            "tpu": tpu,
                            "target": target,
                            "residual": residual,
                            "residual_norm": residual_norm,
                            "relative_residual": relative,
                            "positive_spanning": spanning,
                        }
            if best is None:
                return {
                    "ok": False,
                    "reason": "no mandatory-cable subset could be solved",
                    "allowed": allowed,
                    "tensions_N": {},
                    "tpu_N": tpu_map,
                    "residual": None,
                }
            tensions = best["tensions"]
            selected_ids = best["selected_ids"]
            tpu = best["tpu"]
            target = best["target"]
            residual = best["residual"]
            residual_norm = best["residual_norm"]
            relative = best["relative_residual"]
            positive_spanning = best["positive_spanning"]
            result = {
                "ok": bool(
                    np.all(np.isfinite(tensions))
                    and np.max(tensions, initial=0.0) < maximum - 1e-6
                    and relative < 0.5
                    and positive_spanning["pass"]
                ),
                "allowed": allowed,
                "selected": selected_ids,
                "mandatory": mandatory,
                "positive_spanning": positive_spanning,
                "tensions_N": {
                    route_id: float(tensions[selected_ids.index(route_id)])
                    for route_id in selected_ids
                },
                "tpu_N": {
                    route_id: float(tpu[selected_ids.index(route_id)])
                    for route_id in selected_ids
                },
                "all_tensions_N": {
                    route_id: float(tensions[index])
                    for index, route_id in enumerate(selected_ids)
                },
                "target_generalized_force": target.tolist(),
                "residual": residual.tolist(),
                "residual_norm": residual_norm,
                "relative_residual": relative,
                "max_tension_N": float(np.max(tensions, initial=0.0)),
                "score": best["score"],
            }
            return result
        best = None
        for size in (2, 3, 4):
            for subset in itertools.combinations(range(len(allowed)), size):
                moment = full_moment[:, subset]
                selected_ids = [allowed[index] for index in subset]
                tpu = np.asarray([tpu_map[route_id] for route_id in selected_ids])
                target = bias - moment @ tpu
                solution = lsq_linear(
                    moment,
                    target,
                    bounds=bounds,
                    lsmr_tol="auto",
                    max_iter=2000,
                )
                tensions = np.asarray(solution.x)
                residual = moment @ tensions - target
                residual_norm = float(np.linalg.norm(residual))
                relative = residual_norm / (float(np.linalg.norm(target)) + 1e-9)
                selected_routes = [
                    self._route_geometry(self.route_by_id[route_id], data)
                    for route_id in selected_ids
                ]
                max_contraction = max(
                    route["contraction_percent"] for route in selected_routes
                )
                min_clearance = min(
                    route["thorax_clearance_mm"] for route in selected_routes
                )
                score = (
                    relative
                    + 0.05 * float(np.max(tensions, initial=0.0)) / maximum
                    + 0.03
                    * max(
                        0.0,
                        max_contraction
                        - float(
                            self.config["acceptance"][
                                "cable_contraction_percent_max"
                            ]
                        ),
                    )
                    + 0.05
                    * max(
                        0.0,
                        float(self.config["acceptance"]["cable_clearance_mm_min"])
                        - min_clearance,
                    )
                )
                if max_contraction > float(
                    self.config["acceptance"]["cable_contraction_percent_max"]
                ):
                    score += 100.0
                if best is None or score < best["score"]:
                    best = {
                        "score": score,
                        "subset": subset,
                        "selected_ids": selected_ids,
                        "tensions": tensions,
                        "tpu": tpu,
                        "target": target,
                        "residual": residual,
                        "residual_norm": residual_norm,
                        "relative_residual": relative,
                    }
        if best is None:
            return {
                "ok": False,
                "reason": "no 2-4 cable subset could be solved",
                "allowed": allowed,
                "tensions_N": {},
                "tpu_N": tpu_map,
                "residual": None,
            }
        reduced_tensions = best["tensions"]
        result = {
            "ok": bool(
                np.all(np.isfinite(reduced_tensions))
                and np.max(reduced_tensions, initial=0.0) < maximum - 1e-6
            ),
            "allowed": allowed,
            "selected": best["selected_ids"],
            "tensions_N": {
                route_id: float(reduced_tensions[position])
                for position, route_id in enumerate(best["selected_ids"])
            },
            "tpu_N": {
                route_id: float(best["tpu"][position])
                for position, route_id in enumerate(best["selected_ids"])
            },
            "target_generalized_force": best["target"].tolist(),
            "residual": best["residual"].tolist(),
            "residual_norm": best["residual_norm"],
            "relative_residual": best["relative_residual"],
            "max_tension_N": float(np.max(reduced_tensions, initial=0.0)),
            "positive_spanning": positive_spanning,
        }
        result["ok"] = result["ok"] and result["relative_residual"] < 0.5
        return result

    @staticmethod
    def _positive_spanning(moment: np.ndarray, allowed: list[str]) -> dict:
        if moment.shape[1] < moment.shape[0] + 1:
            return {
                "pass": False,
                "reason": "fewer than n+1 columns",
                "minimum_singular_value": None,
            }
        directions = moment / np.maximum(
            np.linalg.norm(moment, axis=0, keepdims=True), 1e-9
        )
        singular = np.linalg.svd(directions, compute_uv=False)
        equality = np.vstack(
            (
                np.ones(directions.shape[1]),
                directions,
            )
        )
        rhs = np.zeros(directions.shape[0] + 1)
        rhs[0] = 1.0
        solution = linprog(
            np.zeros(directions.shape[1]),
            A_eq=equality,
            b_eq=rhs,
            bounds=(0.01, 1.0),
            method="highs",
        )
        return {
            "pass": bool(solution.success),
            "rank": int(np.linalg.matrix_rank(directions)),
            "minimum_singular_value": float(singular[-1]),
            "columns": allowed,
            "witness_coefficients": (
                solution.x.tolist() if solution.success else None
            ),
        }

    def _route_geometry(self, route: dict, data: mujoco.MjData) -> dict:
        tendon_id = self.tendon_ids[route["id"]]
        address = int(data.ten_wrapadr[tendon_id])
        count = int(data.ten_wrapnum[tendon_id])
        flat = data.wrap_xpos.reshape(-1, 3)
        points_mj = flat[address : address + count]
        points = np.asarray([mj_to_design(point) for point in points_mj])
        neutral = self.neutral_lengths[route["id"]]
        length = float(data.ten_length[tendon_id])
        guides = [
            element for element in route["path"] if element["kind"] == "guide"
        ]
        sampled = self._sample_polyline(points, spacing_mm=2.0)
        interior = sampled[2:-2] if len(sampled) > 5 else sampled
        cable_clearance = (
            float(np.min(self._thorax_clearance(interior)))
            if len(interior)
            else float("inf")
        )
        return {
            "id": route["id"],
            "length_mm": length * 1000.0,
            "neutral_length_mm": neutral * 1000.0,
            "contraction_percent": (neutral - length) / neutral * 100.0,
            "path_points_mm": points.tolist(),
            "guide_count": len(guides),
            "guide_radius_mm": [float(g["radius_mm"]) for g in guides],
            "declared_wrap_deg": [
                float(g.get("wrap_angle_deg", 0.0)) for g in guides
            ],
            "actual_wrap_points": int(count),
            "thorax_clearance_mm": cable_clearance,
            "sampled_path_mm": sampled.tolist(),
        }

    @staticmethod
    def _sample_polyline(points: np.ndarray, spacing_mm: float) -> np.ndarray:
        result = []
        for start, end in zip(points[:-1], points[1:]):
            distance = float(np.linalg.norm(end - start))
            count = max(2, int(math.ceil(distance / spacing_mm)) + 1)
            for fraction in np.linspace(0.0, 1.0, count):
                result.append(start + fraction * (end - start))
        return np.asarray(result) if result else np.empty((0, 3))

    def _thorax_clearance(self, points_design: np.ndarray) -> np.ndarray:
        points = np.asarray(points_design, dtype=float)
        if len(points) == 0:
            return np.empty(0)
        x, y, z = points.T
        radius_pad = (
            np.linalg.norm(points[:, [0, 2]] - self.band_centre, axis=1)
            - self.band_radius
        )
        in_band = (
            (x >= self.band_x[0])
            & (x <= self.band_x[1])
            & (z >= 0.0)
        )
        ellipse_radius = np.sqrt((x / 140.0) ** 2 + (z / 106.0) ** 2)
        clearance = np.where(in_band, radius_pad, (ellipse_radius - 1.0) * 106.0)
        clearance = np.where((y >= 0.0) & (y <= 280.0), clearance, 20.0)
        return clearance

    def _contact_gap(self, data: mujoco.MjData) -> float:
        rotation = data.geom_xmat[self.contact_pad_geom].reshape(3, 3)
        position = data.geom_xpos[self.contact_pad_geom]
        world = (rotation @ self.contact_pad_vertices.T).T + position
        design = np.asarray([mj_to_design(point) for point in world])
        valid = (
            (design[:, 0] >= self.band_x[0])
            & (design[:, 0] <= self.band_x[1])
            & (design[:, 1] >= self.band_y[0])
            & (design[:, 1] <= self.band_y[1])
        )
        if not np.any(valid):
            return float("inf")
        radial = (
            np.linalg.norm(design[valid][:, [0, 2]] - self.band_centre, axis=1)
            - self.band_radius
        )
        return float(np.min(radial))

    def _visual_clearance(self, data: mujoco.MjData) -> float:
        rotation = data.geom_xmat[self.scapula_visual_geom].reshape(3, 3)
        position = data.geom_xpos[self.scapula_visual_geom]
        world = (rotation @ self.scapula_visual_vertices.T).T + position
        design = np.asarray([mj_to_design(point) for point in world])
        return float(np.min(self._thorax_clearance(design)))

    def _cross_cable_minimum(self, routes: list[dict]) -> dict:
        best = None
        for first_index in range(len(routes)):
            first = routes[first_index]
            if len(first["sampled_path_mm"]) == 0:
                continue
            first_tree = cKDTree(np.asarray(first["sampled_path_mm"], dtype=float))
            for second_index in range(first_index + 1, len(routes)):
                second = routes[second_index]
                second_points = np.asarray(
                    second["sampled_path_mm"], dtype=float
                )
                if len(second_points) == 0:
                    continue
                distance, _ = first_tree.query(second_points, k=1)
                current = float(np.min(distance))
                if best is None or current < best[0]:
                    best = (current, first["id"], second["id"])
        if best is None:
            return {"minimum_mm": float("inf"), "between": []}
        return {"minimum_mm": best[0], "between": [best[1], best[2]]}

    def evaluate_action(self, action: dict, fraction: float) -> dict:
        data = self._action_pose(action, fraction)
        achievement = self._action_achievement(action, data)
        tension_solution = self._solve_tensions(action, data)
        routes = [
            self._route_geometry(self.route_by_id[route_id], data)
            for route_id in action["allowed_cables"]
        ]
        cross_cable = self._cross_cable_minimum(routes)
        for route in routes:
            route.pop("sampled_path_mm", None)
        return {
            "action": action["id"],
            "label": action["label"],
            "fraction": float(fraction),
            "qpos_deg": {
                name: float(np.degrees(data.qpos[self.qpos_indices[name]]))
                for name in JOINT_NAMES
            },
            "achievement": achievement,
            "tension_solution": tension_solution,
            "routes": routes,
            "contacts": int(data.ncon),
            "contact_gap_mm": self._contact_gap(data),
            "visual_clearance_mm": self._visual_clearance(data),
            "cross_cable": cross_cable,
        }

    def evaluate_all(self) -> dict:
        fractions = [
            float(value)
            for value in self.config["acceptance"]["action_pose_fractions"]
        ]
        action_reports = []
        for action in self.actions["actions"]:
            samples = [
                self.evaluate_action(action, fraction) for fraction in fractions
            ]
            final = samples[-1]
            selected = final["tension_solution"].get("selected", [])
            selected_routes = [
                route for route in final["routes"] if route["id"] in selected
            ]
            max_contraction = max(
                (
                    route["contraction_percent"]
                    for route in selected_routes
                ),
                default=float("-inf"),
            )
            min_cable_clearance = min(
                (
                    route["thorax_clearance_mm"]
                    for route in selected_routes
                ),
                default=float("inf"),
            )
            guide_count_ok = all(
                route["guide_count"] <= int(self.config["guides"]["maximum_per_cable"])
                and all(
                    radius >= float(self.config["guides"]["minimum_radius_mm"])
                    for radius in route["guide_radius_mm"]
                )
                and (
                    route["guide_count"] == 0
                    or route["actual_wrap_points"] >= 4
                )
                for route in selected_routes
            )
            contraction_ok = (
                max_contraction
                <= float(self.config["acceptance"]["cable_contraction_percent_max"])
            )
            cable_clearance_ok = (
                min_cable_clearance
                >= float(self.config["acceptance"]["cable_clearance_mm_min"])
            )
            visual_ok = (
                final["visual_clearance_mm"]
                >= -float(
                    self.config["acceptance"][
                        "unintended_mesh_penetration_mm_max"
                    ]
                )
            )
            action_reports.append(
                {
                    "id": action["id"],
                    "label": action["label"],
                    "group": action["group"],
                    "samples": samples,
                    "final": final,
                    "selected_cables": selected,
                    "cooperating_cable_count": len(selected),
                    "pass_minimum_visible": final["achievement"][
                        "meets_minimum_visible"
                    ],
                    "pass_80_percent": final["achievement"]["meets_80_percent"],
                    "pass_tension": bool(
                        final["tension_solution"].get("ok", False)
                    ),
                    "max_selected_contraction_percent": max_contraction,
                    "pass_contraction": contraction_ok,
                    "min_selected_cable_clearance_mm": min_cable_clearance,
                    "pass_cable_clearance": cable_clearance_ok,
                    "pass_guide_geometry": guide_count_ok,
                    "visual_clearance_mm": final["visual_clearance_mm"],
                    "pass_visual_clearance": visual_ok,
                    "pass_minimum_functional": bool(
                        final["achievement"]["meets_minimum_visible"]
                        and final["tension_solution"].get("ok", False)
                        and contraction_ok
                        and cable_clearance_ok
                        and guide_count_ok
                        and visual_ok
                    ),
                }
            )
        return {
            "model": str(self.model_path),
            "neutral_lengths_mm": {
                route_id: length * 1000.0
                for route_id, length in self.neutral_lengths.items()
            },
            "actions": action_reports,
        }


def export_action_matrix(path: Path, evaluation: dict, base_matrix: dict) -> dict:
    by_id = {action["id"]: action for action in evaluation["actions"]}
    output = json.loads(json.dumps(base_matrix))
    for action in output["actions"]:
        result = by_id[action["id"]]
        action["selected_cables"] = result["selected_cables"]
        action["cooperating_cable_count"] = result["cooperating_cable_count"]
        action["target_tensions_N"] = result["final"]["tension_solution"].get(
            "tensions_N", {}
        )
        action["pass_minimum_visible"] = result["pass_minimum_visible"]
        action["pass_80_percent"] = result["pass_80_percent"]
        action["pass_tension"] = result["pass_tension"]
        action["positive_spanning"] = result["final"][
            "tension_solution"
        ].get("positive_spanning")
        action["max_selected_contraction_percent"] = result[
            "max_selected_contraction_percent"
        ]
        action["pass_contraction"] = result["pass_contraction"]
        action["min_selected_cable_clearance_mm"] = result[
            "min_selected_cable_clearance_mm"
        ]
        action["pass_cable_clearance"] = result["pass_cable_clearance"]
        action["pass_guide_geometry"] = result["pass_guide_geometry"]
        action["pass_visual_clearance"] = result["pass_visual_clearance"]
        action["pass_minimum_functional"] = result["pass_minimum_functional"]
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    return output


def export_routes_csv(path: Path, config: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "cable_id",
        "group",
        "role",
        "index",
        "kind",
        "body",
        "world_x_mm",
        "world_y_mm",
        "world_z_mm",
        "local_x_mm",
        "local_y_mm",
        "local_z_mm",
        "guide_radius_mm",
        "axis_x",
        "axis_y",
        "axis_z",
        "wrap_angle_deg",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for route in config["routes"]:
            for index, element in enumerate(route["path"]):
                world = element["world_mm"]
                local = element["local_mm"]
                axis = element.get("axis_design", [0.0, 0.0, 0.0])
                writer.writerow(
                    {
                        "cable_id": route["id"],
                        "group": route["group"],
                        "role": route.get("role", ""),
                        "index": index,
                        "kind": element["kind"],
                        "body": element["body"],
                        "world_x_mm": world[0],
                        "world_y_mm": world[1],
                        "world_z_mm": world[2],
                        "local_x_mm": local[0] * 1000.0,
                        "local_y_mm": local[1] * 1000.0,
                        "local_z_mm": local[2] * 1000.0,
                        "guide_radius_mm": element.get("radius_mm", ""),
                        "axis_x": axis[0],
                        "axis_y": axis[1],
                        "axis_z": axis[2],
                        "wrap_angle_deg": element.get("wrap_angle_deg", ""),
                    }
                )


def export_solidworks_csv(path: Path, config: dict) -> None:
    fields = [
        "name",
        "part",
        "kind",
        "x_mm",
        "y_mm",
        "z_mm",
        "normal_x",
        "normal_y",
        "normal_z",
        "radius_mm",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for route in config["routes"]:
            for index, element in enumerate(route["path"]):
                normal = element.get(
                    "side_design", element.get("installation_normal", [0.0, 0.0, 1.0])
                )
                world = element["world_mm"]
                writer.writerow(
                    {
                        "name": f"{route['id']}_{element['kind']}_{index:02d}",
                        "part": element["body"],
                        "kind": element["kind"],
                        "x_mm": world[0],
                        "y_mm": world[1],
                        "z_mm": world[2],
                        "normal_x": normal[0],
                        "normal_y": normal[1],
                        "normal_z": normal[2],
                        "radius_mm": element.get("radius_mm", ""),
                    }
                )


def export_markdown(path: Path, evaluation: dict, config: dict) -> None:
    lines = [
        "# Cable anchors V2",
        "",
        "Generated from the force-aware V2 simulation. This file is provisional "
        "until every staged action passes.",
        "",
        "| Action | Selected cables | 80% | Minimum motion | Tension solution |",
        "| --- | --- | --- | --- | --- |",
    ]
    for action in evaluation["actions"]:
        lines.append(
            f"| {action['id']} | {', '.join(action['selected_cables'])} | "
            f"{'PASS' if action['pass_80_percent'] else 'FAIL'} | "
            f"{'PASS' if action['pass_minimum_visible'] else 'FAIL'} | "
            f"{'PASS' if action['pass_tension'] else 'FAIL'} |"
        )
    lines.extend(
        [
            "",
            "## Route coordinates",
            "",
            "| Cable | Role | Element | Part | World (x,y,z) mm | Radius mm | Wrap deg |",
            "| --- | --- | --- | --- | --- | ---: | ---: |",
        ]
    )
    for route in config["routes"]:
        for index, element in enumerate(route["path"]):
            world = element["world_mm"]
            lines.append(
                f"| {route['id']} | {route.get('role', '')} | "
                f"{index}:{element['kind']} | "
                f"{element['body']} | ({world[0]:.2f}, {world[1]:.2f}, "
                f"{world[2]:.2f}) | {element.get('radius_mm', '')} | "
                f"{element.get('wrap_angle_deg', '')} |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_action_markdown(path: Path, matrix: dict) -> None:
    lines = [
        "# Action-to-cable responsibility matrix V2",
        "",
        "| Action | Allowed | Selected | 80% | Minimum motion | Tension |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for action in matrix["actions"]:
        lines.append(
            f"| {action['id']} | {', '.join(action.get('allowed_cables', []))} | "
            f"{', '.join(action.get('selected_cables', []))} | "
            f"{'PASS' if action.get('pass_80_percent') else 'FAIL'} | "
            f"{'PASS' if action.get('pass_minimum_visible') else 'FAIL'} | "
            f"{'PASS' if action.get('pass_tension') else 'FAIL'} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_tension_report(path: Path, evaluation: dict) -> None:
    actions = []
    for action in evaluation["actions"]:
        final = action["final"]
        result = final["tension_solution"]
        actions.append(
            {
                "id": action["id"],
                "label": action["label"],
                "selected_cables": action["selected_cables"],
                "target_tensions_N": result.get("tensions_N", {}),
                "max_tension_N": result.get("max_tension_N"),
                "relative_residual": result.get("relative_residual"),
                "pass_tension": action["pass_tension"],
                "pass_80_percent": action["pass_80_percent"],
                "pass_minimum_visible": action["pass_minimum_visible"],
                "pass_minimum_functional": action["pass_minimum_functional"],
                "max_contraction_percent": action[
                    "max_selected_contraction_percent"
                ],
                "min_cable_clearance_mm": action[
                    "min_selected_cable_clearance_mm"
                ],
                "visual_clearance_mm": action["visual_clearance_mm"],
                "cross_cable": final["cross_cable"],
                "routes": [
                    {
                        "id": route["id"],
                        "length_mm": route["length_mm"],
                        "neutral_length_mm": route["neutral_length_mm"],
                        "contraction_percent": route["contraction_percent"],
                        "guide_count": route["guide_count"],
                        "actual_wrap_points": route["actual_wrap_points"],
                        "thorax_clearance_mm": route["thorax_clearance_mm"],
                    }
                    for route in final["routes"]
                ],
            }
        )
    path.write_text(
        json.dumps(
            {
                "model": evaluation["model"],
                "neutral_lengths_mm": evaluation["neutral_lengths_mm"],
                "actions": actions,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def export_route_validation(path: Path, config: dict, model_path: Path) -> None:
    model = _load_project_model(str(model_path))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    rows = []
    for route in config["routes"]:
        guides = [
            element for element in route["path"] if element["kind"] == "guide"
        ]
        tendon_id = scalar_object(
            model, mujoco.mjtObj.mjOBJ_TENDON, f"cable_{route['id']}"
        )
        rows.append(
            {
                "id": route["id"],
                "group": route["group"],
                "elements": len(route["path"]),
                "guide_count": len(guides),
                "guide_radii_mm": [float(g["radius_mm"]) for g in guides],
                "declared_wrap_deg": [
                    float(g.get("wrap_angle_deg", 0.0)) for g in guides
                ],
                "initial_length_mm": float(data.ten_length[tendon_id] * 1000.0),
                "actual_wrap_points": int(data.ten_wrapnum[tendon_id]),
                "low_wrap_pass_throughs": sum(
                    element["kind"] == "pass_through"
                    for element in route["path"]
                ),
            }
        )
    path.write_text(
        json.dumps({"routes": rows}, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument("--matrix", type=Path, default=ACTION_MATRIX_PATH)
    parser.add_argument("--out", type=Path, default=EXPORT_ROOT / "evaluation_v2.json")
    args = parser.parse_args()
    config = load_json(args.config)
    actions = load_json(args.matrix)
    evaluator = V2Evaluator(args.model, config, actions)
    evaluation = evaluator.evaluate_all()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(evaluation, indent=2) + "\n", encoding="utf-8")
    matrix = export_action_matrix(ACTION_MATRIX_PATH, evaluation, actions)
    export_action_markdown(
        PROJECT_ROOT / "design" / "ACTION_CABLE_MATRIX_V2.md",
        matrix,
    )
    export_routes_csv(EXPORT_ROOT / "cable_routes_v2.csv", config)
    export_solidworks_csv(EXPORT_ROOT / "SolidWorks_hardpoints_v2.csv", config)
    export_tension_report(EXPORT_ROOT / "tension_report.json", evaluation)
    export_route_validation(
        EXPORT_ROOT / "route_validation_v2.json",
        config,
        args.model,
    )
    export_markdown(
        PROJECT_ROOT / "design" / "CABLE_ANCHORS_V2.md",
        evaluation,
        config,
    )
    print(f"wrote {args.out}")
    for action in evaluation["actions"]:
        print(
            f"{action['id']:<10} selected={','.join(action['selected_cables']):<20} "
            f"80={action['pass_80_percent']} min={action['pass_minimum_visible']} "
            f"tension={action['pass_tension']}"
        )


if __name__ == "__main__":
    main()
