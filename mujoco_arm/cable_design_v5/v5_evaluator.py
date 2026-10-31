"""V5 evaluator with hard tension and contact constraints."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import itertools
import math
from project_runtime import WorkspacePath as Path
import sys

import mujoco
import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path.insert(0, str(V2_ROOT))

from evaluate_v2 import (  # noqa: E402
    JOINT_NAMES,
    V2Evaluator,
    dense_tendon_jacobian,
)
from v2_common import (  # noqa: E402
    POINTS_MM,
    design_to_mj,
    mj_to_design,
)


class V5Evaluator(V2Evaluator):
    """Force-aware evaluator with V5 pose and solver semantics."""

    def __init__(
        self,
        model_path: Path,
        config: dict,
        actions: dict,
        tpu_scale: float = 1.0,
        guide_friction: float | None = None,
        model: mujoco.MjModel | None = None,
    ):
        self.zero_pose_offsets = {
            name: math.radians(float(value))
            for name, value in config.get("zero_pose_offsets_deg", {}).items()
        }
        self.path_corrections = config.get("path_corrections", {})
        super().__init__(
            model_path,
            config,
            actions,
            tpu_scale=tpu_scale,
            guide_friction=guide_friction,
            model=model,
        )
        self._reset_zero_reference()

    def _point_on_body(
        self, data: mujoco.MjData, body_id: int, design_world_mm: np.ndarray
    ) -> np.ndarray:
        if body_id == self.scapula_body:
            local = design_to_mj(
                np.asarray(design_world_mm, dtype=float) - POINTS_MM["AC"]
            )
            rotation = data.xmat[body_id].reshape(3, 3)
            return rotation @ local + data.xpos[body_id]
        return design_to_mj(design_world_mm)

    def _zero_data(self) -> mujoco.MjData:
        data = mujoco.MjData(self.model)
        data.qpos[:] = 0.0
        for joint_name, value in self.zero_pose_offsets.items():
            data.qpos[self.qpos_indices[joint_name]] = value
        mujoco.mj_forward(self.model, data)
        return data

    def _reset_zero_reference(self) -> None:
        zero = self._zero_data()
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

    def _correction_at(self, action_id: str, fraction: float) -> np.ndarray:
        corrections = self.path_corrections.get(action_id)
        if not corrections:
            return np.zeros(3)
        knots = sorted(float(value) for value in corrections)
        values = np.asarray(
            [corrections[str(knot)] if str(knot) in corrections else corrections[knot]
             for knot in knots],
            dtype=float,
        )
        return np.asarray(
            [
                np.interp(float(fraction), knots, values[:, index])
                for index in range(3)
            ]
        )

    def _action_pose(self, action: dict, fraction: float) -> mujoco.MjData:
        data = mujoco.MjData(self.model)
        data.qpos[:] = 0.0
        for joint_name, value in self.zero_pose_offsets.items():
            data.qpos[self.qpos_indices[joint_name]] = value
        for joint_name, degrees in action.get("q_target_deg", {}).items():
            data.qpos[self.qpos_indices[joint_name]] += math.radians(
                float(degrees) * fraction
            )
        correction = self._correction_at(action["id"], fraction)
        for joint_name, delta_deg in zip(
            ("sc_elevation", "sc_protraction", "ac_ext"), correction
        ):
            data.qpos[self.qpos_indices[joint_name]] += math.radians(delta_deg)
        mujoco.mj_forward(self.model, data)
        return data

    def _route_geometry(self, route: dict, data: mujoco.MjData) -> dict:
        result = super()._route_geometry(route, data)
        guides = []
        factor = 1.0
        mu = float(self.config.get("friction_model", {}).get("mu", self.guide_friction))
        for element in route["path"]:
            if element["kind"] != "guide":
                continue
            wrap_deg = float(element.get("wrap_angle_deg", 0.0))
            local_factor = math.exp(-mu * math.radians(wrap_deg))
            factor *= local_factor
            guides.append(
                {
                    "name": element.get("name", "guide"),
                    "radius_mm": float(element.get("radius_mm", 0.0)),
                    "wrap_deg": wrap_deg,
                    "incoming_factor": local_factor,
                    "cumulative_factor": factor,
                }
            )
        result["capstan_guides"] = guides
        result["capstan_efficiency"] = factor
        return result

    def _guide_efficiency(self, route_id: str) -> float:
        route = self.route_by_id[route_id]
        mu = float(
            self.config.get("friction_model", {}).get(
                "mu", self.guide_friction
            )
        )
        total_wrap = sum(
            math.radians(float(element.get("wrap_angle_deg", 0.0)))
            for element in route["path"]
            if element["kind"] == "guide"
        )
        return math.exp(-mu * total_wrap)

    @staticmethod
    def _cone_violation(moment: np.ndarray) -> float:
        if moment.shape[1] < moment.shape[0] + 1:
            return 10.0
        directions = moment / np.maximum(
            np.linalg.norm(moment, axis=0, keepdims=True), 1e-12
        )
        dofs = directions.shape[0]
        columns = directions.shape[1]
        variables = columns + 2 * dofs
        objective = np.zeros(variables)
        objective[columns:] = 1.0
        equality = np.zeros((dofs + 1, variables))
        equality[:dofs, :columns] = directions
        equality[:dofs, columns : columns + dofs] = np.eye(dofs)
        equality[:dofs, columns + dofs :] = -np.eye(dofs)
        equality[dofs, :columns] = 1.0
        rhs = np.zeros(dofs + 1)
        rhs[-1] = 1.0
        bounds = [
            *[(0.0, 1.0) for _ in range(columns)],
            *[(0.0, 10.0) for _ in range(2 * dofs)],
        ]
        solution = linprog(
            objective,
            A_eq=equality,
            b_eq=rhs,
            bounds=bounds,
            method="highs",
        )
        if not solution.success:
            return 10.0
        return float(np.sum(solution.x[columns:]))

    @staticmethod
    def _positive_spanning(moment: np.ndarray, allowed: list[str]) -> dict:
        dofs = int(moment.shape[0])
        columns = int(moment.shape[1])
        if columns < dofs + 1:
            return {
                "pass": False,
                "reason": "fewer_than_dofs_plus_one",
                "rank": int(np.linalg.matrix_rank(moment)),
                "minimum_singular_value": None,
                "columns": allowed,
            }
        singular = np.linalg.svd(moment, compute_uv=False)
        rank = int(np.linalg.matrix_rank(moment))
        if rank < dofs:
            return {
                "pass": False,
                "reason": "rank_deficient",
                "rank": rank,
                "minimum_singular_value": float(singular[-1]),
                "columns": allowed,
            }
        directions = moment / np.maximum(
            np.linalg.norm(moment, axis=0, keepdims=True), 1e-12
        )
        equality = np.vstack((np.ones(columns), directions))
        rhs = np.zeros(dofs + 1)
        rhs[0] = 1.0
        solution = linprog(
            np.zeros(columns),
            A_eq=equality,
            b_eq=rhs,
            bounds=(1e-6, 1.0),
            method="highs",
        )
        return {
            "pass": bool(solution.success),
            "reason": "ok" if solution.success else "cone_intersection_empty",
            "rank": rank,
            "minimum_singular_value": float(singular[-1]),
            "columns": allowed,
            "witness_coefficients": (
                solution.x.tolist() if solution.success else None
            ),
        }

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
        maximum = float(self.config["acceptance"]["tension_max_N"])
        mandatory = list(action.get("mandatory_cables", []))
        minimum_mandatory = float(
            action.get(
                "minimum_mandatory_tension_N",
                self.config.get("tension_solver", {}).get(
                    "minimum_mandatory_tension_N", 2.0
                ),
            )
        )
        residual_tolerance = float(
            self.config.get("tension_solver", {}).get(
                "exact_equilibrium_relative_tolerance", 0.05
            )
        )
        if any(route_id not in allowed for route_id in mandatory):
            return {
                "ok": False,
                "reason": "mandatory_cable_not_allowed",
                "allowed": allowed,
                "selected": [],
                "tension_solution": "hard_constraint_failed",
            }

        def attempt(indices: tuple[int, ...]) -> tuple[dict | None, str]:
            moment = full_moment[:, indices]
            selected_ids = [allowed[index] for index in indices]
            spanning = self._positive_spanning(moment, selected_ids)
            if not spanning["pass"]:
                return None, f"{','.join(selected_ids)}:{spanning['reason']}"
            tpu = np.asarray([tpu_map[route_id] for route_id in selected_ids])
            target = bias - moment @ tpu
            lower = np.zeros(len(indices))
            for route_id in mandatory:
                lower[selected_ids.index(route_id)] = minimum_mandatory
            count = len(indices)
            objective = np.ones(count)
            solution = linprog(
                objective,
                A_eq=moment,
                b_eq=target,
                bounds=[(float(lower[index]), maximum) for index in range(count)],
                method="highs",
            )
            if not solution.success:
                return None, f"{','.join(selected_ids)}:no_bounded_tension"
            tensions = np.asarray(solution.x)
            residual = moment @ tensions - target
            relative = float(np.linalg.norm(residual)) / (
                float(np.linalg.norm(target)) + 1e-9
            )
            if relative > residual_tolerance:
                return None, f"{','.join(selected_ids)}:equilibrium_residual"
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
            if max_contraction > float(
                self.config["acceptance"]["cable_contraction_percent_max"]
            ):
                return None, f"{','.join(selected_ids)}:contraction"
            if min_clearance < float(
                self.config["acceptance"]["cable_clearance_mm_min"]
            ):
                return None, f"{','.join(selected_ids)}:clearance"
            score = (
                float(np.sum(tensions)) / maximum
                + 0.2 * float(np.max(tensions)) / maximum
                + 0.01 * count
            )
            return (
                {
                    "score": score,
                    "selected": selected_ids,
                    "tensions": tensions,
                    "tpu": tpu,
                    "target": target,
                    "residual": residual,
                    "relative_residual": relative,
                    "positive_spanning": spanning,
                    "max_contraction_percent": max_contraction,
                    "min_cable_clearance_mm": min_clearance,
                },
                "",
            )

        failures = []
        mandatory_indices = [allowed.index(route_id) for route_id in mandatory]
        preferred = [
            route_id
            for route_id in action.get("preferred_cables", allowed)
            if route_id in allowed
        ]
        preferred_indices = [allowed.index(route_id) for route_id in preferred]
        preferred_optional = [
            index
            for index in preferred_indices
            if index not in mandatory_indices
        ]
        optional_indices = [
            index for index in range(len(allowed)) if index not in mandatory_indices
        ]
        minimum_subset = len(dofs) + 1
        maximum_optional = int(
            self.config.get("tension_solver", {}).get(
                "maximum_optional_subset", 3
            )
        )
        maximum_subset = min(len(allowed), len(dofs) + maximum_optional)
        minimum_optional = max(0, minimum_subset - len(mandatory_indices))
        maximum_optional_size = min(
            len(preferred_optional), maximum_subset - len(mandatory_indices)
        )
        best = None
        for optional_size in range(
            minimum_optional, maximum_optional_size + 1
        ):
            for optional in itertools.combinations(
                preferred_optional, optional_size
            ):
                subset = tuple(mandatory_indices + list(optional))
                if len(subset) < minimum_subset:
                    continue
                candidate, failure = attempt(subset)
                if candidate is None:
                    failures.append(failure)
                    continue
                if best is None or candidate["score"] < best["score"]:
                    best = candidate
        if best is None:
            full_indices = tuple(range(len(allowed)))
            full_result, full_failure = attempt(full_indices)
            failures.append(full_failure)
            if full_result is not None:
                active_indices = tuple(
                    index
                    for index, route_id in enumerate(allowed)
                    if route_id in mandatory
                    or full_result["tensions"][index] >= 1.0
                )
                pruned_result = None
                if (
                    len(active_indices) < len(full_indices)
                    and len(active_indices) >= minimum_subset
                ):
                    pruned_result, _ = attempt(active_indices)
                best = pruned_result or full_result
        if best is None:
            return {
                "ok": False,
                "reason": "no_hard_constraint_tension_subset",
                "allowed": allowed,
                "selected": [],
                "failed_candidates": failures[:80],
            }
        selected = best["selected"]
        tensions = best["tensions"]
        return {
            "ok": True,
            "allowed": allowed,
            "selected": selected,
            "mandatory": mandatory,
            "solver": "linprog_exact_minimum_total_tension",
            "positive_spanning": best["positive_spanning"],
            "tensions_N": {
                route_id: float(tensions[index])
                for index, route_id in enumerate(selected)
            },
            "tpu_N": {
                route_id: float(best["tpu"][index])
                for index, route_id in enumerate(selected)
            },
            "target_generalized_force": best["target"].tolist(),
            "residual": best["residual"].tolist(),
            "relative_residual": best["relative_residual"],
            "max_tension_N": float(np.max(tensions)),
            "minimum_tension_margin_N": float(tensions.min()),
            "capstan_efficiency": {
                route_id: self._guide_efficiency(route_id)
                for route_id in selected
            },
            "max_contraction_percent": best["max_contraction_percent"],
            "min_cable_clearance_mm": best["min_cable_clearance_mm"],
        }

    def evaluate_pose_sample(
        self,
        action: dict,
        fraction: float,
        solve_tension: bool,
        route_ids: list[str] | None = None,
    ) -> dict:
        data = self._action_pose(action, fraction)
        achievement = self._action_achievement(action, data)
        tension = (
            self._solve_tensions(action, data)
            if solve_tension
            else {"ok": True, "selected": [], "skipped": True}
        )
        ids = route_ids if route_ids is not None else action["allowed_cables"]
        routes = [
            self._route_geometry(self.route_by_id[route_id], data)
            for route_id in ids
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
            "tension_solution": tension,
            "routes": routes,
            "contacts": int(data.ncon),
            "contact_gap_mm": self._contact_gap(data),
            "visual_clearance_mm": self._visual_clearance(data),
            "cross_cable": cross_cable,
        }

    def evaluate_all(self) -> dict:
        tension_fractions = [
            float(value)
            for value in self.config["acceptance"]["action_pose_fractions"]
        ]
        geometry_fractions = [
            float(value)
            for value in self.config["acceptance"]["final_pose_fractions"]
        ]
        action_reports = []
        for action in self.actions["actions"]:
            tension_samples = [
                self.evaluate_pose_sample(action, fraction, True)
                for fraction in tension_fractions
            ]
            selected_union = sorted(
                {
                    route_id
                    for sample in tension_samples
                    for route_id in sample["tension_solution"].get(
                        "selected", []
                    )
                }
            )
            samples = [
                self.evaluate_pose_sample(
                    action,
                    fraction,
                    False,
                    selected_union or action["allowed_cables"],
                )
                for fraction in geometry_fractions
            ]
            final_tension = tension_samples[-1]
            final = samples[-1]
            selected = final_tension["tension_solution"].get("selected", [])
            all_selected_routes = [
                route for sample in samples for route in sample["routes"]
            ]
            max_contraction = max(
                (route["contraction_percent"] for route in all_selected_routes),
                default=float("-inf"),
            )
            min_clearance = min(
                (route["thorax_clearance_mm"] for route in all_selected_routes),
                default=float("inf"),
            )
            min_visual = min(
                sample["visual_clearance_mm"] for sample in samples
            )
            min_contact = min(sample["contact_gap_mm"] for sample in samples)
            tension_ok = all(
                sample["tension_solution"].get("ok", False)
                for sample in tension_samples
            )
            spanning_ok = all(
                sample["tension_solution"]
                .get("positive_spanning", {})
                .get("pass", False)
                for sample in tension_samples
            )
            guide_ok = True
            for sample in samples:
                for route in sample["routes"]:
                    guide_ok &= route["guide_count"] <= 3
                    guide_ok &= all(
                        radius >= float(
                            self.config["guides"]["minimum_radius_mm"]
                        )
                        for radius in route["guide_radius_mm"]
                    )
                    guide_ok &= (
                        route["guide_count"] == 0
                        or route["actual_wrap_points"] >= 4
                    )
            contraction_ok = max_contraction <= float(
                self.config["acceptance"]["cable_contraction_percent_max"]
            )
            cable_clearance_ok = min_clearance >= float(
                self.config["acceptance"]["cable_clearance_mm_min"]
            )
            visual_ok = min_visual >= float(
                self.config["acceptance"]["full_mesh_clearance_mm_min"]
            )
            contact_ok = min_contact >= float(
                self.config["acceptance"]["contact_penetration_mm_min"]
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
                    "pass_minimum_visible": final_tension["achievement"][
                        "meets_minimum_visible"
                    ],
                    "pass_80_percent": final_tension["achievement"][
                        "meets_80_percent"
                    ],
                    "pass_tension": tension_ok,
                    "pass_positive_spanning": spanning_ok,
                    "max_selected_contraction_percent": max_contraction,
                    "pass_contraction": contraction_ok,
                    "min_selected_cable_clearance_mm": min_clearance,
                    "pass_cable_clearance": cable_clearance_ok,
                    "pass_guide_geometry": bool(guide_ok),
                    "visual_clearance_mm": min_visual,
                    "pass_visual_clearance": visual_ok,
                    "contact_gap_min_mm": min_contact,
                    "pass_contact": contact_ok,
                    "pass_minimum_functional": bool(
                        final_tension["achievement"]["meets_minimum_visible"]
                        and tension_ok
                        and spanning_ok
                        and contraction_ok
                        and cable_clearance_ok
                        and guide_ok
                        and visual_ok
                        and contact_ok
                    ),
                }
            )
        return {
            "model": str(self.model_path),
            "schema_version": 5,
            "zero_pose_offsets_deg": self.config.get(
                "zero_pose_offsets_deg", {}
            ),
            "contact_pad_offset_mm": self.config.get(
                "contact_pad_offset_mm", 0.0
            ),
            "neutral_lengths_mm": {
                route_id: length * 1000.0
                for route_id, length in self.neutral_lengths.items()
            },
            "actions": action_reports,
        }
