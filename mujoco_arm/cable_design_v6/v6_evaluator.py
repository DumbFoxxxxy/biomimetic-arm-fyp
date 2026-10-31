"""V6 evaluator with task-conditioned tension feasibility."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import math
from project_runtime import WorkspacePath as Path
import sys

import mujoco
import numpy as np
from scipy.optimize import lsq_linear


ROOT = Path(__file__).resolve().parent
V5_ROOT = ROOT.parent / "cable_design_v5"
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(V5_ROOT), str(V2_ROOT)]

from evaluate_v2 import dense_tendon_jacobian, design_vector_to_mj  # noqa: E402
from v5_evaluator import V5Evaluator  # noqa: E402


class V6Evaluator(V5Evaluator):
    """Relax global wrench closure while retaining local force feasibility."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.roller_site = mujoco.mj_name2id(
            self.model,
            mujoco.mjtObj.mjOBJ_SITE,
            "V6_roller",
        )
        self.wrist_site = mujoco.mj_name2id(
            self.model,
            mujoco.mjtObj.mjOBJ_SITE,
            "V6_W",
        )
        self.elbow_site = mujoco.mj_name2id(
            self.model,
            mujoco.mjtObj.mjOBJ_SITE,
            "V6_E",
        )
        guide = self.config.get("shoulder_guide", {})
        self.preload_n = float(guide.get("nominal_preload_N", 5.0))
        outward = np.asarray([0.6556, 0.0, 0.7552], dtype=float)
        self.preload_direction_mj = -design_vector_to_mj(outward)

    def _dofs_for_action(self, action: dict) -> list[int]:
        if action.get("solve_all_dofs"):
            return [int(self.dof_indices[name]) for name in self.dof_indices]
        return super()._dofs_for_action(action)

    def _preload_generalized_force(self, data: mujoco.MjData) -> np.ndarray:
        result = np.zeros(self.model.nv)
        if self.roller_site < 0:
            return result
        jacobian = np.zeros((3, self.model.nv))
        mujoco.mj_jacSite(
            self.model,
            data,
            jacobian,
            None,
            self.roller_site,
        )
        force = self.preload_direction_mj * self.preload_n
        return jacobian.T @ force

    def _solve_tensions(self, action: dict, data: mujoco.MjData) -> dict:
        allowed = list(
            action.get(
                "active_cables",
                action.get(
                    "preferred_cables",
                    action["allowed_cables"],
                ),
            )
        )
        for route_id in action.get("mandatory_cables", []):
            if route_id not in allowed:
                allowed.append(route_id)
        allowed = [route_id for route_id in allowed if route_id in self.tendon_indices]
        if not allowed:
            return {
                "ok": False,
                "reason": "no_available_cables",
                "selected": [],
            }
        tendon_columns = [self.tendon_indices[route_id] for route_id in allowed]
        dofs = self._dofs_for_action(action)
        jacobian = dense_tendon_jacobian(self.model, data)
        moment = -jacobian[:, dofs][tendon_columns].T
        efficiencies = np.asarray(
            [self._guide_efficiency(route_id) for route_id in allowed]
        )
        moment = moment * efficiencies[None, :]
        tpu_map = self._tpu_forces(
            {"allowed_cables": allowed},
            data,
        )
        tpu = np.asarray([tpu_map[route_id] for route_id in allowed])
        bias = np.asarray(data.qfrc_bias)[dofs]
        preload = self._preload_generalized_force(data)[dofs]
        tracking = np.asarray(
            action.get("torque_tracking", np.zeros(len(dofs))),
            dtype=float,
        )
        if tracking.shape != (len(dofs),):
            raise ValueError(
                f"torque_tracking must have shape ({len(dofs)},), got {tracking.shape}"
            )
        target = bias - moment @ tpu - preload + tracking
        maximum = float(self.config["acceptance"]["tension_max_N"])
        minimum_active = float(
            self.config.get("tension_solver", {}).get(
                "minimum_active_tension_N", 0.25
            )
        )
        mandatory = [
            route_id
            for route_id in action.get("mandatory_cables", [])
            if route_id in allowed
        ]
        minimum_mandatory = float(
            self.config.get("tension_solver", {}).get(
                "minimum_mandatory_tension_N", 1.0
            )
        )
        spanning = self._positive_spanning(
            moment,
            allowed,
        )

        def solve(indices: list[int]) -> tuple[np.ndarray, float]:
            sub_moment = moment[:, indices]
            selected = [allowed[index] for index in indices]
            lower = np.zeros(len(indices))
            for route_id in mandatory:
                lower[selected.index(route_id)] = minimum_mandatory
            result = lsq_linear(
                sub_moment,
                target,
                bounds=(lower, np.full(len(indices), maximum)),
                lsmr_tol="auto",
                max_iter=4000,
            )
            tensions = np.asarray(result.x)
            residual = sub_moment @ tensions - target
            relative = float(np.linalg.norm(residual)) / (
                float(np.linalg.norm(target)) + 1e-9
            )
            return tensions, relative

        all_indices = list(range(len(allowed)))
        tensions, relative = solve(all_indices)
        active_indices = [
            index
            for index in all_indices
            if tensions[index] >= minimum_active
            or allowed[index] in mandatory
        ]
        if active_indices and len(active_indices) < len(all_indices):
            refined, refined_relative = solve(active_indices)
            if refined_relative <= relative + 0.01:
                tensions = np.zeros(len(allowed))
                tensions[active_indices] = refined
                relative = refined_relative
        selected = [
            route_id
            for index, route_id in enumerate(allowed)
            if tensions[index] >= minimum_active or route_id in mandatory
        ]
        selected_indices = [allowed.index(route_id) for route_id in selected]
        selected_routes = [
            self._route_geometry(self.route_by_id[route_id], data)
            for route_id in selected
        ]
        max_contraction = max(
            (route["contraction_percent"] for route in selected_routes),
            default=float("-inf"),
        )
        min_clearance = min(
            (route["thorax_clearance_mm"] for route in selected_routes),
            default=float("inf"),
        )
        tolerance = float(
            self.config.get("tension_solver", {}).get(
                "exact_equilibrium_relative_tolerance", 0.08
            )
        )
        contraction_limit = float(
            self.config["acceptance"]["cable_contraction_percent_max"]
        )
        clearance_limit = float(
            self.config["acceptance"]["cable_clearance_mm_min"]
        )
        selected_tensions = np.asarray(
            [tensions[index] for index in selected_indices], dtype=float
        )
        ok = bool(
            relative <= tolerance
            and max_contraction <= contraction_limit
            and min_clearance >= clearance_limit
            and selected_tensions.size > 0
            and np.all(np.isfinite(selected_tensions))
            and np.all(selected_tensions <= maximum)
        )
        return {
            "ok": ok,
            "allowed": allowed,
            "selected": selected,
            "mandatory": mandatory,
            "solver": "bounded_least_squares_task_conditioned",
            "positive_spanning": spanning,
            "positive_spanning_hard_gate": False,
            "tensions_N": {
                route_id: float(tensions[allowed.index(route_id)])
                for route_id in selected
            },
            "tpu_N": {
                route_id: float(tpu_map[route_id]) for route_id in selected
            },
            "preload_N": float(self.preload_n),
            "target_generalized_force": target.tolist(),
            "relative_residual": relative,
            "max_tension_N": (
                float(np.max(selected_tensions))
                if selected_tensions.size
                else None
            ),
            "minimum_tension_margin_N": (
                float(np.min(selected_tensions))
                if selected_tensions.size
                else None
            ),
            "capstan_efficiency": {
                route_id: self._guide_efficiency(route_id)
                for route_id in selected
            },
            "max_contraction_percent": max_contraction,
            "min_cable_clearance_mm": min_clearance,
        }

    def evaluate_all(self) -> dict:
        evaluation = super().evaluate_all()
        evaluation["schema_version"] = 6
        evaluation["positive_spanning_gate"] = False
        for action in evaluation["actions"]:
            action["pass_positive_spanning_diagnostic"] = action.get(
                "pass_positive_spanning", False
            )
            action["positive_spanning_is_diagnostic"] = True
            action["pass_minimum_functional"] = bool(
                action["pass_minimum_visible"]
                and action["pass_tension"]
                and action["pass_contraction"]
                and action["pass_cable_clearance"]
                and action["pass_guide_geometry"]
                and action["pass_visual_clearance"]
                and action["pass_contact"]
            )
        return evaluation

