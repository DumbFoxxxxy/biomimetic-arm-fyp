"""Optimize V5 zero pose and five-knot shoulder-path corrections."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import json
from project_runtime import WorkspacePath as Path
import sys

import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(ROOT))

from v2_common import load_json  # noqa: E402
from evaluate_v2 import dense_tendon_jacobian  # noqa: E402
from v5_evaluator import V5Evaluator  # noqa: E402


CONFIG = ROOT / "cable_design_v5.json"
MATRIX = ROOT / "action_cable_matrix_v5.json"
MODEL = ROOT / "arm_cable_design_v5.xml"
OUT = ROOT.parent / "exports" / "cable_design_v5" / "pose_optimization_v5.json"
JOINT_NAMES = ("sc_elevation", "sc_protraction", "ac_ext")
FRACTIONS = tuple(round(index / 20.0, 2) for index in range(21))
KNOTS = (0.0, 0.25, 0.5, 0.75, 1.0)
CORRECTION_ACTIONS = ("K6", "K7", "SLIDE_BACK", "SCAP_INT")
VISUAL_LIMIT = -0.5


def set_zero(config: dict, values: np.ndarray) -> None:
    config["zero_pose_offsets_deg"] = {
        name: float(value) for name, value in zip(JOINT_NAMES, values)
    }


def set_correction(config: dict, action_id: str, values: np.ndarray) -> None:
    values = np.asarray(values, dtype=float)
    if values.size == 3:
        node_values = np.tile(values, (len(KNOTS), 1))
    else:
        node_values = values.reshape(len(KNOTS), 3)
    config["path_corrections"][action_id] = {
        str(fraction): node_values[index].tolist()
        for index, fraction in enumerate(KNOTS)
    }


def optimize_zero(config, matrix, model, evaluator) -> dict:
    history = []

    def objective(values):
        set_zero(config, values)
        evaluator.zero_pose_offsets = {
            name: np.radians(float(value))
            for name, value in config["zero_pose_offsets_deg"].items()
        }
        evaluator._reset_zero_reference()
        data = evaluator._zero_data()
        clearance = evaluator._visual_clearance(data)
        violation = max(0.0, VISUAL_LIMIT - clearance) ** 2
        penalty = 1e-5 * float(np.dot(values, values))
        score = 10.0 * violation + penalty
        history.append(
            {
                "values": values.tolist(),
                "clearance_mm": clearance,
                "score": score,
            }
        )
        return score

    result = differential_evolution(
        objective,
        [(-10.0, 10.0)] * 3,
        seed=71,
        maxiter=18,
        popsize=6,
        polish=False,
        workers=1,
        updating="immediate",
    )
    best = min(history, key=lambda item: item["score"])
    set_zero(config, np.asarray(best["values"]))
    return {
        "optimizer_success": bool(result.success),
        "best": best,
        "evaluation_count": len(history),
    }


def optimize_action(config, action, evaluator) -> dict:
    history = []

    def objective(values):
        set_correction(config, action["id"], values)
        evaluator.path_corrections = config["path_corrections"]
        min_clearance = float("inf")
        spanning_failures = 0
        cone_violation = 0.0
        for fraction in FRACTIONS:
            data = evaluator._action_pose(action, fraction)
            min_clearance = min(
                min_clearance, evaluator._visual_clearance(data)
            )
            dofs = evaluator._dofs_for_action(action)
            columns = [
                evaluator.tendon_indices[route_id]
                for route_id in action["allowed_cables"]
            ]
            jacobian = dense_tendon_jacobian(evaluator.model, data)
            moment = -jacobian[:, dofs][columns].T
            moment = moment * np.asarray(
                [
                    evaluator._guide_efficiency(route_id)
                    for route_id in action["allowed_cables"]
                ]
            )[None, :]
            if not evaluator._positive_spanning(
                moment, action["allowed_cables"]
            )["pass"]:
                spanning_failures += 1
            cone_violation += evaluator._cone_violation(moment)
        final = evaluator._action_pose(action, 1.0)
        achievement = evaluator._action_achievement(action, final)
        metric_error = max(
            0.0,
            float(achievement["acceptance_80"])
            - float(achievement["achieved"]),
        )
        violation = max(0.0, VISUAL_LIMIT - min_clearance) ** 2
        score = (
            20.0 * violation
            + 5.0 * spanning_failures
            + 50.0 * cone_violation
            + 0.2 * metric_error**2
            + 1e-5 * float(np.dot(values, values))
        )
        history.append(
            {
                "values": values.tolist(),
                "min_clearance_mm": min_clearance,
                "spanning_failures": spanning_failures,
                "cone_violation": cone_violation,
                "achieved": achievement["achieved"],
                "acceptance": achievement["acceptance_80"],
                "score": score,
            }
        )
        return score

    node_action = action["id"] in ("K6", "SLIDE_BACK")
    variable_count = 15 if node_action else 3
    result = differential_evolution(
        objective,
        [(-10.0, 10.0)] * variable_count,
        seed=71 + int(action["id"][-1]) if action["id"][-1].isdigit() else 89,
        maxiter=16,
        popsize=5,
        polish=False,
        workers=1,
        updating="immediate",
    )
    best = min(history, key=lambda item: item["score"])
    set_correction(config, action["id"], np.asarray(best["values"]))
    return {
        "optimizer_success": bool(result.success),
        "best": best,
        "evaluation_count": len(history),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--skip-zero", action="store_true")
    args = parser.parse_args()
    config = load_json(args.config)
    matrix = load_json(args.matrix)
    evaluator = V5Evaluator(args.model, config, matrix)
    zero = (
        {"skipped": True, "best": {"values": list(config["zero_pose_offsets_deg"].values())}}
        if args.skip_zero
        else optimize_zero(config, matrix, args.model, evaluator)
    )
    evaluator.zero_pose_offsets = {
        name: np.radians(float(value))
        for name, value in config["zero_pose_offsets_deg"].items()
    }
    evaluator._reset_zero_reference()
    actions = []
    for action_id in CORRECTION_ACTIONS:
        action = next(item for item in matrix["actions"] if item["id"] == action_id)
        actions.append({"id": action_id, **optimize_action(config, action, evaluator)})
    args.config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    output = {
        "zero_pose": zero,
        "actions": actions,
        "visual_limit_mm": VISUAL_LIMIT,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
