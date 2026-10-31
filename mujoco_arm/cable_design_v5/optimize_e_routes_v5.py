"""Optimize E01/E02 anchor and pass-through positions for K1 antagonism."""

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

import mujoco
import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from evaluate_v2 import dense_tendon_jacobian  # noqa: E402
from v2_common import design_to_mj, load_json, mj_to_design  # noqa: E402
from v5_evaluator import V5Evaluator  # noqa: E402


CONFIG = ROOT / "cable_design_v5.json"
MATRIX = ROOT / "action_cable_matrix_v5.json"
MODEL = ROOT / "arm_cable_design_v5.xml"
OUT = ROOT.parent / "exports" / "cable_design_v5" / "e_route_optimization_v5.json"
FRACTIONS = np.linspace(0.0, 1.0, 21)

SITE_NAMES = {
    "E01": ("V2_A_E01", "V2_pass_E01_1", "V2_B_E01"),
    "E02": ("V2_A_E02", "V2_pass_E02_1", "V2_B_E02"),
}


def site_ids(model: mujoco.MjModel) -> dict[str, dict[str, int]]:
    return {
        route_id: {
            "A": mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, names[0]),
            "P": mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, names[1]),
            "B": mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, names[2]),
        }
        for route_id, names in SITE_NAMES.items()
    }


def apply_offsets(model, ids, base_positions, values: np.ndarray) -> None:
    cursor = 0
    for route_id in ("E01", "E02"):
        for element in ("A", "P", "B"):
            offset_design = np.asarray(values[cursor : cursor + 3], dtype=float)
            cursor += 3
            model.site_pos[ids[route_id][element]] = (
                base_positions[route_id][element]
                + design_to_mj(offset_design) * 0.001
            )


def refresh_neutral(evaluator: V5Evaluator) -> None:
    zero = evaluator._zero_data()
    evaluator.neutral_lengths = {
        route_id: float(zero.ten_length[tendon_id])
        for route_id, tendon_id in evaluator.tendon_ids.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--maxiter", type=int, default=30)
    parser.add_argument("--popsize", type=int, default=8)
    args = parser.parse_args()
    config = load_json(args.config)
    matrix = load_json(args.matrix)
    evaluator = V5Evaluator(args.model, config, matrix)
    action = next(item for item in matrix["actions"] if item["id"] == "K1")
    ids = site_ids(evaluator.model)
    base_positions = {
        route_id: {
            element: evaluator.model.site_pos[site_id].copy()
            for element, site_id in elements.items()
        }
        for route_id, elements in ids.items()
    }
    history = []

    def objective(values):
        apply_offsets(evaluator.model, ids, base_positions, values)
        refresh_neutral(evaluator)
        direction_penalty = 0.0
        for fraction in FRACTIONS:
            data = evaluator._action_pose(action, float(fraction))
            dofs = evaluator._dofs_for_action(action)
            tendon_columns = [
                evaluator.tendon_indices[route_id] for route_id in ("E01", "E02")
            ]
            jacobian = dense_tendon_jacobian(evaluator.model, data)
            moment = -jacobian[:, dofs][tendon_columns].T
            d0 = moment[:, 0] / max(np.linalg.norm(moment[:, 0]), 1e-12)
            d1 = moment[:, 1] / max(np.linalg.norm(moment[:, 1]), 1e-12)
            opposite = float(np.dot(d0, d1))
            direction_penalty += (opposite + 1.0) ** 2
        score = (
            direction_penalty
            + 1e-5 * float(np.dot(values, values))
        )
        history.append(
            {
                "values": values.tolist(),
                "score": score,
                "direction_penalty": direction_penalty,
            }
        )
        return score

    bounds = []
    for _ in ("E01", "E02"):
        bounds.extend([(-30.0, 30.0), (-30.0, 30.0), (-50.0, 50.0)])
        bounds.extend([(-30.0, 30.0), (-30.0, 30.0), (-50.0, 50.0)])
        bounds.extend([(-20.0, 20.0), (-20.0, 20.0), (-30.0, 30.0)])
    result = differential_evolution(
        objective,
        bounds,
        seed=83,
        maxiter=args.maxiter,
        popsize=args.popsize,
        polish=False,
        workers=1,
        updating="immediate",
    )
    best = min(history, key=lambda item: item["score"])
    apply_offsets(
        evaluator.model, ids, base_positions, np.asarray(best["values"])
    )
    refresh_neutral(evaluator)
    verification = []
    for fraction in FRACTIONS:
        data = evaluator._action_pose(action, float(fraction))
        solution = evaluator._solve_tensions(action, data)
        verification.append(
            {
                "fraction": float(fraction),
                "ok": bool(solution.get("ok", False)),
                "reason": solution.get("reason", ""),
            }
        )
    zero = evaluator._zero_data()
    for route_id, names in SITE_NAMES.items():
        route = next(item for item in config["routes"] if item["id"] == route_id)
        for element, name in zip(route["path"], names):
            site_id = mujoco.mj_name2id(
                evaluator.model, mujoco.mjtObj.mjOBJ_SITE, name
            )
            element["world_mm"] = mj_to_design(
                zero.site_xpos[site_id]
            ).tolist()
            element["local_mm"] = evaluator.model.site_pos[site_id].tolist()
    args.config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    output = {
        "optimizer_success": bool(result.success),
        "best": best,
        "evaluation_count": len(history),
        "verification_failures": sum(
            not item["ok"] for item in verification
        ),
        "verification": verification,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
