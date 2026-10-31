"""Optimize S01-S04 anchor zones toward a four-muscle bladebone topology."""

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
import math
from project_runtime import WorkspacePath as Path

import mujoco
import numpy as np
from scipy.optimize import differential_evolution

from evaluate_v2 import V2Evaluator, design_vector_to_mj, r400_center
from v2_common import (
    ACTION_MATRIX_PATH,
    CONFIG_PATH,
    EXPORT_ROOT,
    MODEL_PATH,
    POINTS_MM,
    design_to_mj,
    load_json,
    normalize,
)


# Chest-anchor zones: (x_min, x_max), (y_min, y_max), fixed outward offset.
CHEST_ZONES = {
    "S01": ((55.0, 120.0), (100.0, 220.0), 8.0),
    "S02": ((70.0, 140.0), (100.0, 220.0), 8.0),
    "S03": ((55.0, 120.0), (220.0, 278.0), 8.0),
    "S04": ((70.0, 140.0), (220.0, 278.0), 8.0),
}

# Bounds are scapula-plane (u, v) coordinates in millimetres.
SCAPULA_ZONES = {
    "S01": ((40.0, 100.0), (-20.0, 40.0)),
    "S02": ((30.0, 100.0), (40.0, 110.0)),
    "S03": ((-110.0, -45.0), (-20.0, 60.0)),
    "S04": ((-110.0, -45.0), (70.0, 130.0)),
}


def scapula_basis(config: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    del config
    p = POINTS_MM["P"]
    gh = POINTS_MM["GH"]
    sup = POINTS_MM["SUP"]
    inf = POINTS_MM["INF"]
    u_axis = normalize(sup - inf)
    lateral = gh - p
    v_axis = normalize(lateral - np.dot(lateral, u_axis) * u_axis)
    w_axis = normalize(np.asarray([0.6556, 0.0, 0.7552]))
    return u_axis, v_axis, w_axis


def chest_anchor(
    config: dict, route_id: str, x: float, y: float
) -> np.ndarray:
    offset = CHEST_ZONES[route_id][2]
    centre = r400_center(config)
    radius = float(config["contact"]["thorax_band"]["radius_mm"])
    z = centre[1] + math.sqrt(max(0.0, radius**2 - (x - centre[0]) ** 2))
    normal = normalize(np.array([x - centre[0], 0.0, z - centre[1]]))
    point = np.array([x, y, z], dtype=float)
    return point + offset * normal


def scapula_anchor(
    config: dict, u: float, v: float, w: float = 0.0
) -> np.ndarray:
    p = POINTS_MM["P"]
    u_axis, v_axis, w_axis = scapula_basis(config)
    return p + u * u_axis + v * v_axis + w * w_axis


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument("--matrix", type=Path, default=ACTION_MATRIX_PATH)
    parser.add_argument(
        "--out",
        type=Path,
        default=EXPORT_ROOT / "four_muscle_optimization.json",
    )
    parser.add_argument("--maxiter", type=int, default=20)
    parser.add_argument("--popsize", type=int, default=5)
    args = parser.parse_args()
    config = load_json(args.config)
    actions = load_json(args.matrix)
    evaluator = V2Evaluator(args.model, config, actions)
    action = next(item for item in actions["actions"] if item["id"] == "K6")

    site_ids = {}
    for route_id in CHEST_ZONES:
        site_ids[route_id] = {
            "A": mujoco.mj_name2id(
                evaluator.model,
                mujoco.mjtObj.mjOBJ_SITE,
                f"V2_A_{route_id}",
            ),
            "B": mujoco.mj_name2id(
                evaluator.model,
                mujoco.mjtObj.mjOBJ_SITE,
                f"V2_B_{route_id}",
            ),
        }

    route_ids = ["S01", "S02", "S03", "S04"]
    bounds = []
    for route_id in route_ids:
        (x_min, x_max), (y_min, y_max), _ = CHEST_ZONES[route_id]
        bounds.extend(((x_min, x_max), (y_min, y_max)))
        bounds.extend(SCAPULA_ZONES[route_id])

    evaluations = []

    def apply(values: np.ndarray) -> float:
        for index, route_id in enumerate(route_ids):
            base = 4 * index
            a_x = float(values[base])
            a_y = float(values[base + 1])
            u = float(values[base + 2])
            v = float(values[base + 3])
            a_design = chest_anchor(config, route_id, a_x, a_y)
            b_design = scapula_anchor(config, u, v)
            evaluator.model.site_pos[site_ids[route_id]["A"]] = design_to_mj(
                a_design
            )
            evaluator.model.site_pos[site_ids[route_id]["B"]] = design_to_mj(
                b_design - POINTS_MM["AC"]
            )
        zero = mujoco.MjData(evaluator.model)
        zero.qpos[:] = 0.0
        mujoco.mj_forward(evaluator.model, zero)
        evaluator.neutral_lengths = {
            route_id: float(zero.ten_length[tendon_id])
            for route_id, tendon_id in evaluator.tendon_ids.items()
        }
        result = evaluator.evaluate_action(action, 1.0)
        solution = result["tension_solution"]
        positive = solution.get("positive_spanning", {})
        residual = float(solution.get("relative_residual", 10.0))
        max_tension = float(solution.get("max_tension_N", 150.0))
        contraction = max(
            (
                route["contraction_percent"]
                for route in result["routes"]
                if route["id"] in route_ids
            ),
            default=100.0,
        )
        clearance = min(
            (
                route["thorax_clearance_mm"]
                for route in result["routes"]
                if route["id"] in route_ids
            ),
            default=-100.0,
        )
        score = residual
        score += 0.0 if positive.get("pass", False) else 2.0
        score += 0.1 * max_tension / 150.0
        score += 0.20 * max(0.0, contraction - 18.0)
        score += 0.10 * max(0.0, 2.0 - clearance)
        evaluations.append(
            {
                "values": values.tolist(),
                "score": score,
                "relative_residual": residual,
                "positive_spanning": positive,
                "max_tension_N": max_tension,
                "max_contraction_percent": contraction,
                "min_clearance_mm": clearance,
            }
        )
        return score

    result = differential_evolution(
        apply,
        bounds,
        maxiter=args.maxiter,
        popsize=args.popsize,
        seed=41,
        polish=False,
        workers=1,
        updating="immediate",
    )
    best = min(evaluations, key=lambda item: item["score"])
    routes = {}
    for index, route_id in enumerate(route_ids):
        base = 4 * index
        a_x = float(best["values"][base])
        a_y = float(best["values"][base + 1])
        u = float(best["values"][base + 2])
        v = float(best["values"][base + 3])
        routes[route_id] = {
            "A_mm": chest_anchor(config, route_id, a_x, a_y).tolist(),
            "B_local_uvw_mm": [u, v, 0.0],
            "B_mm": scapula_anchor(config, u, v).tolist(),
        }
    output = {
        "objective": result.fun,
        "success": bool(result.success),
        "best": best,
        "routes": routes,
        "evaluation_count": len(evaluations),
        "note": (
            "Candidate only. The main V2 config is not overwritten unless the "
            "four-muscle topology passes positive-spanning and force closure."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
