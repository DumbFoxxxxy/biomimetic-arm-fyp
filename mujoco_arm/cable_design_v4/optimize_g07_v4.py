"""Optimize the dedicated K7 GH cable G07."""

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
import sys

import mujoco
import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path.insert(0, str(V2_ROOT))

from build_v2_model import build as build_model  # noqa: E402
from evaluate_v2 import V2Evaluator  # noqa: E402
from v2_common import (  # noqa: E402
    POINTS_MM,
    design_to_mj,
    load_json,
    mj_to_design,
    normalize,
    quat_from_z_axis,
)


V4_CONFIG = ROOT / "cable_design_v4.json"
V4_MATRIX = ROOT / "action_cable_matrix_v4.json"
V4_MODEL = ROOT / "arm_cable_design_v4.xml"
DEFAULT_OUT = ROOT.parent / "exports" / "cable_design_v4" / "g07_optimization_v4.json"


def scapula_basis() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u_axis = normalize(POINTS_MM["SUP"] - POINTS_MM["INF"])
    lateral = POINTS_MM["GH"] - POINTS_MM["P"]
    v_axis = normalize(lateral - np.dot(lateral, u_axis) * u_axis)
    w_axis = normalize(np.array([0.6556, 0.0, 0.7552]))
    return u_axis, v_axis, w_axis


def local_design_to_mj_vector(values_design_mm: np.ndarray) -> np.ndarray:
    values_design_mm = np.asarray(values_design_mm, dtype=float)
    return np.array(
        [
            values_design_mm[0],
            -values_design_mm[2],
            values_design_mm[1],
        ]
    ) * 0.001


def candidate_geometry(
    config: dict,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    route = next(item for item in config["routes"] if item["id"] == "G07")
    a_old = np.asarray(route["path"][0]["world_mm"])
    u_axis, v_axis, w_axis = scapula_basis()
    a_rel = a_old - POINTS_MM["P"]
    a_uvw = np.array(
        [np.dot(a_rel, u_axis), np.dot(a_rel, v_axis), np.dot(a_rel, w_axis)]
    )
    a = (
        POINTS_MM["P"]
        + (a_uvw[0] + values[0]) * u_axis
        + (a_uvw[1] + values[1]) * v_axis
        + np.clip(a_uvw[2] + values[2], 0.0, 10.0) * w_axis
    )
    b_old_local = design_to_mj(
        np.asarray(route["path"][-1]["world_mm"]) - POINTS_MM["GH"]
    )
    b_local = b_old_local + np.array(values[3:6]) * 0.001
    b = POINTS_MM["GH"] + np.array(
        [b_local[0], b_local[2], -b_local[1]]
    ) * 1000.0
    g_old_local = design_to_mj(
        np.asarray(route["path"][1]["world_mm"]) - POINTS_MM["AC"]
    )
    g_local = g_old_local + np.array(values[6:9]) * 0.001
    g = POINTS_MM["AC"] + np.array(
        [g_local[0], g_local[2], -g_local[1]]
    ) * 1000.0
    incoming = g - a
    outgoing = b - g
    axis = normalize(np.cross(incoming, outgoing))
    side = normalize(
        np.cross(axis, normalize(incoming) + normalize(outgoing))
    )
    bend_angle = math.acos(
        float(
            np.clip(
                np.dot(normalize(-incoming), normalize(outgoing)),
                -1.0,
                1.0,
            )
        )
    )
    wrap_angle = math.degrees(math.pi - bend_angle)
    return a, b, g, axis, side, wrap_angle


def segment_distance_mm(start: np.ndarray, end: np.ndarray, point: np.ndarray) -> float:
    segment = end - start
    segment_length_sq = float(np.dot(segment, segment))
    if segment_length_sq <= 1e-12:
        return float(np.linalg.norm(point - start))
    fraction = float(np.clip(np.dot(point - start, segment) / segment_length_sq, 0.0, 1.0))
    closest = start + fraction * segment
    return float(np.linalg.norm(point - closest))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=V4_CONFIG)
    parser.add_argument("--matrix", type=Path, default=V4_MATRIX)
    parser.add_argument("--model", type=Path, default=V4_MODEL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--maxiter", type=int, default=16)
    parser.add_argument("--popsize", type=int, default=4)
    args = parser.parse_args()
    config = load_json(args.config)
    actions = load_json(args.matrix)
    model = _load_project_model(str(args.model))
    evaluator = V2Evaluator(None, config, actions, model=model)
    k7 = next(action for action in actions["actions"] if action["id"] == "K7")
    a_id = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_SITE, "V2_A_G07"
    )
    b_id = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_SITE, "V2_B_G07"
    )
    g_id = mujoco.mj_name2id(
        evaluator.model,
        mujoco.mjtObj.mjOBJ_GEOM,
        "V2_guide_G07_guide_G07_v4",
    )
    side_id = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_SITE, "V2_side_G07"
    )
    if min(a_id, b_id, g_id, side_id) < 0:
        raise RuntimeError("G07 template entities are missing")
    history = []

    def objective(values):
        a, b, g, axis, side, wrap_angle = candidate_geometry(config, values)
        evaluator.model.site_pos[a_id] = design_to_mj(a - POINTS_MM["AC"])
        evaluator.model.site_pos[b_id] = design_to_mj(b - POINTS_MM["GH"])
        evaluator.model.geom_pos[g_id] = design_to_mj(g - POINTS_MM["AC"])
        axis_mj = normalize(np.array([axis[0], -axis[2], axis[1]]))
        side_mj = normalize(np.array([side[0], -side[2], side[1]]))
        evaluator.model.geom_quat[g_id] = quat_from_z_axis(axis_mj)
        evaluator.model.site_pos[side_id] = (
            design_to_mj(g - POINTS_MM["AC"]) + side_mj * 0.012
        )
        route = evaluator.route_by_id["G07"]
        route["path"][0]["world_mm"] = a.tolist()
        route["path"][1]["world_mm"] = g.tolist()
        route["path"][1]["local_mm"] = design_to_mj(g - POINTS_MM["AC"]).tolist()
        route["path"][1]["axis_design"] = axis.tolist()
        route["path"][1]["side_design"] = side.tolist()
        route["path"][1]["wrap_angle_deg"] = wrap_angle
        route["path"][-1]["world_mm"] = b.tolist()

        zero = mujoco.MjData(evaluator.model)
        zero.qpos[:] = 0.0
        mujoco.mj_forward(evaluator.model, zero)
        evaluator.neutral_lengths = {
            route_id: float(zero.ten_length[tendon_id])
            for route_id, tendon_id in evaluator.tendon_ids.items()
        }
        checks = []
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            data = evaluator._action_pose(k7, fraction)
            mujoco.mj_forward(evaluator.model, data)
            route_data = evaluator._route_geometry(route, data)
            line_distance_mm = segment_distance_mm(
                data.site_xpos[a_id],
                data.site_xpos[b_id],
                data.geom_xpos[g_id],
            ) * 1000.0
            checks.append(
                {
                    "fraction": fraction,
                    "wrap_points": route_data["actual_wrap_points"],
                    "line_distance_mm": line_distance_mm,
                    "clearance_mm": route_data["thorax_clearance_mm"],
                    "contraction_percent": route_data["contraction_percent"],
                }
            )
        final = evaluator.evaluate_action(k7, 1.0)
        solution = final["tension_solution"]
        inactive = sum(check["wrap_points"] < 4 for check in checks)
        engagement_error = sum(
            max(0.0, check["line_distance_mm"] - 7.5) ** 2 for check in checks
        )
        max_contraction = max(check["contraction_percent"] for check in checks)
        min_clearance = min(check["clearance_mm"] for check in checks)
        score = 50.0 * inactive + 0.5 * engagement_error
        score += 5.0 * max(0.0, 2.0 - min_clearance)
        score += 0.2 * max(0.0, max_contraction - 18.0)
        score += 2.0 * max(0.0, float(solution.get("relative_residual", 10.0)) - 0.1)
        score += 0.0 if solution.get("positive_spanning", {}).get("pass", False) else 2.0
        score += 0.1 * max(0.0, float(solution.get("max_tension_N", 0.0)) - 149.0)
        history.append(
            {
                "values": values.tolist(),
                "score": score,
                "wrap_checks": checks,
                "min_clearance_mm": min_clearance,
                "max_contraction_percent": max_contraction,
                "tension_ok": bool(solution.get("ok", False)),
                "positive_spanning": solution.get("positive_spanning", {}),
                "relative_residual": solution.get("relative_residual"),
                "max_tension_N": solution.get("max_tension_N"),
                "selected_cables": final["tension_solution"].get("selected", []),
                "a_mm": a.tolist(),
                "b_mm": b.tolist(),
                "g_mm": g.tolist(),
                "axis_design": axis.tolist(),
                "side_design": side.tolist(),
                "wrap_angle_deg": wrap_angle,
            }
        )
        return score

    bounds = [
        (-20.0, 20.0),
        (-20.0, 20.0),
        (0.0, 10.0),
        (-20.0, 20.0),
        (-20.0, 20.0),
        (-30.0, 30.0),
        (-30.0, 30.0),
        (-30.0, 30.0),
        (-30.0, 30.0),
    ]
    result = differential_evolution(
        objective,
        bounds,
        maxiter=args.maxiter,
        popsize=args.popsize,
        seed=71,
        polish=False,
        workers=1,
        updating="immediate",
    )
    best = min(history, key=lambda item: item["score"])
    output = {
        "optimizer_success": bool(result.success),
        "optimizer_objective": float(result.fun),
        "evaluation_count": len(history),
        "best": best,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
