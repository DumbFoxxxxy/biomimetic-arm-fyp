"""Optimize S01-S04 anchors for scapular wrench closure."""

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
OUT = ROOT.parent / "exports" / "cable_design_v5" / "s_route_optimization_v5.json"
ROUTES = ("S01", "S02", "S03", "S04")
ACTION_IDS = ("K6", "K7", "SLIDE_FWD", "SLIDE_BACK")
FRACTIONS = (0.0, 0.25, 0.5, 0.75, 1.0)


def site_names(route_id: str) -> tuple[str, str]:
    return f"V2_A_{route_id}", f"V2_B_{route_id}"


def prepare_matrix(matrix: dict) -> None:
    for action_id in ("SLIDE_FWD", "SLIDE_BACK"):
        action = next(item for item in matrix["actions"] if item["id"] == action_id)
        action["allowed_cables"] = [
            "S01",
            "S02",
            "S03",
            "S04",
            "S05",
            "S06",
        ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--maxiter", type=int, default=35)
    parser.add_argument("--popsize", type=int, default=6)
    args = parser.parse_args()
    config = load_json(args.config)
    matrix = load_json(args.matrix)
    prepare_matrix(matrix)
    evaluator = V5Evaluator(args.model, config, matrix)
    actions = {
        action_id: next(
            item for item in matrix["actions"] if item["id"] == action_id
        )
        for action_id in ACTION_IDS
    }
    ids = {}
    base = {}
    for route_id in ROUTES:
        a_name, b_name = site_names(route_id)
        ids[route_id] = {
            "A": mujoco.mj_name2id(evaluator.model, mujoco.mjtObj.mjOBJ_SITE, a_name),
            "B": mujoco.mj_name2id(evaluator.model, mujoco.mjtObj.mjOBJ_SITE, b_name),
        }
        base[route_id] = {
            element: evaluator.model.site_pos[site_id].copy()
            for element, site_id in ids[route_id].items()
        }
    history = []

    def apply_offsets(values):
        cursor = 0
        for route_id in ROUTES:
            for element in ("A", "B"):
                offset = np.asarray(values[cursor : cursor + 3], dtype=float)
                cursor += 3
                evaluator.model.site_pos[ids[route_id][element]] = (
                    base[route_id][element] + design_to_mj(offset) * 0.001
                )

    def refresh_neutral():
        zero = evaluator._zero_data()
        evaluator.neutral_lengths = {
            route_id: float(zero.ten_length[tendon_id])
            for route_id, tendon_id in evaluator.tendon_ids.items()
        }

    def cone_for(action, fraction):
        data = evaluator._action_pose(action, fraction)
        dofs = evaluator._dofs_for_action(action)
        columns = [
            evaluator.tendon_indices[route_id]
            for route_id in action["allowed_cables"]
        ]
        jacobian = dense_tendon_jacobian(evaluator.model, data)
        moment = -jacobian[:, dofs][columns].T
        efficiency = np.asarray(
            [
                evaluator._guide_efficiency(route_id)
                for route_id in action["allowed_cables"]
            ]
        )
        moment = moment * efficiency[None, :]
        return moment

    def objective(values):
        apply_offsets(values)
        refresh_neutral()
        violations = []
        for action_id in ACTION_IDS:
            action = actions[action_id]
            for fraction in FRACTIONS:
                moment = cone_for(action, fraction)
                violations.append(evaluator._cone_violation(moment))
                if action_id == "K6":
                    indices = [
                        action["allowed_cables"].index(route_id)
                        for route_id in ROUTES
                    ]
                    violations.append(
                        evaluator._cone_violation(moment[:, indices])
                    )
        score = (
            100.0 * float(np.sum(violations))
            + 1e-6 * float(np.dot(values, values))
        )
        history.append(
            {
                "values": values.tolist(),
                "score": score,
                "sum_cone_violation": float(np.sum(violations)),
            }
        )
        return score

    result = differential_evolution(
        objective,
        [(-20.0, 20.0)] * (len(ROUTES) * 6),
        seed=97,
        maxiter=args.maxiter,
        popsize=args.popsize,
        polish=False,
        workers=1,
        updating="immediate",
    )
    best = min(history, key=lambda item: item["score"])
    apply_offsets(np.asarray(best["values"]))
    refresh_neutral()
    zero = evaluator._zero_data()
    for route_id in ROUTES:
        route = next(item for item in config["routes"] if item["id"] == route_id)
        for element, name in zip(route["path"], site_names(route_id)):
            site_id = mujoco.mj_name2id(
                evaluator.model, mujoco.mjtObj.mjOBJ_SITE, name
            )
            element["world_mm"] = mj_to_design(
                zero.site_xpos[site_id]
            ).tolist()
            element["local_mm"] = evaluator.model.site_pos[site_id].tolist()
        route["status"] = "v5_wrench_closure_search"
    args.config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    args.matrix.write_text(json.dumps(matrix, indent=2) + "\n", encoding="utf-8")
    output = {
        "optimizer_success": bool(result.success),
        "best": best,
        "evaluation_count": len(history),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
