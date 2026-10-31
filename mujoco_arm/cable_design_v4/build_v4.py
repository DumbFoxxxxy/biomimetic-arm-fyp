"""Build and evaluate V4 with the dedicated K7 cable G07."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import copy
import json
import math
from project_runtime import WorkspacePath as Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path.insert(0, str(V2_ROOT))

from build_v2_model import build as build_model  # noqa: E402
from evaluate_v2 import (  # noqa: E402
    V2Evaluator,
    export_action_markdown,
    export_action_matrix,
    export_markdown,
    export_route_validation,
    export_routes_csv,
    export_solidworks_csv,
    export_tension_report,
)
from v2_common import (  # noqa: E402
    POINTS_MM,
    PROJECT_ROOT,
    design_to_mj,
    load_json,
    normalize,
)


V4_CONFIG = ROOT / "cable_design_v4.json"
V4_MATRIX = ROOT / "action_cable_matrix_v4.json"
V4_MODEL = ROOT / "arm_cable_design_v4.xml"
V4_EXPORT = ROOT.parent / "exports" / "cable_design_v4"


def anchor(owner: str, body: str, world: np.ndarray, normal: np.ndarray) -> dict:
    origin = {
        "S": POINTS_MM["AC"],
        "H": POINTS_MM["GH"],
    }[owner]
    return {
        "kind": "anchor",
        "owner": owner,
        "body": body,
        "world_mm": np.asarray(world).tolist(),
        "local_mm": design_to_mj(world - origin).tolist(),
        "installation_normal": normalize(normal).tolist(),
    }


def apply_g07(config: dict, result: dict) -> None:
    best = result["best"]
    route = next(item for item in config["routes"] if item["id"] == "G07")
    a = np.asarray(best["a_mm"], dtype=float)
    b = np.asarray(best["b_mm"], dtype=float)
    g = np.asarray(best["g_mm"], dtype=float)
    axis = normalize(np.asarray(best["axis_design"], dtype=float))
    side = normalize(np.asarray(best["side_design"], dtype=float))
    scapula_normal = normalize(np.array([0.6556, 0.0, 0.7552]))
    route["path"] = [
        anchor("S", "scapula", a, scapula_normal),
        {
            "kind": "guide",
            "name": "guide_G07_v4",
            "owner": "S",
            "body": "scapula",
            "world_mm": g.tolist(),
            "local_mm": design_to_mj(g - POINTS_MM["AC"]).tolist(),
            "radius_mm": 8.0,
            "axis_design": axis.tolist(),
            "side_design": side.tolist(),
            "wrap_angle_deg": float(best["wrap_angle_deg"]),
        },
        anchor("H", "humerus", b, side),
    ]
    route["status"] = "v4_optimized"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=V4_CONFIG)
    parser.add_argument("--matrix", type=Path, default=V4_MATRIX)
    parser.add_argument("--model", type=Path, default=V4_MODEL)
    parser.add_argument(
        "--g07-result",
        type=Path,
        default=V4_EXPORT / "g07_optimization_v4.json",
    )
    args, _ = parser.parse_known_args()
    args.out_dir = V4_EXPORT
    config = load_json(args.config)
    actions = load_json(args.matrix)
    if args.g07_result.exists():
        apply_g07(config, load_json(args.g07_result))
    V4_CONFIG.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    V4_MATRIX.write_text(json.dumps(actions, indent=2) + "\n", encoding="utf-8")

    tree, _ = build_model(config)
    tree.write(args.model, encoding="utf-8", xml_declaration=True)
    evaluator = V2Evaluator(args.model, config, actions)
    evaluation = evaluator.evaluate_all()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "evaluation_v4.json").write_text(
        json.dumps(evaluation, indent=2) + "\n", encoding="utf-8"
    )
    matrix = export_action_matrix(args.matrix, evaluation, actions)
    export_action_markdown(
        PROJECT_ROOT / "design" / "ACTION_CABLE_MATRIX_V4.md", matrix
    )
    export_markdown(
        PROJECT_ROOT / "design" / "CABLE_ANCHORS_V4.md", evaluation, config
    )
    export_routes_csv(args.out_dir / "cable_routes_v4.csv", config)
    export_solidworks_csv(
        args.out_dir / "SolidWorks_hardpoints_v4.csv", config
    )
    export_tension_report(
        args.out_dir / "tension_report_v4.json", evaluation
    )
    export_route_validation(
        args.out_dir / "route_validation_v4.json", config, args.model
    )
    k7 = next(action for action in evaluation["actions"] if action["id"] == "K7")
    g07 = next(
        route for route in k7["final"]["routes"] if route["id"] == "G07"
    )
    repair = {
        "G07": {
            "route_count": len(config["routes"]),
            "selected_in_k7": "G07" in k7["selected_cables"],
            "k7_relative_residual": k7["final"]["tension_solution"].get(
                "relative_residual"
            ),
            "k7_positive_spanning": k7["final"]["tension_solution"].get(
                "positive_spanning", {}
            ).get("pass"),
            "k7_max_tension_N": k7["final"]["tension_solution"].get(
                "max_tension_N"
            ),
            "k7_contraction_percent": k7[
                "max_selected_contraction_percent"
            ],
            "k7_clearance_mm": k7["min_selected_cable_clearance_mm"],
            "g07_contraction_percent": g07["contraction_percent"],
            "g07_clearance_mm": g07["thorax_clearance_mm"],
            "g07_wrap_points": g07["actual_wrap_points"],
        },
        "v3_fallback_preserved": {
            "scap_int_deg": -5.0,
            "full_target_deg": -20.0,
        },
    }
    (args.out_dir / "route_repair_report_v4.json").write_text(
        json.dumps(repair, indent=2) + "\n", encoding="utf-8"
    )
    for action in evaluation["actions"]:
        print(
            f"{action['id']:<10} selected={','.join(action['selected_cables']):<32} "
            f"80={action['pass_80_percent']} min={action['pass_minimum_visible']} "
            f"tension={action['pass_tension']} contraction={action['pass_contraction']} "
            f"clearance={action['pass_cable_clearance']}"
        )
    print(f"wrote {args.model}")
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
