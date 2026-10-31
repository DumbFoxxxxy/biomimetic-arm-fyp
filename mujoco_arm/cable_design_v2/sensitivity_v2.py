"""TPU stiffness sensitivity for the V2 cable allocation."""

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

from evaluate_v2 import V2Evaluator
from v2_common import (
    ACTION_MATRIX_PATH,
    CONFIG_PATH,
    EXPORT_ROOT,
    MODEL_PATH,
    load_json,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument("--matrix", type=Path, default=ACTION_MATRIX_PATH)
    parser.add_argument(
        "--out", type=Path, default=EXPORT_ROOT / "sensitivity_v2.json"
    )
    args = parser.parse_args()
    config = load_json(args.config)
    actions = load_json(args.matrix)
    scales = [0.5, 1.0, 2.0]
    guide_frictions = [0.0, 0.05, 0.10]
    reports = {}
    for scale in scales:
        evaluator = V2Evaluator(
            args.model,
            config,
            actions,
            tpu_scale=scale,
        )
        evaluation = evaluator.evaluate_all()
        reports[str(scale)] = {
            action["id"]: {
                "selected_cables": action["selected_cables"],
                "pass_minimum_visible": action["pass_minimum_visible"],
                "pass_tension": action["pass_tension"],
                "pass_minimum_functional": action["pass_minimum_functional"],
            }
            for action in evaluation["actions"]
        }
    guide_reports = {}
    for friction in guide_frictions:
        evaluator = V2Evaluator(
            args.model,
            config,
            actions,
            tpu_scale=1.0,
            guide_friction=friction,
        )
        evaluation = evaluator.evaluate_all()
        guide_reports[str(friction)] = {
            action["id"]: {
                "selected_cables": action["selected_cables"],
                "pass_tension": action["pass_tension"],
                "pass_minimum_functional": action["pass_minimum_functional"],
            }
            for action in evaluation["actions"]
        }
    robust_sets = {}
    for action in actions["actions"]:
        action_id = action["id"]
        sets = [
            set(reports[str(scale)][action_id]["selected_cables"])
            for scale in scales
        ]
        robust_sets[action_id] = {
            "intersection": sorted(set.intersection(*sets)),
            "union": sorted(set.union(*sets)),
            "stable": all(current == sets[0] for current in sets[1:]),
        }
    output = {
        "tpu_scales": scales,
        "by_scale": reports,
        "guide_friction_values": guide_frictions,
        "by_guide_friction": guide_reports,
        "robust_cable_sets": robust_sets,
        "contact_friction_status": (
            "not applied in the current contact-free equilibrium QP; "
            "requires a dynamic contact-force model"
        ),
        "guide_friction_status": (
            "approximated by exp(-mu * total_declared_wrap); a segment-level "
            "capstan model remains future work"
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out}")
    for action_id, data in robust_sets.items():
        print(
            f"{action_id:<10} stable={data['stable']} "
            f"intersection={','.join(data['intersection'])}"
        )


if __name__ == "__main__":
    main()
