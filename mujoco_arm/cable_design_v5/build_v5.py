"""Build and evaluate the V5 cable-design model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(ROOT))

from build_v2_model import build as build_model  # noqa: E402
from evaluate_v2 import (  # noqa: E402
    export_action_markdown,
    export_action_matrix,
    export_routes_csv,
    export_solidworks_csv,
    export_tension_report,
)
from v2_common import PROJECT_ROOT, design_to_mj, load_json, normalize  # noqa: E402
from v5_evaluator import V5Evaluator  # noqa: E402


V5_CONFIG = ROOT / "cable_design_v5.json"
V5_MATRIX = ROOT / "action_cable_matrix_v5.json"
V5_MODEL = ROOT / "arm_cable_design_v5.xml"
V5_EXPORT = ROOT.parent / "exports" / "cable_design_v5"


def xml_vec(values: np.ndarray) -> str:
    return " ".join(f"{float(value):.9f}" for value in values)


def apply_pad_offset(tree, config: dict) -> None:
    offset_mm = float(config.get("contact_pad_offset_mm", 0.0))
    if abs(offset_mm) < 1e-12:
        return
    pad = tree.find(".//geom[@name='v2_scapula_contact_pad']")
    if pad is None:
        raise RuntimeError("scapula contact pad geom is missing")
    current = np.fromstring(pad.attrib["pos"], sep=" ")
    outward = normalize(np.array([0.6556, 0.0, 0.7552]))
    pad.attrib["pos"] = xml_vec(
        current + design_to_mj(outward * offset_mm) * 0.001
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=V5_CONFIG)
    parser.add_argument("--matrix", type=Path, default=V5_MATRIX)
    parser.add_argument("--model", type=Path, default=V5_MODEL)
    parser.add_argument("--out-dir", type=Path, default=V5_EXPORT)
    args = parser.parse_args()

    config = load_json(args.config)
    actions = load_json(args.matrix)
    tree, _ = build_model(config)
    apply_pad_offset(tree, config)
    tree.write(args.model, encoding="utf-8", xml_declaration=True)
    evaluator = V5Evaluator(args.model, config, actions)
    evaluation = evaluator.evaluate_all()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "evaluation_v5.json").write_text(
        json.dumps(evaluation, indent=2) + "\n", encoding="utf-8"
    )
    matrix = export_action_matrix(args.matrix, evaluation, actions)
    args.matrix.write_text(
        json.dumps(matrix, indent=2) + "\n", encoding="utf-8"
    )
    export_action_markdown(
        PROJECT_ROOT / "design" / "ACTION_CABLE_MATRIX_V5.md", matrix
    )
    export_routes_csv(args.out_dir / "cable_routes_v5.csv", config)
    export_solidworks_csv(
        args.out_dir / "SolidWorks_hardpoints_v5.csv", config
    )
    export_tension_report(
        args.out_dir / "tension_report_v5.json", evaluation
    )
    report = {
        "model": str(args.model),
        "route_count": len(config["routes"]),
        "joint_count": int(evaluator.model.nq),
        "zero_pose_offsets_deg": config.get("zero_pose_offsets_deg", {}),
        "contact_pad_offset_mm": config.get("contact_pad_offset_mm", 0.0),
        "actions": {
            action["id"]: {
                "pass_80_percent": action["pass_80_percent"],
                "pass_tension": action["pass_tension"],
                "pass_positive_spanning": action["pass_positive_spanning"],
                "pass_contraction": action["pass_contraction"],
                "pass_cable_clearance": action["pass_cable_clearance"],
                "pass_guide_geometry": action["pass_guide_geometry"],
                "pass_visual_clearance": action["pass_visual_clearance"],
                "pass_contact": action["pass_contact"],
                "visual_clearance_mm": action["visual_clearance_mm"],
                "max_contraction_percent": action[
                    "max_selected_contraction_percent"
                ],
                "min_cable_clearance_mm": action[
                    "min_selected_cable_clearance_mm"
                ],
                "contact_gap_min_mm": action["contact_gap_min_mm"],
                "selected_cables": action["selected_cables"],
            }
            for action in evaluation["actions"]
        },
    }
    (args.out_dir / "summary_v5.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"model={args.model} routes={len(config['routes'])} "
        f"joints={evaluator.model.nq} tendons={evaluator.model.ntendon}"
    )
    for action in evaluation["actions"]:
        print(
            f"{action['id']:<10} 80={action['pass_80_percent']} "
            f"tension={action['pass_tension']} "
            f"span={action['pass_positive_spanning']} "
            f"contract={action['pass_contraction']} "
            f"visual={action['pass_visual_clearance']} "
            f"contact={action['pass_contact']}"
        )


if __name__ == "__main__":
    main()
