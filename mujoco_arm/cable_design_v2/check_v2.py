"""Structural and functional checks for the V2 cable-design model."""

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

import mujoco

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
        "--evaluation",
        type=Path,
        default=EXPORT_ROOT / "evaluation_v2.json",
    )
    parser.add_argument(
        "--out", type=Path, default=EXPORT_ROOT / "check_report_v2.json"
    )
    parser.add_argument("--expected-routes", type=int, default=16)
    parser.add_argument("--expected-tendons", type=int, default=16)
    args = parser.parse_args()
    config = load_json(args.config)
    matrix = load_json(args.matrix)
    evaluation = load_json(args.evaluation)
    model = _load_project_model(str(args.model))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)

    failures = []
    if len(config["routes"]) != args.expected_routes:
        failures.append(
            f"configured route count is not {args.expected_routes}"
        )
    if model.ntendon != args.expected_tendons:
        failures.append(
            f"model tendon count is not {args.expected_tendons}"
        )
    if model.nq != 7:
        failures.append("model DOF count is not 7")
    for route in config["routes"]:
        guides = [
            element for element in route["path"] if element["kind"] == "guide"
        ]
        if len(guides) > int(config["guides"]["maximum_per_cable"]):
            failures.append(f"{route['id']} exceeds guide-count limit")
        for guide in guides:
            if float(guide["radius_mm"]) < float(
                config["guides"]["minimum_radius_mm"]
            ):
                failures.append(f"{route['id']} has an undersized guide")
            if float(guide.get("wrap_angle_deg", 0.0)) < 10.0:
                failures.append(f"{route['id']} has a fake low-wrap guide")

    action_results = evaluation["actions"]
    minimum_motion_failures = [
        action["id"]
        for action in action_results
        if not action["pass_minimum_visible"]
    ]
    tension_failures = [
        action["id"]
        for action in action_results
        if not action["pass_tension"]
    ]
    contraction_failures = [
        action["id"]
        for action in action_results
        if not action["pass_contraction"]
    ]
    clearance_failures = [
        action["id"]
        for action in action_results
        if not action["pass_cable_clearance"]
    ]
    visual_failures = [
        action["id"]
        for action in action_results
        if not action["pass_visual_clearance"]
    ]
    functional_failures = [
        action["id"]
        for action in action_results
        if not action["pass_minimum_functional"]
    ]
    structural_pass = not failures
    minimum_motion_pass = not minimum_motion_failures
    full_pass = bool(
        structural_pass
        and minimum_motion_pass
        and not tension_failures
        and not contraction_failures
        and not clearance_failures
        and not visual_failures
    )
    report = {
        "model": str(args.model),
        "nq": int(model.nq),
        "nu": int(model.nu),
        "ntendon": int(model.ntendon),
        "nwrap": int(model.nwrap),
        "route_count": len(config["routes"]),
        "structural_pass": structural_pass,
        "minimum_motion_pass": minimum_motion_pass,
        "full_pass": full_pass,
        "structural_failures": failures,
        "minimum_motion_failures": minimum_motion_failures,
        "tension_failures": tension_failures,
        "contraction_failures": contraction_failures,
        "cable_clearance_failures": clearance_failures,
        "visual_clearance_failures": visual_failures,
        "minimum_functional_failures": functional_failures,
        "actions_in_matrix": len(matrix["actions"]),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        f"structural={structural_pass} minimum_motion={minimum_motion_pass} "
        f"full={full_pass}"
    )
    print(f"tension failures: {', '.join(tension_failures) or 'none'}")
    print(f"contraction failures: {', '.join(contraction_failures) or 'none'}")
    print(f"cable clearance failures: {', '.join(clearance_failures) or 'none'}")
    print(f"visual clearance failures: {', '.join(visual_failures) or 'none'}")
    print(f"report: {args.out}")
    if not structural_pass or not minimum_motion_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
