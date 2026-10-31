"""Audit V6 static poses, elbow antagonism and continuous trajectories."""

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


ROOT = Path(__file__).resolve().parent
V5_ROOT = ROOT.parent / "cable_design_v5"
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v6"
sys.path[:0] = [str(ROOT), str(V5_ROOT), str(V2_ROOT)]

from evaluate_v2 import dense_tendon_jacobian  # noqa: E402
from trajectory_v6 import TrajectoryRunner  # noqa: E402
from v2_common import load_json  # noqa: E402
from v6_evaluator import V6Evaluator  # noqa: E402


CONFIG = ROOT / "cable_design_v6.json"
MATRIX = ROOT / "action_cable_matrix_v6.json"
MODEL = ROOT / "arm_cable_design_v6.xml"
LIBRARY = ROOT / "trajectory_library_v6.json"


def elbow_audit(
    evaluator: V6Evaluator,
    matrix: dict,
) -> dict:
    action = next(item for item in matrix["actions"] if item["id"] == "K1")
    fractions = np.linspace(0.0, 1.0, 21)
    samples = []
    for fraction in fractions:
        data = evaluator._action_pose(action, float(fraction))
        dofs = evaluator._dofs_for_action(action)
        columns = [
            evaluator.tendon_indices[route_id] for route_id in ("E01", "E02")
        ]
        jacobian = dense_tendon_jacobian(evaluator.model, data)
        moments = -jacobian[:, dofs][columns].T
        moment_values = moments.ravel()
        routes = [
            evaluator._route_geometry(evaluator.route_by_id[route_id], data)
            for route_id in ("E01", "E02")
        ]
        samples.append(
            {
                "fraction": float(fraction),
                "moment_m_per_n": moment_values.tolist(),
                "effective_moment_arm_mm": [
                    abs(float(moment_values[0])) * 1000.0,
                    abs(float(moment_values[1])) * 1000.0,
                ],
                "opposite": bool(moment_values[0] * moment_values[1] < 0.0),
                "wrap_points": [
                    int(route["actual_wrap_points"]) for route in routes
                ],
                "contraction_percent": [
                    float(route["contraction_percent"]) for route in routes
                ],
            }
        )
    minimum_arm = min(
        min(sample["effective_moment_arm_mm"]) for sample in samples
    )
    opposite_all = all(sample["opposite"] for sample in samples)
    wrap_all = all(
        all(count >= 4 for count in sample["wrap_points"])
        for sample in samples
    )
    contraction = max(
        max(sample["contraction_percent"]) for sample in samples
    )
    required_arm = float(
        evaluator.config["acceptance"][
            "elbow_effective_moment_arm_mm_min"
        ]
    )
    return {
        "pass": bool(
            opposite_all
            and wrap_all
            and minimum_arm >= required_arm
            and contraction
            <= float(
                evaluator.config["acceptance"][
                    "cable_contraction_percent_max"
                ]
            )
        ),
        "pass_opposite": opposite_all,
        "pass_wrap": wrap_all,
        "pass_moment_arm": minimum_arm >= required_arm,
        "pass_contraction": contraction
        <= float(
            evaluator.config["acceptance"][
                "cable_contraction_percent_max"
            ]
        ),
        "required_arm_mm": required_arm,
        "minimum_arm_mm": float(minimum_arm),
        "max_contraction_percent": float(contraction),
        "samples": samples,
    }


def trajectory_pass(report: dict, config: dict) -> bool:
    continuity_limit_deg = float(
        config["acceptance"].get("max_joint_step_per_sample_deg", 8.0)
    )
    continuity_limit_mm = float(
        config["acceptance"].get("max_wrist_step_per_sample_mm", 30.0)
    )
    return bool(
        report["pass_path"]
        and report["pass_clearance"]
        and report["all_tensions_ok"]
        and report["max_joint_step_deg"] <= continuity_limit_deg
        and report["max_wrist_step_mm"] <= continuity_limit_mm
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--library", type=Path, default=LIBRARY)
    parser.add_argument("--out-dir", type=Path, default=EXPORT)
    args = parser.parse_args()
    config = load_json(args.config)
    matrix = load_json(args.matrix)
    library = load_json(args.library)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    evaluator = V6Evaluator(args.model, config, matrix)
    static = evaluator.evaluate_all()
    (args.out_dir / "static_evaluation_v6.json").write_text(
        json.dumps(static, indent=2) + "\n",
        encoding="utf-8",
    )

    elbow = elbow_audit(evaluator, matrix)
    runner = TrajectoryRunner(
        args.model,
        config,
        matrix,
        library,
    )
    for trajectory in library["trajectories"]:
        for waypoint in trajectory["waypoints"]:
            qpos, _ = runner.q_for_trajectory(
                trajectory,
                float(waypoint["time_s"]),
            )
            reference = runner.reference_state(qpos)
            waypoint["wrist_target_mm"] = runner.wrist_design_mm(
                reference
            ).tolist()
    args.library.write_text(
        json.dumps(library, indent=2) + "\n",
        encoding="utf-8",
    )
    trajectories = []
    for trajectory in library["trajectories"]:
        report = runner.run(
            trajectory["id"],
            speed=1.0,
            include_samples=True,
        )
        report["pass_continuity"] = bool(
            report["max_joint_step_deg"]
            <= float(
                config["acceptance"].get(
                    "max_joint_step_per_sample_deg", 8.0
                )
            )
            and report["max_wrist_step_mm"]
            <= float(
                config["acceptance"].get(
                    "max_wrist_step_per_sample_mm", 30.0
                )
            )
        )
        report["pass_full"] = trajectory_pass(report, config)
        trajectories.append(report)

    preload_checks = []
    for preload in config["shoulder_guide"]["preload_sensitivity_N"]:
        runner.evaluator.preload_n = float(preload)
        result = runner.run(
            "FRONT_REACH",
            speed=1.0,
            include_samples=False,
        )
        preload_checks.append(
            {
                "preload_N": float(preload),
                "rms_error_mm": result["rms_error_mm"],
                "p95_error_mm": result["p95_error_mm"],
                "min_visual_clearance_mm": result[
                    "min_visual_clearance_mm"
                ],
            }
        )
    runner.evaluator.preload_n = float(
        config["shoulder_guide"]["nominal_preload_N"]
    )

    body_elbow = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_BODY, "ulna"
    )
    body_gh = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_BODY, "gh_flex"
    )
    zero = evaluator._zero_data()
    bone_lengths = {
        "GH_E_mm": float(
            np.linalg.norm(zero.xpos[body_elbow] - zero.xpos[body_gh]) * 1000.0
        ),
        "E_W_design_mm": 237.0,
    }
    structural = bool(
        evaluator.model.nq == 7
        and evaluator.model.ntendon == evaluator.model.nu
        and int(config["channel_plan"]["minimum_count"])
        <= evaluator.model.ntendon
        <= int(config["channel_plan"]["maximum_count"])
        and abs(bone_lengths["GH_E_mm"] - 290.0) <= 0.05
    )
    static_actions = {action["id"]: action for action in static["actions"]}
    static_pass = all(
        action["pass_minimum_functional"] for action in static_actions.values()
    )
    trajectories_pass = all(item["pass_full"] for item in trajectories)
    preload_stable = bool(
        max(item["p95_error_mm"] for item in preload_checks)
        - min(item["p95_error_mm"] for item in preload_checks)
        <= 15.0
    )
    full_pass = bool(
        structural
        and static_pass
        and elbow["pass"]
        and trajectories_pass
        and preload_stable
    )
    failures = []
    if not structural:
        failures.append("structural_model")
    if not static_pass:
        failures.extend(
            f"static:{action_id}"
            for action_id, action in static_actions.items()
            if not action["pass_minimum_functional"]
        )
    if not elbow["pass"]:
        if not elbow["pass_opposite"]:
            failures.append("elbow:not_opposite")
        if not elbow["pass_wrap"]:
            failures.append("elbow:wrap")
        if not elbow["pass_moment_arm"]:
            failures.append("elbow:moment_arm")
        if not elbow["pass_contraction"]:
            failures.append("elbow:contraction")
    failures.extend(
        f"trajectory:{item['id']}"
        for item in trajectories
        if not item["pass_full"]
    )
    if not preload_stable:
        failures.append("shoulder_preload_sensitivity")
    report = {
        "schema_version": 6,
        "model": str(args.model),
        "structural_pass": structural,
        "joint_count": int(evaluator.model.nq),
        "channel_count": int(evaluator.model.ntendon),
        "actuator_count": int(evaluator.model.nu),
        "static_actions_pass": static_pass,
        "elbow_antagonism": elbow,
        "trajectories": trajectories,
        "shoulder_preload_sensitivity": preload_checks,
        "bone_lengths": bone_lengths,
        "positive_spanning_hard_gate": False,
        "full_pass": full_pass,
        "frozen_hardpoints": full_pass,
        "failure_reasons": failures,
    }
    (args.out_dir / "trajectory_report_v6.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    check = {
        key: report[key]
        for key in (
            "schema_version",
            "model",
            "structural_pass",
            "joint_count",
            "channel_count",
            "actuator_count",
            "static_actions_pass",
            "positive_spanning_hard_gate",
            "full_pass",
            "frozen_hardpoints",
            "failure_reasons",
        )
    }
    check["elbow"] = {
        key: elbow[key]
        for key in (
            "pass",
            "pass_opposite",
            "pass_wrap",
            "pass_moment_arm",
            "pass_contraction",
            "minimum_arm_mm",
            "max_contraction_percent",
        )
    }
    check["trajectory_summary"] = {
        item["id"]: {
            "pass_full": item["pass_full"],
            "rms_error_mm": item["rms_error_mm"],
            "p95_error_mm": item["p95_error_mm"],
            "min_visual_clearance_mm": item["min_visual_clearance_mm"],
            "max_relative_residual": item["max_relative_residual"],
        }
        for item in trajectories
    }
    (args.out_dir / "check_report_v6.json").write_text(
        json.dumps(check, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(check, indent=2))


if __name__ == "__main__":
    main()
