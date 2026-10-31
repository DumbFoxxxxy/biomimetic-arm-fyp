"""Audit V5 and export capstan, contact, ellipsoid and SolidWorks data."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import csv
import json
from project_runtime import WorkspacePath as Path
import sys

import mujoco
import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from v2_common import load_json, mj_to_design, normalize  # noqa: E402
from v5_evaluator import V5Evaluator  # noqa: E402


CONFIG = ROOT / "cable_design_v5.json"
MATRIX = ROOT / "action_cable_matrix_v5.json"
MODEL = ROOT / "arm_cable_design_v5.xml"
EXPORT = ROOT.parent / "exports" / "cable_design_v5"


def fit_ellipsoid(points: np.ndarray) -> dict:
    center_guess = np.mean(points, axis=0)
    scale = max(float(np.ptp(points, axis=0).max()), 20.0)
    initial = np.r_[center_guess, scale, scale, scale]
    bounds_low = np.r_[center_guess - 500.0, 10.0, 10.0, 10.0]
    bounds_high = np.r_[center_guess + 500.0, 1000.0, 1000.0, 1000.0]

    def residual(values):
        center = values[:3]
        radii = values[3:]
        normalized = (points - center) / radii
        return np.linalg.norm(normalized, axis=1) - 1.0

    result = least_squares(
        residual,
        initial,
        bounds=(bounds_low, bounds_high),
        max_nfev=5000,
        x_scale="jac",
    )
    radii = result.x[3:]
    normalized = (points - result.x[:3]) / radii
    radial_distance = np.abs(np.linalg.norm(normalized, axis=1) - 1.0) * np.mean(
        radii
    )
    return {
        "success": bool(result.success),
        "center_mm": result.x[:3].tolist(),
        "radii_mm": radii.tolist(),
        "median_distance_mm": float(np.median(radial_distance)),
        "p95_distance_mm": float(np.percentile(radial_distance, 95)),
        "max_distance_mm": float(np.max(radial_distance)),
    }


def route_world_points(evaluator: V5Evaluator, config: dict, zero) -> list[dict]:
    rows = []
    for route in config["routes"]:
        guides = [
            element for element in route["path"] if element["kind"] == "guide"
        ]
        guide_index = 0
        for index, element in enumerate(route["path"]):
            if element["kind"] == "anchor":
                site_name = (
                    f"V2_A_{route['id']}"
                    if index == 0
                    else f"V2_B_{route['id']}"
                )
                entity_id = mujoco.mj_name2id(
                    evaluator.model, mujoco.mjtObj.mjOBJ_SITE, site_name
                )
                world = mj_to_design(zero.site_xpos[entity_id])
                axis = [0.0, 0.0, 0.0]
                radius = 0.0
            elif element["kind"] == "pass_through":
                site_name = f"V2_pass_{route['id']}_{index}"
                entity_id = mujoco.mj_name2id(
                    evaluator.model, mujoco.mjtObj.mjOBJ_SITE, site_name
                )
                world = mj_to_design(zero.site_xpos[entity_id])
                axis = [0.0, 0.0, 0.0]
                radius = 0.0
            else:
                geom_name = f"V2_guide_{route['id']}_{element['name']}"
                entity_id = mujoco.mj_name2id(
                    evaluator.model, mujoco.mjtObj.mjOBJ_GEOM, geom_name
                )
                world = mj_to_design(zero.geom_xpos[entity_id])
                axis = normalize(
                    np.asarray(element.get("axis_design", [0.0, 0.0, 1.0]))
                )
                radius = float(element.get("radius_mm", 0.0))
                guide_index += 1
            normal = normalize(
                np.asarray(
                    element.get(
                        "side_design",
                        element.get("installation_normal", [0.0, 0.0, 1.0]),
                    )
                )
            )
            rows.append(
                {
                    "channel": route["id"],
                    "group": route["group"],
                    "kind": element["kind"],
                    "body": element["body"],
                    "name": f"{route['id']}_{index:02d}",
                    "world_x_mm": float(world[0]),
                    "world_y_mm": float(world[1]),
                    "world_z_mm": float(world[2]),
                    "local_x_m": float(element["local_mm"][0]),
                    "local_y_m": float(element["local_mm"][1]),
                    "local_z_m": float(element["local_mm"][2]),
                    "normal_x": float(normal[0]),
                    "normal_y": float(normal[1]),
                    "normal_z": float(normal[2]),
                    "axis_x": float(axis[0]),
                    "axis_y": float(axis[1]),
                    "axis_z": float(axis[2]),
                    "radius_mm": radius,
                    "pad_offset_mm": (
                        float(config.get("contact_pad_offset_mm", 0.0))
                        if element["body"] == "scapula"
                        else 0.0
                    ),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--evaluation", type=Path, default=EXPORT / "evaluation_v5.json")
    parser.add_argument("--out-dir", type=Path, default=EXPORT)
    args = parser.parse_args()
    config = load_json(args.config)
    matrix = load_json(args.matrix)
    evaluation = load_json(args.evaluation)
    evaluator = V5Evaluator(args.model, config, matrix)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    zero = evaluator._zero_data()
    body_sc = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_BODY, "clavicle_protraction"
    )
    body_ac = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_BODY, "scapula"
    )
    body_gh = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_BODY, "gh_flex"
    )
    body_elbow = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_BODY, "ulna"
    )
    bone_lengths = {
        "SC_AC_mm": float(
            np.linalg.norm(zero.xpos[body_ac] - zero.xpos[body_sc]) * 1000.0
        ),
        "GH_E_mm": float(
            np.linalg.norm(zero.xpos[body_elbow] - zero.xpos[body_gh]) * 1000.0
        ),
        "E_W_design_mm": 237.0,
    }

    capstan = {"mu": float(config["friction_model"]["mu"]), "actions": {}}
    contact = {"actions": {}, "zero": {}}
    ellipsoids = {"actions": {}}
    for action in matrix["actions"]:
        data = evaluator._action_pose(action, 1.0)
        final = evaluator.evaluate_pose_sample(action, 1.0, True)
        solution = final["tension_solution"]
        selected = solution.get("selected", [])
        route_reports = {}
        for route_id in selected:
            route = evaluator.route_by_id[route_id]
            route_data = evaluator._route_geometry(route, data)
            tension = float(solution["tensions_N"][route_id])
            guides = []
            previous = tension
            for guide in route_data["capstan_guides"]:
                outgoing = previous * guide["incoming_factor"]
                guides.append(
                    {
                        **guide,
                        "incoming_tension_N": previous,
                        "outgoing_tension_N": outgoing,
                    }
                )
                previous = outgoing
            route_reports[route_id] = {
                "actuator_tension_N": tension,
                "distal_tension_N": previous,
                "guides": guides,
            }
        capstan["actions"][action["id"]] = route_reports
        contact["actions"][action["id"]] = {
            "min_gap_mm": final["contact_gap_mm"],
            "visual_clearance_mm": final["visual_clearance_mm"],
        }
        pad = evaluator.contact_pad_geom
        points = []
        for fraction in config["acceptance"]["final_pose_fractions"]:
            pose = evaluator._action_pose(action, float(fraction))
            points.append(mj_to_design(pose.geom_xpos[pad]))
        ellipsoids["actions"][action["id"]] = fit_ellipsoid(np.asarray(points))
    contact["zero"] = {
        "visual_clearance_mm": evaluator._visual_clearance(zero),
        "contact_gap_mm": evaluator._contact_gap(zero),
        "pad_offset_mm": config.get("contact_pad_offset_mm", 0.0),
    }
    (args.out_dir / "capstan_report_v5.json").write_text(
        json.dumps(capstan, indent=2) + "\n", encoding="utf-8"
    )
    (args.out_dir / "contact_audit_v5.json").write_text(
        json.dumps(contact, indent=2) + "\n", encoding="utf-8"
    )
    (args.out_dir / "ellipsoid_audit_v5.json").write_text(
        json.dumps(ellipsoids, indent=2) + "\n", encoding="utf-8"
    )
    hardpoints = route_world_points(evaluator, config, zero)
    hardpoint_path = args.out_dir / "SolidWorks_hardpoints_v5.csv"
    with hardpoint_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(hardpoints[0]))
        writer.writeheader()
        writer.writerows(hardpoints)

    failures = {
        action["id"]: {
            key: value
            for key, value in action.items()
            if key.startswith("pass_") and not value
        }
        for action in evaluation["actions"]
    }
    structural = bool(
        evaluator.model.nq == 7
        and evaluator.model.ntendon == evaluator.model.nu
        and int(config["channel_plan"]["minimum_count"])
        <= evaluator.model.ntendon
        <= int(config["channel_plan"]["maximum_count"])
        and abs(bone_lengths["SC_AC_mm"] - 150.0) <= 0.05
        and abs(bone_lengths["GH_E_mm"] - 290.0) <= 0.05
    )
    full_pass = structural and all(not value for value in failures.values())
    report = {
        "model": str(args.model),
        "structural_pass": structural,
        "joint_count": int(evaluator.model.nq),
        "channel_count": int(evaluator.model.ntendon),
        "actuator_count": int(evaluator.model.nu),
        "full_pass": full_pass,
        "frozen_hardpoints": full_pass,
        "failure_actions": {
            action_id: values
            for action_id, values in failures.items()
            if values
        },
        "zero_pose_offsets_deg": config.get("zero_pose_offsets_deg", {}),
        "contact_pad_offset_mm": config.get("contact_pad_offset_mm", 0.0),
        "bone_lengths": bone_lengths,
    }
    (args.out_dir / "check_report_v5.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
