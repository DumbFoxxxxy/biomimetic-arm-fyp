"""Optimize G02 guide topology across K2 and K7 for V3."""

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
import xml.etree.ElementTree as ET

import mujoco
import numpy as np
from scipy.optimize import differential_evolution

from evaluate_v2 import V2Evaluator, r400_center
from v2_common import (
    ACTION_MATRIX_PATH,
    CONFIG_PATH,
    EXPORT_ROOT,
    MODEL_PATH,
    POINTS_MM,
    design_to_mj,
    load_json,
    mj_to_design,
    normalize,
    quat_from_z_axis,
    vec,
)


def scapula_basis() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u_axis = normalize(POINTS_MM["SUP"] - POINTS_MM["INF"])
    lateral = POINTS_MM["GH"] - POINTS_MM["P"]
    v_axis = normalize(lateral - np.dot(lateral, u_axis) * u_axis)
    w_axis = normalize(np.array([0.6556, 0.0, 0.7552]))
    return u_axis, v_axis, w_axis


def find_body(root: ET.Element, name: str) -> ET.Element:
    for body in root.iter("body"):
        if body.get("name") == name:
            return body
    raise KeyError(name)


def build_template(owner: str, model_path: Path) -> tuple[mujoco.MjModel, dict]:
    tree = ET.parse(model_path)
    root = tree.getroot()
    thorax = find_body(root, "thorax")
    scapula = find_body(root, "scapula")
    humerus = find_body(root, "humerus")
    # Remove the current G02 guide from S and from the tendon path.
    for body in (scapula, humerus):
        for geom in list(body.findall("geom")):
            if (geom.get("name") or "").startswith("V2_guide_G02_"):
                body.remove(geom)
        for site in list(body.findall("site")):
            if site.get("name") == "V2_side_G02":
                body.remove(site)

    parent = scapula if owner == "S" else humerus
    guide_name = f"v3_G02_guide_{owner}"
    side_name = f"v3_G02_side_{owner}"
    ET.SubElement(
        parent,
        "geom",
        {
            "name": guide_name,
            "type": "cylinder",
            "size": "0.008 0.008",
            "pos": "0 0 0",
            "contype": "0",
            "conaffinity": "0",
            "density": "0",
        },
    )
    ET.SubElement(
        parent,
        "site",
        {
            "name": side_name,
            "pos": "0 0.012 0",
            "size": "0.002",
            "rgba": "1 0 0 1",
        },
    )
    for spatial in root.iter("spatial"):
        if spatial.get("name") != "cable_G02":
            continue
        for child in list(spatial):
            spatial.remove(child)
        ET.SubElement(spatial, "site", {"site": "V2_A_G02"})
        ET.SubElement(
            spatial,
            "geom",
            {"geom": guide_name, "sidesite": side_name},
        )
        ET.SubElement(spatial, "site", {"site": "V2_B_G02"})
    model = _load_project_xml(ET.tostring(root, encoding="unicode"))
    return model, {
        "owner": owner,
        "guide_geom": guide_name,
        "guide_site": side_name,
    }


def make_alt_actions(actions: dict) -> dict:
    result = copy.deepcopy(actions)
    k7 = next(item for item in result["actions"] if item["id"] == "K7")
    k7["allowed_cables"] = [
        route_id for route_id in k7["allowed_cables"] if route_id != "G01"
    ]
    k7["mandatory_cables"] = ["S01", "S02", "S03", "S04", "G02"]
    return result


def evaluate_candidate(
    evaluator: V2Evaluator,
    owner: str,
    template: dict,
    action_matrix: dict,
    a_base_uvw: np.ndarray,
    b_base_local: np.ndarray,
    guide_base_local: np.ndarray,
    values: np.ndarray,
) -> dict:
    u = float(a_base_uvw[0] + values[0])
    v = float(a_base_uvw[1] + values[1])
    w = float(np.clip(a_base_uvw[2] + values[2], 0.0, 10.0))
    b_dx = float(values[3])
    b_dy = float(values[4])
    g_dx = float(values[5])
    g_dy = float(values[6])
    g_dz = float(values[7])

    u_axis, v_axis, w_axis = scapula_basis()
    a_design = POINTS_MM["P"] + u * u_axis + v * v_axis + w * w_axis
    b_local = b_base_local.copy()
    b_local[0] += b_dx * 0.001
    b_local[1] += b_dy * 0.001
    guide_local = guide_base_local + np.array([g_dx, g_dy, g_dz]) * 0.001

    aid = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_SITE, "V2_A_G02"
    )
    bid = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_SITE, "V2_B_G02"
    )
    gid = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_GEOM, template["guide_geom"]
    )
    sid = mujoco.mj_name2id(
        evaluator.model, mujoco.mjtObj.mjOBJ_SITE, template["guide_site"]
    )
    evaluator.model.site_pos[aid] = design_to_mj(a_design - POINTS_MM["AC"])
    evaluator.model.site_pos[bid] = b_local
    evaluator.model.geom_pos[gid] = guide_local

    if owner == "S":
        g_design = POINTS_MM["AC"] + np.array(
            [guide_local[0], -guide_local[2], guide_local[1]]
        ) * 1000.0
        b_design = POINTS_MM["GH"] + np.array(
            [b_local[0], -b_local[2], b_local[1]]
        ) * 1000.0
        body_name = "scapula"
    else:
        g_design = POINTS_MM["GH"] + np.array(
            [guide_local[0], -guide_local[2], guide_local[1]]
        ) * 1000.0
        b_design = POINTS_MM["GH"] + np.array(
            [b_local[0], -b_local[2], b_local[1]]
        ) * 1000.0
        body_name = "humerus"
    incoming = g_design - a_design
    outgoing = b_design - g_design
    axis = normalize(np.cross(incoming, outgoing))
    side = normalize(
        np.cross(axis, normalize(incoming) + normalize(outgoing))
    )
    axis_mj = normalize(np.array([axis[0], -axis[2], axis[1]]))
    side_mj = normalize(np.array([side[0], -side[2], side[1]]))
    evaluator.model.geom_quat[gid] = quat_from_z_axis(axis_mj)
    evaluator.model.site_pos[sid] = guide_local + side_mj * 0.012
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
    route = evaluator.route_by_id["G02"]
    route["path"][1]["world_mm"] = g_design.tolist()
    route["path"][1]["local_mm"] = guide_local.tolist()
    route["path"][1]["axis_design"] = axis.tolist()
    route["path"][1]["side_design"] = side.tolist()
    route["path"][1]["wrap_angle_deg"] = wrap_angle

    zero = mujoco.MjData(evaluator.model)
    zero.qpos[:] = 0.0
    mujoco.mj_forward(evaluator.model, zero)
    evaluator.neutral_lengths = {
        route_id: float(zero.ten_length[tendon_id])
        for route_id, tendon_id in evaluator.tendon_ids.items()
    }

    action_k2 = next(
        item for item in action_matrix["actions"] if item["id"] == "K2"
    )
    action_k7 = next(
        item for item in action_matrix["actions"] if item["id"] == "K7"
    )
    route_checks = []
    for action in (action_k2, action_k7):
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            data = evaluator._action_pose(action, fraction)
            mujoco.mj_forward(evaluator.model, data)
            route_data = evaluator._route_geometry(route, data)
            route_checks.append(
                {
                    "action": action["id"],
                    "fraction": fraction,
                    "wrap_points": route_data["actual_wrap_points"],
                    "clearance_mm": route_data["thorax_clearance_mm"],
                    "contraction_percent": route_data["contraction_percent"],
                }
            )
    final_data = evaluator._action_pose(action_k7, 1.0)
    tension = evaluator._solve_tensions(action_k7, final_data)
    max_contraction = max(item["contraction_percent"] for item in route_checks)
    min_clearance = min(item["clearance_mm"] for item in route_checks)
    inactive = sum(item["wrap_points"] < 4 for item in route_checks)
    score = 0.0
    score += 5.0 * inactive
    score += 5.0 * max(0.0, 2.0 - min_clearance)
    score += 0.20 * max(0.0, max_contraction - 18.0)
    score += 2.0 * max(0.0, float(tension.get("relative_residual", 10.0)) - 0.1)
    score += 0.0 if tension.get("positive_spanning", {}).get("pass", False) else 2.0
    score += 0.1 * max(0.0, float(tension.get("max_tension_N", 0.0)) - 149.0)
    return {
        "score": score,
        "values": values.tolist(),
        "owner": owner,
        "min_clearance_mm": min_clearance,
        "max_contraction_percent": max_contraction,
        "inactive_wrap_count": inactive,
        "tension_ok": bool(tension.get("ok", False)),
        "positive_spanning": tension.get("positive_spanning", {}),
        "relative_residual": tension.get("relative_residual"),
        "max_tension_N": tension.get("max_tension_N"),
        "checks": route_checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument("--matrix", type=Path, default=ACTION_MATRIX_PATH)
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument(
        "--out", type=Path, default=EXPORT_ROOT / "g02_optimization_v3.json"
    )
    parser.add_argument("--maxiter", type=int, default=14)
    parser.add_argument("--popsize", type=int, default=4)
    args = parser.parse_args()
    config = load_json(args.config)
    actions = load_json(args.matrix)
    alt_actions = make_alt_actions(actions)
    g02 = next(route for route in config["routes"] if route["id"] == "G02")
    a_world = np.asarray(g02["path"][0]["world_mm"], dtype=float)
    q = a_world - POINTS_MM["P"]
    u_axis, v_axis, w_axis = scapula_basis()
    a_base_uvw = np.array(
        [np.dot(q, u_axis), np.dot(q, v_axis), np.dot(q, w_axis)]
    )
    b_world = np.asarray(g02["path"][-1]["world_mm"], dtype=float)
    b_local = design_to_mj(b_world - POINTS_MM["GH"])
    results = []
    bounds = [
        (-20.0, 20.0),
        (-20.0, 20.0),
        (-10.0, 10.0),
        (-20.0, 20.0),
        (-20.0, 20.0),
        (-30.0, 30.0),
        (-30.0, 30.0),
        (-30.0, 30.0),
    ]
    for owner in ("S", "H"):
        model, template = build_template(owner, args.model)
        evaluator = V2Evaluator(None, config, alt_actions, model=model)
        if owner == "S":
            guide_base = design_to_mj(
                np.asarray(g02["path"][1]["world_mm"]) - POINTS_MM["AC"]
            )
        else:
            guide_base = design_to_mj(
                np.asarray(g02["path"][1]["world_mm"]) - POINTS_MM["GH"]
            )
        history = []

        def objective(values):
            result = evaluate_candidate(
                evaluator,
                owner,
                template,
                alt_actions,
                a_base_uvw,
                b_local,
                guide_base,
                values,
            )
            history.append(result)
            return result["score"]

        optimized = differential_evolution(
            objective,
            bounds,
            maxiter=args.maxiter,
            popsize=args.popsize,
            seed=53 if owner == "S" else 59,
            polish=False,
            workers=1,
            updating="immediate",
        )
        best = min(history, key=lambda item: item["score"])
        best["result_x"] = optimized.x.tolist()
        best["success"] = bool(optimized.success)
        results.append(best)
    results.sort(key=lambda item: item["score"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps({"best": results[0], "all": results}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(results[0], indent=2))


if __name__ == "__main__":
    main()
