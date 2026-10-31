"""Build V3 route repairs for G02, S03, S05 and S06."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

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
    r400_center,
)
from v2_common import (  # noqa: E402
    ACTION_MATRIX_PATH,
    CONFIG_PATH,
    MODEL_PATH,
    POINTS_MM,
    PROJECT_ROOT,
    design_to_mj,
    load_json,
    normalize,
    quat_from_z_axis,
)


V3_CONFIG = ROOT / "cable_design_v3.json"
V3_MATRIX = ROOT / "action_cable_matrix_v3.json"
V3_MODEL = ROOT / "arm_cable_design_v3.xml"
V3_EXPORT = ROOT.parent / "exports" / "cable_design_v3"


def scapula_basis() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u_axis = normalize(POINTS_MM["SUP"] - POINTS_MM["INF"])
    lateral = POINTS_MM["GH"] - POINTS_MM["P"]
    v_axis = normalize(lateral - np.dot(lateral, u_axis) * u_axis)
    w_axis = normalize(np.array([0.6556, 0.0, 0.7552]))
    return u_axis, v_axis, w_axis


def surface_point(config: dict, x: float, y: float, offset: float) -> tuple[np.ndarray, np.ndarray]:
    centre = r400_center(config)
    radius = float(config["contact"]["thorax_band"]["radius_mm"])
    z = centre[1] + math.sqrt(max(0.0, radius**2 - (x - centre[0]) ** 2))
    normal = normalize(np.array([x - centre[0], 0.0, z - centre[1]]))
    return np.array([x, y, z]) + offset * normal, normal


def scapula_point(u: float, v: float, w: float) -> np.ndarray:
    u_axis, v_axis, w_axis = scapula_basis()
    return POINTS_MM["P"] + u * u_axis + v * v_axis + w * w_axis


def anchor(
    owner: str,
    body: str,
    world: np.ndarray,
    normal: np.ndarray,
) -> dict:
    origin = {
        "T": np.zeros(3),
        "S": POINTS_MM["AC"],
        "H": POINTS_MM["GH"],
        "X": np.zeros(3),
    }[owner]
    return {
        "kind": "anchor",
        "owner": owner,
        "body": body,
        "world_mm": np.asarray(world).tolist(),
        "local_mm": design_to_mj(world - origin).tolist(),
        "installation_normal": normalize(np.asarray(normal)).tolist(),
    }


G02_STAGE1_VALUES = np.array(
    [
        4.0839983010990055,
        19.34474912695933,
        -1.1717758115214627,
        -11.47618779366657,
        9.17585566364302,
        14.920370215092651,
        4.31739000019018,
        -8.474390398454254,
    ]
)

G02_STAGE2_VALUES = np.array(
    [
        -18.31364184081789,
        17.159877636954725,
        7.220746643595312,
        -16.05133314308909,
        -12.349902496127578,
        -21.919720835581675,
        13.26409639742141,
        -21.17419698003337,
    ]
)


def apply_g02_values(config: dict, values: np.ndarray) -> dict:
    route = next(item for item in config["routes"] if item["id"] == "G02")
    u_axis, v_axis, w_axis = scapula_basis()
    a_old = np.asarray(route["path"][0]["world_mm"])
    a_q = a_old - POINTS_MM["P"]
    base_uvw = np.array(
        [np.dot(a_q, u_axis), np.dot(a_q, v_axis), np.dot(a_q, w_axis)]
    )
    a = scapula_point(
        base_uvw[0] + values[0],
        base_uvw[1] + values[1],
        np.clip(base_uvw[2] + values[2], 0.0, 10.0),
    )
    b_base = design_to_mj(
        np.asarray(route["path"][-1]["world_mm"]) - POINTS_MM["GH"]
    )
    b_local = b_base + np.array([values[3], values[4], 0.0]) * 0.001
    b = POINTS_MM["GH"] + np.array(
        [b_local[0], -b_local[2], b_local[1]]
    ) * 1000.0
    g_base = design_to_mj(
        np.asarray(route["path"][1]["world_mm"]) - POINTS_MM["AC"]
    )
    g_local = g_base + values[5:8] * 0.001
    g = POINTS_MM["AC"] + np.array(
        [g_local[0], -g_local[2], g_local[1]]
    ) * 1000.0
    axis = normalize(
        np.array([0.96732426, -0.05155616, 0.24824532])
    )
    side = normalize(
        np.array([0.24456068, -0.06856683, -0.96720663])
    )
    route["path"] = [
        anchor("S", "scapula", a, w_axis),
        {
            "kind": "guide",
            "name": "guide_G02_v3",
            "owner": "S",
            "body": "scapula",
            "world_mm": g.tolist(),
            "local_mm": g_local.tolist(),
            "radius_mm": 8.0,
            "axis_design": axis.tolist(),
            "side_design": side.tolist(),
            "wrap_angle_deg": 108.41975568589456,
        },
        anchor("H", "humerus", b, side),
    ]
    route["status"] = "v3_optimized"
    return route


def apply_g02(config: dict) -> dict:
    route = next(item for item in config["routes"] if item["id"] == "G02")
    a = np.array([149.63981799, 228.48691091, 26.16143230])
    b = np.array([145.47739806, 346.93115200, 13.82575816])
    g = np.array([145.38624668, 268.94507977, -60.24791529])
    axis = normalize(
        np.array([0.96732426, -0.05155616, 0.24824532])
    )
    side = normalize(
        np.array([0.24456068, -0.06856683, -0.96720663])
    )
    _, _, w_axis = scapula_basis()
    route["path"] = [
        anchor("S", "scapula", a, w_axis),
        {
            "kind": "guide",
            "name": "guide_G02_v3",
            "owner": "S",
            "body": "scapula",
            "world_mm": g.tolist(),
            "local_mm": design_to_mj(g - POINTS_MM["AC"]).tolist(),
            "radius_mm": 8.0,
            "axis_design": axis.tolist(),
            "side_design": side.tolist(),
            "wrap_angle_deg": 108.41975568589456,
        },
        anchor("H", "humerus", b, side),
    ]
    route["status"] = "v3_route_and_contraction_repaired"
    return route


def apply_s03(config: dict) -> dict:
    route = next(item for item in config["routes"] if item["id"] == "S03")
    _, _, w_axis = scapula_basis()
    a, normal = surface_point(config, 55.392766171423666, 267.85462023479374, 8.0)
    b = scapula_point(-97.3761939103724, 3.512533599830494, 6.125560083627154)
    route["path"] = [
        anchor("T", "thorax", a, normal),
        anchor("S", "scapula", b, w_axis),
    ]
    route["status"] = "v3_clearance_optimized"
    return route


def apply_s05_s06(config: dict) -> None:
    u_axis, v_axis, w_axis = scapula_basis()
    base = {
        "S05": (np.array([160.0, 220.0, -30.0]), 4.0, 62.0, 10.0),
        "S06": (np.array([180.0, 220.0, 70.0]), -4.0, 58.0, 10.0),
    }
    values = [
        7.30198049,
        -5.00137951,
        10.81880644,
        2.84061719,
        16.66614982,
        -5.19004654,
        6.54307522,
        -5.65618608,
        9.79963212,
        10.63126389,
    ]
    for index, route_id in enumerate(("S05", "S06")):
        route = next(item for item in config["routes"] if item["id"] == route_id)
        a0, u0, v0, w0 = base[route_id]
        offset = 5 * index
        a = a0 + np.array(values[offset : offset + 3])
        u = u0 + values[offset + 3]
        v = v0 + values[offset + 4]
        b = POINTS_MM["P"] + u * u_axis + v * v_axis + w0 * w_axis
        route["path"] = [
            anchor("X", "world", a, np.array([0.0, 1.0, 0.0])),
            anchor("S", "scapula", b, w_axis),
        ]
        route["status"] = "v3_scap_int_fallback"


def main() -> None:
    config = load_json(CONFIG_PATH)
    actions = load_json(ACTION_MATRIX_PATH)
    apply_g02(config)
    apply_s03(config)
    apply_s05_s06(config)

    k7 = next(item for item in actions["actions"] if item["id"] == "K7")
    k7["mandatory_cables"] = ["S01", "S02", "S03", "S04", "G02"]
    scap_int = next(
        item for item in actions["actions"] if item["id"] == "SCAP_INT"
    )
    scap_int["intended_q_target_deg"] = -20.0
    scap_int["q_target_deg"] = {"ac_ext": -5.0}
    scap_int["degraded"] = True
    scap_int["degradation_reason"] = (
        "At -20 deg the S05/S06 scapula anchors are inside the thorax; "
        "the -5 deg minimum-visible fallback is evaluated separately."
    )

    V3_CONFIG.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    V3_MATRIX.write_text(json.dumps(actions, indent=2) + "\n", encoding="utf-8")

    tree, _ = build_model(config)
    tree.write(V3_MODEL, encoding="utf-8", xml_declaration=True)
    evaluator = V2Evaluator(V3_MODEL, config, actions)
    evaluation = evaluator.evaluate_all()
    by_action = {action["id"]: action for action in evaluation["actions"]}
    k7 = by_action["K7"]
    k7_g02 = next(
        route for route in k7["final"]["routes"] if route["id"] == "G02"
    )
    g02_wrap_checks = []
    for action_id in ("K2", "K7"):
        for sample in by_action[action_id]["samples"]:
            route = next(
                item for item in sample["routes"] if item["id"] == "G02"
            )
            g02_wrap_checks.append(
                {
                    "action": action_id,
                    "fraction": sample["fraction"],
                    "wrap_points": route["actual_wrap_points"],
                    "clearance_mm": route["thorax_clearance_mm"],
                    "contraction_percent": route["contraction_percent"],
                }
            )
    full_scap_int_action = copy.deepcopy(
        next(item for item in actions["actions"] if item["id"] == "SCAP_INT")
    )
    full_scap_int_action["q_target_deg"] = {"ac_ext": -20.0}
    full_scap_int = evaluator.evaluate_action(full_scap_int_action, 1.0)
    full_scap_routes = [
        route
        for route in full_scap_int["routes"]
        if route["id"] in ("S05", "S06")
    ]
    V3_EXPORT.mkdir(parents=True, exist_ok=True)
    (V3_EXPORT / "evaluation_v3.json").write_text(
        json.dumps(evaluation, indent=2) + "\n", encoding="utf-8"
    )
    matrix = export_action_matrix(V3_MATRIX, evaluation, actions)
    export_action_markdown(
        PROJECT_ROOT / "design" / "ACTION_CABLE_MATRIX_V3.md", matrix
    )
    export_markdown(
        PROJECT_ROOT / "design" / "CABLE_ANCHORS_V3.md", evaluation, config
    )
    export_routes_csv(V3_EXPORT / "cable_routes_v3.csv", config)
    export_solidworks_csv(
        V3_EXPORT / "SolidWorks_hardpoints_v3.csv", config
    )
    export_tension_report(V3_EXPORT / "tension_report_v3.json", evaluation)
    export_route_validation(
        V3_EXPORT / "route_validation_v3.json", config, V3_MODEL
    )
    repair_report = {
        "changes": {
            "G02": {
                "status": (
                    "pass"
                    if (
                        k7["pass_tension"]
                        and k7["pass_contraction"]
                        and k7["pass_cable_clearance"]
                    )
                    else "route/contraction repaired but K7 closure incomplete"
                ),
                "k7_mandatory": True,
                "k7_contraction_percent": k7_g02["contraction_percent"],
                "k7_clearance_mm": k7_g02["thorax_clearance_mm"],
                "k7_pass_tension": k7["pass_tension"],
                "k7_pass_positive_spanning": k7["final"][
                    "tension_solution"
                ].get("positive_spanning", {}).get("pass"),
                "wrap_checks": g02_wrap_checks,
            },
            "S03": {
                "status": (
                    "pass"
                    if by_action["SLIDE_BACK"]["pass_cable_clearance"]
                    else "fail"
                ),
                "slide_back_clearance_mm": by_action["SLIDE_BACK"][
                    "min_selected_cable_clearance_mm"
                ],
            },
            "S05/S06": {
                "status": "minimum-visible fallback",
                "intended_scap_int_deg": -20.0,
                "evaluated_scap_int_deg": -5.0,
                "full_target_min_clearance_mm": min(
                    route["thorax_clearance_mm"] for route in full_scap_routes
                ),
                "fallback_min_clearance_mm": by_action["SCAP_INT"][
                    "min_selected_cable_clearance_mm"
                ],
                "reason": (
                    "At -20 deg both scapula-side anchors are inside the "
                    "thorax envelope, so no guide can recover 2 mm clearance."
                ),
            },
        },
        "evaluation": str(V3_EXPORT / "evaluation_v3.json"),
    }
    (V3_EXPORT / "route_repair_report_v3.json").write_text(
        json.dumps(repair_report, indent=2) + "\n", encoding="utf-8"
    )
    for action in evaluation["actions"]:
        print(
            f"{action['id']:<10} selected={','.join(action['selected_cables']):<28} "
            f"80={action['pass_80_percent']} min={action['pass_minimum_visible']} "
            f"tension={action['pass_tension']} contraction={action['pass_contraction']} "
            f"clearance={action['pass_cable_clearance']}"
        )
    print(f"wrote {V3_CONFIG}")
    print(f"wrote {V3_MODEL}")
    print(f"wrote {V3_EXPORT}")


if __name__ == "__main__":
    main()
