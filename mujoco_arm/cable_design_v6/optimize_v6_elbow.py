"""Search a real R10 antagonist route with opposing elbow moments."""

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
from project_runtime import WorkspacePath as Path
import sys
import xml.etree.ElementTree as ET

import mujoco
import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v6"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from build_v2_model import build  # noqa: E402
from evaluate_v2 import V2Evaluator, dense_tendon_jacobian  # noqa: E402
from v2_common import (  # noqa: E402
    POINTS_MM,
    design_to_mj,
    load_json,
    normalize,
)


CONFIG = ROOT / "cable_design_v6.json"
MATRIX = ROOT / "action_cable_matrix_v6.json"
FRACTIONS = np.linspace(0.0, 1.0, 21)
ELBOW_AXIS = normalize(np.array([0.9962, 0.0872, 0.0]))
SIDE_01 = normalize(np.array([-1.0, 0.0, 0.0]))
SIDE_02 = normalize(np.array([-1.0, 0.0, 0.0]))


def route_from_candidate(values: np.ndarray, mirrored: bool) -> dict:
    values = np.asarray(values, dtype=float)
    a = values[0:3]
    g = values[3:6]
    b = values[6:9]
    centre_z = float(POINTS_MM["E"][2])
    if mirrored:
        a = a.copy()
        g = g.copy()
        b = b.copy()
        a[2] = 2.0 * centre_z - a[2]
        g[2] = 2.0 * centre_z - g[2]
        b[2] = 2.0 * centre_z - b[2]
    side = SIDE_01 if not mirrored else SIDE_02
    route_id = "E02" if mirrored else "E01"
    origin_f = POINTS_MM["E"]
    local_g = np.array(
        [g[0] - origin_f[0], g[1] - origin_f[1], g[2] - origin_f[2]]
    )
    local_b = np.array(
        [b[0] - origin_f[0], b[1] - origin_f[1], b[2] - origin_f[2]]
    )
    return {
        "id": route_id,
        "group": "elbow",
        "role": "elbow_extension" if mirrored else "elbow_flexion",
        "diameter_mm": 1.0,
        "max_tension_N": 150.0,
        "max_guides": 1,
        "path": [
            {
                "kind": "anchor",
                "owner": "X",
                "body": "world",
                "world_mm": a.tolist(),
                "local_mm": design_to_mj(a).tolist(),
                "installation_normal": [0.0, 1.0, 0.0],
            },
            {
                "kind": "guide",
                "name": f"guide_{route_id}_v6",
                "owner": "F",
                "body": "ulna",
                "world_mm": g.tolist(),
                "local_mm": [
                    local_g[0] * 0.001,
                    -local_g[2] * 0.001,
                    local_g[1] * 0.001,
                ],
                "radius_mm": 10.0,
                "axis_design": ELBOW_AXIS.tolist(),
                "side_design": side.tolist(),
                "wrap_angle_deg": 90.0,
                "owner_part": "ulna",
            },
            {
                "kind": "anchor",
                "owner": "F",
                "body": "ulna",
                "world_mm": b.tolist(),
                "local_mm": [
                    local_b[0] * 0.001,
                    -local_b[2] * 0.001,
                    local_b[1] * 0.001,
                ],
                "installation_normal": side.tolist(),
            },
        ],
        "tpu": {
            "enabled": True,
            "reference_strain": 0.1,
            "reference_force_N": 20.0,
            "compression_force_N": 0.0,
        },
        "status": "v6_antagonist_search",
        "notes": [
            "Two external anchors around an R10 pulley fixed to the ulna."
        ],
    }


def install_candidate(config: dict, values: np.ndarray) -> None:
    routes = [
        route
        for route in config["routes"]
        if route["id"] not in ("E01", "E02")
    ]
    config["routes"] = [
        route_from_candidate(values, False),
        route_from_candidate(values, True),
        *routes,
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=35)
    parser.add_argument("--popsize", type=int, default=8)
    parser.add_argument("--seed", type=int, default=106)
    args = parser.parse_args()

    base_config = load_json(CONFIG)
    matrix = load_json(MATRIX)
    action = next(item for item in matrix["actions"] if item["id"] == "K1")

    initial_values = np.array(
        [
            184.34,
            285.0,
            -28.0,
            194.34,
            -81.35,
            -23.0,
            203.27,
            -191.35,
            -40.0,
        ],
        dtype=float,
    )
    install_candidate(base_config, initial_values)
    tree, _ = build(base_config)
    model = _load_project_xml(
        ET.tostring(tree.getroot(), encoding="unicode")
    )
    evaluator = V2Evaluator(
        ROOT / "in_memory_v6_elbow.xml",
        base_config,
        matrix,
        model=model,
    )
    names = {
        "E01A": "V2_A_E01",
        "E01B": "V2_B_E01",
        "E01G": "V2_guide_E01_guide_E01_v6",
        "E01S": "V2_side_E01",
        "E02A": "V2_A_E02",
        "E02B": "V2_B_E02",
        "E02G": "V2_guide_E02_guide_E02_v6",
        "E02S": "V2_side_E02",
    }
    ids = {
        key: int(
            mujoco.mj_name2id(
                model,
                (
                    mujoco.mjtObj.mjOBJ_GEOM
                    if key.endswith("G")
                    else mujoco.mjtObj.mjOBJ_SITE
                ),
                name,
            )
        )
        for key, name in names.items()
    }
    history: list[dict] = []

    def set_positions(values: np.ndarray) -> None:
        states = {"E01": route_from_candidate(values, False), "E02": route_from_candidate(values, True)}
        for route_id, state in states.items():
            path = state["path"]
            model.site_pos[ids[f"{route_id}A"]] = np.asarray(path[0]["local_mm"])
            model.site_pos[ids[f"{route_id}B"]] = np.asarray(path[2]["local_mm"])
            guide_local = np.asarray(path[1]["local_mm"])
            model.geom_pos[ids[f"{route_id}G"]] = guide_local
            side = normalize(np.asarray(path[1]["side_design"]))
            side_mj = normalize(
                np.array([side[0], -side[2], side[1]], dtype=float)
            )
            model.site_pos[ids[f"{route_id}S"]] = (
                guide_local + side_mj * 0.010 * 1.5
            )

    def objective(values: np.ndarray) -> float:
        set_positions(values)
        penalty = 0.0
        signs: list[float] = []
        wraps: list[int] = []
        arms: list[float] = []
        states = states_for_values(values)
        try:
            for fraction in FRACTIONS:
                data = evaluator._action_pose(action, float(fraction))
                dofs = evaluator._dofs_for_action(action)
                columns = [
                    evaluator.tendon_indices[route_id]
                    for route_id in ("E01", "E02")
                ]
                jacobian = dense_tendon_jacobian(evaluator.model, data)
                moments = -jacobian[:, dofs][columns].T
                values_now = moments.ravel()
                product = float(values_now[0] * values_now[1])
                signs.append(product)
                arms.extend([abs(float(values_now[0])), abs(float(values_now[1]))])
                if product >= 0.0:
                    penalty += 100.0 + min(20.0, abs(product) * 1000.0)
                for route in states:
                    wrap_count = int(
                        data.ten_wrapnum[evaluator.tendon_ids[route["id"]]]
                    )
                    wraps.append(wrap_count)
                    if wrap_count < 4:
                        penalty += 20.0
            if arms:
                penalty += sum(max(0.0, 0.008 - arm) * 5000.0 for arm in arms)
        except Exception:  # noqa: BLE001
            return 1e9
        score = penalty
        history.append(
            {
                "values": values.tolist(),
                "score": score,
                "moment_products": signs,
                "min_abs_arm_m": min(arms) if arms else None,
                "min_wrap_points": min(wraps) if wraps else None,
            }
        )
        return score

    bounds = [
        (165.0, 205.0),
        (220.0, 330.0),
        (-50.0, -15.0),
        (175.0, 205.0),
        (-150.0, -60.0),
        (-60.0, -15.0),
        (175.0, 205.0),
        (-230.0, -140.0),
        (-80.0, -20.0),
    ]
    result = differential_evolution(
        objective,
        bounds,
        seed=args.seed,
        maxiter=args.maxiter,
        popsize=args.popsize,
        polish=False,
        workers=1,
        updating="immediate",
    )
    best = min(history, key=lambda item: item["score"])
    install_candidate(base_config, np.asarray(best["values"]))
    CONFIG.write_text(
        json.dumps(base_config, indent=2) + "\n", encoding="utf-8"
    )
    output = {
        "optimizer_success": bool(result.success),
        "best": best,
        "evaluations": len(history),
        "candidate_found": bool(
            best["score"] < 1e-6
            and best["moment_products"]
            and max(best["moment_products"]) < 0.0
            and best["min_abs_arm_m"] >= 0.008
            and best["min_wrap_points"] >= 4
        ),
    }
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "elbow_optimization_v6.json").write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, indent=2))


def states_for_values(values: np.ndarray) -> list[dict]:
    return [
        route_from_candidate(values, False),
        route_from_candidate(values, True),
    ]


if __name__ == "__main__":
    main()
