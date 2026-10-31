"""Diagnose real E01/E02 guide wrap and moment signs for V6."""

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
from project_runtime import WorkspacePath as Path
import sys
import xml.etree.ElementTree as ET

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v6"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from build_v2_model import build  # noqa: E402
from evaluate_v2 import V2Evaluator, dense_tendon_jacobian  # noqa: E402
from v2_common import load_json  # noqa: E402


CONFIG = ROOT / "cable_design_v6.json"
MATRIX = ROOT / "action_cable_matrix_v6.json"
SIDES = {
    "+x": [1.0, 0.0, 0.0],
    "-x": [-1.0, 0.0, 0.0],
    "+y": [0.0, 1.0, 0.0],
    "-y": [0.0, -1.0, 0.0],
    "+z": [0.0, 0.0, 1.0],
    "-z": [0.0, 0.0, -1.0],
}


def evaluate_variant(config: dict, matrix: dict, action: dict) -> dict:
    tree, _ = build(config)
    model = _load_project_xml(
        ET.tostring(tree.getroot(), encoding="unicode")
    )
    evaluator = V2Evaluator(
        ROOT / "in_memory_v6.xml",
        config,
        matrix,
        model=model,
    )
    samples = []
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        data = evaluator._action_pose(action, fraction)
        dofs = evaluator._dofs_for_action(action)
        columns = [
            evaluator.tendon_indices[route_id] for route_id in ("E01", "E02")
        ]
        jacobian = dense_tendon_jacobian(evaluator.model, data)
        moments = -jacobian[:, dofs][columns].T
        wrap = [
            evaluator._route_geometry(evaluator.route_by_id[route_id], data)[
                "actual_wrap_points"
            ]
            for route_id in ("E01", "E02")
        ]
        samples.append(
            {
                "fraction": fraction,
                "moments_m_per_n": moments.ravel().tolist(),
                "wrap_points": wrap,
            }
        )
    opposite = all(
        sample["moments_m_per_n"][0] * sample["moments_m_per_n"][1] < 0.0
        for sample in samples
    )
    arm_ok = all(
        min(abs(value) for value in sample["moments_m_per_n"]) >= 0.008
        for sample in samples
    )
    wrap_ok = all(
        all(count >= 4 for count in sample["wrap_points"])
        for sample in samples
    )
    return {
        "opposite": opposite,
        "arm_ok": arm_ok,
        "wrap_ok": wrap_ok,
        "pass": opposite and arm_ok and wrap_ok,
        "samples": samples,
    }


def main() -> None:
    config = load_json(CONFIG)
    matrix = load_json(MATRIX)
    action = next(item for item in matrix["actions"] if item["id"] == "K1")
    results = {}
    for e01_side in SIDES:
        for e02_side in SIDES:
            trial = copy.deepcopy(config)
            for route_id, side in (("E01", e01_side), ("E02", e02_side)):
                route = next(
                    item for item in trial["routes"] if item["id"] == route_id
                )
                guide = next(
                    item for item in route["path"] if item["kind"] == "guide"
                )
                guide["side_design"] = SIDES[side]
            name = f"{e01_side}/{e02_side}"
            try:
                results[name] = evaluate_variant(trial, matrix, action)
            except Exception as error:  # noqa: BLE001
                results[name] = {"pass": False, "error": str(error)}
    best = sorted(
        results.items(),
        key=lambda item: (
            not item[1].get("pass", False),
            not item[1].get("opposite", False),
            not item[1].get("arm_ok", False),
            not item[1].get("wrap_ok", False),
        ),
    )
    output = {
        "best": best[0][0] if best else None,
        "passing": [name for name, value in results.items() if value.get("pass")],
        "variants": results,
    }
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "elbow_side_diagnostic_v6.json").write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
