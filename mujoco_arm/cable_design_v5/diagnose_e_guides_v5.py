"""Diagnose which radial sidesite keeps E01/E02 antagonistic."""

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

import numpy as np


ROOT = Path(__file__).resolve().parent
EXPORT = ROOT.parent / "exports" / "cable_design_v5"
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from build_v2_model import build  # noqa: E402
from evaluate_v2 import dense_tendon_jacobian  # noqa: E402
from v2_common import load_json  # noqa: E402
from v5_evaluator import V5Evaluator  # noqa: E402


def main() -> None:
    config = load_json(ROOT / "cable_design_v5.json")
    matrix = load_json(ROOT / "action_cable_matrix_v5.json")
    action = next(item for item in matrix["actions"] if item["id"] == "K1")
    sides = {
        "+x": [1.0, 0.0, 0.0],
        "-x": [-1.0, 0.0, 0.0],
        "+y": [0.0, 1.0, 0.0],
        "-y": [0.0, -1.0, 0.0],
        "+z": [0.0, 0.0, 1.0],
        "-z": [0.0, 0.0, -1.0],
    }
    output = {}
    EXPORT.mkdir(parents=True, exist_ok=True)
    temp_model = EXPORT / "e_side_diagnostic.xml"
    for e01_side in sides:
        for e02_side in sides:
            trial = copy.deepcopy(config)
            for route_id, side in (("E01", e01_side), ("E02", e02_side)):
                route = next(
                    item for item in trial["routes"] if item["id"] == route_id
                )
                guide = next(
                    item for item in route["path"] if item["kind"] == "guide"
                )
                guide["side_design"] = sides[side]
            tree, _ = build(trial)
            tree.write(temp_model, encoding="utf-8", xml_declaration=True)
            evaluator = V5Evaluator(temp_model, trial, matrix)
            rows = []
            for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
                data = evaluator._action_pose(action, fraction)
                dofs = evaluator._dofs_for_action(action)
                columns = [
                    evaluator.tendon_indices[route_id]
                    for route_id in ("E01", "E02")
                ]
                jacobian = dense_tendon_jacobian(evaluator.model, data)
                moment = -jacobian[:, dofs][columns].T
                rows.append(
                    {
                        "fraction": fraction,
                        "moments": moment.ravel().tolist(),
                        "wrap": [
                            evaluator._route_geometry(
                                evaluator.route_by_id[route_id], data
                            )["actual_wrap_points"]
                            for route_id in ("E01", "E02")
                        ],
                    }
                )
            output[f"{e01_side}/{e02_side}"] = rows
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
