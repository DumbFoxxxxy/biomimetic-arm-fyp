"""Create the V4 configuration with the dedicated K7 cable G07."""

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


ROOT = Path(__file__).resolve().parent
V3_ROOT = ROOT.parent / "cable_design_v3"
V4_CONFIG = ROOT / "cable_design_v4.json"
V4_MATRIX = ROOT / "action_cable_matrix_v4.json"


def main() -> None:
    config = json.loads(
        (V3_ROOT / "cable_design_v3.json").read_text(encoding="utf-8")
    )
    actions = json.loads(
        (V3_ROOT / "action_cable_matrix_v3.json").read_text(encoding="utf-8")
    )
    if any(route["id"] == "G07" for route in config["routes"]):
        raise RuntimeError("V4 config already contains G07")
    g01 = next(route for route in config["routes"] if route["id"] == "G01")
    g07 = copy.deepcopy(g01)
    g07["id"] = "G07"
    g07["role"] = "k7_gh_auxiliary"
    g07["status"] = "v4_search_initial"
    for element in g07["path"]:
        if element["kind"] == "guide":
            element["name"] = "guide_G07_v4"
    config["routes"].append(g07)

    k7 = next(action for action in actions["actions"] if action["id"] == "K7")
    if "G07" not in k7["allowed_cables"]:
        k7["allowed_cables"].append("G07")
    if "G07" not in k7["mandatory_cables"]:
        k7["mandatory_cables"].append("G07")

    V4_CONFIG.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    V4_MATRIX.write_text(json.dumps(actions, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {V4_CONFIG}")
    print(f"wrote {V4_MATRIX}")


if __name__ == "__main__":
    main()
