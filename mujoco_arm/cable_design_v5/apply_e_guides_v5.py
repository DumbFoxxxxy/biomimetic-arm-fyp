"""Convert E01/E02 pass-through points into antagonistic elbow guides."""

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


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "cable_design_v5.json"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    for route_id, side in (("E01", [0.0, 0.0, -1.0]), ("E02", [0.0, 0.0, -1.0])):
        route = next(item for item in config["routes"] if item["id"] == route_id)
        path = route["path"]
        pass_index = next(
            index
            for index, item in enumerate(path)
            if item["kind"] in ("pass_through", "guide")
        )
        point = path[pass_index]
        path[pass_index] = {
            "kind": "guide",
            "name": f"guide_{route_id}_v5",
            "owner": "H",
            "body": "humerus",
            "world_mm": point["world_mm"],
            "local_mm": point["local_mm"],
            "radius_mm": 8.0,
            "axis_design": [0.9962, 0.0872, 0.0],
            "side_design": side,
            "wrap_angle_deg": 60.0,
        }
        route["max_guides"] = 2
        route["status"] = "v5_true_elbow_guide"
    args.config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"updated {args.config}")


if __name__ == "__main__":
    main()
