"""Sweep V6 channel counts in the fixed order S07, G08, S08."""

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


ROOT = Path(__file__).resolve().parent
V5_ROOT = ROOT.parent / "cable_design_v5"
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v6"
sys.path[:0] = [str(ROOT), str(V5_ROOT), str(V2_ROOT)]

from build_v2_model import build as build_base  # noqa: E402
from build_v6 import add_v6_markers_and_roller  # noqa: E402
from trajectory_v6 import TrajectoryRunner, load_inputs  # noqa: E402
from v6_evaluator import V6Evaluator  # noqa: E402


SOURCE_ROUTE = {
    "S07": "S06",
    "G08": "G07",
    "S08": "S04",
}


def clone_route(config: dict, add_id: str) -> dict:
    source_id = SOURCE_ROUTE[add_id]
    route = copy.deepcopy(
        next(item for item in config["routes"] if item["id"] == source_id)
    )
    route["id"] = add_id
    route["role"] = config["optional_channels"][add_id]["role"]
    route["status"] = "v6_channel_sweep_candidate"
    for element in route["path"]:
        if element["kind"] == "guide":
            element["name"] = f"guide_{add_id}_v6"
    return route


def sync_lists(config: dict, matrix: dict, library: dict) -> None:
    channels = [route["id"] for route in config["routes"]]
    for action in matrix["actions"]:
        action["allowed_cables"] = list(channels)
        action["preferred_cables"] = [
            route_id
            for route_id in action.get("preferred_cables", [])
            if route_id in channels
        ]
        action["mandatory_cables"] = [
            route_id
            for route_id in action.get("mandatory_cables", [])
            if route_id in channels
        ]
    for trajectory in library["trajectories"]:
        trajectory["active_cables"] = [
            route_id
            for route_id in trajectory["active_cables"]
            if route_id in channels
        ]


def evaluate_variant(
    name: str,
    config: dict,
    matrix: dict,
    library: dict,
) -> dict:
    tree, _ = build_base(config)
    add_v6_markers_and_roller(tree, config)
    xml = ET.tostring(tree.getroot(), encoding="unicode")
    model = _load_project_xml(xml)
    evaluator = V6Evaluator(
        ROOT / f"{name}.xml",
        config,
        matrix,
        model=model,
    )
    static = evaluator.evaluate_all()
    runner = TrajectoryRunner(
        ROOT / f"{name}.xml",
        config,
        matrix,
        library,
        model=model,
    )
    trajectories = []
    for trajectory in library["trajectories"]:
        report = runner.run(trajectory["id"], 1.0, False)
        report["pass_full"] = bool(
            report["pass_path"]
            and report["pass_clearance"]
            and report["all_tensions_ok"]
        )
        trajectories.append(report)
    static_pass = all(
        action["pass_minimum_functional"]
        for action in static["actions"]
    )
    full_pass = static_pass and all(
        trajectory["pass_full"] for trajectory in trajectories
    )
    return {
        "name": name,
        "channel_count": len(config["routes"]),
        "static_pass": static_pass,
        "trajectories_pass": all(
            item["pass_full"] for item in trajectories
        ),
        "full_pass": full_pass,
        "trajectory_summary": {
            item["id"]: {
                "pass_full": item["pass_full"],
                "rms_error_mm": item["rms_error_mm"],
                "p95_error_mm": item["p95_error_mm"],
                "min_visual_clearance_mm": item[
                    "min_visual_clearance_mm"
                ],
            }
            for item in trajectories
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=EXPORT)
    args = parser.parse_args()
    base_config, base_matrix, base_library = load_inputs()
    variants = []
    config = copy.deepcopy(base_config)
    matrix = copy.deepcopy(base_matrix)
    library = copy.deepcopy(base_library)
    sync_lists(config, matrix, library)
    variants.append(
        evaluate_variant("channels17", config, matrix, library)
    )
    for add_id in ("S07", "G08", "S08"):
        config["routes"].append(clone_route(config, add_id))
        for action in matrix["actions"]:
            if action["group"] in ("scapula", "combined", "gh"):
                action.setdefault("preferred_cables", []).append(add_id)
        for trajectory in library["trajectories"]:
            if trajectory["id"] != "ALL_ACTIONS":
                trajectory["active_cables"].append(add_id)
        sync_lists(config, matrix, library)
        name = f"channels{len(config['routes'])}"
        variants.append(
            evaluate_variant(name, config, matrix, library)
        )
    passing = [
        item for item in variants if item["full_pass"]
    ]
    selected = (
        min(passing, key=lambda item: item["channel_count"])
        if passing
        else variants[0]
    )
    output = {
        "candidate_order": ["S07", "G08", "S08"],
        "variants": variants,
        "selected": selected["name"],
        "selected_channel_count": selected["channel_count"],
        "selected_full_pass": selected["full_pass"],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "channel_sweep_v6.json").write_text(
        json.dumps(output, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
