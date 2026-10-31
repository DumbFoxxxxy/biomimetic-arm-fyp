"""Group V2 action results into the agreed staged optimization order."""

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

from v2_common import EXPORT_ROOT


STAGES = {
    "A_scapula_up": ["K6"],
    "B_scapula_aux": ["SLIDE_FWD", "SLIDE_BACK", "SCAP_EXT", "SCAP_INT"],
    "C_gh": ["K2", "K3", "K4", "K5"],
    "D_elbow": ["K1"],
    "E_combined": ["K7"],
    "F_pectoral": ["PECTORAL"],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--evaluation",
        type=Path,
        default=EXPORT_ROOT / "evaluation_v2.json",
    )
    parser.add_argument(
        "--out", type=Path, default=EXPORT_ROOT / "staged_report.json"
    )
    args = parser.parse_args()
    evaluation = json.loads(args.evaluation.read_text(encoding="utf-8"))
    by_id = {action["id"]: action for action in evaluation["actions"]}
    stages = {}
    for name, action_ids in STAGES.items():
        actions = [by_id[action_id] for action_id in action_ids]
        stages[name] = {
            "actions": action_ids,
            "pass_minimum_visible": all(
                action["pass_minimum_visible"] for action in actions
            ),
            "pass_80_percent": all(
                action["pass_80_percent"] for action in actions
            ),
            "pass_tension": all(
                action["pass_tension"] for action in actions
            ),
            "pass_contraction": all(
                action["pass_contraction"] for action in actions
            ),
            "pass_cable_clearance": all(
                action["pass_cable_clearance"] for action in actions
            ),
            "pass_guide_geometry": all(
                action["pass_guide_geometry"] for action in actions
            ),
            "pass_visual_clearance": all(
                action["pass_visual_clearance"] for action in actions
            ),
            "pass_minimum_functional": all(
                action["pass_minimum_functional"] for action in actions
            ),
            "selected_cables": {
                action["id"]: action["selected_cables"] for action in actions
            },
        }
    args.out.write_text(
        json.dumps({"stages": stages}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.out}")
    for name, stage in stages.items():
        print(
            f"{name:<18} min={stage['pass_minimum_visible']} "
            f"80={stage['pass_80_percent']} tension={stage['pass_tension']} "
            f"functional={stage['pass_minimum_functional']}"
        )


if __name__ == "__main__":
    main()
