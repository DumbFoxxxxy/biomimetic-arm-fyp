"""Create the independent V5 configuration from the V4 baseline."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
V4_ROOT = ROOT.parent / "cable_design_v4"
V5_CONFIG = ROOT / "cable_design_v5.json"
V5_MATRIX = ROOT / "action_cable_matrix_v5.json"

CLEAN_ACTION_KEYS = (
    "selected_cables",
    "cooperating_cable_count",
    "target_tensions_N",
    "pass_minimum_visible",
    "pass_80_percent",
    "pass_tension",
    "positive_spanning",
    "max_selected_contraction_percent",
    "pass_contraction",
    "min_selected_cable_clearance_mm",
    "pass_cable_clearance",
    "pass_guide_geometry",
    "pass_visual_clearance",
    "pass_minimum_functional",
)


def configure_actions(matrix: dict, all_cables: list[str]) -> None:
    design = {
        "K1": (145.0, 116.0),
        "K2": (120.0, 96.0),
        "K3": (120.0, 96.0),
        "K4": (70.0, 56.0),
        "K5": (90.0, 72.0),
        "K6": (30.0, 24.0),
        "K7": (126.0, 120.0),
        "SCAP_EXT": (20.0, 16.0),
        "SCAP_INT": (5.0, 5.0),
        "PECTORAL": (20.0, 16.0),
    }
    for action in matrix["actions"]:
        for key in CLEAN_ACTION_KEYS:
            action.pop(key, None)
        action_id = action["id"]
        action["validation_pose"] = True
        if action_id in design:
            target, acceptance = design[action_id]
            action["target_deg"] = target
            action["acceptance_deg"] = acceptance
        if action_id == "K7":
            action["legacy_target_deg"] = 150.0
        if action_id == "SCAP_INT":
            action["legacy_target_deg"] = 20.0
            action["limitation"] = (
                "Full -20 degree internal rotation is infeasible with the "
                "current thorax/scapula geometry; -5 degrees is the formal V5 target."
            )
        if action_id in ("SLIDE_FWD", "SLIDE_BACK"):
            action["target_mm"] = 15.0
            action["acceptance_mm"] = 12.0
            action["minimum_visible_mm"] = 5.0
        action["allowed_cables"] = list(all_cables)


def main() -> None:
    config = json.loads(
        (V4_ROOT / "cable_design_v4.json").read_text(encoding="utf-8")
    )
    matrix = json.loads(
        (V4_ROOT / "action_cable_matrix_v4.json").read_text(encoding="utf-8")
    )
    config = copy.deepcopy(config)
    matrix = copy.deepcopy(matrix)
    config["model_name"] = "arm_cable_design_v5"
    config["schema_version"] = 5
    config["zero_pose_offsets_deg"] = {
        "sc_elevation": 0.0,
        "sc_protraction": 0.0,
        "ac_ext": 0.0,
    }
    config["zero_pose_limits_deg"] = 10.0
    config["contact_pad_offset_mm"] = 0.0
    config["contact_pad_offset_limit_mm"] = 10.0
    config["max_channels"] = 18
    config["channel_plan"] = {
        "preferred_count": 17,
        "minimum_count": 16,
        "maximum_count": 18,
        "optional_eighteenth": ["G08", "S07"],
    }
    config["friction_model"] = {
        "type": "per_guide_capstan",
        "mu": 0.05,
        "sensitivity": [0.0, 0.05, 0.1],
        "static_direction": "worst_case",
        "dahl_dynamic": "not_applicable_to_quasistatic",
    }
    config["positive_spanning"] = {
        "mode": "hard",
        "minimum_columns": "dofs_plus_one",
    }
    config["tension_solver"] = {
        "minimum_mandatory_tension_N": 2.0,
        "maximum_optional_subset": 3,
        "exact_equilibrium_relative_tolerance": 0.05,
        "nullspace_objective": "maximize_minimum_tension_margin",
    }
    config["acceptance"]["action_pose_fractions"] = [
        0.0,
        0.25,
        0.5,
        0.75,
        1.0,
    ]
    config["acceptance"]["final_pose_fractions"] = [
        round(index / 20.0, 2) for index in range(21)
    ]
    config["acceptance"]["full_mesh_clearance_mm_min"] = -0.5
    config["acceptance"]["contact_penetration_mm_min"] = -0.5
    config["acceptance"]["positive_spanning_required"] = True
    config["acceptance"]["minimum_channel_count"] = 16
    config["acceptance"]["maximum_channel_count"] = 18
    config["path_corrections"] = {
        action_id: {
            fraction: [0.0, 0.0, 0.0]
            for fraction in (0.0, 0.25, 0.5, 0.75, 1.0)
        }
        for action_id in ("K6", "K7", "SLIDE_BACK", "SCAP_INT")
    }
    config["path_correction_limits_deg"] = 10.0
    config["route_search"] = {
        "seed": 71,
        "G": {
            "tangent_mm": 20.0,
            "outward_mm": [0.0, 10.0],
            "guides": [0, 1, 2],
        },
        "S": {
            "tangent_mm": 15.0,
            "outward_mm": [0.0, 10.0],
            "guides": [0, 1, 2],
        },
        "P": {
            "tangent_mm": 15.0,
            "outward_mm": [0.0, 10.0],
            "guides": [0, 1],
        },
    }
    config["optional_channels"] = {
        "G08": {
            "group": "gh",
            "role": "v5_gh_auxiliary",
            "allowed_actions": ["K2", "K3", "K7"],
            "max_guides": 2,
        },
        "S07": {
            "group": "scapula",
            "role": "v5_scapula_stabilizer",
            "allowed_actions": ["K6", "K7", "SLIDE_BACK", "SCAP_INT"],
            "max_guides": 2,
        },
    }
    configure_actions(matrix, [route["id"] for route in config["routes"]])
    (ROOT / "cable_design_v5.json").write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8"
    )
    (ROOT / "action_cable_matrix_v5.json").write_text(
        json.dumps(matrix, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {ROOT / 'cable_design_v5.json'}")
    print(f"wrote {ROOT / 'action_cable_matrix_v5.json'}")


if __name__ == "__main__":
    main()
