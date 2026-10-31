"""Create the isolated V6 configuration and trajectory library from V5."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
V5_ROOT = ROOT.parent / "cable_design_v5"
V6_CONFIG = ROOT / "cable_design_v6.json"
V6_MATRIX = ROOT / "action_cable_matrix_v6.json"
V6_TRAJECTORIES = ROOT / "trajectory_library_v6.json"

RESULT_KEYS = (
    "selected_cables",
    "cooperating_cable_count",
    "target_tensions_N",
    "pass_minimum_visible",
    "pass_80_percent",
    "pass_tension",
    "pass_positive_spanning",
    "positive_spanning",
    "max_selected_contraction_percent",
    "pass_contraction",
    "min_selected_cable_clearance_mm",
    "pass_cable_clearance",
    "pass_guide_geometry",
    "pass_visual_clearance",
    "pass_contact",
    "pass_minimum_functional",
)


def world_to_local_mj(
    world_mm: list[float],
    origin_mm: list[float],
) -> list[float]:
    delta = [
        float(world_mm[0]) - float(origin_mm[0]),
        float(world_mm[1]) - float(origin_mm[1]),
        float(world_mm[2]) - float(origin_mm[2]),
    ]
    return [
        delta[0] * 0.001,
        -delta[2] * 0.001,
        delta[1] * 0.001,
    ]


def mirror_z(world_mm: list[float], centre_z: float) -> list[float]:
    return [
        float(world_mm[0]),
        float(world_mm[1]),
        2.0 * centre_z - float(world_mm[2]),
    ]


def initial_elbow_route(route_id: str, mirrored: bool) -> dict:
    """A conservative starting geometry, replaced by the route search."""
    origin_h = [159.06, 187.54, 7.30]
    origin_f = [184.34, -101.35, 7.30]
    elbow_axis = [0.9962, 0.0872, 0.0]
    centre_z = 7.30
    a = [165.50, 35.00, -22.00]
    g = [194.00, -120.00, -22.00]
    b = [194.00, -210.00, 7.30]
    if mirrored:
        a = mirror_z(a, centre_z)
        g = mirror_z(g, centre_z)
        b = mirror_z(b, centre_z)
    side = [0.0, 0.0, -1.0 if not mirrored else 1.0]
    return {
        "id": route_id,
        "group": "elbow",
        "role": "elbow_flexion" if not mirrored else "elbow_extension",
        "diameter_mm": 1.0,
        "max_tension_N": 150.0,
        "max_guides": 1,
        "path": [
            {
                "kind": "anchor",
                "owner": "H",
                "body": "humerus",
                "world_mm": a,
                "local_mm": world_to_local_mj(a, origin_h),
                "installation_normal": side,
            },
            {
                "kind": "guide",
                "name": f"guide_{route_id}_v6",
                "owner": "F",
                "body": "ulna",
                "world_mm": g,
                "local_mm": world_to_local_mj(g, origin_f),
                "radius_mm": 10.0,
                "axis_design": elbow_axis,
                "side_design": side,
                "wrap_angle_deg": 90.0,
                "owner_part": "ulna",
            },
            {
                "kind": "anchor",
                "owner": "F",
                "body": "ulna",
                "world_mm": b,
                "local_mm": world_to_local_mj(b, origin_f),
                "installation_normal": side,
            },
        ],
        "tpu": {
            "enabled": True,
            "reference_strain": 0.1,
            "reference_force_N": 20.0,
            "compression_force_N": 0.0,
        },
        "status": "v6_antagonistic_elbow_initial",
        "notes": [
            "R10 guide on the ulna; route is replaced by optimize_v6_elbow.py."
        ],
    }


def configure_config(config: dict) -> None:
    config["model_name"] = "arm_cable_design_v6"
    config["schema_version"] = 6
    config["max_channels"] = 20
    config["channel_plan"] = {
        "preferred_count": 17,
        "minimum_count": 17,
        "maximum_count": 20,
        "candidate_order": ["S07", "G08", "S08"],
    }
    config["positive_spanning"] = {
        "mode": "diagnostic",
        "used_as_pass_gate": False,
    }
    config["tension_solver"] = {
        "mode": "task_conditioned_nonnegative",
        "positive_spanning_hard_gate": False,
        "exact_equilibrium_relative_tolerance": 0.08,
        "minimum_mandatory_tension_N": 1.0,
        "minimum_active_tension_N": 0.25,
        "nullspace_objective": "minimum_norm",
    }
    config["acceptance"]["positive_spanning_required"] = False
    config["acceptance"]["minimum_channel_count"] = 17
    config["acceptance"]["maximum_channel_count"] = 20
    config["acceptance"]["elbow_effective_moment_arm_mm_min"] = 8.0
    config["acceptance"]["elbow_wrap_deg"] = [30.0, 170.0]
    config["acceptance"]["wrist_path_rms_mm_max"] = 15.0
    config["acceptance"]["wrist_path_p95_mm_max"] = 25.0
    config["acceptance"]["trajectory_sample_period_s"] = 0.02
    config["simulation"]["control_period_s"] = 0.02
    config["simulation"]["actuator_type"] = "position"
    config["simulation"]["position_kp_n_per_rad"] = 3.0
    config["simulation"]["position_kd_nms_per_rad"] = 0.18
    config["simulation"]["position_pretension_mm"] = 0.1
    config["simulation"]["position_jacobian_gain"] = 5.0
    config["simulation"]["position_mode"] = "jacobian"
    config["simulation"]["playback_speed_min"] = 0.25
    config["simulation"]["playback_speed_max"] = 1.0
    config["elbow_antagonism"] = {
        "guide_radius_mm": 10.0,
        "guide_owner": "F",
        "moment_sign_opposite_required": True,
        "effective_moment_arm_mm_min": 8.0,
        "sample_fractions": [0.0, 0.25, 0.5, 0.75, 1.0],
        "search_seed": 106,
    }
    config["shoulder_guide"] = {
        "type": "r400_single_sided_crowned_roller",
        "roller_radius_mm": 8.0,
        "roller_width_mm": 16.0,
        "roller_material": "POM",
        "nominal_preload_N": 5.0,
        "preload_sensitivity_N": [3.0, 5.0, 10.0],
        "equivalent_stiffness_N_per_mm": 2.0,
        "normal_travel_mm": 4.0,
        "hard_slot": False,
        "contact_owner": "scapula",
        "contact_band_owner": "thorax",
    }
    config["trajectory_policy"] = {
        "primary_metric": "wrist_path",
        "joint_angles_are_diagnostic": True,
        "speed_and_acceleration_pass_gate": False,
        "all_action_sequence_required": True,
        "shoulder_rhythm_ratio": 2.0,
        "scapular_up_saturation_deg": 30.0,
    }
    config["optional_channels"] = {
        "S07": {
            "group": "scapula",
            "role": "v6_chest_wall_stabilizer",
            "allowed_actions": [
                "K6",
                "K7",
                "SLIDE_BACK",
                "SCAP_INT",
                "SCAP_EXT",
            ],
            "max_guides": 2,
        },
        "G08": {
            "group": "gh",
            "role": "v6_gh_auxiliary",
            "allowed_actions": ["K2", "K3", "K7"],
            "max_guides": 2,
        },
        "S08": {
            "group": "scapula",
            "role": "v6_lower_scapula_antagonist",
            "allowed_actions": [
                "K6",
                "K7",
                "SLIDE_BACK",
                "SCAP_INT",
                "SCAP_EXT",
            ],
            "max_guides": 2,
        },
    }
    routes = [
        route for route in config["routes"] if route["id"] not in ("E01", "E02")
    ]
    config["routes"] = [
        initial_elbow_route("E01", False),
        initial_elbow_route("E02", True),
        *routes,
    ]


def configure_matrix(matrix: dict, channel_ids: list[str]) -> None:
    preferred_by_group = {
        "elbow": ["E01", "E02"],
        "gh": ["G01", "G02", "G03", "G04", "G05", "G06", "G07"],
        "scapula": ["S01", "S02", "S03", "S04", "S05", "S06"],
        "combined": [
            "S01",
            "S02",
            "S03",
            "S04",
            "G01",
            "G02",
            "G07",
        ],
        "pectoral": ["P01", "P02", "G05", "G06"],
    }
    for action in matrix["actions"]:
        for key in RESULT_KEYS:
            action.pop(key, None)
        action["allowed_cables"] = list(channel_ids)
        action["preferred_cables"] = [
            route_id
            for route_id in preferred_by_group.get(action["group"], channel_ids)
            if route_id in channel_ids
        ]
        action["mandatory_cables"] = [
            route_id
            for route_id in action.get("mandatory_cables", [])
            if route_id in channel_ids
        ]
        if action["id"] == "K1":
            action["preferred_cables"] = ["E01", "E02"]
        if action["id"] == "K7":
            action["target_deg"] = 126.0
            action["acceptance_deg"] = 120.0
            action["legacy_target_deg"] = 150.0
        if action["id"] == "SCAP_INT":
            action["target_deg"] = 5.0
            action["acceptance_deg"] = 5.0
            action["legacy_target_deg"] = 20.0
    matrix["policy"]["positive_spanning_required"] = False
    matrix["policy"]["tension_mode"] = "task_conditioned"


def q_target_from_action(matrix: dict, action_id: str) -> dict[str, float]:
    action = next(item for item in matrix["actions"] if item["id"] == action_id)
    return {
        str(name): float(value)
        for name, value in action.get("q_target_deg", {}).items()
    }


def combine_q(*parts: dict[str, float]) -> dict[str, float]:
    result: dict[str, float] = {}
    for part in parts:
        result.update(part)
    return result


def build_trajectory_library(matrix: dict) -> dict:
    def scaled(values: dict[str, float], factor: float) -> dict[str, float]:
        return {
            name: float(value) * factor
            for name, value in values.items()
        }

    k6 = q_target_from_action(matrix, "K6")
    trajectories = [
        {
            "id": "FRONT_REACH",
            "label": "前伸取物",
            "duration_s": 1.6,
            "return_to_neutral": True,
            "waypoints": [
                {"time_s": 0.0, "q_target_deg": {}},
                {
                    "time_s": 0.8,
                    "q_target_deg": combine_q(
                        scaled(k6, 1.0),
                        {"gh_flex": 66.0, "elbow": 15.0},
                    ),
                },
                {"time_s": 1.6, "q_target_deg": {}},
            ],
            "active_cables": [
                "E01",
                "E02",
                "G01",
                "G02",
                "G03",
                "G04",
                "S01",
                "S02",
                "S03",
                "S04",
                "S05",
                "S06",
            ],
        },
        {
            "id": "HAND_TO_MOUTH",
            "label": "手到口鼻",
            "duration_s": 1.4,
            "return_to_neutral": False,
            "waypoints": [
                {"time_s": 0.0, "q_target_deg": {}},
                {
                    "time_s": 1.4,
                    "q_target_deg": combine_q(
                        scaled(k6, 25.0 / 30.0),
                        {"gh_flex": 25.0, "gh_rot": -10.0, "elbow": 116.0},
                    ),
                },
            ],
            "active_cables": [
                "E01",
                "E02",
                "G01",
                "G02",
                "G05",
                "G06",
                "S01",
                "S02",
                "S03",
                "S04",
                "S05",
                "S06",
            ],
        },
        {
            "id": "LATERAL_RAISE",
            "label": "侧向上举",
            "duration_s": 1.8,
            "return_to_neutral": True,
            "waypoints": [
                {"time_s": 0.0, "q_target_deg": {}},
                {
                    "time_s": 0.9,
                    "q_target_deg": combine_q(
                        scaled(k6, 1.0),
                        {"gh_abd": 66.0, "elbow": 15.0},
                    ),
                },
                {"time_s": 1.8, "q_target_deg": {}},
            ],
            "active_cables": [
                "E01",
                "E02",
                "G03",
                "G04",
                "G01",
                "G02",
                "S01",
                "S02",
                "S03",
                "S04",
                "S05",
                "S06",
            ],
        },
        {
            "id": "ALL_ACTIONS",
            "label": "全部动作连续演示",
            "duration_s": 20.0,
            "return_to_neutral": True,
            "waypoints": [
                {"time_s": 0.0, "q_target_deg": {}},
                {"time_s": 1.2, "q_target_deg": q_target_from_action(matrix, "K1")},
                {"time_s": 2.2, "q_target_deg": {}},
                {"time_s": 3.4, "q_target_deg": q_target_from_action(matrix, "K2")},
                {"time_s": 4.4, "q_target_deg": {}},
                {"time_s": 5.7, "q_target_deg": q_target_from_action(matrix, "K3")},
                {"time_s": 6.7, "q_target_deg": {}},
                {"time_s": 7.8, "q_target_deg": q_target_from_action(matrix, "K4")},
                {"time_s": 8.7, "q_target_deg": {}},
                {"time_s": 9.8, "q_target_deg": q_target_from_action(matrix, "K5")},
                {"time_s": 10.7, "q_target_deg": {}},
                {"time_s": 11.8, "q_target_deg": q_target_from_action(matrix, "K6")},
                {"time_s": 12.7, "q_target_deg": {}},
                {
                    "time_s": 13.6,
                    "q_target_deg": q_target_from_action(matrix, "SLIDE_FWD"),
                },
                {
                    "time_s": 14.4,
                    "q_target_deg": q_target_from_action(matrix, "SLIDE_BACK"),
                },
                {"time_s": 15.1, "q_target_deg": {}},
                {
                    "time_s": 16.0,
                    "q_target_deg": q_target_from_action(matrix, "SCAP_EXT"),
                },
                {
                    "time_s": 16.8,
                    "q_target_deg": q_target_from_action(matrix, "SCAP_INT"),
                },
                {
                    "time_s": 17.7,
                    "q_target_deg": q_target_from_action(matrix, "PECTORAL"),
                },
                {"time_s": 18.4, "q_target_deg": {}},
                {"time_s": 19.2, "q_target_deg": q_target_from_action(matrix, "K7")},
                {"time_s": 20.0, "q_target_deg": {}},
            ],
            "active_cables": [
                route["id"] for route in matrix.get("all_cables", [])
            ]
            or [
                "E01",
                "E02",
                "G01",
                "G02",
                "G03",
                "G04",
                "G05",
                "G06",
                "G07",
                "S01",
                "S02",
                "S03",
                "S04",
                "S05",
                "S06",
                "P01",
                "P02",
            ],
        },
    ]
    return {
        "schema_version": 6,
        "frame_rate_hz": 50.0,
        "wrist_marker": "V6_W",
        "elbow_marker": "V6_E",
        "trajectories": trajectories,
    }


def main() -> None:
    config = copy.deepcopy(
        json.loads(
            (V5_ROOT / "cable_design_v5.json").read_text(encoding="utf-8")
        )
    )
    matrix = copy.deepcopy(
        json.loads(
            (V5_ROOT / "action_cable_matrix_v5.json").read_text(
                encoding="utf-8"
            )
        )
    )
    configure_config(config)
    configure_matrix(matrix, [route["id"] for route in config["routes"]])
    trajectories = build_trajectory_library(matrix)
    V6_CONFIG.write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8"
    )
    V6_MATRIX.write_text(
        json.dumps(matrix, indent=2) + "\n", encoding="utf-8"
    )
    V6_TRAJECTORIES.write_text(
        json.dumps(trajectories, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {V6_CONFIG}")
    print(f"wrote {V6_MATRIX}")
    print(f"wrote {V6_TRAJECTORIES}")


if __name__ == "__main__":
    main()
