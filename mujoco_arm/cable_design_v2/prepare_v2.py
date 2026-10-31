"""Generate the V2 cable and action configuration from the accepted 16 routes."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from v2_common import (
    ACTION_MATRIX_PATH,
    BODY_NAME,
    BODY_ORIGIN_MM,
    CONFIG_PATH,
    GROUP_BY_CABLE,
    MUJOCO_ROOT,
    POINTS_MM,
    PROJECT_ROOT,
    ROOT,
    action_definitions,
    design_to_mj,
    guide_geometry_from_bend,
)


CURRENT_CONFIG = MUJOCO_ROOT / "config" / "tendon_demo_16.json"
ENGINEERING_CONFIG = MUJOCO_ROOT / "engineering" / "k6_baseline.json"
FOUR_MUSCLE_LAYOUT = ROOT / "four_muscle_layout.json"
ROUTE_ROLE = {
    "E01": "elbow_flexion",
    "E02": "elbow_extension",
    "G01": "gh_flexion",
    "G02": "gh_extension",
    "G03": "gh_abduction",
    "G04": "gh_adduction",
    "G05": "gh_internal_rotation",
    "G06": "gh_external_rotation",
    "S01": "upper_scapula_pull_down_medial",
    "S02": "upper_scapula_pull_down_lateral",
    "S03": "lower_scapula_pull_up_medial",
    "S04": "lower_scapula_pull_up_lateral",
    "S05": "auxiliary_scapula_rotation_positive",
    "S06": "antagonistic_scapula_rotation",
    "P01": "pectoral_clavicular_head",
    "P02": "pectoral_sternocostal_head",
}


def make_route(route: dict) -> dict:
    route_id = route["id"]
    owners = route["owners"]
    points = [np.asarray(point, dtype=float) for point in route["points_mm"]]
    path = []
    guide_note = None
    if route["has_guide"]:
        a, g, b = points
        radius = 8.0
        centre, axis, side, wrap_angle = guide_geometry_from_bend(
            a, g, b, radius
        )
        if wrap_angle < 10.0:
            guide_note = (
                f"guide omitted: bend point produces only {wrap_angle:.2f} deg "
                "of wrap for an 8 mm guide"
            )
            route["has_guide"] = False
            path.append(
                {
                    "kind": "anchor",
                    "owner": owners[0],
                    "body": BODY_NAME[owners[0]],
                    "world_mm": a.tolist(),
                    "local_mm": design_to_mj(
                        a, BODY_ORIGIN_MM[owners[0]]
                    ).tolist(),
                    "installation_normal": side.tolist(),
                }
            )
            path.append(
                {
                    "kind": "pass_through",
                    "owner": owners[1],
                    "body": BODY_NAME[owners[1]],
                    "world_mm": g.tolist(),
                    "local_mm": design_to_mj(
                        g, BODY_ORIGIN_MM[owners[1]]
                    ).tolist(),
                    "installation_normal": side.tolist(),
                }
            )
            path.append(
                {
                    "kind": "anchor",
                    "owner": owners[2],
                    "body": BODY_NAME[owners[2]],
                    "world_mm": b.tolist(),
                    "local_mm": design_to_mj(
                        b, BODY_ORIGIN_MM[owners[2]]
                    ).tolist(),
                    "installation_normal": side.tolist(),
                }
            )
        else:
            path.append(
                {
                    "kind": "anchor",
                    "owner": owners[0],
                    "body": BODY_NAME[owners[0]],
                    "world_mm": a.tolist(),
                    "local_mm": design_to_mj(
                        a, BODY_ORIGIN_MM[owners[0]]
                    ).tolist(),
                    "installation_normal": side.tolist(),
                }
            )
            path.append(
                {
                    "kind": "guide",
                    "name": f"guide_{route_id}_01",
                    "owner": owners[1],
                    "body": BODY_NAME[owners[1]],
                    "bend_point_mm": g.tolist(),
                    "world_mm": centre.tolist(),
                    "local_mm": design_to_mj(
                        centre, BODY_ORIGIN_MM[owners[1]]
                    ).tolist(),
                    "radius_mm": radius,
                    "axis_design": axis.tolist(),
                    "side_design": side.tolist(),
                    "wrap_angle_deg": wrap_angle,
                }
            )
            path.append(
                {
                    "kind": "anchor",
                    "owner": owners[2],
                    "body": BODY_NAME[owners[2]],
                    "world_mm": b.tolist(),
                    "local_mm": design_to_mj(
                        b, BODY_ORIGIN_MM[owners[2]]
                    ).tolist(),
                    "installation_normal": side.tolist(),
                }
            )
    else:
        a, b = points
        path.append(
            {
                "kind": "anchor",
                "owner": owners[0],
                "body": BODY_NAME[owners[0]],
                "world_mm": a.tolist(),
                "local_mm": design_to_mj(
                    a, BODY_ORIGIN_MM[owners[0]]
                ).tolist(),
                "installation_normal": [0.0, 1.0, 0.0],
            }
        )
        path.append(
            {
                "kind": "anchor",
                "owner": owners[-1],
                "body": BODY_NAME[owners[-1]],
                "world_mm": b.tolist(),
                "local_mm": design_to_mj(
                    b, BODY_ORIGIN_MM[owners[-1]]
                ).tolist(),
                "installation_normal": [0.0, 1.0, 0.0],
            }
        )
    return {
        "id": route_id,
        "group": GROUP_BY_CABLE[route_id],
        "role": ROUTE_ROLE[route_id],
        "diameter_mm": 1.0,
        "max_tension_N": 150.0,
        "max_guides": 3,
        "path": path,
        "tpu": {
            "enabled": True,
            "reference_strain": 0.10,
            "reference_force_N": 20.0,
            "compression_force_N": 0.0,
        },
        "status": "initial",
        "notes": [guide_note] if guide_note else [],
    }


def apply_four_muscle_layout(route: dict, layout: dict, config: dict) -> dict:
    route_id = route["id"]
    if route_id not in layout["routes"] or len(route["path"]) != 2:
        return route
    item = layout["routes"][route_id]
    a = np.asarray(item["A_mm"], dtype=float)
    b = np.asarray(item["B_mm"], dtype=float)
    band = config["contact"]["thorax_band"]
    radius = float(band["radius_mm"])
    x0 = float(band["x_min_mm"])
    x1 = float(band["x_max_mm"])
    semi_z = 106.0
    z0 = semi_z * np.sqrt(max(0.0, 1.0 - (x0 / 140.0) ** 2))
    z1 = semi_z * np.sqrt(max(0.0, 1.0 - (x1 / 140.0) ** 2))
    p0 = np.array([x0, z0])
    p1 = np.array([x1, z1])
    chord = p1 - p0
    midpoint = 0.5 * (p0 + p1)
    perpendicular = np.array([-chord[1], chord[0]]) / np.linalg.norm(chord)
    height = np.sqrt(radius**2 - (np.linalg.norm(chord) / 2.0) ** 2)
    c0 = midpoint + height * perpendicular
    c1 = midpoint - height * perpendicular
    centre = c0 if np.linalg.norm(c0) < np.linalg.norm(c1) else c1
    normal_xz = np.array([a[0] - centre[0], a[2] - centre[1]])
    normal_xz /= np.linalg.norm(normal_xz)
    thorax_normal = [float(normal_xz[0]), 0.0, float(normal_xz[1])]
    scapula_normal = [0.6556, 0.0, 0.7552]
    route["path"][0] = {
        "kind": "anchor",
        "owner": "T",
        "body": "thorax",
        "world_mm": a.tolist(),
        "local_mm": design_to_mj(a).tolist(),
        "installation_normal": thorax_normal,
    }
    route["path"][1] = {
        "kind": "anchor",
        "owner": "S",
        "body": "scapula",
        "world_mm": b.tolist(),
        "local_mm": design_to_mj(b - np.asarray(POINTS_MM["AC"])).tolist(),
        "installation_normal": scapula_normal,
    }
    route["status"] = "four_muscle_optimized"
    route["layout_source"] = layout["source"]
    return route


def main() -> None:
    current = json.loads(CURRENT_CONFIG.read_text(encoding="utf-8"))
    engineering = json.loads(ENGINEERING_CONFIG.read_text(encoding="utf-8"))
    routes = [make_route(route) for route in current["routes"]]
    config = {
        "model_name": "arm_cable_design_v2",
        "coordinate_frame": {
            "origin": "thorax bottom centre",
            "x": "right",
            "y": "up",
            "z": "back",
            "units": "mm",
        },
        "materials": {
            "polymer_density_kg_m3": 1100.0,
            "tpu_density_kg_m3": 1200.0,
            "cable_diameter_mm": 1.0,
        },
        "load_case": {
            "gravity_m_s2": 9.81,
            "external_payload_kg": 0.0,
        },
        "contact": {
            "band_mesh": str(
                (MUJOCO_ROOT / "engineering" / "assets" / "thorax_contact_band.stl")
                .resolve()
            ),
            "pad_mesh": str(
                (MUJOCO_ROOT / "engineering" / "assets" / "scapula_contact_pad.stl")
                .resolve()
            ),
            "thorax_band": engineering["contact"]["thorax_band"],
            "scapula_pad": engineering["contact"]["scapula_pad"],
            "friction": [0.2, 0.005, 0.0001],
            "friction_sensitivity": [0.1, 0.2, 0.3],
            "condim": 3,
            "solref": [0.002, 1.0],
            "solimp": [0.9, 0.95, 0.001],
        },
        "guides": {
            "minimum_radius_mm": 8.0,
            "maximum_per_cable": 3,
            "maximum_wrap_deg": 180.0,
            "friction_coefficient": 0.05,
            "friction_sensitivity": [0.0, 0.05, 0.10],
        },
        "joint_limits_deg": {
            "sc_elevation": [-10.0, 55.0],
            "sc_protraction": [-25.0, 25.0],
            "ac_ext": [-35.0, 35.0],
            "gh_flex": [-120.0, 120.0],
            "gh_abd": [-120.0, 120.0],
            "gh_rot": [-110.0, 110.0],
            "elbow": [0.0, 145.0],
        },
        "acceptance": {
            "minimum_visible_deg": 5.0,
            "minimum_visible_mm": 5.0,
            "cable_contraction_percent_max": 18.0,
            "cable_clearance_mm_min": 2.0,
            "unintended_mesh_penetration_mm_max": 0.5,
            "tension_max_N": 150.0,
            "action_pose_fractions": [0.0, 0.25, 0.5, 0.75, 1.0],
        },
        "simulation": {
            "timestep_s": 0.0005,
            "settle_steps": 1200,
            "pose_steps": 900,
        },
        "routes": routes,
    }
    if FOUR_MUSCLE_LAYOUT.exists():
        layout = json.loads(FOUR_MUSCLE_LAYOUT.read_text(encoding="utf-8"))
        config["routes"] = [
            apply_four_muscle_layout(route, layout, config)
            for route in config["routes"]
        ]
    CONFIG_PATH.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    actions = action_definitions()
    matrix = {
        "policy": {
            "action_first": True,
            "allowed_cable_groups": True,
            "minimum_visible_deg": 5.0,
            "minimum_visible_mm": 5.0,
        },
        "actions": actions,
        "initial_assignment": {
            action["id"]: action["allowed_cables"] for action in actions
        },
    }
    ACTION_MATRIX_PATH.write_text(json.dumps(matrix, indent=2) + "\n", encoding="utf-8")

    md_lines = [
        "# Action-to-cable responsibility matrix",
        "",
        "The optimizer may choose 2-4 cooperating cables from each allowed set.",
        "",
        "| Action | Goal | 80% / minimum | Allowed cables |",
        "| --- | --- | --- | --- |",
    ]
    for action in actions:
        if "target_deg" in action:
            goal = f"{action['target_deg']:.0f} deg"
            minimum = action.get("acceptance_deg", action["target_deg"] * 0.8)
            minimum_text = f"{minimum:.0f} deg / {action['minimum_visible_deg']:.0f} deg"
        else:
            goal = f"{action['target_mm']:.0f} mm"
            minimum_text = (
                f"{action['target_mm'] * 0.8:.0f} mm / "
                f"{action['minimum_visible_mm']:.0f} mm"
            )
        md_lines.append(
            f"| {action['id']} | {goal} | {minimum_text} | "
            f"{', '.join(action['allowed_cables'])} |"
        )
    (PROJECT_ROOT / "design" / "ACTION_CABLE_MATRIX_V2.md").write_text(
        "\n".join(md_lines) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {CONFIG_PATH}")
    print(f"wrote {ACTION_MATRIX_PATH}")
    print(f"wrote {PROJECT_ROOT / 'design' / 'ACTION_CABLE_MATRIX_V2.md'}")


if __name__ == "__main__":
    main()
