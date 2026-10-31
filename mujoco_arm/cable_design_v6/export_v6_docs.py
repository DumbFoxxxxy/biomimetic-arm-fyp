"""Export V6 reports, route tables and SolidWorks hardpoints."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
V5_ROOT = ROOT.parent / "cable_design_v5"
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v6"
DESIGN = ROOT.parent.parent / "design"
sys.path[:0] = [str(ROOT), str(V5_ROOT), str(V2_ROOT)]

from evaluate_v2 import (  # noqa: E402
    export_action_markdown,
    export_action_matrix,
    export_routes_csv,
    export_tension_report,
)
from v2_common import load_json, mj_to_design, normalize  # noqa: E402
from v6_evaluator import V6Evaluator  # noqa: E402


CONFIG = ROOT / "cable_design_v6.json"
MATRIX = ROOT / "action_cable_matrix_v6.json"
MODEL = ROOT / "arm_cable_design_v6.xml"
LIBRARY = ROOT / "trajectory_library_v6.json"


def route_rows(evaluator: V6Evaluator, config: dict) -> list[dict]:
    zero = evaluator._zero_data()
    rows = []
    for route in config["routes"]:
        for index, element in enumerate(route["path"]):
            if element["kind"] == "anchor":
                site_name = f"V2_{'AB'[sum(e['kind'] == 'anchor' for e in route['path'][:index])]}_{route['id']}"
                site_id = mujoco.mj_name2id(
                    evaluator.model,
                    mujoco.mjtObj.mjOBJ_SITE,
                    site_name,
                )
                world = mj_to_design(zero.site_xpos[site_id])
                axis = np.zeros(3)
                radius = 0.0
            else:
                geom_name = (
                    f"V2_guide_{route['id']}_{element['name']}"
                )
                geom_id = mujoco.mj_name2id(
                    evaluator.model,
                    mujoco.mjtObj.mjOBJ_GEOM,
                    geom_name,
                )
                world = mj_to_design(zero.geom_xpos[geom_id])
                axis = normalize(
                    np.asarray(
                        element.get("axis_design", [0.0, 0.0, 1.0])
                    )
                )
                radius = float(element.get("radius_mm", 0.0))
            normal = normalize(
                np.asarray(
                    element.get(
                        "installation_normal",
                        element.get("side_design", [0.0, 0.0, 1.0]),
                    )
                )
            )
            local = np.asarray(element["local_mm"], dtype=float)
            rows.append(
                {
                    "channel": route["id"],
                    "group": route["group"],
                    "role": route.get("role", ""),
                    "sequence": index,
                    "kind": element["kind"],
                    "part": element.get("body", ""),
                    "world_x_mm": float(world[0]),
                    "world_y_mm": float(world[1]),
                    "world_z_mm": float(world[2]),
                    "local_x_m": float(local[0]),
                    "local_y_m": float(local[1]),
                    "local_z_m": float(local[2]),
                    "normal_x": float(normal[0]),
                    "normal_y": float(normal[1]),
                    "normal_z": float(normal[2]),
                    "axis_x": float(axis[0]),
                    "axis_y": float(axis[1]),
                    "axis_z": float(axis[2]),
                    "radius_mm": radius,
                    "wrap_angle_deg": float(
                        element.get("wrap_angle_deg", 0.0)
                    ),
                    "preload_N": "",
                    "stiffness_N_per_mm": "",
                    "travel_mm": "",
                }
            )
    roller_site = mujoco.mj_name2id(
        evaluator.model,
        mujoco.mjtObj.mjOBJ_SITE,
        "V6_roller",
    )
    roller_world = mj_to_design(zero.site_xpos[roller_site])
    guide = config["shoulder_guide"]
    rows.append(
        {
            "channel": "SHOULDER_GUIDE",
            "group": "scapula_contact",
            "role": guide["type"],
            "sequence": 0,
            "kind": "roller",
            "part": "scapula",
            "world_x_mm": float(roller_world[0]),
            "world_y_mm": float(roller_world[1]),
            "world_z_mm": float(roller_world[2]),
            "local_x_m": "",
            "local_y_m": "",
            "local_z_m": "",
            "normal_x": 0.6556,
            "normal_y": 0.0,
            "normal_z": 0.7552,
            "axis_x": 0.0,
            "axis_y": 1.0,
            "axis_z": 0.0,
            "radius_mm": float(guide["roller_radius_mm"]),
            "wrap_angle_deg": "",
            "preload_N": float(guide["nominal_preload_N"]),
            "stiffness_N_per_mm": float(
                guide["equivalent_stiffness_N_per_mm"]
            ),
            "travel_mm": float(guide["normal_travel_mm"]),
        }
    )
    return rows


def write_anchors_markdown(path: Path, rows: list[dict], report: dict) -> None:
    lines = [
        "# V6 钢丝与肩胸导向硬点",
        "",
        f"- Full pass: `{report['full_pass']}`",
        f"- Frozen hardpoints: `{report['frozen_hardpoints']}`",
        f"- Positive spanning hard gate: `False`",
        "",
        "| 通道 | 类型 | 零件 | 顺序 | 世界坐标 (mm) | 导轮半径 | 安装法线 |",
        "| --- | --- | --- | ---: | --- | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['channel']} | {row['kind']} | {row['part']} | "
            f"{row['sequence']} | ({row['world_x_mm']:.2f}, "
            f"{row['world_y_mm']:.2f}, {row['world_z_mm']:.2f}) | "
            f"{row['radius_mm']:.2f} | "
            f"({row['normal_x']:.4f}, {row['normal_y']:.4f}, "
            f"{row['normal_z']:.4f}) |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_decisions(
    path: Path,
    config: dict,
    report: dict,
    library: dict,
) -> None:
    lines = [
        "# V6 design decisions",
        "",
        "## Architecture",
        "",
        "- Seven controlled DOFs, 17-20 independent cable channels.",
        "- Global wrench closure is diagnostic only; local nonnegative tension is checked.",
        "- Two active elbow antagonist cables with R10 guides.",
        "- R400 single-sided crowned POM roller on the scapula with compliant preload.",
        "- Three ADL trajectories plus one all-action sequence.",
        "",
        "## Current audit",
        "",
        f"- Structural pass: `{report['structural_pass']}`",
        f"- Static actions pass: `{report['static_actions_pass']}`",
        f"- Elbow antagonism pass: `{report['elbow_antagonism']['pass']}`",
        f"- Full trajectory pass: `{report['full_pass']}`",
        f"- Frozen hardpoints: `{report['frozen_hardpoints']}`",
        "",
        "## Trajectory library",
        "",
    ]
    for trajectory in library["trajectories"]:
        lines.append(
            f"- `{trajectory['id']}` {trajectory['label']}: "
            f"{trajectory['duration_s']:.1f} s"
        )
    lines.extend(["", "## Failure reasons", ""])
    if report["failure_reasons"]:
        lines.extend(f"- `{reason}`" for reason in report["failure_reasons"])
    else:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Shoulder preload",
            "",
            f"- Nominal: {config['shoulder_guide']['nominal_preload_N']:.1f} N",
            f"- Stiffness: {config['shoulder_guide']['equivalent_stiffness_N_per_mm']:.1f} N/mm",
            f"- Travel: {config['shoulder_guide']['normal_travel_mm']:.1f} mm",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    config = load_json(CONFIG)
    matrix = load_json(MATRIX)
    library = load_json(LIBRARY)
    evaluator = V6Evaluator(MODEL, config, matrix)
    static = load_json(EXPORT / "static_evaluation_v6.json")
    report = load_json(EXPORT / "trajectory_report_v6.json")
    rows = route_rows(evaluator, config)
    hardpoint_path = EXPORT / "SolidWorks_hardpoints_v6.csv"
    with hardpoint_path.open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    export_routes_csv(EXPORT / "cable_routes_v6.csv", config)
    export_solidworks_legacy = EXPORT / "cable_routes_v6.csv"
    export_tension_report(EXPORT / "tension_report_v6.json", static)
    updated_matrix = export_action_matrix(
        MATRIX,
        static,
        matrix,
    )
    MATRIX.write_text(
        json.dumps(updated_matrix, indent=2) + "\n",
        encoding="utf-8",
    )
    DESIGN.mkdir(parents=True, exist_ok=True)
    export_action_markdown(
        DESIGN / "ACTION_CABLE_MATRIX_V6.md",
        updated_matrix,
    )
    write_anchors_markdown(
        DESIGN / "CABLE_ANCHORS_V6.md",
        rows,
        report,
    )
    write_decisions(
        DESIGN / "V6_DESIGN_DECISIONS.md",
        config,
        report,
        library,
    )
    print(f"wrote {hardpoint_path}")
    print(f"wrote {export_solidworks_legacy}")


if __name__ == "__main__":
    main()
