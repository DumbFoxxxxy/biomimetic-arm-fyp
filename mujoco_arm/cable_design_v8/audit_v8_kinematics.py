"""Audit V8 five-DOF kinematics, IK residuals and full-mesh clearance."""

from __future__ import annotations

import json
import math

import mujoco
import numpy as np

from v8_shoulder import (
    EXPORT_ROOT,
    OPEN_CHAIN_PATH,
    TARGETS_PATH,
    actual_body_geometry,
    compile_tree,
    build_shoulder_tree,
    direct_mesh_clearance,
    geom_mesh_vertices,
    load_config,
    mesh_clearance,
    points,
    pose_error,
    set_joint_state,
    solve_ik,
    surface_target,
)


def main() -> None:
    config = load_config()
    target_data = json.loads(TARGETS_PATH.read_text(encoding="utf-8"))
    targets = target_data["targets"]
    model = compile_tree(build_shoulder_tree(config, roller_count=0))
    data = mujoco.MjData(model)
    scapula_vertices = geom_mesh_vertices(model, "scapula_mesh")
    pad_vertices = geom_mesh_vertices(model, "v8_scapula_contact_pad")
    point = points(config)
    bone_length_expected = float(
        config["geometry_mm"]["expected_sc_ac_mm"]
    )
    previous_by_action: dict[str, np.ndarray] = {}
    rows = []
    for target in targets:
        action = target["action"]
        initial = previous_by_action.get(action)
        solution = solve_ik(config, target, initial)
        q = np.asarray(solution["q_rad"], dtype=float)
        set_joint_state(model, data, q)
        mujoco.mj_forward(model, data)
        actual = actual_body_geometry(model, data, config)
        _, pose_metrics = pose_error(config, q, target)
        orientation_error = pose_metrics["orientation_error_deg"]
        position_error = float(
            np.linalg.norm(
                actual["contact_mm"]
                - np.asarray(target["contact_mm"], dtype=float)
            )
        )
        sc_ac = float(np.linalg.norm(actual["ac_mm"] - point["SC"]))
        mesh_min = mesh_clearance(
            model,
            data,
            "scapula_mesh",
            scapula_vertices,
            config,
        )
        pad_min = mesh_clearance(
            model,
            data,
            "v8_scapula_contact_pad",
            pad_vertices,
            config,
        )
        sc_deg = np.degrees(q[:2])
        ac_deg = np.degrees(q[2:])
        previous = previous_by_action.get(action)
        sc_step = None
        ac_step = None
        if previous is not None:
            sc_step = float(np.linalg.norm(sc_deg - np.degrees(previous[:2])))
            ac_step = float(
                np.linalg.norm(ac_deg - np.degrees(previous[2:]))
            )
        previous_by_action[action] = q
        rows.append(
            {
                "id": target["id"],
                "action": action,
                "fraction": float(target["fraction"]),
                "surface": {
                    "s_mm": float(target["s_mm"]),
                    "t_mm": float(target["t_mm"]),
                    "phi_deg": float(target["phi_deg"]),
                },
                "q_deg": {
                    "sc_elevation": float(sc_deg[0]),
                    "sc_protraction": float(sc_deg[1]),
                    "ac_u": float(ac_deg[0]),
                    "ac_v": float(ac_deg[1]),
                    "ac_w": float(ac_deg[2]),
                },
                "ik_success": bool(solution["success"]),
                "ik_nfev": int(solution["nfev"]),
                "position_residual_mm": position_error,
                "orientation_residual_deg": orientation_error,
                "perpendicular_tilt_deg": pose_metrics[
                    "perpendicular_tilt_deg"
                ],
                "sc_ac_mm": sc_ac,
                "bone_length_error_mm": abs(sc_ac - bone_length_expected),
                "full_mesh_clearance_mm": mesh_min,
                "contact_pad_clearance_mm": pad_min,
                "sc_step_deg": sc_step,
                "ac_step_deg": ac_step,
                "finite": bool(
                    np.all(np.isfinite(q))
                    and np.all(np.isfinite(data.qpos))
                    and np.all(np.isfinite(data.qvel))
                ),
            }
        )

    acceptance = config["acceptance"]
    max_position = max(row["position_residual_mm"] for row in rows)
    max_orientation = max(row["orientation_residual_deg"] for row in rows)
    max_bone_error = max(row["bone_length_error_mm"] for row in rows)
    min_mesh = min(row["full_mesh_clearance_mm"] for row in rows)
    min_pad = min(row["contact_pad_clearance_mm"] for row in rows)
    max_sc_step = max(
        (row["sc_step_deg"] for row in rows if row["sc_step_deg"] is not None),
        default=0.0,
    )
    max_ac_step = max(
        (row["ac_step_deg"] for row in rows if row["ac_step_deg"] is not None),
        default=0.0,
    )
    failure_diagnostics = {}
    worst_row = min(rows, key=lambda row: row["full_mesh_clearance_mm"])
    if worst_row["full_mesh_clearance_mm"] < float(
        acceptance["full_mesh_clearance_mm_min"]
    ):
        best_nearby = None
        target_row = next(
            target for target in targets if target["id"] == worst_row["id"]
        )
        for t_value in np.linspace(
            target_row["t_mm"] - 30.0,
            target_row["t_mm"] + 30.0,
            25,
        ):
            candidate = surface_target(
                config,
                float(target_row["s_mm"]),
                float(t_value),
                float(target_row["phi_deg"]),
            )
            candidate_solution = solve_ik(config, candidate, None)
            candidate_clearance = direct_mesh_clearance(
                config,
                np.asarray(candidate_solution["q_rad"], dtype=float),
            )
            if (
                best_nearby is None
                or candidate_clearance > best_nearby["clearance_mm"]
            ):
                best_nearby = {
                    "t_mm": float(t_value),
                    "clearance_mm": candidate_clearance,
                    "position_residual_mm": candidate_solution[
                        "position_error_mm"
                    ],
                    "orientation_residual_deg": candidate_solution[
                        "orientation_error_deg"
                    ],
                }
        feasible_s = None
        if worst_row["action"] == "SLIDE_BACK":
            for magnitude in np.arange(0.0, 15.01, 0.5):
                candidate = surface_target(
                    config,
                    -float(magnitude),
                    float(best_nearby["t_mm"]),
                    float(target_row["phi_deg"]),
                )
                candidate_solution = solve_ik(config, candidate, None)
                candidate_clearance = direct_mesh_clearance(
                    config,
                    np.asarray(candidate_solution["q_rad"], dtype=float),
                )
                if candidate_clearance >= float(
                    acceptance["full_mesh_clearance_mm_min"]
                ):
                    feasible_s = {
                        "s_mm": -float(magnitude),
                        "clearance_mm": candidate_clearance,
                    }
        failure_diagnostics = {
            "worst_target": worst_row["id"],
            "measured_clearance_mm": worst_row[
                "full_mesh_clearance_mm"
            ],
            "allowed_clearance_mm": float(
                acceptance["full_mesh_clearance_mm_min"]
            ),
            "best_axial_correction": best_nearby,
            "maximum_clearance_feasible_slide_back": feasible_s,
        }
    sc_ranges = config["joints"]["sc_elevation"]["range_deg"]
    sp_ranges = config["joints"]["sc_protraction"]["range_deg"]
    ac_limit = float(config["joints"]["ac_ball"]["component_limit_deg"])
    sc_limits_ok = all(
        float(sc_ranges[0]) - 1e-9
        <= row["q_deg"]["sc_elevation"]
        <= float(sc_ranges[1]) + 1e-9
        and float(sp_ranges[0]) - 1e-9
        <= row["q_deg"]["sc_protraction"]
        <= float(sp_ranges[1]) + 1e-9
        for row in rows
    )
    ac_limits_ok = all(
        abs(row["q_deg"][name]) <= ac_limit + 1e-9
        for row in rows
        for name in ("ac_u", "ac_v", "ac_w")
    )
    finite_ok = all(row["finite"] for row in rows)
    joint_objects = int(model.njnt)
    shoulder_dofs = int(model.nv)
    checks = {
        "model_loads": True,
        "joint_objects_exact": joint_objects == int(
            acceptance["joint_object_count"]
        ),
        "scalar_dof_exact": shoulder_dofs == int(
            acceptance["scalar_dof_count"]
        ),
        "finite": finite_ok,
        "position_residual": max_position
        <= float(acceptance["position_residual_mm_max"]),
        "orientation_residual": max_orientation
        <= float(acceptance["orientation_residual_deg_max"]),
        "bone_length": max_bone_error
        <= float(acceptance["bone_length_error_mm_max"]),
        "full_mesh_clearance": min_mesh
        >= float(acceptance["full_mesh_clearance_mm_min"]),
        "contact_pad_clearance": min_pad
        >= float(acceptance["contact_pad_penetration_mm_min"]),
        "sc_limits": sc_limits_ok,
        "ac_limits": ac_limits_ok,
        "sc_continuity": max_sc_step
        <= float(config["ik"]["max_adjacent_sc_deg"]),
        "ac_continuity": max_ac_step
        <= float(config["ik"]["max_adjacent_ac_deg"]),
    }
    failure_reasons = [
        name for name, passed in checks.items() if not bool(passed)
    ]
    report = {
        "version": 8,
        "stage": "kinematics",
        "model": str(OPEN_CHAIN_PATH),
        "joint_objects": joint_objects,
        "shoulder_scalar_dofs": shoulder_dofs,
        "target_count": len(rows),
        "checks": checks,
        "max_position_residual_mm": max_position,
        "max_orientation_residual_deg": max_orientation,
        "max_bone_length_error_mm": max_bone_error,
        "min_full_mesh_clearance_mm": min_mesh,
        "min_contact_pad_clearance_mm": min_pad,
        "max_sc_step_deg": max_sc_step,
        "max_ac_step_deg": max_ac_step,
        "pass": not failure_reasons,
        "failure_reasons": failure_reasons,
        "failure_diagnostics": failure_diagnostics,
        "poses": rows,
    }
    EXPORT_ROOT.mkdir(parents=True, exist_ok=True)
    output = EXPORT_ROOT / "kinematics_report_v8.json"
    output.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "pass": report["pass"],
                "target_count": len(rows),
                "joint_objects": joint_objects,
                "shoulder_scalar_dofs": shoulder_dofs,
                "max_position_residual_mm": max_position,
                "max_orientation_residual_deg": max_orientation,
                "min_full_mesh_clearance_mm": min_mesh,
                "min_contact_pad_clearance_mm": min_pad,
                "failure_reasons": failure_reasons,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

