"""Apply the V8 contact gate only after the kinematic gate passes."""

from __future__ import annotations

import json
import math

import mujoco
import numpy as np

from v8_shoulder import (
    CONTACT_PATH,
    EXPORT_ROOT,
    TARGETS_PATH,
    build_shoulder_tree,
    compile_tree,
    contact_force,
    direct_mesh_clearance,
    load_config,
    roller_gap_mm,
    set_joint_state,
    solve_ik,
    write_model,
)


def evaluate_candidate(config: dict, targets: list[dict], count: int, preload: float) -> dict:
    model = compile_tree(build_shoulder_tree(config, roller_count=count))
    data = mujoco.MjData(model)
    stiffness = float(config["contact"]["roller"]["stiffness_N_per_m"])
    roller_joints = [
        joint
        for joint in range(model.njnt)
        if (
            mujoco.mj_id2name(
                model,
                mujoco.mjtObj.mjOBJ_JOINT,
                joint,
            )
            or ""
        ).startswith("roller_normal")
    ]
    for joint in roller_joints:
        model.jnt_stiffness[joint] = stiffness
        qpos = int(model.jnt_qposadr[joint])
        model.qpos_spring[qpos] = -preload / stiffness
    marker_names = [
        f"v8_roller_marker{suffix}"
        for suffix in (("",) if count == 1 else ("_1", "_2"))
    ]
    roller_names = [
        f"v8_shoulder_roller{suffix}"
        for suffix in (("",) if count == 1 else ("_1", "_2"))
    ]
    rows = []
    for target in targets:
        solution = solve_ik(config, target, None)
        q = np.asarray(solution["q_rad"], dtype=float)
        set_joint_state(model, data, q)
        mujoco.mj_forward(model, data)
        for _ in range(1200):
            mujoco.mj_step(model, data)
            set_joint_state(model, data, q)
        mujoco.mj_forward(model, data)
        gaps = [
            roller_gap_mm(model, data, marker, config)
            for marker in marker_names
        ]
        force, contact_count = contact_force(model, data, roller_names)
        rows.append(
            {
                "id": target["id"],
                "action": target["action"],
                "roller_gap_mm": gaps,
                "roller_contact_count": contact_count,
                "roller_contact_force_N": force,
                "full_mesh_clearance_mm": direct_mesh_clearance(config, q),
                "finite": bool(
                    np.all(np.isfinite(data.qpos))
                    and np.all(np.isfinite(data.qvel))
                ),
            }
        )
    max_detachment = max(
        max(0.0, max(row["roller_gap_mm"])) for row in rows
    )
    min_penetration = min(
        min(0.0, min(row["roller_gap_mm"])) for row in rows
    )
    min_mesh = min(row["full_mesh_clearance_mm"] for row in rows)
    max_force = max(row["roller_contact_force_N"] for row in rows)
    return {
        "roller_count": count,
        "preload_N": preload,
        "max_roller_detachment_mm": max_detachment,
        "min_roller_penetration_mm": min_penetration,
        "min_full_mesh_clearance_mm": min_mesh,
        "max_contact_force_N": max_force,
        "finite": all(row["finite"] for row in rows),
        "poses": rows,
    }


def candidate_pass(config: dict, candidate: dict) -> tuple[bool, list[str]]:
    acceptance = config["acceptance"]
    force_min, force_max = acceptance["contact_force_N"]
    checks = {
        "finite": candidate["finite"],
        "roller_detachment": candidate["max_roller_detachment_mm"]
        <= float(acceptance["roller_detachment_mm_max"]),
        "roller_penetration": candidate["min_roller_penetration_mm"]
        >= float(acceptance["roller_penetration_mm_min"]),
        "full_mesh_clearance": candidate["min_full_mesh_clearance_mm"]
        >= float(acceptance["full_mesh_clearance_mm_min"]),
        "contact_force": float(force_min)
        <= candidate["max_contact_force_N"]
        <= float(force_max),
    }
    return all(checks.values()), [
        name for name, passed in checks.items() if not passed
    ]


def main() -> None:
    config = load_config()
    kinematics_path = EXPORT_ROOT / "kinematics_report_v8.json"
    kinematics = json.loads(kinematics_path.read_text(encoding="utf-8"))
    report = {
        "version": 8,
        "stage": "contact",
        "prerequisite_kinematics_pass": bool(kinematics["pass"]),
        "pass": False,
        "selected_roller_count": None,
        "frozen_hardpoints": False,
        "nominal": [],
        "preload_sensitivity": [],
        "failure_reasons": [],
    }
    if not kinematics["pass"]:
        report["failure_reasons"] = [
            "blocked_by_kinematics",
            *[
                f"kinematics:{reason}"
                for reason in kinematics["failure_reasons"]
            ],
        ]
        EXPORT_ROOT.mkdir(parents=True, exist_ok=True)
        (EXPORT_ROOT / "contact_report_v8.json").write_text(
            json.dumps(report, indent=2) + "\n",
            encoding="utf-8",
        )
        print(
            json.dumps(
                {
                    "pass": False,
                    "failure_reasons": report["failure_reasons"],
                },
                indent=2,
            )
        )
        return

    targets = json.loads(
        TARGETS_PATH.read_text(encoding="utf-8")
    )["targets"]
    nominal = []
    for count in (1, 2):
        candidate = evaluate_candidate(
            config,
            targets,
            count,
            float(config["contact"]["roller"]["preload_N"]),
        )
        passed, reasons = candidate_pass(config, candidate)
        candidate["pass"] = passed
        candidate["failure_reasons"] = reasons
        nominal.append(candidate)
    passing = [candidate for candidate in nominal if candidate["pass"]]
    selected = passing[0] if passing else min(
        nominal,
        key=lambda item: (
            item["max_roller_detachment_mm"],
            max(0.0, -item["min_full_mesh_clearance_mm"]),
        ),
    )
    count = int(selected["roller_count"])
    sensitivity = []
    for preload in config["contact"]["roller"]["preload_sensitivity_N"]:
        candidate = evaluate_candidate(
            config,
            targets,
            count,
            float(preload),
        )
        passed, reasons = candidate_pass(config, candidate)
        candidate["pass"] = passed
        candidate["failure_reasons"] = reasons
        sensitivity.append(candidate)
    preload_stable = all(item["pass"] for item in sensitivity)
    reasons = []
    if not selected["pass"]:
        reasons.extend(selected["failure_reasons"])
    if not preload_stable:
        reasons.append("preload_sensitivity")
    report.update(
        {
            "pass": bool(selected["pass"] and preload_stable),
            "selected_roller_count": count,
            "nominal": nominal,
            "preload_sensitivity": sensitivity,
            "failure_reasons": reasons,
        }
    )
    if report["pass"]:
        write_model(build_shoulder_tree(config, roller_count=count), CONTACT_PATH)
        report["frozen_hardpoints"] = True
    EXPORT_ROOT.mkdir(parents=True, exist_ok=True)
    (EXPORT_ROOT / "contact_report_v8.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "pass": report["pass"],
                "selected_roller_count": count,
                "max_roller_detachment_mm": selected[
                    "max_roller_detachment_mm"
                ],
                "min_full_mesh_clearance_mm": selected[
                    "min_full_mesh_clearance_mm"
                ],
                "max_contact_force_N": selected["max_contact_force_N"],
                "failure_reasons": reasons,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

