"""Run scripted slider motions and report the joint angles reached."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import mujoco
import numpy as np

from cable_drive_common import (
    apply_one_sided_tension,
    cable_assist_pairs,
    DEFAULT_MANIFEST,
    DEFAULT_OUTPUT_XML,
    DEFAULT_REPORT_DIR,
    load_cable_manifest,
)
from check_cable_drive import cable_channels, sensitivity_matrix
from shoulder_coupling import apply_shoulder_rhythm_from_state


def actuator_id(model: mujoco.MjModel, name: str) -> int:
    value = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, name)
    if value < 0:
        raise RuntimeError(f"Actuator missing: {name}")
    return int(value)


def tendon_id(model: mujoco.MjModel, name: str) -> int:
    value = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_TENDON, name)
    if value < 0:
        raise RuntimeError(f"Tendon missing: {name}")
    return int(value)


def joint_id(model: mujoco.MjModel, name: str) -> int:
    value = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
    if value < 0:
        raise RuntimeError(f"Joint missing: {name}")
    return int(value)


def tendon_force(model: mujoco.MjModel, data: mujoco.MjData, tid: int) -> float:
    stiffness = float(model.tendon_stiffness[tid])
    springlength = float(model.tendon_lengthspring[tid, 0])
    damping = float(model.tendon_damping[tid])
    return (
        stiffness * (float(data.ten_length[tid]) - springlength)
        + damping * float(data.ten_velocity[tid])
    )


def select_channels(
    channels: list[dict],
    sensitivities: dict[str, dict[str, float]],
    joint_name: str,
    ids: tuple[str, ...],
) -> list[str]:
    candidates = []
    for channel in channels:
        if ids and not any(token in channel["muscle_id"] for token in ids):
            continue
        value = sensitivities[joint_name][channel["id"]]
        candidates.append((abs(value), channel["id"]))
    candidates.sort(reverse=True)
    return [
        channel_id
        for magnitude, channel_id in candidates[:3]
        if magnitude > 1e-6
    ]


def run_motion(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    channels: list[dict],
    selected: list[str],
    joint_name: str,
    gravity: tuple[float, float, float],
    output_csv: Path,
    release: float = 0.08,
    ramp_seconds: float = 2.0,
    hold_seconds: float = 0.5,
) -> dict:
    model.opt.gravity[:] = gravity
    data.qpos[:] = 0.0
    data.qvel[:] = 0.0
    data.ctrl[:] = 0.0
    mujoco.mj_forward(model, data)

    actuator_ids = {
        channel_id: actuator_id(model, f"pull_{channel_id}")
        for channel_id in selected
    }
    assist_pairs = cable_assist_pairs(model)
    slack_actuators = {
        channel["id"]: actuator_id(model, f"pull_{channel['id']}")
        for channel in channels
        if channel["id"] not in selected
    }
    travel = {
        channel["id"]: float(channel["travel_m"])
        for channel in channels
        if channel["id"] in selected
    }
    target_joint = joint_id(model, joint_name)
    address = int(model.jnt_qposadr[target_joint])

    rows = []
    total_steps = int((ramp_seconds + hold_seconds) / model.opt.timestep)
    initial_angle = float(data.qpos[address])
    maximum_angle = initial_angle
    minimum_angle = initial_angle

    for step in range(total_steps):
        seconds = step * model.opt.timestep
        ratio = min(1.0, seconds / ramp_seconds)
        for channel_id, actuator in actuator_ids.items():
            data.ctrl[actuator] = travel[channel_id] * ratio
        # Slacken every cable that is not doing the pulling, so the idle wires
        # do not act as springs fighting the motion.
        for actuator in slack_actuators.values():
            data.ctrl[actuator] = -release
        apply_one_sided_tension(model, data, assist_pairs)
        apply_shoulder_rhythm_from_state(model, data)
        mujoco.mj_step(model, data)

        angle = float(data.qpos[address])
        maximum_angle = max(maximum_angle, angle)
        minimum_angle = min(minimum_angle, angle)
        row = {
            "time_s": seconds,
            "joint": joint_name,
            "angle_deg": np.degrees(angle),
        }
        for channel_id, actuator in actuator_ids.items():
            row[f"pull_{channel_id}_m"] = float(data.ctrl[actuator])
        for channel in channels:
            if channel["id"] not in selected:
                continue
            tid = tendon_id(model, channel["cable"])
            row[f"tension_{channel['id']}_n"] = tendon_force(model, data, tid)
        rows.append(row)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    return {
        "joint": joint_name,
        "selected_cables": selected,
        "gravity": list(gravity),
        "initial_angle_deg": float(np.degrees(initial_angle)),
        "minimum_angle_deg": float(np.degrees(minimum_angle)),
        "maximum_angle_deg": float(np.degrees(maximum_angle)),
        "relative_positive_range_deg": float(
            np.degrees(maximum_angle - initial_angle)
        ),
        "relative_negative_range_deg": float(
            np.degrees(initial_angle - minimum_angle)
        ),
        "csv": str(output_csv),
    }


def build_cases(
    channels: list[dict],
    sensitivities: dict[str, dict[str, float]],
) -> list[dict]:
    cases = []
    definitions = (
        (
            "elbow_flexion",
            ("elbow_flexors",),
            90.0,
            60.0,
        ),
        (
            "shoulder_flexion",
            ("deltoid_anterior",),
            90.0,
            60.0,
        ),
        (
            "shoulder_abduction",
            ("deltoid_middle",),
            90.0,
            60.0,
        ),
        (
            "shoulder_rotation",
            ("infraspinatus",),
            45.0,
            30.0,
        ),
    )
    for joint_name, ids, no_gravity_target, gravity_target in definitions:
        selected = select_channels(channels, sensitivities, joint_name, ids)
        cases.append(
            {
                "joint": joint_name,
                "selected": selected,
                "no_gravity_target_deg": no_gravity_target,
                "gravity_target_deg": gravity_target,
            }
        )
    return cases


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_OUTPUT_XML)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_REPORT_DIR / "cable_drive_motion_test.json",
    )
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument(
        "--ramp-seconds",
        type=float,
        default=6.0,
        help=(
            "Time allowed for the slider command to ramp up to its full "
            "travel.  Scaling a bone up lengthens and thickens the limb, so a "
            "slow ramp tests the reachable range instead of the joint speed."
        ),
    )
    parser.add_argument(
        "--hold-seconds",
        type=float,
        default=1.0,
        help="Time held at full travel before the final angle is reported.",
    )
    args = parser.parse_args()

    manifest = load_cable_manifest(args.manifest)
    channels = cable_channels(manifest)
    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    for channel in channels:
        actuator = actuator_id(model, f"pull_{channel['id']}")
        channel["travel_m"] = float(model.actuator_ctrlrange[actuator, 1])
    sensitivities = sensitivity_matrix(model, channels)
    cases = build_cases(channels, sensitivities)

    results = []
    failures = []
    for case in cases:
        if not case["selected"]:
            failures.append(
                {
                    "case": case["joint"],
                    "reason": "No cable in the intended muscle group pulls in "
                    "the required direction.",
                }
            )
            continue
        for gravity, target, label in (
            ((0.0, 0.0, 0.0), case["no_gravity_target_deg"], "no_gravity"),
            ((0.0, 0.0, -9.81), case["gravity_target_deg"], "gravity"),
        ):
            result = run_motion(
                model,
                data,
                channels,
                case["selected"],
                case["joint"],
                gravity,
                args.csv_dir / f"{case['joint']}_{label}.csv",
                float(manifest.get("defaults", {}).get("release_m", 0.08)),
                args.ramp_seconds,
                args.hold_seconds,
            )
            result["condition"] = label
            result["target_deg"] = target
            result["range_deg"] = max(
                result["relative_positive_range_deg"],
                result["relative_negative_range_deg"],
            )
            result["passed"] = result["range_deg"] >= target
            results.append(result)
            if not result["passed"]:
                failures.append(
                    {
                        "case": case["joint"],
                        "condition": label,
                        "reason": (
                            f"reached {result['range_deg']:.1f} "
                            f"deg, target {target:.1f} deg"
                        ),
                    }
                )

    report = {
        "model": str(args.model),
        "results": results,
        "failures": failures,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    for result in results:
        print(
            f"{result['joint']:20s} {result['condition']:10s} "
            f"range={result['range_deg']:7.2f} deg "
            f"target={result['target_deg']:.1f} "
            f"{'OK' if result['passed'] else 'FAIL'}"
        )
    for failure in failures:
        print(f"[FAIL] {failure}")
    print(f"report: {args.report}")

    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
