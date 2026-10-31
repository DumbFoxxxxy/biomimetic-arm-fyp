"""Automated checks for the interactive tendon demo model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_MODEL = ROOT / "arm_tendon_demo.xml"
DEFAULT_CONFIG = ROOT / "config" / "tendon_demo_16.json"
DEFAULT_REPORT = ROOT / "exports" / "tendon_demo" / "check_report.json"

JOINT_NAMES = (
    "sc_elevation",
    "sc_protraction",
    "ac_ext",
    "gh_flex",
    "gh_abd",
    "gh_rot",
    "elbow",
)


def actuator_map(model: mujoco.MjModel) -> dict[str, int]:
    result = {}
    for index in range(model.nu):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, index) or ""
        if name.startswith("pull_"):
            result[name[5:]] = index
    return result


def neutral_controls(model, data, neutral, actuators) -> None:
    data.ctrl[:] = 0.0
    for cable_id, actuator in actuators.items():
        data.ctrl[actuator] = neutral[cable_id]


def joint_values(model, data) -> dict[str, float]:
    values = {}
    for name in JOINT_NAMES:
        joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
        address = int(model.jnt_qposadr[joint])
        value = float(data.qpos[address])
        values[name] = float(np.degrees(value))
    return values


def simulate(model, neutral, actuators, pulls, seconds=5.0):
    data = mujoco.MjData(model)
    data.eq_active[:] = 0
    neutral_controls(model, data, neutral, actuators)
    for cable_id, pull_mm in pulls.items():
        data.ctrl[actuators[cable_id]] = neutral[cable_id] - pull_mm * 0.001
    steps = int(round(seconds / model.opt.timestep))
    for _ in range(steps):
        mujoco.mj_step(model, data)
    return data, joint_values(model, data)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    model = mujoco.MjModel.from_xml_path(str(args.model))
    config = json.loads(args.config.read_text(encoding="utf-8"))
    neutral = {
        route["id"]: float(route["neutral_length_m"])
        for route in config["routes"]
    }
    actuators = actuator_map(model)

    rest_data = mujoco.MjData(model)
    rest_data.eq_active[:] = 0
    neutral_controls(model, rest_data, neutral, actuators)
    for _ in range(1000):
        mujoco.mj_step(model, rest_data)
    rest = joint_values(model, rest_data)

    tests = {
        "K1_E01": {"E01": 40.0},
        "K2_G01": {"G01": 30.0},
        "K3_G03": {"G03": 40.0},
        "K4_G05": {"G05": 30.0},
        "K5_G06": {"G06": 30.0},
        "K6_S01": {"S01": 40.0},
        "K7_G01_S01": {"G01": 30.0, "S01": 40.0},
        "SLIDE_S03": {"S03": 15.0},
        "SCAP_EXT_S05": {"S05": 20.0},
        "PECTORALIS_P01": {"P01": 30.0},
        "PECTORALIS_P02": {"P02": 30.0},
        "PECTORALIS_BOTH": {"P01": 30.0, "P02": 40.0},
    }

    print("model:", args.model)
    print("nq:", model.nq, "nu:", model.nu, "ntendon:", model.ntendon)
    print("rest:", " ".join(f"{key}={value:+.3f}" for key, value in rest.items()))

    failures = []
    results = {}
    for name, pulls in tests.items():
        data, values = simulate(model, neutral, actuators, pulls)
        finite = bool(np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all())
        delta = {key: values[key] - rest[key] for key in JOINT_NAMES}
        results[name] = {
            "pulls_mm": pulls,
            "pose": values,
            "delta": delta,
            "finite": finite,
            "min_contact_mm": float(
                min((contact.dist for contact in data.contact), default=0.0)
                * 1000.0
            ),
        }
        print(f"\n{name}: pulls={pulls}")
        print("  pose:", " ".join(f"{key}={value:+.2f}" for key, value in values.items()))
        print("  delta:", " ".join(f"{key}={value:+.2f}" for key, value in delta.items()))
        if not finite:
            failures.append(f"{name}: non-finite state")
        if max(abs(value) for value in delta.values()) < 0.5:
            failures.append(f"{name}: no visible motion")

    single_cable = {}
    for cable_id in actuators:
        data, values = simulate(model, neutral, actuators, {cable_id: 30.0})
        delta = {key: values[key] - rest[key] for key in JOINT_NAMES}
        dominant = max(delta, key=lambda key: abs(delta[key]))
        single_cable[cable_id] = {
            "dominant_joint": dominant,
            "delta": delta,
        }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(
            {
                "model": str(args.model),
                "rest_pose": rest,
                "tests": results,
                "single_cable_30mm": single_cable,
                "failures": failures,
                "pass": not failures,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"\nreport: {args.report}")

    if failures:
        print("\nFAIL")
        for failure in failures:
            print(" ", failure)
        raise SystemExit(1)
    print("\nPASS: all demo groups produce finite, visible motion.")


if __name__ == "__main__":
    main()
