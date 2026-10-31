"""Check that every slider can move the arm on its own.

The objective is a simulation that starts and that the arm can be driven by
moving sliders, so every ``pull_<cable>`` channel has to have authority: pull
it alone and at least one skeleton joint must move.  This does not judge the
range of the motion, only that the channel does something, in a stable way.

    D:\\MujocoEnv\\Scripts\\python.exe test_slider_authority.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mujoco
import numpy as np

from cable_drive_common import (
    DEFAULT_MANIFEST,
    DEFAULT_OUTPUT_XML,
    DEFAULT_REPORT_DIR,
    SKELETON_JOINT_NAMES,
    apply_one_sided_tension,
    cable_assist_pairs,
    load_cable_manifest,
)
from check_cable_drive import cable_channels
from shoulder_coupling import apply_shoulder_rhythm_from_state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_OUTPUT_XML)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_REPORT_DIR / "slider_authority.json",
    )
    parser.add_argument("--seconds", type=float, default=4.0)
    parser.add_argument(
        "--gravity",
        action="store_true",
        help="Keep gravity on; the default runs the sweep without it.",
    )
    parser.add_argument(
        "--min-deg",
        type=float,
        default=1.0,
        help="Motion a channel must produce somewhere to count as connected.",
    )
    args = parser.parse_args()

    manifest = load_cable_manifest(args.manifest)
    channels = cable_channels(manifest)
    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    if not args.gravity:
        model.opt.gravity[:] = 0.0

    addresses = {
        name: int(
            model.jnt_qposadr[
                mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
            ]
        )
        for name in SKELETON_JOINT_NAMES
    }
    actuators = {
        channel["id"]: int(
            mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"pull_{channel['id']}"
            )
        )
        for channel in channels
    }
    travel = {
        channel_id: float(model.actuator_ctrlrange[actuator, 1])
        for channel_id, actuator in actuators.items()
    }
    assist_pairs = cable_assist_pairs(model)
    steps = int(args.seconds / model.opt.timestep)

    results = []
    for channel in channels:
        channel_id = channel["id"]
        data.qpos[:] = 0.0
        data.qvel[:] = 0.0
        data.ctrl[:] = 0.0
        mujoco.mj_forward(model, data)
        start = {name: float(data.qpos[addresses[name]]) for name in addresses}
        unstable = False
        for step in range(steps):
            for other_id, actuator in actuators.items():
                data.ctrl[actuator] = (
                    travel[channel_id] * (step / steps)
                    if other_id == channel_id
                    else 0.0
                )
            apply_one_sided_tension(model, data, assist_pairs)
            apply_shoulder_rhythm_from_state(model, data)
            mujoco.mj_step(model, data)
            if not np.isfinite(data.qpos).all():
                unstable = True
                break
        moved = {
            name: float(np.degrees(data.qpos[addresses[name]] - start[name]))
            for name in addresses
        }
        dominant = max(moved, key=lambda name: abs(moved[name]))
        entry = {
            "channel": channel_id,
            "travel_m": travel[channel_id],
            "dominant_joint": dominant,
            "dominant_delta_deg": round(moved[dominant], 2),
            "max_abs_delta_deg": round(max(abs(v) for v in moved.values()), 2),
            "deltas_deg": {k: round(v, 2) for k, v in moved.items()},
            "unstable": unstable,
            "connected": (
                not unstable
                and max(abs(v) for v in moved.values()) >= args.min_deg
            ),
        }
        results.append(entry)
        flag = "ok  " if entry["connected"] else "DEAD"
        print(
            f"{flag} {channel_id:22s} {entry['max_abs_delta_deg']:7.2f} deg "
            f"on {dominant:24s} travel={travel[channel_id]:.3f} m"
        )

    dead = [entry["channel"] for entry in results if not entry["connected"]]
    print(f"\n{len(results) - len(dead)}/{len(results)} sliders move the arm")
    if dead:
        print("no authority: " + ", ".join(dead))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"report: {args.report}")


if __name__ == "__main__":
    main()
