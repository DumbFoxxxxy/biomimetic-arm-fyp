"""Trace one cable pull and report where every skeleton joint goes.

The motion test only reports the joint a cable is supposed to drive.  When a
pull produces almost no motion in that joint the reason is usually that the
wire is loading a different joint of the chain, so this tool logs every
skeleton angle together with the slider command and the tendon tension.

    D:\\MujocoEnv\\Scripts\\python.exe diagnose_pull.py --cable deltoid_anterior
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mujoco
import numpy as np

from cable_drive_common import (
    DEFAULT_MANIFEST,
    DEFAULT_OUTPUT_XML,
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
    parser.add_argument("--cable", type=str, required=True)
    parser.add_argument(
        "--co-cable",
        action="append",
        default=[],
        metavar="NAME=FRACTION",
        help=(
            "Also pull a channel, as a fraction of its travel.  Real motions "
            "are produced by groups of muscles, so a single cable rarely acts "
            "alone.  May be repeated."
        ),
    )
    parser.add_argument("--travel-fraction", type=float, default=1.0)
    parser.add_argument(
        "--pretension-fraction",
        type=float,
        default=0.0,
        help=(
            "Command given to every cable that is not under test, as a "
            "fraction of its travel.  A real arm always has resting tension in "
            "the antagonist muscles; releasing them completely leaves the "
            "humeral twist unopposed."
        ),
    )
    parser.add_argument("--seconds", type=float, default=6.0)
    parser.add_argument("--gravity", action="store_true")
    parser.add_argument("--rows", type=int, default=12)
    args = parser.parse_args()

    manifest = load_cable_manifest(args.manifest)
    channels = cable_channels(manifest)
    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    if not args.gravity:
        model.opt.gravity[:] = 0.0

    def actuator_id(name: str) -> int:
        value = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, name)
        if value < 0:
            raise ValueError(f"Unknown actuator: {name}")
        return int(value)

    selected = [
        channel for channel in channels if channel["id"].startswith(args.cable)
    ]
    if not selected:
        raise SystemExit(f"No channel matches {args.cable!r}")
    print("channels: " + ", ".join(channel["id"] for channel in selected))

    pulls = {
        channel["id"]: (
            actuator_id(f"pull_{channel['id']}"),
            float(model.actuator_ctrlrange[actuator_id(f"pull_{channel['id']}"), 1])
            * args.travel_fraction,
        )
        for channel in selected
    }
    for item in args.co_cable:
        name, _, fraction = item.partition("=")
        matching = [
            channel for channel in channels if channel["id"].startswith(name)
        ]
        if not matching:
            raise SystemExit(f"No channel matches co-cable {name!r}")
        for channel in matching:
            actuator = actuator_id(f"pull_{channel['id']}")
            pulls[channel["id"]] = (
                actuator,
                float(model.actuator_ctrlrange[actuator, 1])
                * float(fraction or 1.0),
            )
    idle = [
        actuator_id(f"pull_{channel['id']}")
        for channel in channels
        if channel["id"] not in pulls
    ]
    pretension = {
        actuator: args.pretension_fraction
        * float(model.actuator_ctrlrange[actuator, 1])
        for actuator in idle
    }
    assist_pairs = cable_assist_pairs(model)
    addresses = {
        name: int(
            model.jnt_qposadr[
                mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
            ]
        )
        for name in SKELETON_JOINT_NAMES
    }

    steps = int(args.seconds / model.opt.timestep)
    every = max(1, steps // args.rows)
    print(
        f"{'t s':>6s} " + " ".join(f"{name[:11]:>11s}" for name in SKELETON_JOINT_NAMES)
    )
    for step in range(steps + 1):
        if step:
            for actuator, travel in pulls.values():
                data.ctrl[actuator] = travel * min(
                    1.0, step * model.opt.timestep / args.seconds
                )
            for actuator, value in pretension.items():
                data.ctrl[actuator] = value
            apply_one_sided_tension(model, data, assist_pairs)
            apply_shoulder_rhythm_from_state(model, data)
            mujoco.mj_step(model, data)
        else:
            mujoco.mj_forward(model, data)
        if step % every == 0:
            angles = " ".join(
                f"{np.degrees(data.qpos[addresses[name]]):11.2f}"
                for name in SKELETON_JOINT_NAMES
            )
            print(f"{step * model.opt.timestep:6.2f} {angles}")

    print("\nfinal slider and tension")
    for channel_id, (actuator, travel) in pulls.items():
        tendon_id = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_TENDON, f"cable_{channel_id}"
        )
        print(
            f"  {channel_id:22s} command={data.ctrl[actuator]:.4f} m  "
            f"tendon={data.ten_length[tendon_id]:.4f} m  "
            f"tension={data.actuator_force[actuator]:.1f} N"
        )


if __name__ == "__main__":
    main()
