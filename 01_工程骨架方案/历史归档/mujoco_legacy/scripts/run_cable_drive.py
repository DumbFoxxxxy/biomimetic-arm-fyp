"""Open the cable-driven arm model with live slider controls."""

from __future__ import annotations

import argparse
from pathlib import Path
import time

import mujoco
import mujoco.viewer

from cable_drive_common import (
    DEFAULT_OUTPUT_XML,
    apply_one_sided_tension,
    cable_assist_pairs,
)
from shoulder_coupling import apply_shoulder_rhythm_from_state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_OUTPUT_XML)
    parser.add_argument("--no-rhythm", action="store_true")
    args = parser.parse_args()

    if not args.model.exists():
        raise FileNotFoundError(
            f"Missing {args.model}. Run build_cable_drive_model.py first."
        )

    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    assist_pairs = cable_assist_pairs(model)

    with mujoco.viewer.launch_passive(model, data) as viewer:
        while viewer.is_running():
            if not args.no_rhythm:
                apply_shoulder_rhythm_from_state(model, data)
            # Wires only pull: cancel the compressive half of every tendon
            # spring so a slack cable releases instead of pushing.
            apply_one_sided_tension(model, data, assist_pairs)
            mujoco.mj_step(model, data)
            viewer.sync()
            time.sleep(model.opt.timestep)


if __name__ == "__main__":
    main()
