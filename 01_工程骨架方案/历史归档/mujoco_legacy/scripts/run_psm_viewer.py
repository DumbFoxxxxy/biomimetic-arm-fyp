"""Open the generated flexible-muscle model in the MuJoCo viewer."""

import argparse
from pathlib import Path

import mujoco
import mujoco.viewer


MODEL_PATH = Path(__file__).with_name("arm_psm_generated.xml")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    args = parser.parse_args()

    if not args.model.exists():
        raise FileNotFoundError(
            "Run build_arm_model.py before run_psm_viewer.py"
        )
    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    mujoco.viewer.launch(model, data)


if __name__ == "__main__":
    main()
