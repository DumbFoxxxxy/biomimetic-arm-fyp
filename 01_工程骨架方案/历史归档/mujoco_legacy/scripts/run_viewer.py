"""Open the reference arm in the MuJoCo interactive viewer."""

from pathlib import Path

import mujoco
import mujoco.viewer


MODEL_PATH = Path(__file__).with_name("arm_reference.xml")


def main() -> None:
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
    data = mujoco.MjData(model)
    mujoco.viewer.launch(model, data)


if __name__ == "__main__":
    main()
