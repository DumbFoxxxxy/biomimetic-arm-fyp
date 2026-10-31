"""Run the arm with simplified scapulohumeral rhythm coupling."""

from pathlib import Path
import time

import mujoco
import mujoco.viewer

from shoulder_coupling import apply_shoulder_rhythm


MODEL_PATH = Path(__file__).with_name("arm_reference.xml")


def main() -> None:
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
    data = mujoco.MjData(model)

    with mujoco.viewer.launch_passive(model, data) as viewer:
        while viewer.is_running():
            # The shoulder sliders remain user-controlled. The clavicle and
            # scapula actuator values are overwritten by the rhythm model.
            apply_shoulder_rhythm(model, data)
            mujoco.mj_step(model, data)
            viewer.sync()
            time.sleep(model.opt.timestep)


if __name__ == "__main__":
    main()
