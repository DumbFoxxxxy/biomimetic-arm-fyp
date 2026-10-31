"""Contact-friction sensitivity under gravity-only settling."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import copy
import json
from project_runtime import WorkspacePath as Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from v2_common import CONFIG_PATH, EXPORT_ROOT, MODEL_PATH, load_json


def set_contact_friction(root: ET.Element, value: float) -> None:
    friction = f"{value:.4f} 0.005 0.0001"
    for geom in root.iter("geom"):
        if geom.get("name") in (
            "v2_thorax_contact_band",
            "v2_scapula_contact_pad",
        ):
            geom.set("friction", friction)


def settle(xml: str, steps: int) -> dict:
    model = _load_project_xml(xml)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    for actuator in range(model.nu):
        tendon = int(model.actuator_trnid[actuator, 0])
        data.ctrl[actuator] = float(data.ten_length[tendon])
    mujoco.mj_forward(model, data)
    for _ in range(steps):
        mujoco.mj_step(model, data)
    joints = {}
    for name in (
        "sc_elevation",
        "sc_protraction",
        "ac_ext",
        "gh_flex",
        "gh_abd",
        "gh_rot",
        "elbow",
    ):
        joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
        joints[name] = float(np.degrees(data.qpos[int(model.jnt_qposadr[joint])]))
    return {
        "joints_deg": joints,
        "contacts": int(data.ncon),
        "min_contact_mm": float(
            min((contact.dist for contact in data.contact), default=0.0) * 1000.0
        ),
        "finite": bool(np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=MODEL_PATH)
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument(
        "--out",
        type=Path,
        default=EXPORT_ROOT / "friction_sensitivity_v2.json",
    )
    args = parser.parse_args()
    config = load_json(args.config)
    base = ET.parse(args.model)
    values = [
        float(value) for value in config["contact"]["friction_sensitivity"]
    ]
    results = {}
    for friction in values:
        tree = copy.deepcopy(base)
        set_contact_friction(tree.getroot(), friction)
        xml = ET.tostring(tree.getroot(), encoding="unicode")
        results[f"{friction:.3f}"] = settle(
            xml,
            int(config["simulation"]["settle_steps"]),
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(
            {
                "friction_values": values,
                "results": results,
                "guide_friction_status": (
                    "requires segment-level capstan tension; not included "
                    "in the current single-tendon static solver"
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.out}")
    for friction, result in results.items():
        print(
            f"mu={friction} contacts={result['contacts']} "
            f"min={result['min_contact_mm']:.3f} mm "
            f"finite={result['finite']}"
        )


if __name__ == "__main__":
    main()
