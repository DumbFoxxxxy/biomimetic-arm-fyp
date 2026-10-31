"""Apply V8 gates and execute historical regressions without modifying them."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import csv
import hashlib
import json
from project_runtime import WorkspacePath as Path
import subprocess
import sys

import mujoco
import numpy as np

from v8_shoulder import (
    CONTACT_PATH,
    EXPORT_ROOT,
    OPEN_CHAIN_PATH,
    PROJECT_ROOT,
    build_shoulder_tree,
    compile_tree,
    load_config,
)


ROOT = Path(__file__).resolve().parent
DEMO_CHECK = ROOT.parent / "check_tendon_demo.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    config = load_config()
    kinematics_path = EXPORT_ROOT / "kinematics_report_v8.json"
    contact_path = EXPORT_ROOT / "contact_report_v8.json"
    kinematics = json.loads(kinematics_path.read_text(encoding="utf-8"))
    contact = json.loads(contact_path.read_text(encoding="utf-8"))
    model_checks = {}
    model_loads_pass = True
    for name, path, count in (
        ("open_chain", OPEN_CHAIN_PATH, 0),
        ("contact", CONTACT_PATH, 1),
    ):
        try:
            model = compile_tree(build_shoulder_tree(config, roller_count=count))
            data = mujoco.MjData(model)
            mujoco.mj_forward(model, data)
            finite = bool(
                np.all(np.isfinite(data.qpos))
                and np.all(np.isfinite(data.qvel))
            )
            model_checks[name] = {
                "path": str(path),
                "loads": True,
                "njoint": int(model.njnt),
                "nq": int(model.nq),
                "nv": int(model.nv),
                "finite": finite,
            }
            model_loads_pass = model_loads_pass and finite
        except Exception as exc:  # pragma: no cover - audit failure path
            model_checks[name] = {
                "path": str(path),
                "loads": False,
                "error": repr(exc),
            }
            model_loads_pass = False
    v7_elbow_path = PROJECT_ROOT / config["regression"]["v7_elbow_report"]
    v7_elbow = json.loads(v7_elbow_path.read_text(encoding="utf-8"))
    demo = subprocess.run(
        [sys.executable, str(DEMO_CHECK)],
        cwd=DEMO_CHECK.parent,
        capture_output=True,
        text=True,
        check=False,
    )
    demo_pass = demo.returncode == 0 and "PASS" in demo.stdout
    gate_pass = bool(
        model_loads_pass
        and
        kinematics["pass"]
        and contact["pass"]
        and v7_elbow["pass"]
        and demo_pass
    )
    failures = []
    if not model_loads_pass:
        failures.append("model_load")
    if not kinematics["pass"]:
        failures.extend(
            f"kinematics:{reason}" for reason in kinematics["failure_reasons"]
        )
    if not contact["pass"]:
        failures.extend(
            f"contact:{reason}" for reason in contact["failure_reasons"]
        )
    if not v7_elbow["pass"]:
        failures.append("v7_elbow_regression")
    if not demo_pass:
        failures.append("sixteen_cable_demo_regression")
    inherited_files = [
        ROOT.parent / "cable_design_v5" / "arm_cable_design_v5.xml",
        ROOT.parent / "cable_design_v5" / "cable_design_v5.json",
        ROOT.parent / "cable_design_v6" / "arm_cable_design_v6.xml",
        ROOT.parent / "cable_design_v6" / "cable_design_v6.json",
        ROOT.parent / "cable_design_v7" / "v7_elbow_standalone.xml",
        ROOT.parent / "cable_design_v7" / "v7_shoulder_standalone.xml",
    ]
    report = {
        "version": 8,
        "stage": "integration",
        "model_checks": model_checks,
        "model_loads_pass": model_loads_pass,
        "kinematics_gate_pass": bool(kinematics["pass"]),
        "contact_gate_pass": bool(contact["pass"]),
        "v7_elbow_regression_pass": bool(v7_elbow["pass"]),
        "sixteen_cable_demo_pass": bool(demo_pass),
        "demo_stdout": demo.stdout.strip(),
        "demo_stderr": demo.stderr.strip(),
        "gate_pass": gate_pass,
        "status": "ready_for_integration" if gate_pass else "blocked",
        "frozen_hardpoints": bool(gate_pass),
        "integrated_model_written": False,
        "failure_reasons": failures,
        "inherited_file_sha256": {
            str(path.relative_to(PROJECT_ROOT)): sha256(path)
            for path in inherited_files
            if path.exists()
        },
    }
    if gate_pass:
        hardpoints = EXPORT_ROOT / "shoulder_hardpoints_v8.csv"
        with hardpoints.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "name",
                    "part",
                    "world_x_mm",
                    "world_y_mm",
                    "world_z_mm",
                    "description",
                ]
            )
            values = config["geometry_mm"]
            for name in ("SC", "AC", "P", "SUP", "INF"):
                writer.writerow(
                    [
                        name,
                        "thorax" if name == "P" else "scapula",
                        *values[name],
                        "V8 reference geometry",
                    ]
                )
    EXPORT_ROOT.mkdir(parents=True, exist_ok=True)
    (EXPORT_ROOT / "integration_report_v8.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "gate_pass": gate_pass,
                "v7_elbow_regression_pass": report[
                    "v7_elbow_regression_pass"
                ],
                "sixteen_cable_demo_pass": report[
                    "sixteen_cable_demo_pass"
                ],
                "failure_reasons": failures,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

