"""Run the V7 subsystem gates in dependency order."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import subprocess
import sys
from project_runtime import WorkspacePath as Path


ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable


def run(script: str) -> None:
    subprocess.run(
        [PYTHON, str(ROOT / script)],
        check=True,
        cwd=ROOT,
    )


def main() -> None:
    run("build_elbow_standalone.py")
    run("audit_elbow_v7.py")
    run("build_shoulder_standalone.py")
    run("audit_shoulder_v7.py")
    run("integration_v7.py")


if __name__ == "__main__":
    main()
