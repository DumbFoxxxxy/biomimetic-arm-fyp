"""Run the V7 subsystem gates in dependency order."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


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
