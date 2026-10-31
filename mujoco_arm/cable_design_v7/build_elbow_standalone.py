"""Write the isolated R18 V7 elbow model."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import json
from project_runtime import WorkspacePath as Path
import sys


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v7"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from elbow_spool_v7 import DEFAULT_PARAMETERS, build_elbow_model  # noqa: E402


MODEL = ROOT / "v7_elbow_standalone.xml"
PARAMETERS = ROOT / "elbow_parameters_v7.json"


def main() -> None:
    parameters = dict(DEFAULT_PARAMETERS)
    parameters["radius_mm"] = 18.0
    tree = build_elbow_model(parameters)
    tree.write(MODEL, encoding="utf-8", xml_declaration=True)
    PARAMETERS.write_text(
        json.dumps(parameters, indent=2) + "\n",
        encoding="utf-8",
    )
    EXPORT.mkdir(parents=True, exist_ok=True)
    print(f"wrote {MODEL}")
    print(f"wrote {PARAMETERS}")


if __name__ == "__main__":
    main()
