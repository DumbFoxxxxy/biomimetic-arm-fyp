"""Write the isolated R18 V7 elbow model."""

from __future__ import annotations

import json
from pathlib import Path
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
