"""Write the initial V7 isolated shoulder model."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v7"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from shoulder_roller_v7 import build_shoulder_model  # noqa: E402


MODEL = ROOT / "v7_shoulder_standalone.xml"


def main() -> None:
    tree = build_shoulder_model(1)
    tree.write(MODEL, encoding="utf-8", xml_declaration=True)
    EXPORT.mkdir(parents=True, exist_ok=True)
    print(f"wrote {MODEL}")


if __name__ == "__main__":
    main()
