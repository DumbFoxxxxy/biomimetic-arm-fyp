"""Deterministic coarse search for a valid R18 sector-spool cable layout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import mujoco


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v7"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from elbow_spool_v7 import build_elbow_model, evaluate_elbow  # noqa: E402


ANGLES = [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]


def second_angles(anchor: float, entry_side: float, sector_side: float, mode: int):
    if mode == 0:
        return 180.0 - anchor, 180.0 - entry_side, 180.0 - sector_side
    if mode == 1:
        return -anchor, -entry_side, -sector_side
    if mode == 2:
        return anchor + 180.0, entry_side + 180.0, sector_side + 180.0
    return anchor, entry_side, sector_side


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", type=float, default=18.0)
    parser.add_argument("--limit", type=int, default=40)
    args = parser.parse_args()
    candidates = []
    evaluations = 0
    for anchor in ANGLES:
        for sector_side in ANGLES:
            for entry_side in ANGLES:
                for mode in range(4):
                    anchor2, entry2, sector2 = second_angles(
                        anchor,
                        entry_side,
                        sector_side,
                        mode,
                    )
                    parameters = {
                        "radius_mm": args.radius,
                        "anchor_angle_E01_deg": anchor,
                        "anchor_angle_E02_deg": anchor2 % 360.0,
                        "entry_side_angle_E01_deg": entry_side,
                        "entry_side_angle_E02_deg": entry2 % 360.0,
                        "sector_side_angle_E01_deg": sector_side,
                        "sector_side_angle_E02_deg": sector2 % 360.0,
                    }
                    tree = build_elbow_model(parameters)
                    xml = ET.tostring(tree.getroot(), encoding="unicode")
                    try:
                        model = mujoco.MjModel.from_xml_string(xml)
                        result = evaluate_elbow(model, args.radius)
                    except Exception as error:  # noqa: BLE001
                        evaluations += 1
                        candidates.append(
                            {
                                "parameters": parameters,
                                "pass": False,
                                "error": str(error),
                            }
                        )
                        continue
                    evaluations += 1
                    score = (
                        (0.0 if result["opposite"] else 100.0)
                        + (0.0 if result["wrap_ok"] else 50.0)
                        + max(0.0, 8.0 - result["minimum_arm_mm"]) * 20.0
                        + max(
                            0.0,
                            result["contraction_percent"]["max"] - 18.0,
                        )
                        * 5.0
                    )
                    candidates.append(
                        {
                            "parameters": parameters,
                            "score": score,
                            "radius_mm": args.radius,
                            "opposite": result["opposite"],
                            "wrap_ok": result["wrap_ok"],
                            "minimum_arm_mm": result["minimum_arm_mm"],
                            "contraction_percent": result[
                                "contraction_percent"
                            ],
                            "pass": result["pass"],
                        }
                    )
    candidates.sort(
        key=lambda item: (
            not item.get("pass", False),
            float(item.get("score", 1e9)),
        )
    )
    output = {
        "radius_mm": args.radius,
        "evaluations": evaluations,
        "passing_count": sum(item.get("pass", False) for item in candidates),
        "best": candidates[: args.limit],
    }
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / f"elbow_spool_search_r{args.radius:g}.json").write_text(
        json.dumps(output, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "radius_mm": args.radius,
                "evaluations": evaluations,
                "passing_count": output["passing_count"],
                "best": candidates[0] if candidates else None,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
