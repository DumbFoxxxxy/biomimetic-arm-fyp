"""Allow every installed channel as a candidate while preserving mandatory sets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "cable_design_v5.json"
MATRIX = ROOT / "action_cable_matrix_v5.json"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    matrix = json.loads(args.matrix.read_text(encoding="utf-8"))
    all_cables = [route["id"] for route in config["routes"]]
    preferred = {
        "K1": ["E01", "E02"],
        "K2": ["G01", "G02", "G03", "G04", "G05", "G06", "G07", "P01", "P02"],
        "K3": ["G01", "G02", "G03", "G04", "G05", "G06", "G07"],
        "K4": ["G01", "G02", "G03", "G04", "G05", "G06", "P01", "P02"],
        "K5": ["G01", "G02", "G03", "G04", "G05", "G06", "P01", "P02"],
        "K6": ["S01", "S02", "S03", "S04", "S05", "S06"],
        "K7": [
            "G01",
            "G02",
            "G03",
            "G07",
            "S01",
            "S02",
            "S03",
            "S04",
            "S05",
            "S06",
            "P01",
            "P02",
        ],
        "SLIDE_FWD": ["S01", "S02", "S03", "S04", "S05", "S06"],
        "SLIDE_BACK": ["S01", "S02", "S03", "S04", "S05", "S06"],
        "SCAP_EXT": ["S01", "S02", "S03", "S04", "S05", "S06"],
        "SCAP_INT": ["S01", "S02", "S03", "S04", "S05", "S06"],
        "PECTORAL": ["P01", "P02", "G05", "G06"],
    }
    for action in matrix["actions"]:
        action["allowed_cables"] = list(all_cables)
        action["preferred_cables"] = [
            route_id
            for route_id in preferred[action["id"]]
            if route_id in all_cables
        ]
    args.matrix.write_text(json.dumps(matrix, indent=2) + "\n", encoding="utf-8")
    print(f"updated {len(matrix['actions'])} actions with {len(all_cables)} channels")


if __name__ == "__main__":
    main()
