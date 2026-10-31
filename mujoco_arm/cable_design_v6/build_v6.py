"""Build the isolated V6 MuJoCo model with shoulder roller markers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from build_v2_model import build as build_base_model  # noqa: E402
from build_v2_model import find_body  # noqa: E402
from v2_common import POINTS_MM, load_json, normalize, vec  # noqa: E402


CONFIG = ROOT / "cable_design_v6.json"
MODEL = ROOT / "arm_cable_design_v6.xml"
EXPORT = ROOT.parent / "exports" / "cable_design_v6"


def design_world_to_body_mj(
    world_mm: np.ndarray,
    origin_mm: np.ndarray,
) -> np.ndarray:
    delta = np.asarray(world_mm, dtype=float) - np.asarray(origin_mm, dtype=float)
    return np.array([delta[0], -delta[2], delta[1]], dtype=float) * 0.001


def add_v6_markers_and_roller(tree: ET.ElementTree, config: dict) -> None:
    root = tree.getroot()
    asset = root.find("asset")
    if asset is None:
        raise RuntimeError("missing asset section")
    ET.SubElement(
        asset,
        "material",
        {"name": "v6_roller_mat", "rgba": "0.18 0.65 0.82 1"},
    )
    ulna = find_body(root, "ulna")
    ET.SubElement(
        ulna,
        "site",
        {
            "name": "V6_E",
            "pos": "0 0 0",
            "size": "0.004",
            "rgba": "0.15 0.9 0.45 1",
            "group": "3",
        },
    )
    ET.SubElement(
        ulna,
        "site",
        {
            "name": "V6_W",
            "pos": "0 0 -0.237",
            "size": "0.005",
            "rgba": "0.15 0.55 1 1",
            "group": "3",
        },
    )
    scapula = find_body(root, "scapula")
    outward = normalize(np.array([0.6556, 0.0, 0.7552]))
    guide = config["shoulder_guide"]
    radius_mm = float(guide["roller_radius_mm"])
    clearance_mm = 0.35
    centre_world = POINTS_MM["P"] + outward * (radius_mm + clearance_mm)
    local = design_world_to_body_mj(centre_world, POINTS_MM["AC"])
    ET.SubElement(
        scapula,
        "geom",
        {
            "name": "v6_shoulder_roller",
            "type": "cylinder",
            "size": f"{radius_mm * 0.001:.9f} {float(guide['roller_width_mm']) * 0.0005:.9f}",
            "pos": vec(local),
            "material": "v6_roller_mat",
            "contype": "2",
            "conaffinity": "1",
            "condim": "3",
            "friction": "0.10 0.005 0.0001",
            "solref": "0.002 1",
            "solimp": "0.9 0.95 0.001",
            "density": "1400",
            "group": "1",
        },
    )
    ET.SubElement(
        scapula,
        "site",
        {
            "name": "V6_roller",
            "pos": vec(local),
            "size": "0.003",
            "rgba": "0.1 0.85 1 0.9",
            "group": "3",
        },
    )
    actuator_root = root.find("actuator")
    if actuator_root is None:
        actuator_root = ET.SubElement(root, "actuator")
    for child in list(actuator_root):
        actuator_root.remove(child)
    for route in config["routes"]:
        if config["simulation"].get("actuator_type") == "position":
            ET.SubElement(
                actuator_root,
                "position",
                {
                    "name": f"pull_{route['id']}",
                    "tendon": f"cable_{route['id']}",
                    "kp": str(config["simulation"].get("kp", 1200.0)),
                    "kv": str(config["simulation"].get("kv", 80.0)),
                    "ctrllimited": "true",
                    "ctrlrange": "0 1",
                    "forcelimited": "true",
                    "forcerange": "-160 0",
                    "group": "2",
                },
            )
        else:
            ET.SubElement(
                actuator_root,
                "motor",
                {
                    "name": f"pull_{route['id']}",
                    "tendon": f"cable_{route['id']}",
                    "ctrllimited": "true",
                    "ctrlrange": "0 150",
                    "forcelimited": "true",
                    "forcerange": "0 150",
                    "gear": "1",
                    "group": "2",
                },
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--model", type=Path, default=MODEL)
    args = parser.parse_args()
    config = load_json(args.config)
    tree, metadata = build_base_model(config)
    add_v6_markers_and_roller(tree, config)
    ET.indent(tree, space="  ")
    args.model.parent.mkdir(parents=True, exist_ok=True)
    tree.write(args.model, encoding="utf-8", xml_declaration=True)
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "route_metadata_v6.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.model}")


if __name__ == "__main__":
    main()
