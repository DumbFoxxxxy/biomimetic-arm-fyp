"""Create the cable-drive manifest from the enabled PSM muscle pool."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from cable_drive_common import ROOT, load_json, vec3


DEFAULT_SOURCE = ROOT / "config" / "muscle_manifest.json"
DEFAULT_OUTPUT = ROOT / "config" / "cable_drive_manifest.json"
DEFAULT_REFERENCE_XML = ROOT / "arm_reference.xml"

GUIDE_SPECS = {
    "supraspinatus": (
        "torso",
        "shoulder_abduction",
        -0.060,
        0.025,
        1.0,
    ),
    "infraspinatus": (
        "torso",
        "shoulder_rotation",
        -0.060,
        0.025,
        1.0,
    ),
    "deltoid_anterior": (
        "torso",
        "shoulder_flexion",
        -0.060,
        0.030,
        1.0,
    ),
    "deltoid_middle": (
        "torso",
        "shoulder_abduction",
        -0.060,
        0.030,
        1.0,
    ),
    "deltoid_posterior": (
        "torso",
        "shoulder_flexion",
        -0.060,
        0.030,
        -1.0,
    ),
    "elbow_flexors": (
        "shoulder_abduction",
        "elbow_flexion",
        -0.080,
        0.020,
        1.0,
    ),
    "triceps": (
        "shoulder_abduction",
        "elbow_flexion",
        -0.080,
        0.020,
        -1.0,
    ),
}


def body_local_positions(reference_xml: Path) -> dict[str, np.ndarray]:
    tree = ET.parse(reference_xml)
    result: dict[str, np.ndarray] = {}
    for body in tree.getroot().iter("body"):
        name = body.get("name")
        if not name:
            continue
        result[name] = vec3(
            [float(value) for value in body.get("pos", "0 0 0").split()],
            f"{name}.pos",
        )
    return result


def convert_muscle(
    muscle: dict,
    body_positions: dict[str, np.ndarray],
) -> dict:
    insertion_body = muscle["insertion_body"]
    insertion_local = vec3(muscle["insertion_local"], "insertion_local")

    if insertion_body == "radius":
        radius_offset = body_positions["radius"]
        insertion_body = "forearm_ulna"
        insertion_local = radius_offset + insertion_local

    guide_body, target_joint, guide_offset, _moment_arm, moment_sign = (
        GUIDE_SPECS[muscle["id"]]
    )
    return {
        "id": muscle["id"],
        "group": muscle.get("group", muscle["id"]),
        "enabled": True,
        "insertion_body": insertion_body,
        "insertion_local_m": insertion_local.tolist(),
        "via_body": muscle["origin_body"],
        "via_local_m": vec3(muscle["origin_local"], "origin_local").tolist(),
        "guides": [
            {
                "body": muscle["origin_body"],
                "local_m": vec3(
                    muscle["origin_local"], "origin_local"
                ).tolist(),
                "radius_m": 0.006,
            }
        ],
        "guide_mode": "constant_moment_arm",
        "guide_body": guide_body,
        "target_joint": target_joint,
        "guide_offset_m": abs(guide_offset),
        "guide_moment_sign": moment_sign,
        "channels": int(muscle.get("cables", 1)),
        "channel_offset_m": 0.002,
        "max_force_n": float(muscle["max_force_n"]),
        "shape_class": muscle.get("shape_class", "straight_diamond"),
    }


def build_manifest(source: Path, reference_xml: Path) -> dict:
    source_payload = load_json(source)
    body_positions = body_local_positions(reference_xml)
    enabled = [
        muscle
        for muscle in source_payload["muscles"]
        if muscle.get("enabled", False)
    ]
    cables = [convert_muscle(muscle, body_positions) for muscle in enabled]
    return {
        "schema_version": 1,
        "model_name": "right_arm_cable_drive_v1",
        "description": (
            "World-fixed cable sliders pull the humerus and ulna. "
            "Only one via point is used per cable channel."
        ),
        "units": "meters",
        "defaults": {
            "travel_m": 0.08,
            "rated_stretch_m": 0.02,
            "pretension_n": 0.2,
            "spring_safety_margin_m": 0.01,
            "damping_ratio": 0.7,
            "slider_offset_m": 0.12,
            "slider_radius_m": 0.008,
            "slider_kp": 50000.0,
            "slider_kv": 500.0,
            "mesh_units": "meters",
            "cable_width_m": 0.0015,
        },
        "cables": cables,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--reference-xml", type=Path, default=DEFAULT_REFERENCE_XML)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    payload = build_manifest(args.source, args.reference_xml)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    total_channels = sum(int(cable["channels"]) for cable in payload["cables"])
    print(f"wrote: {args.output}")
    print(f"muscle groups: {len(payload['cables'])}")
    print(f"cable channels: {total_channels}")


if __name__ == "__main__":
    main()
