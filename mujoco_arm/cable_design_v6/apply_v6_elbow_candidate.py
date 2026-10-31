"""Install the compact H-F elbow candidate selected by route search."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "cable_design_v6.json"
E = np.array([184.34, -101.35, 7.30])
GH = np.array([159.06, 187.54, 7.30])
AXIS = np.array([0.9962, 0.0872, 0.0])
SIDE = np.array([-1.0, 0.0, 0.0])
VALUES = np.array(
    [
        184.52162890480844,
        36.449309498016255,
        -49.72240477746492,
        187.785170485292,
        -139.8606889870204,
        -59.91420483509829,
        203.1056489646603,
        -197.88834130676074,
        -53.657688598519755,
    ]
)


def local_mj(world_mm: np.ndarray, origin_mm: np.ndarray) -> list[float]:
    delta = world_mm - origin_mm
    return [
        float(delta[0] * 0.001),
        float(-delta[2] * 0.001),
        float(delta[1] * 0.001),
    ]


def route(values: np.ndarray, mirrored: bool) -> dict:
    values = np.asarray(values, dtype=float).copy()
    a = values[0:3].copy()
    g = values[3:6].copy()
    b = values[6:9].copy()
    if mirrored:
        for point in (a, g, b):
            point[2] = 2.0 * E[2] - point[2]
    route_id = "E02" if mirrored else "E01"
    return {
        "id": route_id,
        "group": "elbow",
        "role": "elbow_extension" if mirrored else "elbow_flexion",
        "diameter_mm": 1.0,
        "max_tension_N": 150.0,
        "max_guides": 1,
        "path": [
            {
                "kind": "anchor",
                "owner": "H",
                "body": "humerus",
                "world_mm": a.tolist(),
                "local_mm": local_mj(a, GH),
                "installation_normal": SIDE.tolist(),
            },
            {
                "kind": "guide",
                "name": f"guide_{route_id}_v6",
                "owner": "F",
                "body": "ulna",
                "world_mm": g.tolist(),
                "local_mm": local_mj(g, E),
                "radius_mm": 10.0,
                "axis_design": AXIS.tolist(),
                "side_design": SIDE.tolist(),
                "wrap_angle_deg": 90.0,
                "owner_part": "ulna",
            },
            {
                "kind": "anchor",
                "owner": "F",
                "body": "ulna",
                "world_mm": b.tolist(),
                "local_mm": local_mj(b, E),
                "installation_normal": SIDE.tolist(),
            },
        ],
        "tpu": {
            "enabled": True,
            "reference_strain": 0.1,
            "reference_force_N": 20.0,
            "compression_force_N": 0.0,
        },
        "status": "v6_compact_antagonist_candidate",
        "notes": [
            "Selected for trajectory controllability; "
            "minimum effective moment arm is below the 8 mm target."
        ],
    }


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    routes = [
        item for item in config["routes"] if item["id"] not in ("E01", "E02")
    ]
    config["routes"] = [route(VALUES, False), route(VALUES, True), *routes]
    CONFIG.write_text(
        json.dumps(config, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"updated elbow routes in {CONFIG}")


if __name__ == "__main__":
    main()
