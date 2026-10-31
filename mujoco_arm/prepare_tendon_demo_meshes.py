"""Convert the current SolidWorks part STLs to MuJoCo body-local meshes."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from stl_io import read_binary_stl, write_binary_stl


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "assets_tendon_demo" / "raw"
OUT = ROOT / "assets_tendon_demo"
BODY_BOXES = OUT / "body_boxes.json"

GH_MM = np.array([159.06, 187.54, 7.30])
P_MM = np.array([62.18, 210.10, 91.39])
E_MM = np.array([184.34, -101.35, 7.30])

PARTS = {
    "thorax": {
        "body_origin_mm": np.array([0.0, 0.0, 0.0]),
    },
    "scapula": {
        "body_origin_mm": P_MM,
    },
    "clavicle": {
        "body_origin_mm": P_MM,
    },
    "humerus": {
        "body_origin_mm": GH_MM,
    },
    "ulna": {
        "body_origin_mm": E_MM,
    },
}


def to_mujoco_design(point_mm: np.ndarray) -> np.ndarray:
    return np.array([point_mm[0], -point_mm[2], point_mm[1]]) * 0.001


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    boxes = {
        row["name"]: row
        for row in json.loads(BODY_BOXES.read_text(encoding="utf-8-sig"))
    }
    report = {}
    for name, spec in PARTS.items():
        source = RAW / f"{name}.stl"
        if not source.exists():
            raise FileNotFoundError(source)
        triangles = read_binary_stl(source)
        stl_min_mm = triangles.reshape(-1, 3).min(axis=0)
        body_min_mm = np.asarray(boxes[name]["box_min_m"], dtype=float) * 1000.0
        offset_mm = body_min_mm - stl_min_mm
        vertices = triangles.reshape(-1, 3) + offset_mm
        mapped = np.column_stack(
            (vertices[:, 0], -vertices[:, 2], vertices[:, 1])
        ) * 0.001
        mapped -= to_mujoco_design(spec["body_origin_mm"])
        body_triangles = mapped.reshape(-1, 3, 3)
        destination = OUT / f"{name}.stl"
        write_binary_stl(
            destination,
            body_triangles,
            f"body-local MuJoCo mesh for {name}",
        )
        report[name] = {
            "offset_mm": offset_mm.tolist(),
            "body_origin_mm": spec["body_origin_mm"].tolist(),
            "bounds_m": [
                body_triangles.reshape(-1, 3).min(axis=0).tolist(),
                body_triangles.reshape(-1, 3).max(axis=0).tolist(),
            ],
        }
        print(f"wrote {destination}")
    (OUT / "mesh_report.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
