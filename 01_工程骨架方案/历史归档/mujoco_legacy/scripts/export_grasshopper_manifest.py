"""Export the enabled muscle interfaces for Rhino/Grasshopper."""

import argparse
import csv
from pathlib import Path
import json

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
ARM_CONFIG = ROOT / "config" / "arm_model.json"
MUSCLE_CONFIG = ROOT / "config" / "muscle_manifest.json"
OUTPUT_JSON = ROOT / "exports" / "grasshopper_muscles.json"
OUTPUT_CSV = ROOT / "exports" / "grasshopper_muscles.csv"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm-config", type=Path, default=ARM_CONFIG)
    parser.add_argument(
        "--muscle-config", type=Path, default=MUSCLE_CONFIG
    )
    parser.add_argument("--output-json", type=Path, default=OUTPUT_JSON)
    parser.add_argument("--output-csv", type=Path, default=OUTPUT_CSV)
    args = parser.parse_args()

    arm = json.loads(args.arm_config.read_text(encoding="utf-8"))
    manifest = json.loads(
        args.muscle_config.read_text(encoding="utf-8")
    )
    base_path = ROOT / arm["base_xml"]

    model = mujoco.MjModel.from_xml_path(str(base_path))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)

    def body_transform(name: str) -> tuple[np.ndarray, np.ndarray]:
        body_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, name)
        if body_id < 0:
            raise ValueError(f"Unknown body: {name}")
        return (
            np.array(data.xpos[body_id]),
            np.array(data.xmat[body_id]).reshape(3, 3),
        )

    def to_world(body: str, local) -> list[float]:
        position, rotation = body_transform(body)
        return (position + rotation @ np.asarray(local, dtype=float)).tolist()

    output = []
    for muscle in manifest["muscles"]:
        if not muscle.get("enabled", False):
            continue
        item = dict(muscle)
        item["origin_world"] = to_world(
            muscle["origin_body"], muscle["origin_local"]
        )
        item["insertion_world"] = to_world(
            muscle["insertion_body"], muscle["insertion_local"]
        )
        item["via_world"] = [
            {
                "body": via["body"],
                "position": to_world(via["body"], via["local"]),
            }
            for via in muscle.get("via", [])
        ]
        item["psm_parameters"] = {
            "x_mm": manifest["defaults"]["lattice_x_mm"],
            "y_mm": manifest["defaults"]["lattice_y_mm"],
            "z_mm": manifest["defaults"]["lattice_z_mm"],
            "p_mm": manifest["defaults"]["branch_diameter_mm"],
        }
        output.append(item)

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(output, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    with args.output_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "id",
                "group",
                "shape_class",
                "cables",
                "origin_body",
                "insertion_body",
                "origin_world",
                "insertion_world",
                "x_mm",
                "y_mm",
                "z_mm",
                "p_mm",
            ],
        )
        writer.writeheader()
        for item in output:
            writer.writerow(
                {
                    "id": item["id"],
                    "group": item["group"],
                    "shape_class": item["shape_class"],
                    "cables": item["cables"],
                    "origin_body": item["origin_body"],
                    "insertion_body": item["insertion_body"],
                    "origin_world": json.dumps(item["origin_world"]),
                    "insertion_world": json.dumps(item["insertion_world"]),
                    "x_mm": item["psm_parameters"]["x_mm"],
                    "y_mm": item["psm_parameters"]["y_mm"],
                    "z_mm": item["psm_parameters"]["z_mm"],
                    "p_mm": item["psm_parameters"]["p_mm"],
                }
            )

    print(f"exported {len(output)} enabled muscles")
    print(args.output_json)
    print(args.output_csv)


if __name__ == "__main__":
    main()
