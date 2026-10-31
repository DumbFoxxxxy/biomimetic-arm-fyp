"""Validate the JSON contracts before generating or optimizing the model."""

import argparse
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parent
ARM_CONFIG = ROOT / "config" / "arm_model.json"
MUSCLE_CONFIG = ROOT / "config" / "muscle_manifest.json"

ALLOWED_SHAPES = {
    "straight_diamond",
    "fan",
    "y_branch",
    "wrap",
}
REQUIRED_MUSCLES = {
    "supraspinatus",
    "infraspinatus",
    "deltoid_anterior",
    "deltoid_middle",
    "deltoid_posterior",
    "elbow_flexors",
    "triceps",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm-config", type=Path, default=ARM_CONFIG)
    parser.add_argument(
        "--muscle-config", type=Path, default=MUSCLE_CONFIG
    )
    args = parser.parse_args()

    arm = json.loads(args.arm_config.read_text(encoding="utf-8"))
    muscles = json.loads(
        args.muscle_config.read_text(encoding="utf-8")
    )

    body_names = set(arm["bodies"].values())
    ids = [muscle["id"] for muscle in muscles["muscles"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Muscle ids must be unique.")

    enabled = {
        muscle["id"]
        for muscle in muscles["muscles"]
        if muscle.get("enabled", False)
    }
    missing = REQUIRED_MUSCLES - enabled
    if missing:
        raise ValueError(f"Required muscles are not enabled: {sorted(missing)}")

    for muscle in muscles["muscles"]:
        muscle_id = muscle["id"]
        if muscle["shape_class"] not in ALLOWED_SHAPES:
            raise ValueError(f"Unsupported shape for {muscle_id}.")
        if int(muscle["cables"]) not in (1, 2):
            raise ValueError(f"Unsupported cable count for {muscle_id}.")
        if muscle["origin_body"] not in body_names:
            raise ValueError(f"Unknown origin body for {muscle_id}.")
        if muscle["insertion_body"] not in body_names:
            raise ValueError(f"Unknown insertion body for {muscle_id}.")
        for via in muscle.get("via", []):
            if via["body"] not in body_names:
                raise ValueError(f"Unknown via body for {muscle_id}.")

    defaults = muscles["defaults"]
    if not 1.5 <= defaults["branch_diameter_mm"] <= 3.0:
        raise ValueError("Default branch diameter must be within 1.5-3 mm.")
    if not 5.0 <= defaults["lattice_x_mm"] <= 20.0:
        raise ValueError("Default lattice x must be within 5-20 mm.")
    if not 5.0 <= defaults["lattice_y_mm"] <= 20.0:
        raise ValueError("Default lattice y must be within 5-20 mm.")
    if not 20.0 <= defaults["lattice_z_mm"] <= 40.0:
        raise ValueError("Default lattice z must be within 20-40 mm.")

    print(f"arm contract: {args.arm_config}")
    print(f"muscle contract: {args.muscle_config}")
    print(f"enabled muscles: {len(enabled)}")
    print(f"candidate muscles: {len(ids)}")
    print("contract validation: OK")


if __name__ == "__main__":
    main()
