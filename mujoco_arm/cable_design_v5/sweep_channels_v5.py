"""Evaluate 16-, 17- and 18-channel V5 variants without modifying the baseline."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from build_v2_model import build  # noqa: E402
from v2_common import load_json  # noqa: E402
from v5_evaluator import V5Evaluator  # noqa: E402


CONFIG = ROOT / "cable_design_v5.json"
MATRIX = ROOT / "action_cable_matrix_v5.json"
EXPORT = ROOT.parent / "exports" / "cable_design_v5"


def sync_channels(config: dict, matrix: dict) -> None:
    channel_ids = [route["id"] for route in config["routes"]]
    for action in matrix["actions"]:
        action["allowed_cables"] = list(channel_ids)
        action["preferred_cables"] = [
            route_id
            for route_id in action.get("preferred_cables", [])
            if route_id in channel_ids
        ]
        action["mandatory_cables"] = [
            route_id
            for route_id in action.get("mandatory_cables", [])
            if route_id in channel_ids
        ]


def variant(
    base_config: dict,
    base_matrix: dict,
    name: str,
    remove: str | None = None,
    add_from: str | None = None,
    add_id: str | None = None,
) -> tuple[dict, dict]:
    config = copy.deepcopy(base_config)
    matrix = copy.deepcopy(base_matrix)
    if remove:
        config["routes"] = [
            route for route in config["routes"] if route["id"] != remove
        ]
    if add_from and add_id:
        route = copy.deepcopy(
            next(item for item in config["routes"] if item["id"] == add_from)
        )
        route["id"] = add_id
        route["role"] = f"v5_sweep_{add_id.lower()}"
        route["status"] = "channel_sweep_candidate"
        for element in route["path"]:
            if element["kind"] == "guide":
                element["name"] = f"guide_{add_id}_v5"
        config["routes"].append(route)
        for action in matrix["actions"]:
            if add_id == "G08" and action["group"] in ("gh", "combined"):
                action.setdefault("preferred_cables", []).append(add_id)
            if add_id == "S07" and action["group"] in ("scapula", "combined"):
                action.setdefault("preferred_cables", []).append(add_id)
    sync_channels(config, matrix)
    return config, matrix


def evaluate_variant(name: str, config: dict, matrix: dict) -> dict:
    tree, _ = build(config)
    EXPORT.mkdir(parents=True, exist_ok=True)
    model_path = EXPORT / f"sweep_{name}.xml"
    tree.write(model_path, encoding="utf-8", xml_declaration=True)
    evaluator = V5Evaluator(model_path, config, matrix)
    evaluation = evaluator.evaluate_all()
    failures = {
        action["id"]: [
            key
            for key in (
                "pass_80_percent",
                "pass_tension",
                "pass_positive_spanning",
                "pass_contraction",
                "pass_cable_clearance",
                "pass_guide_geometry",
                "pass_visual_clearance",
                "pass_contact",
            )
            if not action[key]
        ]
        for action in evaluation["actions"]
    }
    return {
        "name": name,
        "channel_count": len(config["routes"]),
        "full_pass": all(not value for value in failures.values()),
        "failures": {
            action_id: values
            for action_id, values in failures.items()
            if values
        },
    }


def main() -> None:
    config = load_json(CONFIG)
    matrix = load_json(MATRIX)
    variants = [
        ("channels16_no_g07", *variant(
            config, matrix, "channels16_no_g07", remove="G07"
        )),
        ("channels17", config, matrix),
        ("channels18_g08", *variant(
            config,
            matrix,
            "channels18_g08",
            add_from="G07",
            add_id="G08",
        )),
        ("channels18_s07", *variant(
            config,
            matrix,
            "channels18_s07",
            add_from="S06",
            add_id="S07",
        )),
    ]
    results = []
    for name, variant_config, variant_matrix in variants:
        results.append(
            evaluate_variant(name, variant_config, variant_matrix)
        )
    output = {"variants": results}
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "channel_sweep_v5.json").write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
