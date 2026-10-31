"""Search S01 anchor and thorax-guide coordinates for the K6 baseline."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_k6_baseline import build_model
from k6_audit import K6Evaluator
from k6_common import (
    DEFAULT_ASSETS,
    DEFAULT_CONFIG,
    load_config,
    scapula_local_to_design,
    surface_point,
)


COARSE_A = {
    "x_mm": (85.0, 105.0, 120.0),
    "y_mm": (240.0, 270.0),
    "outward_offset_mm": (8.0, 12.0, 15.0),
}
COARSE_B_ACTIVE = {
    "u_mm": (-110.0, -90.0, -70.0),
    "v_mm": (80.0, 100.0, 120.0),
    "w_mm": (0.0, 10.0),
}


def route_from_values(
    a_x: float,
    a_y: float,
    a_offset: float,
    g_x: float,
    g_y: float,
    g_offset: float,
    b_u: float,
    b_v: float,
    b_w: float,
) -> dict:
    return {
        "A": {
            "x_mm": float(a_x),
            "y_mm": float(a_y),
            "outward_offset_mm": float(a_offset),
        },
        "G": {
            "x_mm": float(g_x),
            "y_mm": float(g_y),
            "outward_offset_mm": float(g_offset),
        },
        "B": {
            "u_mm": float(b_u),
            "v_mm": float(b_v),
            "w_mm": float(b_w),
        },
    }


def evaluator_from_route(config: dict, route: dict, steps: int) -> K6Evaluator:
    tree, _ = build_model(config, route)
    xml = ET.tostring(tree.getroot(), encoding="unicode")
    model = mujoco.MjModel.from_xml_string(xml)
    evaluator = K6Evaluator(None, config, model=model)
    return evaluator


def find_target_pull(
    evaluator: K6Evaluator,
    target_deg: float,
    tolerance_deg: float,
    max_pull_mm: float,
    steps: int,
) -> tuple[float | None, dict | None]:
    previous_pull = 0.0
    previous = evaluator.measure(0.0)
    if previous["up_deg"] >= target_deg - tolerance_deg:
        return 0.0, previous
    pull = 5.0
    while pull <= max_pull_mm + 1e-9:
        sample = evaluator.measure(pull)
        if sample["up_deg"] >= target_deg - tolerance_deg:
            low = previous_pull
            high = pull
            best = sample
            for _ in range(8):
                middle = 0.5 * (low + high)
                candidate = evaluator.measure(middle)
                if candidate["up_deg"] >= target_deg - tolerance_deg:
                    high = middle
                    best = candidate
                else:
                    low = middle
            return high, best
        previous_pull = pull
        previous = sample
        pull += 5.0
    return None, None


def evaluate_route(
    config: dict,
    route: dict,
    steps: int,
) -> dict:
    evaluator = evaluator_from_route(config, route, steps)
    pull, pose = find_target_pull(
        evaluator,
        float(config["target"]["up_deg"]),
        float(config["target"]["up_tolerance_deg"]),
        float(config["target"]["max_pull_mm"]),
        steps,
    )
    result = {
        "route": route,
        "target_pull_mm": pull,
        "pose": pose,
    }
    if pull is None or pose is None:
        max_pose = evaluator.measure(float(config["target"]["max_pull_mm"]))
        result["max_pull_pose"] = max_pose
        result["score"] = 1000.0 + max(
            0.0, float(config["target"]["up_deg"]) - max_pose["up_deg"]
        ) * 20.0
        return result

    delta = 0.5
    before = evaluator.measure(max(0.0, pull - delta))
    after = evaluator.measure(min(float(config["target"]["max_pull_mm"]), pull + delta))
    angle_delta = math.radians(after["up_deg"] - before["up_deg"])
    moment_arm = None
    if abs(angle_delta) > 1e-9:
        moment_arm = -(
            after["cable_length_mm"] - before["cable_length_mm"]
        ) / angle_delta
    contraction = (
        evaluator.neutral_length_m * 1000.0 - pose["cable_length_mm"]
    ) / (evaluator.neutral_length_m * 1000.0) * 100.0
    target_arm = float(config["target"]["moment_arm_target_mm"])
    score = 0.0
    if contraction > float(config["target"]["max_contraction_percent"]):
        score += (contraction - float(config["target"]["max_contraction_percent"])) * 20.0
    if moment_arm is None:
        score += 100.0
    else:
        score += abs(moment_arm - target_arm)
    score += max(0.0, pose["visual_clearance_mm"] * -1.0 - 0.5) * 5.0
    score += max(0.0, pose["contact_gap_mm"] * -1.0 - 0.5) * 20.0
    score += max(0.0, pose["contact_gap_mm"] - 5.0) * 0.5
    result.update(
        {
            "contraction_percent": contraction,
            "moment_arm_mm": moment_arm,
            "score": float(score),
        }
    )
    return result


def closest_guide_to_line(
    config: dict,
    a_design: np.ndarray,
    b_design: np.ndarray,
) -> dict | None:
    guide_radius = float(config["s01"]["guide_radius_mm"])
    direction = b_design - a_design
    denominator = float(np.dot(direction, direction))
    best = None
    for g_x in np.arange(95.0, 140.0 + 1e-9, 5.0):
        for g_y in np.arange(220.0, 275.0 + 1e-9, 5.0):
            for offset in (12.0, 16.0, 20.0):
                guide = surface_point(g_x, g_y, offset, config)
                parameter = float(np.dot(guide - a_design, direction) / denominator)
                if not 0.08 <= parameter <= 0.92:
                    continue
                nearest = a_design + parameter * direction
                distance = float(np.linalg.norm(guide - nearest))
                if distance > guide_radius:
                    continue
                candidate = {
                    "distance_mm": distance,
                    "parameter": parameter,
                    "G": {
                        "x_mm": float(g_x),
                        "y_mm": float(g_y),
                        "outward_offset_mm": float(offset),
                    },
                }
                if best is None or distance < best["distance_mm"]:
                    best = candidate
    return best


def coarse_routes(config: dict) -> list[dict]:
    routes = []
    for a_x in COARSE_A["x_mm"]:
        for a_y in COARSE_A["y_mm"]:
            for a_offset in COARSE_A["outward_offset_mm"]:
                a_design = surface_point(a_x, a_y, a_offset, config)
                for b_u in COARSE_B_ACTIVE["u_mm"]:
                    for b_v in COARSE_B_ACTIVE["v_mm"]:
                        for b_w in COARSE_B_ACTIVE["w_mm"]:
                            b_design = scapula_local_to_design(
                                np.asarray([b_u, b_v, b_w]), config
                            )
                            guide = closest_guide_to_line(
                                config, a_design, b_design
                            )
                            if guide is None:
                                continue
                            route = route_from_values(
                                a_x,
                                a_y,
                                a_offset,
                                guide["G"]["x_mm"],
                                guide["G"]["y_mm"],
                                guide["G"]["outward_offset_mm"],
                                b_u,
                                b_v,
                                b_w,
                            )
                            route["guide_line_distance_mm"] = guide[
                                "distance_mm"
                            ]
                            routes.append(route)
    return routes


def refined_routes(
    seeds: list[dict],
    count: int,
    seed: int = 17,
) -> list[dict]:
    rng = np.random.default_rng(seed)
    routes = []
    for _ in range(count):
        base = seeds[int(rng.integers(0, len(seeds)))]
        a = base["A"]
        g = base["G"]
        b = base["B"]
        routes.append(
            route_from_values(
                np.clip(a["x_mm"] + rng.uniform(-25.0, 25.0), 55.0, 120.0),
                np.clip(a["y_mm"] + rng.uniform(-35.0, 35.0), 180.0, 270.0),
                np.clip(
                    a["outward_offset_mm"] + rng.uniform(-3.0, 3.0),
                    8.0,
                    15.0,
                ),
                np.clip(g["x_mm"] + rng.uniform(-12.0, 12.0), 95.0, 140.0),
                np.clip(g["y_mm"] + rng.uniform(-12.0, 12.0), 220.0, 275.0),
                np.clip(
                    g["outward_offset_mm"] + rng.uniform(-4.0, 4.0),
                    12.0,
                    20.0,
                ),
                np.clip(b["u_mm"] + rng.uniform(-12.0, 12.0), -110.0, -60.0),
                np.clip(b["v_mm"] + rng.uniform(-12.0, 12.0), 80.0, 130.0),
                np.clip(b["w_mm"] + rng.uniform(-5.0, 5.0), -10.0, 10.0),
            )
        )
    return routes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--out", type=Path, default=DEFAULT_ASSETS / "route_search.json")
    parser.add_argument("--coarse-steps", type=int, default=700)
    parser.add_argument("--refine-steps", type=int, default=1000)
    parser.add_argument("--refine-count", type=int, default=240)
    args = parser.parse_args()
    config = load_config(args.config)

    coarse_results = []
    coarse = coarse_routes(config)
    for index, route in enumerate(coarse, start=1):
        result = evaluate_route(config, route, args.coarse_steps)
        coarse_results.append(result)
        if index % 12 == 0:
            print(f"coarse {index}/{len(coarse)}")
    coarse_results.sort(key=lambda item: item["score"])
    seeds = [item["route"] for item in coarse_results[:8]]

    refined_results = []
    routes = refined_routes(seeds, args.refine_count)
    for index, route in enumerate(routes, start=1):
        result = evaluate_route(config, route, args.refine_steps)
        refined_results.append(result)
        if index % 12 == 0:
            print(f"refine {index}/{len(routes)}")
    all_results = coarse_results + refined_results
    all_results.sort(key=lambda item: item["score"])
    best = all_results[0]
    output = {
        "best": best,
        "top": all_results[:20],
        "coarse_count": len(coarse_results),
        "refined_count": len(refined_results),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(best, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
