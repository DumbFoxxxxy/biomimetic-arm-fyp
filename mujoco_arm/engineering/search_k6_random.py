"""Random active-wrap search for the S01 guide and anchors."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import json
from project_runtime import WorkspacePath as Path
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
from search_k6_route import evaluate_route, route_from_values


def active_guide_candidates(
    config: dict,
    a_design: np.ndarray,
    b_design: np.ndarray,
) -> list[dict]:
    radius = float(config["s01"]["guide_radius_mm"])
    direction = b_design - a_design
    denominator = float(np.dot(direction, direction))
    candidates = []
    for g_x in np.arange(95.0, 140.0 + 1e-9, 2.5):
        for g_y in np.arange(220.0, 275.0 + 1e-9, 2.5):
            for offset in (12.0, 14.0, 16.0, 18.0, 20.0):
                guide = surface_point(g_x, g_y, offset, config)
                parameter = float(np.dot(guide - a_design, direction) / denominator)
                if not 0.08 <= parameter <= 0.92:
                    continue
                nearest = a_design + parameter * direction
                distance = float(np.linalg.norm(guide - nearest))
                if 5.0 <= distance <= radius - 0.05:
                    candidates.append(
                        {
                            "x_mm": float(g_x),
                            "y_mm": float(g_y),
                            "outward_offset_mm": float(offset),
                            "line_distance_mm": distance,
                        }
                    )
    return candidates


def random_routes(config: dict, count: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    routes = []
    seen = set()
    attempts = 0
    while len(routes) < count and attempts < count * 30:
        attempts += 1
        a = {
            "x_mm": float(rng.uniform(55.0, 120.0)),
            "y_mm": float(rng.uniform(180.0, 270.0)),
            "outward_offset_mm": float(rng.uniform(8.0, 15.0)),
        }
        b = {
            "u_mm": float(rng.uniform(-110.0, -60.0)),
            "v_mm": float(rng.uniform(80.0, 130.0)),
            "w_mm": float(rng.uniform(-10.0, 10.0)),
        }
        a_design = surface_point(
            a["x_mm"], a["y_mm"], a["outward_offset_mm"], config
        )
        b_design = scapula_local_to_design(
            np.asarray([b["u_mm"], b["v_mm"], b["w_mm"]]), config
        )
        candidates = active_guide_candidates(config, a_design, b_design)
        if not candidates:
            continue
        candidate = candidates[int(rng.integers(0, len(candidates)))]
        route = {
            "A": a,
            "G": {
                "x_mm": candidate["x_mm"],
                "y_mm": candidate["y_mm"],
                "outward_offset_mm": candidate["outward_offset_mm"],
            },
            "B": b,
            "guide_line_distance_mm": candidate["line_distance_mm"],
        }
        key = json.dumps(route, sort_keys=True)
        if key in seen:
            continue
        seen.add(key)
        routes.append(route)
    return routes


def screen_route(config: dict, route: dict, steps: int) -> dict:
    tree, _ = build_model(config, route)
    xml = ET.tostring(tree.getroot(), encoding="unicode")
    model = _load_project_xml(xml)
    evaluator = K6Evaluator(None, config, model=model)
    pose = evaluator.measure(float(config["target"]["max_pull_mm"]))
    return {
        "route": route,
        "max_pull_pose": pose,
        "screen_up_deg": pose["up_deg"],
        "screen_gap_mm": pose["contact_gap_mm"],
        "screen_visual_mm": pose["visual_clearance_mm"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--out", type=Path, default=DEFAULT_ASSETS / "route_search_active.json"
    )
    parser.add_argument("--count", type=int, default=450)
    parser.add_argument("--screen-steps", type=int, default=650)
    parser.add_argument("--full-steps", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=29)
    args = parser.parse_args()
    config = load_config(args.config)
    routes = random_routes(config, args.count, args.seed)
    screened = []
    for index, route in enumerate(routes, start=1):
        screened.append(screen_route(config, route, args.screen_steps))
        if index % 25 == 0:
            print(f"screen {index}/{len(routes)}")
    screened.sort(
        key=lambda item: (
            -item["screen_up_deg"],
            -item["screen_visual_mm"],
            abs(item["screen_gap_mm"]),
        )
    )
    full_results = []
    for item in screened[:36]:
        result = evaluate_route(config, item["route"], args.full_steps)
        full_results.append(result)
    full_results.sort(key=lambda item: item["score"])
    output = {
        "best": full_results[0] if full_results else None,
        "top": full_results[:20],
        "best_screen_only": screened[:20],
        "route_count": len(routes),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["best"], indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
