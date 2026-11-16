"""Find the first K6 obstruction and the largest verified cable-driven sweep."""

from itertools import combinations

import mujoco
import numpy as np

from scapular_cable_control import (
    ASM, FOLDER, OUTPUT, ROOT, Geometry, find_hits, pose_metrics, provenance, read, save, simulate,
)
from project_runtime import load_mujoco_model


SCAN_DEG = .5
BOUNDARY_DEG = .025
OLD_SERRATUS = "FJ1459_g3"


def first_collision(probe, extent):
    n = round(abs(extent) / SCAN_DEG)
    clear = probe(0.)
    if clear["new_hits"] or clear["existing_hits"]:
        raise ValueError(f"Neutral collision: {clear}")
    for angle in np.linspace(0., extent, n + 1)[1:]:
        row = probe(angle)
        if not row["new_hits"] and not row["existing_hits"]:
            clear = row
            continue
        blocked = row
        while abs(blocked["input_deg"] - clear["input_deg"]) > BOUNDARY_DEG:
            mid = probe((blocked["input_deg"] + clear["input_deg"]) / 2)
            if mid["new_hits"] or mid["existing_hits"]:
                blocked = mid
            else:
                clear = mid
        return clear, blocked
    return clear, None


def motion_summary(result):
    return {key: value for key, value in result.items() if key != "frames"}


def make_probe(geometry, action, routes):
    """Check new routing and existing pairs whose separation can change on K6."""
    names = [*geometry.parts, *(key for key in geometry.routed if key != OLD_SERRATUS)]
    owners = [{geometry.owners[key]} if key in geometry.parts else set(geometry.routed[key]["owners"]) for key in names]
    groups = [{owner == "thorax" for owner in attached} for attached in owners]
    legacy_pairs = np.array([(i, j) for i, j in combinations(range(len(names)), 2)
                             if len(groups[i] | groups[j]) > 1])

    def probe(angle):
        state = geometry.prescribed(action, angle)
        new_hits = geometry.screen(state, routes)["hits"]
        bones, wires = geometry.posed_solids(state, excluded={OLD_SERRATUS})
        existing_hits = find_hits([*bones.values(), *wires.values()], names, legacy_pairs)
        geometry.cached.clear()
        return dict(input_deg=float(angle), new_hits=new_hits, existing_hits=existing_hits, **pose_metrics(state, geometry))

    return probe


def main():
    source = OUTPUT / "scapular_cable_control.json"
    saved = read(source)
    original = saved["trials"]["upward_rotation"]
    action_file = FOLDER / "whole_arm_anatomical_actions.json"
    action = read(action_file)["scapular"]["upward_rotation"]
    geometry = Geometry()
    model_file = ROOT / original["model"]["file"]
    if provenance(ASM) != saved["assembly"] or provenance(model_file) != original["model"]:
        raise ValueError("The selected assembly or dynamic model has changed")
    if (action["SC_axis"], action["SC_direction"], action["AC_per_SC_same_axis"]) != (1, 1., 0.):
        raise ValueError("Collision-pair filtering is defined for the verified SC-only K6 path")
    routes = [saved["routes"][key] for key in original["routes"]]
    for route in routes:
        if not route["installed"]:
            route["stock"] = geometry.stock(route)
            if route["stock"] is None:
                raise ValueError(f"Could not rebuild {route['id']} support envelope")
    probe = make_probe(geometry, action, routes)
    forward, block = first_collision(probe, 50.)
    print(f"Positive geometry: {forward['input_deg']:.5f} deg; first sampled collision {block}", flush=True)
    reverse, reverse_block = first_collision(probe, -6.)
    print(f"Negative return clearance: {reverse['input_deg']:.5f} deg; first sampled collision {reverse_block}", flush=True)
    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    if model.nq != 1:
        raise ValueError("Expected the same one-coordinate active K6 model")
    driven = original["driven"]
    max_geometry = forward["input_deg"]
    limits = (reverse["input_deg"], max_geometry)
    tests = []
    original_return = simulate(model, geometry, action, routes, driven, original["target_deg"],
                               neutral_preload=True, limits_deg=limits)
    brake = 1.5

    def try_target(target):
        result = simulate(model, geometry, action, routes, driven, target,
                          neutral_preload=True, limits_deg=limits, return_brake_Nms=brake)
        row = dict(target_deg=target, **motion_summary(result))
        tests.append(row)
        print(f"Target {target:.3f} deg: peak {result.get('peak_deg', 0.):.5f}; "
              f"return {result.get('final_deg', 'guard')}", flush=True)
        return result

    best = None
    best_target = None
    last_safe_target, first_unsafe_target = 0., None
    for target in np.arange(.5, max(50., max_geometry * 2.) + .5, .5):
        trial = try_target(float(target))
        if not trial["usable"] and ("reason" in trial or trial.get("peak_deg", 0.) >= 2.):
            first_unsafe_target = float(target)
            break
        last_safe_target = float(target)
        if trial["usable"] and (best is None or trial["peak_deg"] > best["peak_deg"]):
            best, best_target = trial, float(target)
    if first_unsafe_target is not None:
        while first_unsafe_target - last_safe_target > BOUNDARY_DEG:
            target = (last_safe_target + first_unsafe_target) / 2
            trial = try_target(target)
            if trial["usable"]:
                last_safe_target = target
                if best is None or trial["peak_deg"] > best["peak_deg"]:
                    best, best_target = trial, target
            else:
                first_unsafe_target = target
    output = dict(assembly=saved["assembly"], selected_control=provenance(source), model=original["model"],
                  action=provenance(action_file), scan_step_deg=SCAN_DEG, boundary_step_deg=BOUNDARY_DEG,
                  last_positive_clear=forward, first_positive_collision=block, last_negative_clear=reverse,
                  first_negative_collision=reverse_block, active={**saved["settings"], "return_brake_Nms_per_rad": brake},
                  original_return=motion_summary(original_return), target_trials=tests,
                  best_target_deg=best_target, best=best, first_unsafe_target_deg=first_unsafe_target)
    if best is not None:
        angles = np.linspace(best["minimum_angle_deg"], best["peak_deg"],
                             int(np.ceil((best["peak_deg"] - best["minimum_angle_deg"]) / .25)) + 1)
        interval = []
        for angle in angles:
            row = probe(angle)
            if row["new_hits"] or row["existing_hits"]:
                raise ValueError(f"Extra geometry sample found a collision: {row}")
            lengths = []
            for sign in (-1., 1.):
                data.qpos[0] = np.deg2rad(angle) + sign * 1e-5
                mujoco.mj_forward(model, data)
                lengths.append(data.ten_length.copy())
            leverage = 1000. * (lengths[0] - lengths[1]) / 2e-5
            if any(leverage[i] <= 0. for i, route in enumerate(routes) if route["id"] in driven):
                raise ValueError(f"Active tendon reverses pulling direction at {angle}")
            interval.append(dict(input_deg=float(angle), moment_arms_mm=leverage.tolist()))
        output["active_path_check"] = interval
        print(f"Verified active peak {best['peak_deg']:.6f} deg at target {best_target:.6f} deg", flush=True)
    save(OUTPUT / "upward_rotation_limit.json", output)


if __name__ == "__main__":
    main()
