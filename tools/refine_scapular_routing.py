"""Try one supported local guide for the current serratus wire, without a muscle mesh."""

import argparse
import csv
from itertools import product
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from check_scapular_upward_limit import OLD_SERRATUS, first_collision, make_probe, motion_summary
from scapular_cable_control import (
    FOLDER, OUTPUT, ROOT, RADIUS_MM, Geometry, active_model, length, mesh_geom, posed_points,
    provenance, read, rod, route_loads, save, serialize, simulate, transforms, view_trajectory,
)
from study_elbow_series_package import bounded, overlaps
from project_runtime import load_mujoco_model


RIB = "thorax_reference:FJ3344"
TARGETS = (10., 8., 6.)
RESULT = OUTPUT / "serratus_guide_refinement.json"


def support_vertices(geometry, key):
    entry = geometry.sources[key]
    return mesh_geom.load_obj(ROOT / entry["file"])[0] + entry["translation_source_mm"]


def guide_options(geometry, owner, region):
    parts = (RIB, "thorax_reference:FJ3342", "thorax_reference:FJ3346") if owner == "thorax" else ("bone:scapula",)
    axes = ((-1, 0, 0), (-1, -1, 0), (-1, 1, 0), (-1, 0, 1), (-1, 0, -1),
            (0, -1, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1))
    if owner == "scapula":
        axes = (*axes, (1, -1, 0), (1, 0, -1), (1, 0, 0))
    anchors = np.concatenate([support_vertices(geometry, key) for key in parts])
    anchors = anchors[np.all(abs(anchors - region) <= [20., 27., 24.], axis=1)]
    _, unique = np.unique(np.round(anchors / 3.), return_index=True, axis=0)
    tested = set()
    for point, axis, lift in product(anchors[unique], axes, (3., 6., 9.)):
        vector = np.array(axis) / np.linalg.norm(axis)
        point = point + lift * vector
        key = tuple(np.round(point, 2))
        if key in tested:
            continue
        tested.add(key)
        site = geometry.mount(owner, point)
        if site["support"] in parts and 1.5 <= site["standoff_mm"] <= 9.01:
            site.update(role="new_wrap", reference_distance_mm=float(np.linalg.norm(point - region)))
            yield site


def turn_deflection(points):
    directions = np.diff(points, axis=0)
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    return float(np.rad2deg(np.arccos(np.clip(directions[-2] @ directions[-1], -1., 1.))))


def screen_candidate_points(geometry, action, original, selection, angles, regions):
    old_points, frames, solid_targets = [], [], []
    for angle in angles:
        data = geometry.prescribed(action, angle)
        old_points.append(posed_points(data, geometry.origins, original))
        frames.append(transforms(geometry.model, data, geometry.origins))
        solids, _ = geometry.posed_solids(data, excluded={OLD_SERRATUS})
        targets = bounded(solids)
        solid_targets.append((targets, list(targets), np.asarray([value[1] for value in targets.values()])))
    candidates = []
    for owner, centre in regions.items():
        for site in guide_options(geometry, owner, centre):
            port = np.asarray(site["port_mm"])
            guide_points = [port @ frame[owner][:, :3].T + frame[owner][:, 3] for frame in frames]
            paths = [np.vstack([points[:2], guide, points[2:]]) for points, guide in zip(old_points, guide_points, strict=True)]
            long = np.array([length(points) for points in paths])
            leverage = -np.diff(long) / np.deg2rad(np.diff(angles))
            deflection = max(turn_deflection(points) for points in paths)
            if min(leverage) <= 0. or any(np.linalg.norm(np.diff(points, axis=0), axis=1).min() < 5. for points in paths):
                continue
            # Prefer modest bends and short supports around the section that hits the rib.
            cost = deflection + 1.8 * site["standoff_mm"] + .35 * site["reference_distance_mm"]
            candidates.append((float(cost), site, paths, leverage, deflection))
    candidates.sort(key=lambda row: row[0])
    failures = {}
    for cost, site, paths, leverage, deflection in candidates:
        reason = None
        for index in reversed(range(len(angles))):
            path = paths[index]
            targets, names, bounds = solid_targets[index]
            # Fast exact rejection on the two altered spans before checking all interfaces.
            for start, end in zip(path[1:], path[2:]):
                low = np.minimum(start, end) - RADIUS_MM
                high = np.maximum(start, end) + RADIUS_MM
                near = np.flatnonzero(np.all(low <= bounds[:, 3:], axis=1) & np.all(high >= bounds[:, :3], axis=1))
                nearby = {names[i]: targets[names[i]] for i in near}
                hits = overlaps(rod(start, end, RADIUS_MM), nearby)
                if hits:
                    reason = max(hits, key=hits.get)
                    break
            if reason:
                break
        if reason:
            failures[reason] = failures.get(reason, 0) + 1
            continue
        trial = {**original, "sites": [*original["sites"][:2], site, *original["sites"][2:]],
                 "guide_score": cost, "coarse_moment_arm_min_mm": float(min(leverage)),
                 "maximum_added_guide_deflection_deg": deflection}
        trial["stock"] = geometry.stock(trial)
        if trial["stock"] is None:
            failures["support footprint"] = failures.get("support footprint", 0) + 1
            continue
        selected = [trial, *selection]
        for angle in angles:
            hits = geometry.screen(geometry.prescribed(action, angle), selected)["hits"]
            if hits:
                reason = f"{hits[0]['left']} / {hits[0]['right']}"
                failures[reason] = failures.get(reason, 0) + 1
                break
        if reason is None:
            return trial, dict(screened=len(candidates), failures=failures)
    return None, dict(screened=len(candidates), failures=failures)


def angles_for(target):
    return sorted({0., 2., 4., target - 2., target})


def direction_limit(model, routes, driven, maximum):
    data = mujoco.MjData(model)
    indices = [index for index, route in enumerate(routes) if route["id"] in driven]

    def check(angle):
        length_samples = []
        for sign in (-1., 1.):
            data.qpos[0] = np.deg2rad(angle) + sign * 1e-5
            mujoco.mj_forward(model, data)
            length_samples.append(data.ten_length.copy())
        arms = 1000. * (length_samples[0] - length_samples[1]) / 2e-5
        return dict(angle_deg=float(angle), moment_arms_mm=arms.tolist(), positive=bool(min(arms[indices]) > 0.))

    previous = check(0.)
    if not previous["positive"]:
        raise ValueError("Refined route has the wrong direction at neutral")
    for angle in np.linspace(0., maximum, int(np.ceil(maximum / .25)) + 1)[1:]:
        row = check(angle)
        if row["positive"]:
            previous = row
            continue
        while row["angle_deg"] - previous["angle_deg"] > .025:
            mid = check((row["angle_deg"] + previous["angle_deg"]) / 2)
            if mid["positive"]:
                previous = mid
            else:
                row = mid
        return previous, row
    return previous, None


def dynamic_trial(model, geometry, action, routes, limit, return_limit, driven):
    cap = (return_limit, limit)
    trials = []
    best = None
    command = None
    for target in np.arange(limit, limit * 1.7 + .25, .25):
        row = simulate(model, geometry, action, routes, driven, float(target),
                       neutral_preload=True, limits_deg=cap, return_brake_Nms=1.5)
        trials.append(dict(target_deg=float(target), **motion_summary(row)))
        if row["usable"]:
            best, command = row, float(target)
        elif best is not None:
            break
    if best is None:
        return None, trials
    lo, hi = command, command + .25
    for _ in range(4):
        mid = (lo + hi) / 2
        row = simulate(model, geometry, action, routes, driven, mid,
                       neutral_preload=True, limits_deg=cap, return_brake_Nms=1.5)
        trials.append(dict(target_deg=mid, **motion_summary(row)))
        if row["usable"]:
            lo = mid
            if row["peak_deg"] > best["peak_deg"]:
                best = row
        else:
            hi = mid
    return best, trials


def restore_stock(geometry, routes):
    for route in routes:
        if not route["installed"]:
            route["stock"] = geometry.stock(route)
            if not route["stock"]:
                raise ValueError(f"Support could not be rebuilt: {route['id']}")


def main():
    saved_path = OUTPUT / "scapular_cable_control.json"
    saved = read(saved_path)
    original = saved["trials"]["upward_rotation"]
    action = read(FOLDER / "whole_arm_anatomical_actions.json")["scapular"]["upward_rotation"]
    geometry = Geometry()
    routes = [saved["routes"][name] for name in original["routes"]]
    restore_stock(geometry, routes)
    spots = {"thorax": np.array([-98., -28., 1230.]), "scapula": np.array([-81., -13., 1217.])}
    attempts = []
    for target in TARGETS:
        trial, stats = screen_candidate_points(geometry, action, routes[0], routes[1:], angles_for(target), spots)
        attempts.append(dict(coarse_target_deg=target, **stats))
        print(f"{target:.0f}° routing search: {stats}; selected {trial is not None}", flush=True)
        if trial:
            break
    if trial is None:
        save(RESULT, dict(source=provenance(saved_path), candidate=None, attempts=attempts))
        return
    selected = [trial, *routes[1:]]
    probe = make_probe(geometry, action, selected)
    last, first = first_collision(probe, 50.)
    floor, first_negative = first_collision(probe, -6.)
    print(f"Geometric boundary {last['input_deg']:.3f}°, next point {first}", flush=True)
    tree, model = active_model(geometry, action, selected)
    positive, sign_loss = direction_limit(model, selected, original["driven"], last["input_deg"])
    print(f"Both active wires aid K6 up to {positive['angle_deg']:.3f}°; next point {sign_loss}", flush=True)
    best, commands = dynamic_trial(
        model, geometry, action, selected, positive["angle_deg"], floor["input_deg"], original["driven"])
    if best is None:
        raise ValueError(f"No active K6 trial returned safely: {commands}")
    checked = []
    data = mujoco.MjData(model)
    low, high = best["minimum_angle_deg"], best["peak_deg"]
    samples = np.linspace(low, high, int(np.ceil((high - low) / .25)) + 1)
    for angle in samples:
        row = probe(angle)
        if row["new_hits"] or row["existing_hits"]:
            raise ValueError(f"Refined trajectory found earlier interference: {row}")
        lengths = []
        for sign in (-1., 1.):
            data.qpos[0] = np.deg2rad(angle) + sign * 1e-5
            mujoco.mj_forward(model, data)
            lengths.append(data.ten_length.copy())
        sign_test = 1000. * (lengths[0] - lengths[1]) / 2e-5
        if min(sign_test[:2]) <= 0.:
            raise ValueError(f"An active route changes sign at {angle}")
        state = geometry.prescribed(action, angle)
        gap = geometry.screen(state, selected, measure_gap=True)
        geometry.cached.clear()
        checked.append(dict(**row, moment_arms_mm=sign_test.tolist(), **gap))
    for angle in samples[::-1]:
        row = probe(angle)
        if row["new_hits"] or row["existing_hits"]:
            raise ValueError(f"Interference in reverse audit: {row}")
    destination = OUTPUT / "upward_rotation_extra_guide.xml"
    destination.write_bytes(ET.tostring(tree, encoding="utf-8", xml_declaration=True).replace(b"\r\n", b"\n"))
    reactions = route_loads(mujoco.MjData(model), geometry, selected, best, model)
    with (OUTPUT / "serratus_guide_points.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("role", "owner", "port_x_mm", "port_y_mm", "port_z_mm",
                         "base_x_mm", "base_y_mm", "base_z_mm", "support", "peak_reaction_N"))
        for index, site in enumerate(trial["sites"]):
            writer.writerow((site["role"], site["owner"], *site["port_mm"], *site["base_mm"],
                             site["support"], reactions[f"{trial['id']}:{index}"]["resultant_N"]))
    result = dict(source=provenance(saved_path), candidate=serialize(trial), attempts=attempts,
                  model=provenance(destination), last_positive_clear=last, first_positive_collision=first,
                  last_negative_clear=floor, first_negative_collision=first_negative,
                  last_both_producing_upward_moment=positive, first_loss_of_active_direction=sign_loss,
                  best=best, commands=commands, checked_angles=checked, checked_return=len(checked),
                  peak_local_reaction=reactions, scope="One additional fixed guide on current K6 geometry; "
                  "stock clearance, nominal cable line and constrained tension simulation.")
    save(RESULT, result)
    print(f"Active peak {best['peak_deg']:.4f}°, returns to {best['final_deg']:.4f}°", flush=True)


def check_candidate():
    record = read(RESULT)
    source = read(OUTPUT / "scapular_cable_control.json")
    path = ROOT / record["model"]["file"]
    if provenance(path) != record["model"]:
        raise ValueError(f"Replay model changed: {path}")
    model = load_mujoco_model(path)
    data = mujoco.MjData(model)
    for frame in record["best"]["frames"]:
        data.qpos[:] = frame["qpos"]
        mujoco.mj_forward(model, data)
        np.testing.assert_allclose(1000. * data.ten_length, frame["length_mm"], rtol=0., atol=1e-7)
        if min(frame["tension_N"]) < 0.:
            raise ValueError("The cable motor applied negative tension")
    print(f'K6: {len(record["best"]["frames"])} replay frames verified', flush=True)
    geometry = Geometry()
    selected = [record["candidate"], source["routes"]["trapezius_superior"], source["routes"]["FJ1536_g2"]]
    restore_stock(geometry, selected)
    actions = read(FOLDER / "whole_arm_anatomical_actions.json")["scapular"]
    for key in ("protraction", "retraction_alternative"):
        old = source["trials"][key]
        action = actions[key]
        _, model = active_model(geometry, action, selected)
        motion = simulate(model, geometry, action, selected, old["driven"], old["target_deg"], neutral_preload=True)
        if not motion["frames"]:
            print(key, motion.get("reason"), flush=True)
            continue
        frames = motion["frames"]
        data = mujoco.MjData(model)
        hits = []
        for index in sorted({*range(0, len(frames), 10), int(np.argmax([frame["sc_deg"] for frame in frames])),
                             len(frames) - 1}):
            data.qpos[:] = frames[index]["qpos"]
            mujoco.mj_forward(model, data)
            hits = geometry.screen(data, selected, model)["hits"]
            if hits:
                break
        travel = [min(frame["anterior_mm"] for frame in frames), max(frame["anterior_mm"] for frame in frames)]
        print(key, motion["usable"], travel, hits[:1], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--view", action="store_true", help="Replay the selected routing candidate")
    mode.add_argument("--check", action="store_true", help="Check the replay and neighboring motion directions")
    args = parser.parse_args()
    if args.view:
        record = read(RESULT)
        path = ROOT / record["model"]["file"]
        if provenance(path) != record["model"]:
            raise ValueError(f"Replay model changed: {path}")
        view_trajectory(path, record["best"]["frames"])
    elif args.check:
        check_candidate()
    else:
        main()
