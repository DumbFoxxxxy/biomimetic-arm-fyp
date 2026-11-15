"""Verify the elbow on the same fixed assembly used for the whole-arm K7 replay."""

from collections import Counter
import argparse
import csv
from itertools import combinations
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from project_runtime import load_mujoco_model
from study_s1_local_force_geometry import load_reference
from study_elbow_series_package import union


CAP = "elbow:FJ1477_proximal_rigid_cap"


def section(route, frames):
    points = np.asarray(route["points_source_mm"])
    weights = np.asarray(route["distal_weights"])[:, None]
    a, b = (frames[owner] for owner in route["owners"])
    moved = (1. - weights) * (points @ a[:, :3].T + a[:, 3]) + weights * (points @ b[:, :3].T + b[:, 3])
    shape = rod(*moved, .30) if len(moved) == 2 else A.D.outer_envelope(moved, .30)
    return shape, float(np.sum(np.linalg.norm(np.diff(moved, axis=0), axis=1)))


def closest_gap(placed, pairs, cap=3.):
    bounds = np.asarray([solid.bounding_box() for solid in placed])
    left, right = bounds[pairs[:, 0]], bounds[pairs[:, 1]]
    floor = np.linalg.norm(np.maximum(0., np.maximum(left[:, :3] - right[:, 3:], right[:, :3] - left[:, 3:])), axis=1)
    closest = None
    for index in np.argsort(floor):
        if floor[index] >= cap:
            break
        i, j = pairs[index]
        gap = A.I.gap(placed[i], placed[j], cap)
        if gap < cap:
            cap, closest = gap, (i, j)
    return cap, closest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cap-relief", action="store_true", help="Check a guarded notch on the elbow cap")
    args = parser.parse_args()
    folder = DEST.parent / "whole_arm_first_pass"
    combined_ref = provenance(folder / "whole_arm_combined_cuff_k7.json")
    combined = load_reference(combined_ref)
    motion = load_reference(combined["motion"])
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    routed = load_reference(combined["routed"])
    routes = routed["routes"]
    support = load_reference(combined["support"])
    revision = load_reference(load_reference(support["audit"])["humeral_revision"])["selected"]["meshes"]
    base_ref = provenance(folder / "tabletop_base_check.json")
    base = load_reference(base_ref)
    replay_ref = provenance(folder / "whole_arm_passive_elbow_audit.json")
    replay = load_reference(replay_ref)
    if replay["scene"] != motion["scene_check"] or replay["base"] != base_ref:
        raise ValueError("Elbow trajectory belongs to a different reference assembly")
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
               "elbow:FJ1477_bone_adapter": support["source"], "base:T12": base["source"]}
    sources.update({key: entry for key, entry in revision.items() if key.startswith("link:") and key in sources})
    owners = {**transfer["owner_by_source"], "base:T12": "thorax"}
    assembly_file = folder / "whole_arm_latest_neutral.xml"
    assembly_ref = provenance(assembly_file)
    xml = ET.parse(assembly_file)
    xml_assets = {row.get("name"): row.get("file") for row in xml.findall("./asset/mesh")}
    origins = {key: np.asarray(value) for key, value in scene["neutral_body_origins_mm"].items()}

    def asset_placement(file, owner, local_position):
        return file.replace("\\", "/"), owner, tuple(np.round(local_position, 10))

    rigid_files = Counter(asset_placement(xml_assets[geom.get("mesh")], body.get("name"),
                                          np.fromstring(geom.get("pos", "0 0 0"), sep=" "))
                          for body in xml.findall(".//body") for geom in body.findall("geom") if geom.get("group") != "3")
    placements = Counter(asset_placement(entry["file"], owners[key],
                                        (entry["translation_source_mm"] - origins[owners[key]]) * .001)
                         for key, entry in sources.items())
    if rigid_files != placements:
        raise ValueError("Rigid audit list differs from the K7 MJCF assembly")
    model = load_mujoco_model(assembly_file)
    data = mujoco.MjData(model)
    parts = {key: load_component(entry) for key, entry in sources.items()}
    extended = load_reference(replay["extended_rails"])["shafts"]
    long_rails = {key: load_component(entry) for key, entry in extended.items()}
    keys = list(parts) + list(routes)
    part_count = len(parts)
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                        if j >= part_count or owners[keys[i]] != owners[keys[j]]])
    moving_bodies = {"forearm", "FJ1486_distal_pivot", "FJ1477_distal_pivot", "FJ1477_carriage"}
    elbow_lines = {"FJ1486", "FJ1477"}
    moving_indices = {i for i, key in enumerate(keys) if key in elbow_lines or owners.get(key) in moving_bodies}
    dynamic_pairs = np.asarray([(i, j) for i, j in pairs if i in moving_indices or j in moving_indices])
    moving_wire_pairs = np.asarray([(i, j) for i, j in dynamic_pairs if j >= part_count])
    mujoco.mj_forward(model, data)
    neutral_frames = transforms(model, data, origins)
    static_wire = {key: section(route, neutral_frames) for key, route in routes.items() if key not in elbow_lines}

    def k1(qpos, tick):
        data.qpos[:] = qpos
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        placed = [parts[key].transform(frames[owners[key]]) if owners[key] in moving_bodies else parts[key] for key in parts]
        lines = {key: static_wire[key] if key in static_wire else section(route, frames) for key, route in routes.items()}
        placed.extend(shape for shape, _ in lines.values())
        hits = find_hits(placed, keys, pairs if tick == 0 else dynamic_pairs)
        distance, contact = closest_gap(placed, moving_wire_pairs) if not hits else (0., None)
        return dict(trace_frame=tick, qpos=qpos.tolist(), elbow_deg=float(np.rad2deg(qpos[elbow_index])),
                    carriage_mm=float(1000. * qpos[slide_index]), hits=hits,
                    length_mm={key: length for key, (_, length) in lines.items()},
                    moving_wire_interface_gap_capped_3mm=distance,
                    nearest_pair=[keys[i] for i in contact] if contact is not None else None)

    donor = load_reference(replay["donor"])
    donor_trace = ROOT / donor["trace"]["file"]
    if provenance(donor_trace) != donor["trace"]:
        raise ValueError("Elbow donor trace has changed")
    donor_model_path = ROOT / donor["model"]["file"]
    if provenance(donor_model_path) != donor["model"]:
        raise ValueError("Elbow donor model has changed")
    donor_model = load_mujoco_model(donor_model_path)
    donor_joints = ("elbow", "FJ1486_swivel", "FJ1477_swivel", "FJ1477_series_travel")
    shared_indices = [model.joint(joint).qposadr.item() for joint in donor_joints]
    donor_indices = [donor_model.joint(joint).qposadr.item() for joint in donor_joints]
    elbow_index, slide_index = shared_indices[0], shared_indices[-1]
    with donor_trace.open(encoding="utf-8", newline="") as stream:
        trace = list(csv.DictReader(stream))
    qpos = np.tile(model.qpos0, (len(trace), 1))
    qpos[:, shared_indices] = [[float(row[f"qpos_{i}"]) for i in donor_indices] for row in trace]
    for row in replay["samples"]:
        np.testing.assert_allclose(qpos[row["trace_frame"]], row["qpos"], rtol=0., atol=1e-11)
    relief = None
    if args.cap_relief:
        native_tpu = load_reference(provenance(folder / "whole_arm_k7_tpu_moving_cuff_audit.json"))
        if native_tpu["combined"] != combined_ref:
            raise ValueError("TPU audit uses another K7 assembly")
        extensor = load_component(native_tpu["tpu_sources"]["tpu:elbow_extensor"])
        original = parts[CAP]
        carriage = model.joint("FJ1477_series_travel")
        travel_axis = carriage.axis.copy()
        if not np.isclose(np.linalg.norm(travel_axis), 1.):
            raise ValueError("Expected a unit carriage axis")
        # Form the track of the middle-deltoid wire in the translating cap frame.
        # This probe is for a modest first elbow movement of 30 degrees.
        stop = next(i for i, angle in enumerate(np.rad2deg(qpos[:, elbow_index])) if angle >= 30.)
        max_travel = 1000. * np.max(qpos[:stop + 1, slide_index])
        min_travel = min(0., 1000. * np.min(qpos[:stop + 1, slide_index]))
        step = .20
        stroke = max_travel - min_travel
        offsets = np.linspace(min_travel, max_travel, int(np.ceil(stroke / step)) + 1)
        wire = A.D.outer_envelope(np.asarray(routes["FJ1467"]["points_source_mm"]), .52)
        cuts = [original ^ wire.translate(-offset * travel_axis) for offset in offsets]
        removed = union(cut for cut in cuts if cut.volume() > 1e-7)
        candidate = (original - removed).as_original().simplify(1e-7)
        wall = A.I._manifold().Manifold.cube([1.2] * 3, True)
        fixed_positions = {key: part for key, part in parts.items()
                           if owners[key] == "FJ1477_carriage" and key != CAP}
        fixed_positions["TPU:extensor_at_neutral"] = extensor
        removals_at_protection = {}
        for key, part in fixed_positions.items():
            protected = part if key.startswith("TPU:") else part.minkowski_sum(wall)
            removals_at_protection[key] = A.I.intersection_volume(removed, protected)
        for key in long_rails:
            rail = parts[key].minkowski_sum(wall)
            rail_min = rail.translate(-min_travel * travel_axis)
            rail_max = rail.translate(-max_travel * travel_axis)
            removals_at_protection[key] = A.I.intersection_volume(removed, (rail_min + rail_max).hull())
        components = [p.volume() for p in candidate.decompose() if p.volume() > 1e-5]
        accepted_geometry = bool(len(components) == 1 and removed.volume() > .001
                                 and max(removals_at_protection.values()) <= .001)
        relief = dict(original=sources[CAP], selected_source=None, trial_angle_deg=30.,
                      wire_tube_radius_mm=.52, travel_interval_mm=[min_travel, max_travel], swept_samples=len(offsets),
                      removed_mm3=removed.volume(), removed_bounds_source_mm=list(removed.bounding_box()),
                      retained_components_mm3=components, removed_at_protected_interfaces_mm3=removals_at_protection,
                      geometry_accepted=accepted_geometry)
        print("Cap relief trial:", relief, flush=True)
        if not accepted_geometry:
            save(folder / "whole_arm_k1_cap_relief.json", dict(
                baseline=provenance(folder / "whole_arm_k1_shared.json"), trial=relief,
                scope="Contact relief rejected by the connected-stock or interface guard check."))
            return
        origin = np.asarray(sources[CAP]["translation_source_mm"])
        _, relief["selected_source"] = export(candidate, folder / "FJ1477_shared_cap_wire_relief.obj", origin)
        parts[CAP] = load_component(relief["selected_source"])
        relief["added_after_roundtrip_mm3"] = max(0., (parts[CAP] - original).volume())
        relief["removed_after_roundtrip_mm3"] = (original - parts[CAP]).volume()
        if relief["added_after_roundtrip_mm3"] > .001 or not np.isclose(parts[CAP].volume(), candidate.volume(), atol=.001, rtol=0.):
            raise ValueError("Round-trip cap is not a material-only relief of the original")
    accepted, first = [], None
    tick = 0
    while tick <= replay["peak_trace_frame"]:
        row = k1(qpos[tick], tick)
        if row["hits"]:
            # Recheck each source pose back to the last accepted forward sample.
            for i in range(accepted[-1]["trace_frame"] + 1 if accepted else tick, tick):
                fine = k1(qpos[i], i)
                if fine["hits"]:
                    row = fine
                    break
                accepted.append(fine)
            first = row
            break
        accepted.append(row)
        next_tick = min(tick + 4, replay["peak_trace_frame"])
        if next_tick == tick:
            break
        # Bound spacing by elbow motion even during the fast portion of the trace.
        while next_tick > tick + 1 and abs(np.rad2deg(qpos[next_tick, elbow_index] - qpos[tick, elbow_index])) > 1.:
            next_tick -= 1
        tick = next_tick
    if accepted and first:
        left, right = np.asarray(accepted[-1]["qpos"]), np.asarray(first["qpos"])
        first_tick, last_tick = accepted[-1]["trace_frame"], first["trace_frame"]
        subdivisions = max(1, int(np.ceil(abs(first["elbow_deg"] - accepted[-1]["elbow_deg"]) / .1)))
        for fraction in np.arange(1, subdivisions) / subdivisions:
            fine = k1((1. - fraction) * left + fraction * right, float(first_tick + fraction * (last_tick - first_tick)))
            if fine["hits"]:
                first = fine
                break
            accepted.append(fine)
    for row in reversed(accepted):
        repeat = k1(np.asarray(row["qpos"]), row["trace_frame"])
        if repeat != row:
            raise ValueError(f"Return differs at trace frame {row['trace_frame']}")
    segments = {}
    angles = np.asarray([row["elbow_deg"] for row in accepted])
    for key in routes:
        lengths = np.asarray([row["length_mm"][key] for row in accepted])
        da = np.diff(angles)
        bending = da > .05
        gradients = -np.diff(lengths)[bending] / np.deg2rad(da[bending])
        if lengths.size:
            segments[key] = dict(net_mm=float(lengths[-1] - lengths[0]), excursion_mm=float(np.ptp(lengths)),
                                 effective_moment_arm_mm=[float(gradients.min()), float(gradients.max())]
                                 if gradients.size else None)
    print(f"K1 current rails: {len(accepted)} poses; first hit {first['hits'][0] if first else None}", flush=True)

    cuff_ref = provenance(folder / "whole_arm_k7_moving_cuff_geometry.json")
    cuff = load_reference(cuff_ref)
    if cuff["combined"] != combined_ref:
        raise ValueError("K7 cuff poses belong to a different model")
    cuff_parts = {key: load_component(entry) for key, entry in cuff["moving_sources"].items()}
    # Isolate what the long rails add to the existing K7 result.
    rail_added = {key: shape - parts[key] for key, shape in long_rails.items()}
    rail_compatibility = []
    cap_compatibility = []
    for pose, reference, cuff_frame in zip(motion["paths"]["K7_SC1_pos"]["accepted"][:len(combined["combined"])],
                                           combined["combined"], cuff["frames"], strict=True):
        if pose["input_deg"] != reference["input_deg"] or pose["input_deg"] != cuff_frame["input_deg"]:
            raise ValueError("K7 and cuff trajectories are out of phase")
        data.qpos[:] = pose["qpos"]
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        parts_at_pose = [parts[key].transform(frames[owners[key]]) for key in parts]
        parts_at_pose.extend(section(route, frames)[0] for route in routes.values())
        high = next(row for row in reference["steps"] if row["load_N"] == 10.)
        chords = high["central_tendon_tips_world_mm"]
        parts_at_pose.extend(part.transform(np.asarray(cuff_frame["matrices"][key])) for key, part in cuff_parts.items())
        parts_at_pose.extend(rod(*np.asarray(points), .30) for points in chords.values())
        names = [*keys, *cuff_parts, *[f"chord:{key}" for key in chords]]
        original_count = len(names)
        if relief:
            cap_index = names.index(CAP)
            comparison = np.asarray([(i, cap_index) for i, key in enumerate(names) if i != cap_index
                                    and owners.get(key) != owners[CAP]])
            hits = find_hits(parts_at_pose, names, comparison)
            minimum, pair = closest_gap(parts_at_pose, comparison, cap=3.)
            cap_compatibility.append(dict(input_deg=pose["input_deg"], hits=hits, minimum_gap_capped_3mm=minimum,
                                          nearest_pair=[names[i] for i in pair] if pair else None))
        names.extend(f"added:{key}" for key in rail_added)
        parts_at_pose.extend(shape.transform(frames[owners[key]]) for key, shape in rail_added.items())
        rail_pairs = np.asarray([(i, j) for j in range(original_count, len(names)) for i in range(j)
                                 if names[i] not in long_rails])
        hits = find_hits(parts_at_pose, names, rail_pairs)
        rail_compatibility.append(dict(input_deg=pose["input_deg"], humerus_deg=pose["humerus_direction_change_deg"], hits=hits))
    peak = max(accepted, key=lambda row: row["elbow_deg"], default=None)
    target = "whole_arm_k1_cap_relief.json" if args.cap_relief else "whole_arm_k1_shared.json"
    save(folder / target, dict(
        combined=combined_ref, base=base_ref, elbow_replay=replay_ref, assembly=assembly_ref, cuff=cuff_ref,
        baseline=provenance(folder / "whole_arm_k1_shared.json") if args.cap_relief else None,
        generator=provenance(ROOT / "tools/audit_whole_arm_k1_shared.py"),
        accepted=accepted, first_limit=first, checked_return=len(accepted), focal_peak=peak,
        local_lines=segments, K7_additional_rail_material=rail_compatibility, cap_relief=relief,
        K7_relieved_cap=cap_compatibility,
        scope="Prescribed K1 elbow, cap swivels and carriage coordinates with the same fixed mesh list as K7. "
              "Original rails, eight fixed-end wires, no spring, moving cuff or deformed TPU. K1 return is the reverse "
              "of the checked outward path. Source spacing may exceed one degree. The final bracket is interpolated "
              "between neighboring trace rows at no more than 0.1 degree; no interpolated dynamics. "
              "Effective moment arm follows the coupled elbow/carriage trajectory, not the isolated elbow moment. "
              "K7 long-rail trial checks added material against rigid geometry, eight wires and the 10 N cuff pose and chord. "
              "TPU and assembly under load remain open."))


if __name__ == "__main__":
    main()
