"""Transfer the existing local wire polylines to the shared anatomical arm."""

import argparse
import csv
from itertools import combinations

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import weigh
from audit_s1_infraspinatus_distributed_deformation import blend
from audit_s1_series_elbow_integration import find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save, rod
from project_runtime import load_mujoco_model
from study_elbow_extensor_long_stroke import turn_radius
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def nominal_rail_hits(qpos, model, data, origins, owners, nominal_solids, parts, changed):
    data.qpos[:] = qpos
    mujoco.mj_forward(model, data)
    frames = transforms(model, data, origins)
    changed_targets = {key: nominal_solids[key].transform(frames[owners[key]]) for key in sorted(changed)}
    hits = []
    for key, target in changed_targets.items():
        for other in parts:
            if owners[other] == owners[key]:
                continue
            part = nominal_solids[other].transform(frames[owners[other]])
            volume = A.I.intersection_volume(target, part)
            if volume > .001:
                hits.append(dict(left=key, right=other, overlap_mm3=volume))
    return hits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--localize", action="store_true",
                        help="Locate first collisions using the saved route audit")
    mode.add_argument("--passive-elbow", action="store_true",
                      help="Audit historical elbow swivel and carriage coordinates on the shared assembly")
    mode.add_argument("--documented-deltoid", action="store_true",
                      help="Use the deltoid guide motion from related local studies on the shared K7 path")
    mode.add_argument("--posterior-detour", action="store_true",
                      help="Screen a small posterior-bundle offset on the shared K7 path")
    args = parser.parse_args()
    folder = DEST.parent / "whole_arm_first_pass"
    scene_ref = provenance(folder / "whole_arm_shared_scene_check.json")
    motion_ref = provenance(folder / "whole_arm_shared_motion.json")
    scene, motion = map(load_reference, (scene_ref, motion_ref))
    if motion["model"] != scene["model"]:
        raise ValueError("Motion and geometry refer to different shared models")
    transfer = load_reference(scene["transfer"])
    register = load_reference(transfer["s1_elbow_register"])
    inventory = load_reference(register["inventory"])
    items = {row["id"]: row for row in inventory["routes"]}
    model_path = ROOT / scene["model"]["file"]
    if provenance(model_path) != scene["model"]:
        raise ValueError("Shared model changed after motion audit")
    model = load_mujoco_model(model_path)
    data = mujoco.MjData(model)
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"]}
    owners = transfer["owner_by_source"]
    solids = {key: load_component(source) for key, source in sources.items()}
    origins = {key: np.asarray(scene["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    passive_ref = provenance(DEST / "elbow_passive_model/model_check.json")
    elbow = load_reference(load_reference(passive_ref)["assembly"])
    routes, unresolved = {}, {}
    for row in register["channels"]:
        key = row["id"]
        first, last = row["entry_body"], row["terminal_body"]
        if first == "floating TPU terminal":
            unresolved[key] = dict(reason="Floating proximal terminal requires a supported tensioned pose",
                                   endpoints_source_mm=[row["entry_xyz_mm"], row["terminal_xyz_mm"]])
            continue
        if row["region"] == "elbow":
            specification = next(r for r in elbow["muscles"].values() if r["anatomical_id"] == key)
            definition = specification["source_definition"]
            definition = load_reference(definition) if isinstance(definition, dict) else read(ROOT / definition)
            path = np.asarray(specification["wire_path_source_mm"])[1:-1]
            axis = np.asarray(definition["basis_columns_display"])[:, 2]
            weights = weigh(path, path[0], axis / ((path[-1] - path[0]) @ axis))
            source = passive_ref
        else:
            item = items[key]
            source = row["route_reference"]
            candidate = load_reference(source)
            if key == "FJ1468":
                path = np.asarray(candidate["neutral_path_source_mm"])
            elif key == "FJ1467":
                path = np.asarray(candidate["neutral_candidate_mm"])
            else:
                path = np.asarray(item["local_route"]["points_source_mm"])
            chord = path[-1] - path[0]
            length = np.linalg.norm(chord)
            weights, _ = blend((path - path[0]) @ (chord / length), 6., length - 6., .1)
        routes[key] = dict(owners=(first, last), source=source, points=path, weights=weights,
                           type="polyline with prescribed attachment blending")
    if args.documented_deltoid or args.posterior_detour:
        front = routes["FJ1468"]
        guide_motion = load_reference(load_reference(front["source"])["guide_motion"])
        front_points = front["points"]
        axis = front_points[-1] - front_points[0]
        axis /= np.linalg.norm(axis)
        window = guide_motion["deformation_window_mm"]
        for key in ("FJ1468", "FJ1513", "FJ1467"):
            route = routes[key]
            path = route["points"]
            weights, _ = blend((path - front_points[0]) @ axis, *window, .1)
            # The front bundle's guide field does not start at the middle and
            # posterior terminals. Transition from their fixed terminal zones.
            station = np.r_[0., np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
            entrance = np.clip((station - 3.) / 12., 0., 1.)
            weights *= entrance**2 * (3. - 2. * entrance)
            if not np.allclose(weights[[0, -1]], [0., 1.], rtol=0, atol=1e-10):
                raise ValueError(f"Documented deltoid guide motion drifts from the terminal ports: {key}")
            route["weights"] = weights
            route["type"] = "historical deltoid guide field with fixed proximal terminal and tapered transition"
            route["guide_motion_reference"] = load_reference(front["source"])["guide_motion"]
    seats_ref = provenance(folder / "girdle_seat_pair.json")
    seats = load_reference(seats_ref)["seats"]
    for key in ("FJ1459_g3", "FJ1536_g2"):
        ends = [seats[f"{key}_{role}"] for role in ("origin", "insertion")]
        routes[key] = dict(owners=[end["owner"] for end in ends], source=seats_ref,
                           points=np.asarray([end["port_source_mm"] for end in ends]),
                           weights=np.array([0., 1.]), type="straight through two existing seat eyes")

    def pose_wire(row, selected):
        data.qpos[:] = row["qpos"]
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        obstacles = bounded({key: part.transform(frames[owners[key]]) for key, part in solids.items()})
        paths = {}
        for key in selected:
            spec = routes[key]
            a, b = (frames[owner] for owner in spec["owners"])
            path = spec["points"]
            weights = spec["weights"][:, None]
            posed = (1. - weights) * (path @ a[:, :3].T + a[:, 3]) + weights * (path @ b[:, :3].T + b[:, 3])
            paths[key] = posed
        return obstacles, paths

    def probe(row, selected):
        obstacles, paths = pose_wire(row, selected)
        wires, records = {}, {}
        for key, posed in paths.items():
            distance = np.linalg.norm(np.diff(posed, axis=0), axis=1)
            if min(distance) < 1e-8:
                raise ValueError(f"Collapsed wire segment in {key}")
            wire = rod(*posed, .3) if len(posed) == 2 else A.D.outer_envelope(posed, .3)
            wires[key] = wire
            records[key] = dict(length_mm=distance.sum().item(), rigid_hits_mm3=overlaps(wire, obstacles),
                                turn_radius_mm=turn_radius(posed) if len(posed) > 2 else None)
        pair_hits = []
        for (left, a), (right, b) in combinations(wires.items(), 2):
            volume = A.I.intersection_volume(a, b)
            if volume > .001:
                pair_hits.append(dict(left=left, right=right, mm3=volume))
        return dict(input_deg=row["input_deg"], routes=records, neighboring_wire_hits=pair_hits)

    detour = None
    if args.posterior_detour:
        route = routes["FJ1467"]
        original = route["points"].copy()
        first, last = 66, 75
        section = slice(first, last + 1)
        tangent = original[last] - original[first]
        tangent /= np.linalg.norm(tangent)
        side = np.cross(tangent, np.eye(3)[np.argmin(abs(tangent))])
        side /= np.linalg.norm(side)
        up = np.cross(tangent, side)
        station = np.r_[0., np.cumsum(np.linalg.norm(np.diff(original[section], axis=0), axis=1))]
        fraction = station / station[-1]
        taper = 16. * fraction**2 * (1. - fraction)**2
        seeds = []
        for row in motion["paths"]["K7_SC1_pos"]["accepted"][::5]:
            targets, _ = pose_wire(row, ())
            frames = transforms(model, data, origins)
            seeds.append((targets, [frames[owner] for owner in route["owners"]]))
        proposals = []
        for radius in (1., 2., 3., 4.):
            for direction in np.arange(8) * np.pi / 4.:
                offset = radius * (np.cos(direction) * side + np.sin(direction) * up)
                points = original[section] + taper[:, None] * offset
                weight = route["weights"][section, None]
                blocked = []
                for targets, (a, b) in seeds:
                    moved = (1. - weight) * (points @ a[:, :3].T + a[:, 3]) + weight * (points @ b[:, :3].T + b[:, 3])
                    blocked.append(overlaps(A.D.outer_envelope(moved, .3), targets))
                option = dict(offset_source_mm=offset.tolist(), blocked_seed_count=sum(map(bool, blocked)),
                              peak_seed_overlap_mm3=max(map(lambda hits: sum(hits.values()), blocked)))
                proposals.append((option, points))
            if any(not row["blocked_seed_count"] for row, _ in proposals):
                break
        chosen, points = min(proposals, key=lambda item: (item[0]["blocked_seed_count"],
                                                         item[0]["peak_seed_overlap_mm3"]))
        route["points"][section] = points
        detour = dict(route="FJ1467", edited_node_indices=[first, last], selected=chosen,
                      trials=[option for option, _ in proposals],
                      support_status="Local channel or guide redesign required")
        print("Posterior candidate", chosen, flush=True)

    if args.localize:
        audit_ref = provenance(folder / "whole_arm_local_wire_audit.json")
        audit = load_reference(audit_ref)
        if audit["motion"] != motion_ref:
            raise ValueError("Saved wire audit refers to another motion")
        localized = {}
        for name, result in audit["paths"].items():
            entries = []
            for key, summary in result["summary"].items():
                first = summary["first_limit"]
                if first is None or not first["rigid_hits_mm3"]:
                    continue
                row = next(row for row in motion["paths"][name]["accepted"]
                           if row["input_deg"] == first["input_deg"])
                obstacles, paths = pose_wire(row, [key])
                posed = paths[key]
                cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(posed, axis=0), axis=1))]
                mf = A.I._manifold()
                for obstacle in first["rigid_hits_mm3"]:
                    target = {obstacle: obstacles[obstacle]}
                    spans = []
                    for index, (a, b) in enumerate(zip(posed, posed[1:])):
                        shape = rod(a, b, .3) + mf.Manifold.sphere(.3, 32).translate(a)
                        shape += mf.Manifold.sphere(.3, 32).translate(b)
                        hits = overlaps(shape, target, threshold=1e-6)
                        if hits:
                            spans.append(dict(segment=index, arc_start_mm=float(cumulative[index]),
                                              arc_end_mm=float(cumulative[index + 1]), overlap_mm3=hits[obstacle],
                                              distal_weights=routes[key]["weights"][index:index + 2].tolist(),
                                              source_endpoints_mm=routes[key]["points"][index:index + 2].tolist()))
                    entries.append(dict(route=key, input_deg=row["input_deg"], obstacle=obstacle,
                                        obstacle_owner=owners[obstacle], segments=spans))
            localized[name] = entries
        save(folder / "whole_arm_local_wire_limits.json", dict(audit=audit_ref, motion=motion_ref, cases=localized))
        print(f"Localized {sum(map(len, localized.values()))} route–obstacle limits")
        return

    if args.passive_elbow:
        donor_file = DEST.parent / "integration_review/s1_series_elbow/nominal_spring_spool_seating_k5_coupling.json"
        donor_ref = provenance(donor_file)
        donor_report = load_reference(donor_ref)
        donor_xml = ROOT / donor_report["model"]["file"]
        if provenance(donor_xml) != donor_report["model"]:
            raise ValueError("The elbow donor model has changed")
        donor = load_mujoco_model(donor_xml)
        joints = ("elbow", "FJ1486_swivel", "FJ1477_swivel", "FJ1477_series_travel")
        for name in joints:
            np.testing.assert_allclose(model.joint(name).axis, donor.joint(name).axis, atol=1e-10)
        donor_indices = [donor.joint(name).qposadr.item() for name in joints]
        shared_indices = [model.joint(name).qposadr.item() for name in joints]
        trace_path = ROOT / donor_report["trace"]["file"]
        if provenance(trace_path) != donor_report["trace"]:
            raise ValueError("The elbow donor trace has changed")
        with trace_path.open(encoding="utf-8", newline="") as stream:
            trace = list(csv.DictReader(stream))

        base_ref = provenance(folder / "tabletop_base_check.json")
        base = load_reference(base_ref)
        if base["scene"] != scene_ref:
            raise ValueError("The tabletop base belongs to another assembly")
        solids["base:T12"] = load_component(base["source"])
        owners = {**owners, "base:T12": "thorax"}
        rail_ref = provenance(donor_file.parent / "elbow_extended_rail_free_drive.json")
        rails = load_reference(rail_ref)
        if rails["model"] != donor_report["model"]:
            raise ValueError("Extended rails and donated coordinates have different model sources")
        nominal_solids = dict(solids)
        solids.update({key: load_component(spec) for key, spec in rails["shafts"].items()})
        parts = list(solids)
        active = {"forearm", "FJ1486_distal_pivot", "FJ1477_distal_pivot", "FJ1477_carriage"}
        pairs = np.asarray([(a, b) for a, b in combinations(range(len(parts)), 2)
                            if owners[parts[a]] != owners[parts[b]]
                            and (owners[parts[a]] in active or owners[parts[b]] in active)])
        selected = ("FJ1486", "FJ1477")
        peak = max(range(len(trace)), key=lambda i: float(trace[i]["elbow_deg"]))
        carriage_max = max(range(len(trace)), key=lambda i: float(trace[i]["slide_mm"]))
        indices = sorted(set(range(0, len(trace), 4)) | {peak, carriage_max, len(trace) - 1})
        rest = data.qpos.copy()
        samples = []
        changed = set(rails["shafts"])
        neutral_check_pairs = np.asarray([(a, b) for a, b in combinations(range(len(parts)), 2)
                                          if parts[a] in changed or parts[b] in changed])
        static_overlaps = find_hits(list(solids.values()), parts, neutral_check_pairs)
        for index in indices:
            frame = trace[index]
            values = np.array([float(frame[f"qpos_{i}"]) for i in donor_indices])
            qpos = rest.copy()
            qpos[shared_indices] = values
            result = probe(dict(qpos=qpos, input_deg=float(np.rad2deg(values[0]))), selected)
            frames = transforms(model, data, origins)
            placed = [solids[key].transform(frames[owners[key]]) for key in parts]
            result.update(trace_frame=index, time_s=float(frame["time_s"]), qpos=qpos.tolist(),
                          swiveled_deg=np.rad2deg(values[1:3]).tolist(), carriage_mm=values[3] * 1000.,
                          rigid_hits=find_hits(placed, parts, pairs))
            samples.append(result)
        def blocked(row):
            return bool(row["rigid_hits"] or row["neighboring_wire_hits"]
                        or any(record["rigid_hits_mm3"] for record in row["routes"].values()))

        first = next((index for index, row in enumerate(samples) if blocked(row)), len(samples))
        accepted = samples[:first]
        save(folder / "whole_arm_passive_elbow_audit.json", dict(
            scene=scene_ref, donor=donor_ref, base=base_ref, peak_trace_frame=peak,
            extended_rails=rail_ref, nominal_scene_rail_sources={key: sources[key] for key in changed},
            replacement_neutral_overlaps=static_overlaps,
            nominal_rail_first_failed_sample=next((
                dict(trace_frame=row["trace_frame"], input_deg=row["input_deg"], hits=hits)
                for row in samples
                if (hits := nominal_rail_hits(row["qpos"], model, data, origins, owners,
                                             nominal_solids, parts, changed))
            ), None),
            sampled_count=len(samples), checked_rigid_pairs=len(pairs),
            contiguous_samples=len(accepted), first_limit=samples[first] if first < len(samples) else None,
            max_initial_clear_elbow_deg=max((row["input_deg"] for row in accepted), default=None),
            samples=samples,
            scope="Donor elbow, swivel and slider coordinates prescribed on the shared arm with matching extended shafts and "
                  "SC/AC/GH fixed at neutral. Wire paths use geometric interpolation. No new force "
                  "integration, shoulder wires, TPU contact, or spring clearance included. "
                  "The donor phase and mounting offsets do not transfer drive validation."))
        print(f"Elbow transfer: {len(accepted)}/{len(samples)} initial clear samples; "
              f"first limit {samples[first]['input_deg'] if first < len(samples) else None}")
        return

    tests = {
        "K7_SC1_pos": list(routes),
        "K1_elbow_fixed_passive_interfaces": ["FJ1486", "FJ1477"],
        "K2_GH_flexion": ["FJ1468", "FJ1467", "FJ1513", "FJ1507"],
        "K3_GH_abduction": ["FJ1468", "FJ1467", "FJ1513", "FJ1507"],
        "GH_axial_positive": ["FJ1468", "FJ1467", "FJ1513", "FJ1507"],
        "GH_axial_negative": ["FJ1468", "FJ1467", "FJ1513", "FJ1507"],
    }
    if args.documented_deltoid or args.posterior_detour:
        tests = {"K7_SC1_pos": tests["K7_SC1_pos"]}
    runs = {}
    for name, selected in tests.items():
        path = motion["paths"][name]["accepted"]
        if name != "K7_SC1_pos":
            path = [*path[::5], *([] if (len(path) - 1) % 5 == 0 else [path[-1]])]
        samples = [probe(row, selected) for row in path]
        summary = {}
        for key in selected:
            clear = []
            first = None
            for row in samples:
                record = row["routes"][key]
                crossings = [hit for hit in row["neighboring_wire_hits"] if key in (hit["left"], hit["right"])]
                if crossings or record["rigid_hits_mm3"]:
                    first = dict(input_deg=row["input_deg"], **record, neighboring_wire_hits=crossings)
                    break
                clear.append(row)
            lengths = [row["routes"][key]["length_mm"] for row in samples]
            derivative = -np.diff(lengths) / np.diff(np.deg2rad([row["input_deg"] for row in samples]))
            summary[key] = dict(first_limit=first, initial_clear_samples=len(clear),
                                final_clear_input_deg=clear[-1]["input_deg"] if clear else None,
                                nominal_excursion_mm=float(np.ptp(lengths)),
                                nominal_leverage_along_input_mm=[float(derivative.min()), float(derivative.max())]
                                if derivative.size else None)
        runs[name] = dict(source_pose_count=len(path), summary=summary, samples=samples)
        print(name, len(path), "poses", {key: row["final_clear_input_deg"] for key, row in summary.items()}, flush=True)
    filename = ("whole_arm_posterior_detour_audit.json" if args.posterior_detour else
                "whole_arm_documented_deltoid_audit.json" if args.documented_deltoid else
                "whole_arm_local_wire_audit.json")
    save(folder / filename, dict(
        scene=scene_ref, motion=motion_ref, anchor_register=transfer["s1_elbow_register"],
        generator=provenance(ROOT / "tools/audit_whole_arm_local_wires.py"),
        routes={key: dict(**{k: v for k, v in spec.items() if k not in ("points", "weights")},
                          points_source_mm=spec["points"].tolist(), distal_weights=spec["weights"].tolist())
                for key, spec in routes.items()},
        unresolved=unresolved, paths=runs, posterior_detour=detour,
        scope="Geometric transfer of selected local wire paths. 0.60 mm wire envelopes. "
              "The combined path uses all saved samples, single-axis probes use every fifth saved pose "
              "and the endpoint to diagnose limits. Reversed prescribed poses reproduce the geometry. "
              "Nominal excursions and derivatives may include blocked poses and do not establish usable stroke. "
              "No cable equilibrium, TPU contacts, external supply, or force acceptance."))


if __name__ == "__main__":
    main()
