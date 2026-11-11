"""Move FJ1500 fastener engagement to the free side of the scapula."""

from itertools import combinations
from pathlib import Path

import numpy as np

from audit_s1_series_elbow_integration import MOBILE, assert_source, find_hits
from audit_elbow_slider_short_guides import export_lattice
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    parent_path = folder / "assembly_FJ1504_local_interface.json"
    parent = read(parent_path)
    revision = read(assert_source(parent["scapular_mount_revision"]))
    routes = read(assert_source(revision["route_reference"]))["routes"]
    route = next(item for item in routes if item["head"] == "FJ1500")
    holes = next(item["holes"] for item in read(assert_source(revision["fixation_geometry"]))["rows"]
                 if item["head"] == "FJ1500")
    port = route["ports"]["origin"]
    normal = np.asarray(port["normal_source"])
    inner = np.asarray(port["inner_face_source_mm"])
    cap_thickness = float((np.asarray(port["outer_face_source_mm"]) - inner) @ normal)
    points = [np.asarray(hole["axis_point_source_mm"]) for hole in holes]
    sources = parent["sources"]
    owners = parent["owners"]
    prefix = "origin:FJ1500/"
    replaced = {prefix + key for key in ("front_support", "backing_plate", "bolt_0", "bolt_1",
                                         "nut_0", "back_washer_0", "front_nut_1")}
    reference_bone = load_component(sources["bone:scapula"])
    support = load_component(sources[prefix + "front_support"])
    tpu = load_component(parent["flexible_geometry_sources"]["tpu:infraspinatus"])
    fixed = {key: load_component(row) for key, row in sources.items()
             if key not in replaced and key != "bone:scapula"}
    baseline = read(folder / "scapular_thorax_relief_audit.json")
    path = np.asarray(route["full_wire_path_source_mm"])
    wire = A.D.outer_envelope(path, .3)
    fixed_bounds = bounded(fixed)
    trials = []
    chosen = None
    for radius in (3.4, 3., 2.6, 2.2):
        for floor in (-6.7, -7.1, -7.5):
            bosses = union(rod(point + floor * normal, point - 2.5 * normal, radius) for point in points)
            # Threads are nominal envelopes. Drill stops before the existing
            # posterior face; no new fastener or seat projects toward the rib.
            drills = union(rod(point + 3. * normal, point + (floor + .25) * normal, 1.25)
                           for point in points)
            stock = union([reference_bone, support, bosses])
            body = (stock - drills).as_original().simplify(1e-6)
            added = body - reference_bone
            body_hits = overlaps(added, fixed_bounds)
            soft_overlap = float((added ^ tpu).volume())
            wire_overlap = float((added ^ wire).volume())
            components = [float(part.volume()) for part in body.decompose() if part.volume() > 1e-8]
            row = dict(boss_radius_mm=radius, boss_floor_from_cap_mm=floor, fixed_hits_mm3=body_hits,
                       soft_overlap_mm3=soft_overlap, wire_overlap_mm3=wire_overlap, components=components)
            trials.append(row)
            if (len(components) != 1 or body_hits or soft_overlap > .001 or wire_overlap > .001):
                continue
            bolt_tip = floor + .45
            bolts = {}
            for index, point in enumerate(points):
                head_start = cap_thickness + .65
                shaft = rod(point + bolt_tip * normal, point + (head_start + .15) * normal, 1.19)
                head = rod(point + head_start * normal, point + (head_start + 2.) * normal, 2.75)
                bolts[prefix + f"bolt_{index}"] = shaft + head
            bolt_obstacles = {**fixed_bounds, **bounded({"scapula_with_bosses": body, "TPU": tpu, "wire": wire})}
            bolt_hits = {key: hits for key, part in bolts.items() if (hits := overlaps(part, bolt_obstacles))}
            if bolt_hits:
                row["bolt_hits_mm3"] = bolt_hits
                continue
            part_sources = {}
            try:
                for key, part in {"bone:scapula": body, **bolts}.items():
                    writer = export_lattice if key == "bone:scapula" else export
                    _, part_sources[key] = writer(
                        part, folder / f"FJ1500_front_{key.split('/')[-1].replace(':', '_')}.obj", inner)
            except (AssertionError, ValueError) as error:
                row["export_error"] = str(error) or type(error).__name__
                continue
            chosen = dict(**row, sources=part_sources, screw_tip_mm=bolt_tip,
                          nominal_engagement_length_mm=-2.5 - bolt_tip)
            break
        if chosen:
            break
    report = dict(generator=provenance(Path(__file__).resolve()), parent=provenance(parent_path),
                  trials=trials, selected=chosen, superseded_FJ1500_parts=sorted(replaced))
    if not chosen:
        save(folder / "FJ1500_front_mount.json", report)
        print("No front mount passed", trials, flush=True)
        return

    final_sources = {**{key: row for key, row in sources.items() if key not in replaced}, **chosen["sources"]}
    final_owners = {key: owners[key] for key in final_sources}
    solids = {key: load_component(row) for key, row in final_sources.items()}
    keys = list(solids)
    revised = set(chosen["sources"])
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                        if final_owners[keys[i]] != final_owners[keys[j]] or revised.intersection((keys[i], keys[j]))])
    neutral = find_hits(list(solids.values()), keys, pairs)
    centres = read(DEST / "elbow_passive_model/model_check.json")["neutral_origins_m"]
    tests = []
    for label, pivot, moving in (("SC", centres["clavicle"], MOBILE | {"clavicle", "scapula"}),
                                  ("AC", centres["scapula"], MOBILE | {"scapula"})):
        moving_keys = {key for key in keys if final_owners[key] in moving}
        cross = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                            if (keys[i] in moving_keys) != (keys[j] in moving_keys)])
        for axis in np.eye(3):
            for sign in (-1., 1.):
                def hits(angle):
                    matrix = A.I.rotation_about(np.asarray(pivot) * 1000., axis, angle * sign)
                    placed = [solid.transform(matrix) if key in moving_keys else solid for key, solid in solids.items()]
                    return find_hits(placed, keys, cross)

                last, first = 0., None
                for angle in range(1, 11):
                    contact = hits(float(angle))
                    if contact:
                        low, high = last, float(angle)
                        while high - low > .1:
                            middle = (low + high) / 2.
                            probe = hits(middle)
                            if probe:
                                high, contact = middle, probe
                            else:
                                low = middle
                        last = low
                        first = dict(angle_deg=sign * high, pairs=contact)
                        break
                    last = float(angle)
                prior = next(row for row in baseline["shoulder_pose_probes"] if row["joint"] == label
                             and row["direction"] == sign and row["axis_source"] == axis.tolist())
                tests.append(dict(joint=label, axis_source=axis.tolist(), direction=sign, first_hit=first,
                                  last_sample_clear_deg=sign * last,
                                  prior_last_sample_clear_deg=prior["global_last_sample_clear_deg"]))
    scapula = solids["bone:scapula"]
    report.update(neutral_overlaps=neutral, neutral_pairs_checked=len(pairs), shoulder_axis_probes=tests,
                  added_stock_mm3=float((scapula - reference_bone).volume()),
                  removed_reference_bone_mm3=float((reference_bone - scapula).volume()),
                  readback_new_bone_TPU_overlap_mm3=float(((scapula - reference_bone) ^ tpu).volume()),
                  readback_new_bone_wire_overlap_mm3=float(((scapula - reference_bone) ^ wire).volume()),
                  scope="Integrated scapula and FJ1500 front saddle, nominal blind bore and screw core. "
                  "Neutral and individual SC/AC angle probes, not motion acceptance. Screw strength, "
                  "thread retention, assembly and flexible shoulder response remain open.")
    save(folder / "FJ1500_front_mount.json", report)
    if not neutral:
        save(folder / "assembly_FJ1500_front_mount.json", dict(
            **{key: value for key, value in parent.items() if key not in ("sources", "owners", "scope")},
            sources=final_sources, owners=final_owners,
            FJ1500_mount_revision=provenance(folder / "FJ1500_front_mount.json"),
            scope="Neutral assembly candidate with same-side FJ1500 engagement and revised FJ1504 origin."))
    print("neutral hits", neutral, "motion probes", [(r["joint"], r["axis_source"], r["direction"],
                                                     r["last_sample_clear_deg"]) for r in tests], flush=True)


if __name__ == "__main__":
    main()
