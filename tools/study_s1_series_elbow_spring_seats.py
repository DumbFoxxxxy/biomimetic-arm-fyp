"""Test local spring pockets while preserving the transverse pin seats."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_series_elbow_return_spring import spring_path, sweep
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    report_path = folder / "return_spring_geometry.json"
    report = read(report_path)
    size = load_reference(report["size"])
    audit = load_reference(size["geometry"])
    extended = load_reference(audit["extended_drive"])
    assembly = load_reference(load_reference(extended["budget"])["assembly"])
    mount = load_reference(size["spring_mount"])["spring_mount"]
    start = np.array(mount["fixed_eye_source_mm"])
    end = np.array(mount["moving_eye_neutral_source_mm"])
    axis = A.D.unit(end-start)
    x = np.array(mount["pin_axis_source"])
    y = np.cross(axis, x)
    frame = np.c_[x, y, axis, start]
    rest = float(np.linalg.norm(end-start))
    stroke = size["full_rail_stroke_mm"]
    design = report["design"]
    radius = .5 * design["wire_diameter_mm"]
    mean = .5 * (design["coil_OD_mm"] - design["wire_diameter_mm"])
    turns = design["active_turns"] + 2
    mf = A.I._manifold()
    keys = ("elbow:FJ1477_bone_adapter", "elbow:FJ1477_proximal_rigid_cap")
    pins = ("series:fixed_pin", "series:moving_pin")
    parts = {k: load_component(v) for k, v in {**assembly["sources"], **extended["shafts"]}.items()
             if assembly["owners"][k] in ("humerus", "FJ1477_carriage")}
    before = {k: parts[k] for k in keys}
    attachment_guards = {}
    for key, point in zip(keys, (start, end), strict=True):
        # The lug starts at pin-axis station +1.3 mm. Keep its full bearing
        # cylinder plus 0.20 mm beyond the loaded face and outer wall.
        lug = rod(point + 1.10*x, point + 4.2*x, 2.60)
        attachment_guards[key] = parts[key] ^ lug
    pad = mf.Manifold.cube([.30] * 3, True)
    travel = np.linspace(0., stroke, 21)
    springs = [sweep(spring_path(rest+t, turns, radius, mean), radius).transform(frame) for t in travel]
    targets = {key: [] for key in keys}
    # Reverse carriage travel on its spring contacts before cutting its seat.
    for t, spring in zip(travel, springs, strict=True):
        for i, key in enumerate(keys):
            local = spring.translate(-i*t*axis)
            intersection = local ^ before[key]
            targets[key].extend(part.hull().minkowski_sum(pad) for part in intersection.decompose()
                                if part.volume() > 1e-5)
    cuts = {key: union(patches).as_original().simplify(1e-6) for key, patches in targets.items()}
    rest_spring = load_component(report["exports"]["neutral"]["mesh"])
    base_options = {}
    for key in keys:
        base_options[key] = (before[key] - cuts[key]).as_original().simplify(1e-6)
    clearance = []
    for direction, vector in (("negative_pin_axis", -x), ("positive_pin_axis", x),
                              ("negative_side", -y), ("positive_side", y),
                              ("negative_side_diagonal", A.D.unit(-x-y)),
                              ("positive_side_diagonal", A.D.unit(-x+y))):
        obstacles = bounded({**{k: v for k, v in parts.items() if k not in pins}, **base_options})
        hits = []
        for offset in np.linspace(6., 0., 25):
            spring = rest_spring.translate(offset*vector)
            contact = overlaps(spring, obstacles)
            if contact:
                hits.append(dict(offset_mm=float(offset), hits_mm3=contact))
        clearance.append(dict(direction=direction, vector=vector.tolist(), sampled_count=25, hits=hits))
    best = min(clearance, key=lambda row: sum(sum(item["hits_mm3"].values()) for item in row["hits"]))
    direction = np.array(best["vector"])
    for key in keys:
        pocket = base_options[key]
        extra = []
        for offset in np.linspace(0., 6., 49):
            hit = rest_spring.translate(offset*direction) ^ pocket
            extra.extend(part.hull().minkowski_sum(pad) for part in hit.decompose() if part.volume() > 1e-5)
        if extra:
            base_options[key] = (pocket - union(extra)).as_original().simplify(1e-6)
    geometry = {}
    accepted = {}
    for key, component in base_options.items():
        removed = float((before[key] - component).volume())
        protected_loss = float((attachment_guards[key] - component).volume())
        connected = sum(part.volume() > 1e-5 for part in component.decompose()) == 1
        geometry[key] = dict(removed_mm3=removed, pin_seat_protected_loss_mm3=protected_loss, connected=connected)
        if protected_loss < .001 and connected:
            accepted[key], meta = export(component, folder / f"return_spring_pocket_{key.split('_')[-1]}.obj",
                                         start if key == keys[0] else end)
            geometry[key]["candidate_mesh"] = meta
    stroke_results = []
    insertion_hits = []
    pin_insertion = {}
    if len(accepted) == len(keys):
        for t, spring in zip(travel, springs, strict=True):
            placed = {k: v.translate(t*axis) if assembly["owners"][k] == "FJ1477_carriage" else v
                      for k, v in accepted.items()}
            stroke_results.append(dict(travel_mm=float(t), hits_mm3=overlaps(spring, bounded(placed))))
        surroundings = {**parts, **accepted}
        mount_without_pins = bounded({k: v for k, v in surroundings.items() if k not in pins})
        for offset in np.linspace(6., 0., 97):
            hit = overlaps(rest_spring.translate(offset*direction), mount_without_pins)
            if hit:
                insertion_hits.append(dict(offset_mm=float(offset), hits_mm3=hit))
        for key in pins:
            fixed = bounded({**{k: v for k, v in surroundings.items() if k != key}, "spring": rest_spring})
            pin_insertion[key] = []
            for sign in (-1., 1.):
                failures = []
                for offset in np.linspace(6., 0., 49):
                    hit = overlaps(parts[key].translate(sign*offset*x), fixed)
                    if hit:
                        failures.append(dict(offset_mm=float(offset), hits_mm3=hit))
                pin_insertion[key].append(dict(pin_axis_sign=sign, sampled_count=49, hits=failures))
    save(folder / "return_spring_seat_study.json", dict(
        generator=provenance(ROOT / "tools/study_s1_series_elbow_spring_seats.py"), parent=provenance(report_path),
        parts=geometry, reloaded_stroke_screen=stroke_results, preassembly_approaches=clearance,
        refined_approach=dict(direction=best["direction"], source_vector=best["vector"],
                              sampled_count=97 if accepted.keys() == before.keys() else 0, hits=insertion_hits),
        pin_insertion=pin_insertion,
        scope="Local seating pocket study with rigid neighbours in humerus coordinates."
              " Lug bearing material protected. Full preassembly, tool space, loaded hook rotation and pin retention remain open."))
    print("Pocket fit", geometry, "clear approach trials", [row["direction"] for row in clearance if not row["hits"]], flush=True)


if __name__ == "__main__":
    main()
