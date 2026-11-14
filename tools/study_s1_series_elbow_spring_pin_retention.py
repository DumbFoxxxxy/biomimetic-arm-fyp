"""Screen radial dog-point retainers for the two removable spring pins."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_series_elbow_return_spring import spring_path, sweep
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    seat_path = folder / "return_spring_seat_study.json"
    seats = read(seat_path)
    spring_report = load_reference(seats["parent"])
    size = load_reference(spring_report["size"])
    audit = load_reference(size["geometry"])
    extended = load_reference(audit["extended_drive"])
    assembly = load_reference(load_reference(extended["budget"])["assembly"])
    mount = load_reference(size["spring_mount"])["spring_mount"]
    start, end = (np.array(mount[key]) for key in ("fixed_eye_source_mm", "moving_eye_neutral_source_mm"))
    slide_axis = A.D.unit(end-start)
    pin_axis = np.array(mount["pin_axis_source"])
    transverse = np.cross(slide_axis, pin_axis)
    stroke = size["full_rail_stroke_mm"]
    source_meshes = {**assembly["sources"], **extended["shafts"],
                     **{key: row["candidate_mesh"] for key, row in seats["parts"].items()}}
    owners = assembly["owners"]
    parts = {key: load_component(source) for key, source in source_meshes.items()}
    local = {key: solid for key, solid in parts.items() if owners[key] in ("humerus", "FJ1477_carriage")}
    spring = load_component(spring_report["exports"]["neutral"]["mesh"])
    outcomes = {}
    selected_parts = {}
    for label, origin, lug_key in (("fixed", start, "elbow:FJ1477_bone_adapter"),
                                  ("moving", end, "elbow:FJ1477_proximal_rigid_cap")):
        pin_key = f"series:{label}_pin"
        original = parts[lug_key]
        pin = parts[pin_key]
        other = {key: part for key, part in local.items() if key not in (pin_key, lug_key)}
        obstacles = bounded({**other, "spring": spring})
        candidates = []
        selected = None
        for station in (2.65, 2.35, 2.95):
            centre = origin + station*pin_axis
            groove = rod(centre-.50*pin_axis, centre+.50*pin_axis, 1.02) - rod(centre-.55*pin_axis, centre+.55*pin_axis, .62)
            grooved = (pin-groove).as_original().simplify(1e-6)
            bearing_guard = rod(origin+1.25*pin_axis, origin+4.03*pin_axis, .971)
            for angle in range(0, 360, 30):
                radial = np.cos(np.deg2rad(angle))*transverse + np.sin(np.deg2rad(angle))*slide_axis
                point = lambda d: centre+d*radial
                boss = rod(point(1.04), point(3.17), 1.13)
                if (boss ^ original).volume() < .15:
                    continue
                bore = rod(point(.61), point(3.19), .81)
                lug = ((original+boss)-bore-bearing_guard).as_original().simplify(1e-6)
                screw = (rod(point(1.06), point(3.13), .78) + rod(point(.66), point(1.10), .39)).as_original()
                static = overlaps((lug-original)+screw, obstacles)
                if static or sum(part.volume() > 1e-5 for part in lug.decompose()) != 1:
                    candidates.append(dict(station_mm=station, angle_deg=angle, static_hits_mm3=static))
                    continue
                if overlaps(screw, bounded({"lug": lug, "grooved_pin": grooved})):
                    continue
                # Reserve a straight driver envelope behind the screw and
                # check the screw's actual approach with the pin already seated.
                driver = rod(point(3.18), point(11.13), 1.05)
                driver_hits = overlaps(driver, bounded({**other, "spring": spring, "lug": lug}))
                screw_path = bounded({**other, "spring": spring, "lug": lug, "grooved_pin": grooved})
                screw_hits = [float(d) for d in np.linspace(6., 0., 49) if overlaps(screw.translate(d*radial), screw_path)]
                pin_path = bounded({**other, "spring": spring, "lug": lug})
                pin_hits = [float(d) for d in np.linspace(6., 0., 49) if overlaps(grooved.translate(-d*pin_axis), pin_path)]
                witnesses = [float((grooved.translate(d*pin_axis) ^ screw).volume()) for d in (-.25, .25)]
                row = dict(station_mm=station, angle_deg=angle, radial_source=radial.tolist(),
                           driver_hits_mm3=driver_hits, screw_insertion_hit_offsets_mm=screw_hits,
                           pin_insertion_hit_offsets_mm=pin_hits, axial_025mm_stop_witness_mm3=witnesses)
                candidates.append(row)
                if driver_hits or screw_hits or pin_hits or min(witnesses) <= .001:
                    continue
                # Both directions must be blocked before the pin can leave
                # the eye. Record the first sampled solid stop for each.
                axial = {}
                for sign in (-1., 1.):
                    axial[str(sign)] = next(float(d) for d in np.arange(.01, .251, .01)
                                             if (grooved.translate(sign*d*pin_axis) ^ screw).volume() > .001)
                exported = {}
                reloaded = {}
                for name, solid in (("lug", lug), ("pin", grooved), ("dog_screw", screw)):
                    reloaded[name], exported[name] = export(solid, folder / f"spring_{label}_{name}.obj", origin)
                if overlaps(reloaded["dog_screw"], bounded({"lug": reloaded["lug"], "pin": reloaded["pin"]})):
                    raise ValueError("Retention fit changed on readback")
                for sign, distance in axial.items():
                    if (reloaded["pin"].translate(float(sign)*distance*pin_axis) ^ reloaded["dog_screw"]).volume() <= .001:
                        raise ValueError("Retention stop changed on readback")
                selected = dict(**row, first_axial_stop_mm=axial, meshes=exported)
                selected_parts.update({lug_key: reloaded["lug"], pin_key: reloaded["pin"],
                                       f"series:{label}_dog_screw": reloaded["dog_screw"]})
                local.update(selected_parts)
                break
            if selected:
                break
        outcomes[label] = dict(trials=candidates, selected=selected)
    motion = []
    spring_hits = []
    assembled_path = {}
    if all(row["selected"] for row in outcomes.values()):
        added = {key: solid if key.endswith("dog_screw") else solid-parts[key] for key, solid in selected_parts.items()}
        owner = lambda key: ("FJ1477_carriage" if key.startswith("series:moving") else "humerus") if key.endswith("dog_screw") else owners[key]
        design = spring_report["design"]
        radius = .5 * design["wire_diameter_mm"]
        mean = .5 * (design["coil_OD_mm"] - design["wire_diameter_mm"])
        turns = design["active_turns"] + 2
        frame = np.c_[pin_axis, transverse, slide_axis, start]
        for travel in np.linspace(0., stroke, 31):
            posed = {key: solid.translate(travel*slide_axis) if owner(key) == "FJ1477_carriage" else solid
                     for key, solid in local.items()}
            failures = {}
            for key, material in added.items():
                if owner(key) == "FJ1477_carriage":
                    material = material.translate(travel*slide_axis)
                other = {k: solid for k, solid in posed.items() if k != key and owner(k) != owner(key)}
                hit = overlaps(material, bounded(other))
                if hit:
                    failures[key] = hit
            if failures:
                motion.append(dict(travel_mm=float(travel), hits_mm3=failures))
            moving_spring = sweep(spring_path(np.linalg.norm(end-start)+travel, turns, radius, mean), radius).transform(frame)
            hit = overlaps(moving_spring, bounded({key: posed[key] for key in selected_parts}))
            if hit:
                spring_hits.append(dict(travel_mm=float(travel), hits_mm3=hit))
        approach = np.array(seats["refined_approach"]["source_vector"])
        install_order = ["series:fixed_pin", "series:moving_pin", "series:fixed_dog_screw", "series:moving_dog_screw"]
        stage = {key: solid for key, solid in local.items() if key not in install_order}
        assembled_path["spring"] = [float(d) for d in np.linspace(6., 0., 97)
                                    if overlaps(spring.translate(d*approach), bounded(stage))]
        stage["spring"] = spring
        for key in install_order:
            label = "fixed" if ":fixed" in key else "moving"
            radial = np.array(outcomes[label]["selected"]["radial_source"])
            direction = radial if key.endswith("dog_screw") else -pin_axis
            targets = bounded(stage)
            assembled_path[key] = [float(d) for d in np.linspace(6., 0., 49)
                                   if overlaps(local[key].translate(d*direction), targets)]
            stage[key] = local[key]
    save(folder / "return_spring_pin_retention.json", dict(
        generator=provenance(ROOT / "tools/study_s1_series_elbow_spring_pin_retention.py"), parent=provenance(seat_path),
        pin_outer_diameter_mm=1.8, groove_core_diameter_mm=1.24, groove_width_mm=1.,
        nominal_screw_body_diameter_mm=1.56, nominal_dog_point_diameter_mm=.78,
        results=outcomes, cross_body_slider_sample_count=31 if all(row["selected"] for row in outcomes.values()) else 0,
        cross_body_slider_hits=motion, spring_full_stroke_hits=spring_hits,
        combined_preassembly_hit_offsets_mm=assembled_path,
        scope="Solid dog point and groove witnesses. Thread form, tapped metal or insert, axial capacity,"
              " driver grip and the rest of the complete moving arm require additional verification."))
    print("Retainers", {key: (row["selected"] or {}).get("first_axial_stop_mm") for key, row in outcomes.items()},
          "moving hits", motion, flush=True)


if __name__ == "__main__":
    main()
