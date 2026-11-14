"""Recheck the split FJ1500 belly against the actual shoulder fasteners."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def insertion(moving, targets, direction, travel):
    hits = []
    for offset in np.linspace(travel, 0., round(travel * 10) + 1):
        overlap = overlaps(moving.translate(offset * direction), targets)
        if overlap:
            hits.append(dict(offset_mm=float(offset), hits_mm3=overlap))
    return dict(offset_direction_source=direction.tolist(), travel_mm=travel,
                step_mm=.1, hits=hits)


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    split_file = folder / "FJ1500_split_origin_screen.json"
    split = read(split_file)
    relief = load_reference(split["parent"])
    side = load_reference(relief["parent"])
    removable = load_reference(side["source"])
    feed = load_reference(removable["integrated_sleeve"])
    assembly = load_reference(feed["mount"])
    cap = load_component(removable["selected"]["cap"])
    sleeve = load_component(removable["selected"]["sleeve"])
    scene = {key: load_component(value) for key, value in assembly["sources"].items()
             if key != "cap:FJ1500/origin"}
    scene["cap:FJ1500/origin"] = cap
    configurations = {
        "FJ1504_origin_deferred": {key: part for key, part in scene.items()
                                  if not key.startswith("origin:FJ1504/")},
        "all_origin_hardware_installed": scene,
    }
    path = np.asarray(load_reference(feed["route"])["exported"]["trial"]["new_path_source_mm"])
    axis = path[-1] - path[2]
    axis /= np.linalg.norm(axis)
    sleeve_offset_direction = -np.asarray(removable["bore_axis_source"])
    cases = []
    for row in split["cases"]:
        if "belly_mesh" not in row:
            continue
        belly = load_component(row["belly_mesh"])
        centre = np.asarray(row["channel_crossings_source_mm"][0])
        stock = []
        for start, end in ((0., 1.), (1., 2.), (2., 4.)):
            outer = rod(centre + start * axis, centre + end * axis, 3., 64)
            inner = rod(centre + (start - .1) * axis, centre + (end + .1) * axis, .86, 64)
            ring = outer - inner
            material = float((ring ^ belly).volume())
            stock.append(dict(interval_mm=[start, end], stock_mm3=material,
                              annulus_volume_mm3=float(ring.volume()),
                              fraction=material / float(ring.volume())))
        results = {}
        for label, parts in configurations.items():
            targets = bounded(parts)
            directions = [np.asarray(entry["direction_source"]) for entry in row["belly_without_pad"]["fine"]
                          if not entry["hits"]]
            results[label] = dict(
                belly_insertion=[insertion(belly, targets, direction, 8.) for direction in directions],
                seated_hits_mm3=overlaps(belly, targets, threshold=1e-6),
                sleeve_after_belly=insertion(sleeve, bounded({**parts, "FJ1500_split_belly": belly}),
                                             sleeve_offset_direction, 8.),
            )
        cases.append(dict(split_station_mm=row["split_station_mm"], belly_mesh=row["belly_mesh"],
                          stock=stock, results=results))
        print(row["split_station_mm"], [round(item["fraction"], 3) for item in stock],
              {label: ([len(trace["hits"]) for trace in result["belly_insertion"]],
                       len(result["sleeve_after_belly"]["hits"])) for label, result in results.items()},
              flush=True)
    save(folder / "FJ1500_split_installation_audit.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_infraspinatus_split_installation.py"),
        split=provenance(split_file), intersection_threshold_mm3=.001, cases=cases,
        scope="Last 8 mm of fixed-orientation installation at neutral. The installed sleeve offset is "
              "measured outward from its seated position. The planar proximal pad is absent. "
              "Annulus intersections quantify available lattice material, with no strength or retention claim.",
    ))


if __name__ == "__main__":
    main()
