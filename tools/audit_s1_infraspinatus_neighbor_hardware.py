"""Probe FJ1504 reassembly after the split FJ1500 belly is seated."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


PREFIX = "origin:FJ1504/"


def stage(moving, fixed, direction, travel):
    targets = bounded(fixed)
    hits = []
    for offset in np.linspace(travel, 0., round(travel * 10) + 1):
        for name, solid in moving.items():
            found = overlaps(solid.translate(direction * offset), targets)
            if found:
                hits.append(dict(offset_mm=float(offset), moving=name, overlaps_mm3=found))
    return dict(travel_mm=travel, step_mm=.1, outward_direction_source=direction.tolist(), hits=hits)


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_split_installation_audit.json"
    audit = read(source)
    split = load_reference(audit["split"])
    relief = load_reference(split["parent"])
    side = load_reference(relief["parent"])
    removable = load_reference(side["source"])
    feed = load_reference(removable["integrated_sleeve"])
    assembly = load_reference(feed["mount"])
    scene = {key: load_component(value) for key, value in assembly["sources"].items()}
    scene["cap:FJ1500/origin"] = load_component(removable["selected"]["cap"])
    scene["FJ1500/feed_sleeve"] = load_component(removable["selected"]["sleeve"])
    routes_file = A.DESIGN / "s1_port_path_revision/routes_compact_distal_candidate.json"
    route = next(row for row in read(routes_file)["routes"] if row["head"] == "FJ1504")
    normal = np.asarray(route["ports"]["origin"]["normal_source"])
    normal /= np.linalg.norm(normal)
    phases = [
        ("front_support", ["front_support"], normal, 12.),
        ("backing_plate", ["backing_plate"], -normal, 12.),
        ("front_bolt_0", ["bolt_0", "front_washer_0"], normal, 30.),
        ("front_bolt_1", ["bolt_1", "front_washer_1"], normal, 30.),
        ("rear_washer_0", ["back_washer_0"], -normal, 12.),
        ("rear_nut_0", ["nut_0"], -normal, 12.),
        ("rear_washer_1", ["back_washer_1"], -normal, 12.),
        ("rear_nut_1", ["nut_1"], -normal, 12.),
    ]
    results = []
    for case in [None, *audit["cases"]]:
        fixed = {key: solid for key, solid in scene.items() if not key.startswith(PREFIX)}
        if case is not None:
            fixed["FJ1500/split_belly"] = load_component(case["belly_mesh"])
        rows = []
        for name, keys, direction, travel in phases:
            moving = {PREFIX + key: scene[PREFIX + key] for key in keys}
            result = stage(moving, fixed, direction, travel)
            rows.append(dict(stage=name, **result))
            # Later stages are conditional geometry probes even if an earlier
            # stage is blocked. Keep the earlier failure in the same report.
            fixed.update(moving)
        tools = []
        for index in range(2):
            for member, outward, radius in (("bolt", normal, 1.5), ("nut", -normal, 4.0)):
                key = f"{PREFIX}{member}_{index}"
                solid = fixed[key]
                vertices = np.asarray(solid.to_mesh64().vert_properties)[:, :3]
                axial = vertices @ outward
                end_face = vertices[axial >= axial.max() - .001]
                start = (end_face.max(axis=0) + end_face.min(axis=0)) / 2.
                start += (axial.max() + .05 - start @ outward) * outward
                probe = rod(start, start + 20. * outward, radius, 32)
                if member == "nut":
                    # Leave room for the protruding M3 bolt in the socket.
                    probe -= rod(start - outward, start + 21. * outward, 1.6, 32)
                tools.append(dict(hardware=key, start_source_mm=start.tolist(), radius_mm=radius,
                                  hits_mm3=overlaps(probe, bounded(fixed))))
        station = case["split_station_mm"] if case is not None else None
        print(station, [(row["stage"], len(row["hits"])) for row in rows],
              [(row["hardware"], list(row["hits_mm3"])) for row in tools], flush=True)
        sequence_conflicts = []
        if case is not None:
            blocked_late = {
                hit["moving"] for row in rows for hit in row["hits"]
                if "FJ1500/split_belly" in hit["overlaps_mm3"]
            }
            for trace in case["results"]["all_origin_hardware_installed"]["belly_insertion"]:
                blocked_early = {
                    key for hit in trace["hits"] for key in hit["hits_mm3"]
                    if key.startswith(PREFIX)
                }
                sequence_conflicts.append(dict(
                    outward_direction_source=trace["offset_direction_source"],
                    hardware_blocking_belly_when_preinstalled=sorted(blocked_early),
                    hardware_blocked_by_seated_belly=sorted(blocked_late),
                    hardware_conflicting_in_both_orders=sorted(blocked_early & blocked_late),
                ))
        results.append(dict(split_station_mm=station, stages=rows, axial_tool_probes=tools,
                            sequence_conflicts=sequence_conflicts))
    save(folder / "FJ1504_hardware_after_split_belly.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_infraspinatus_neighbor_hardware.py"),
        split_audit=provenance(source), routes=provenance(routes_file),
        overlap_threshold_mm3=.001, results=results,
        scope="Nominal shoulder zero-position translations using the FJ1504 fastener axis. "
              "The FJ1504 origin cap stays in place. Later stages assume earlier parts seated. "
              "Tool probes are coaxial straight clearance envelopes outside the hardware faces; "
              "thread engagement, tool rotation and operator access are not simulated.",
    ))


if __name__ == "__main__":
    main()
