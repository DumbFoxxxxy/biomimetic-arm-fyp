"""Check sequential AC closure after the preassembled scapula is seated."""

import numpy as np

from audit_s1_infraspinatus_split_installation import insertion
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_preassembled_scapula_entry.json"
    entry = read(source)
    cap = load_reference(entry["source"])
    blind = load_reference(cap["source"])
    split = load_reference(load_reference(load_reference(blind["parent"])["split_audit"])["split"])
    side = load_reference(load_reference(split["parent"])["parent"])
    removable = load_reference(side["source"])
    assembly = load_reference(load_reference(removable["integrated_sleeve"])["mount"])
    fixed = {key: load_component(meta) for key, meta in entry["sources"].items()}
    hardware = {key.removeprefix("AC:AC_"): load_component(meta)
                for key, meta in assembly["sources"].items() if key.startswith("AC:")}
    vertices = np.asarray(hardware["bolt_left"].to_mesh64().vert_properties)[:, :3]
    axis = np.linalg.svd(vertices - vertices.mean(axis=0), full_matrices=False)[2][0]
    if axis[0] < 0.:
        axis = -axis
    directions = [axis, -axis, *np.eye(3), *-np.eye(3)]
    for row in entry["cases"]["clavicle_seated_AC_hardware_deferred"]["coarse"]:
        directions.append(np.asarray(row["outward_direction_source"]))

    def cap_entry(part, targets):
        scored = []
        for direction in directions:
            hits = [overlaps(part.translate(offset * direction), targets)
                    for offset in (.2, .5, 1., 2., 4., 8., 12.)]
            scored.append((sum(bool(hit) for hit in hits), sum(sum(hit.values()) for hit in hits), direction))
        best = sorted(scored, key=lambda row: row[:2])[:5]
        return [insertion(part, targets, direction, 12.) for _, _, direction in best]

    cases = []
    for order in (("left", "right"), ("right", "left")):
        world = dict(fixed)
        rows = []
        phases = [(f"cap_{side}", None, 12.) for side in order]
        for side in order:
            phases.extend(((f"top_washer_{side}", axis, 12.), (f"bolt_{side}", axis, 22.),
                           (f"bottom_washer_{side}", -axis, 12.), (f"nut_{side}", -axis, 12.)))
        for name, direction, travel in phases:
            part = hardware[name]
            targets = bounded(world)
            traces = (cap_entry(part, targets) if direction is None else
                      [insertion(part, targets, direction, travel)])
            rows.append(dict(part=name, traces=traces, sampled_clear=any(not row["hits"] for row in traces)))
            world[name] = part
        tools = []
        for side in order:
            for name, direction, radius in ((f"bolt_{side}", axis, 1.5), (f"nut_{side}", -axis, 4.5)):
                verts = np.asarray(world[name].to_mesh64().vert_properties)[:, :3]
                axial = verts @ direction
                face = verts[axial >= axial.max() - .001]
                point = face.mean(axis=0) + .05 * direction
                probe = rod(point, point + 20. * direction, radius, 32)
                if name.startswith("nut_"):
                    probe -= rod(point - direction, point + 21. * direction, 2.1, 32)
                tools.append(dict(part=name, radius_mm=radius, hits_mm3=overlaps(probe, bounded(world))))
        print(order, [(row["part"], min(len(t["hits"]) for t in row["traces"])) for row in rows],
              [(row["part"], list(row["hits_mm3"])) for row in tools], flush=True)
        cases.append(dict(cap_order=order, stages=rows, tool_probes=tools,
                          nominal_closure_clear=all(row["sampled_clear"] for row in rows)))
    scapula_case = entry["cases"]["clavicle_seated_AC_hardware_deferred"]
    seated_world = bounded({key: fixed[key] for key in scapula_case["fixed"]})
    preload_cases = []
    # The original scapula and all parts in its previous clear trajectories
    # have already been replayed. Check the added AC parts on those same paths.
    for names in (["cap_right"], [key for key in hardware if key.endswith("_right")], list(hardware)):
        carried = union(hardware[name] for name in names)
        traces = [insertion(carried, seated_world, np.asarray(row["outward_direction_source"]), 24.)
                  for row in scapula_case["fine"] if not row["hits"]]
        preload_cases.append(dict(carried_parts=names, additional_parts_traces=traces))
        print("scapula preload", names, [len(row["hits"]) for row in traces], flush=True)
    held_on_clavicle = []
    for names in (["cap_right"], [key for key in hardware if key.endswith("_right")]):
        targets = bounded({name: hardware[name] for name in names})
        traces = []
        for row in scapula_case["fine"]:
            if row["hits"]:
                continue
            direction = np.asarray(row["outward_direction_source"])
            traces.append({key: insertion(fixed[key], targets, direction, 24.)
                           for key in scapula_case["moving"]})
        held_on_clavicle.append(dict(held_parts=names, added_pair_traces=traces))
        print("held on clavicle", names, [sum(len(t["hits"]) for t in route.values()) for route in traces],
              flush=True)
    save(folder / "FJ1504_AC_closure_screen.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_ac_closure.py"), source=provenance(source),
        hardware_sources={key: meta for key, meta in assembly["sources"].items() if key.startswith("AC:")},
        axis_source=axis.tolist(), cap_screen_directions=len(directions), cases=cases,
        scapula_preload_cases=preload_cases, clavicle_held_cases=held_on_clavicle,
        scope="Fixed-orientation translations after scapula entry. Later stages are conditional on "
              "prior parts being seated. Tools are straight clearance probes. Thread engagement, "
              "tolerance, final preload and later GH/humerus stages remain unverified.",
    ))


if __name__ == "__main__":
    main()
