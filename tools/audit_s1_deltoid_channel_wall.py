"""Locate anterior deltoid wall deficits and recheck neighbouring wire passages."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1468_proximal_transition_screen.json"
    transition = read(source)
    selected = next(row for row in transition["recuts"] if "mesh" in row)
    screen = load_reference(transition["source"])
    stroke = load_reference(screen["path_stroke"])
    inventory = load_reference(stroke["inventory"])
    route = {row["id"]: row["local_route"]["points_source_mm"] for row in inventory["routes"]}
    mount = load_reference(screen["proximal_mount"])["sources"]
    original = load_component(mount["tpu:deltoid"])
    revised = load_component(selected["mesh"])
    solids = {key: load_component(meta) for key, meta in mount.items()
              if key.startswith(("bone:", "origin:FJ1468/")) or key == "cap:FJ1468/origin"}
    targets = bounded(solids)
    full = np.asarray(selected["centreline_source_mm"])
    start = 3
    stop = len(full) - (len(route["FJ1468"]) - selected["join_index"] - 1)
    bend = full[start:stop]
    lumen = A.D.outer_envelope(full[start - 1:stop + 1], .86)
    wall = A.D.outer_envelope(bend, 1.41) - lumen
    missing = wall - revised
    owners = overlaps(missing, targets, threshold=1e-7)
    blocked = union([missing ^ solids[key] for key in owners]) if owners else missing - missing
    clear_missing = missing - blocked
    cap = solids["cap:FJ1468/origin"]
    combined_guide = revised + cap
    passage = {}
    for radius in (.30, .50, .70, .86):
        swept = A.D.outer_envelope(bend, radius)
        passage[f"{radius:.2f}"] = dict(cap_overlap_mm3=float((swept ^ cap).volume()),
                                        TPU_overlap_mm3=float((swept ^ revised).volume()))
    lengths = np.linalg.norm(np.diff(bend, axis=0), axis=1)
    stations = np.r_[0., np.cumsum(lengths)]
    sections = []
    for index, (a, b) in enumerate(zip(bend, bend[1:])):
        window = rod(a, b, 1.42, 48)
        deficit = window ^ missing
        sections.append(dict(interval_mm=stations[index:index + 2].tolist(),
                             midpoint_source_mm=((a + b) / 2.).tolist(),
                             missing_wall_mm3=float(deficit.volume()),
                             wall_missing_after_cap_mm3=float((window ^ (wall - combined_guide)).volume()),
                             solid_owners_mm3=overlaps(deficit, targets, threshold=1e-7)))
    wires = {}
    for key in ("FJ1468", "FJ1467", "FJ1513"):
        old_wire = A.D.outer_envelope(np.asarray(route[key]), .3)
        wire = A.D.outer_envelope(full, .3) if key == "FJ1468" else old_wire
        wires[key] = dict(original_TPU_overlap_mm3=float((old_wire ^ original).volume()),
                          revised_TPU_overlap_mm3=float((wire ^ revised).volume()),
                          new_stock_overlap_mm3=float((wire ^ (revised - original)).volume()),
                          wall_overlap_with_this_wire_mm3=float((wall ^ wire).volume()))
    save(folder / "FJ1468_channel_wall_audit.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_deltoid_channel_wall.py"),
        source=provenance(source), total_missing_mm3=float(missing.volume()),
        rigid_occupied_missing_mm3=float(blocked.volume()),
        unoccupied_missing_mm3=float(clear_missing.volume()),
        missing_by_solid_mm3=owners, sections=sections, wires=wires,
        assembled_wall_missing_mm3=float((wall - combined_guide).volume()),
        proximal_passage_by_radius_mm=passage,
        scope="Neutral pose. The intervals measure local cylindrical windows; neighbouring windows "
              "may overlap at a bend and their volumes should not be summed. Rigid occupancy is "
              "computed on the full nominal wall solid.",
    ))
    print("wall", missing.volume(), "occupied", blocked.volume(), owners, flush=True)
    print("most affected", sorted(sections, key=lambda r: -r["missing_wall_mm3"])[:3], flush=True)
    print(wires, flush=True)
    print("assembled guide", float((wall - combined_guide).volume()), passage, flush=True)


if __name__ == "__main__":
    main()
