"""Rebuild the FJ1468 proximal lumen while preserving the existing port axis."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from audit_elbow_slider_short_guides import export_lattice
from study_elbow_series_package import bounded, overlaps
from study_s1_infraspinatus_smooth_tendon import centreline
from study_s1_infraspinatus_wire_travel import route_metrics
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1468_negative_channel_screen.json"
    screen = read(source)
    stroke = load_reference(screen["path_stroke"])
    inventory = load_reference(stroke["inventory"])
    route = next(row for row in inventory["routes"] if row["id"] == "FJ1468")
    path = np.array(route["local_route"]["points_source_mm"])
    mount = load_reference(screen["proximal_mount"])["sources"]
    lattice = load_component(mount["tpu:deltoid"])
    fixed = {key: load_component(meta) for key, meta in mount.items()
             if key.startswith(("bone:", "origin:FJ1468/")) or key == "cap:FJ1468/origin"}
    targets = bounded(fixed)
    start_idx = 3  # Last point on the existing rigid proximal axis.
    start = path[start_idx]
    lead = A.D.unit(start - path[start_idx - 1])
    baseline = route_metrics(path)
    trials = []
    ranked = []
    for join, a, b in product(range(6, 16), (.15, .3, .5, .75, 1.), (.15, .3, .5, .75, 1.)):
        end = path[join]
        distance = np.linalg.norm(end - start)
        tangent = A.D.unit(path[join + 1] - end)
        bend = centreline([start, start + a * distance * lead,
                           end - b * distance * tangent, end], 49)
        full = np.vstack([path[:start_idx], bend, path[join + 1:]])
        length, turns, radii = route_metrics(full)
        transition = slice(start_idx - 1, start_idx + len(bend) - 1)
        row = dict(join_index=join, handle_ratios=[a, b],
                   transition_min_radius_mm=float(min(radii[transition])),
                   full_route_min_radius_mm=float(min(radii)), arc_mm=length,
                   largest_turn_deg=float(np.rad2deg(max(turns))))
        trials.append(row)
        ranked.append((row, bend, full))
    checked = []
    # Retain a spread of join locations for local fit, including short blends.
    for join in range(6, 16):
        options = [item for item in ranked if item[0]["join_index"] == join]
        for row, bend, full in sorted(options, key=lambda item: -item[0]["transition_min_radius_mm"])[:3]:
            wire = A.D.outer_envelope(bend, .30)
            row["rigid_hits_mm3"] = overlaps(wire, targets)
            row["existing_TPU_wire_overlap_mm3"] = float((wire ^ lattice).volume())
            checked.append((row, bend, full))
    recuts = []
    clear = [(row, bend, full) for row, bend, full in checked if not row["rigid_hits_mm3"]]
    for row, bend, full in sorted(clear, key=lambda item: -item[0]["transition_min_radius_mm"])[:5]:
        # Carry the bore through both ends of the added outer envelope.
        # Otherwise its end caps can intrude into the retained channel.
        drill_points = np.vstack([path[start_idx - 1], bend, path[row["join_index"] + 1]])
        lumen = A.D.outer_envelope(drill_points, .86)
        outer = A.D.outer_envelope(bend, 1.41)
        # The local guide grows from the existing TPU while clearing adjacent
        # rigid hardware. Report missing wall separately instead of hiding it.
        guard = outer - lumen
        obstacles = [solid for key, solid in fixed.items() if key in overlaps(outer, targets)]
        for part in obstacles:
            outer -= part
        revised = ((lattice + outer) - lumen).as_original().simplify(1e-6)
        result = dict(
            join_index=row["join_index"], handle_ratios=row["handle_ratios"],
            transition_min_radius_mm=row["transition_min_radius_mm"],
            centreline_source_mm=full.tolist(),
            component_volumes_mm3=sorted((float(p.volume()) for p in revised.decompose()
                                         if p.volume() > 1e-6), reverse=True),
            TPU_added_mm3=float((revised - lattice).volume()),
            TPU_removed_mm3=float((lattice - revised).volume()),
            missing_nominal_wall_mm3=float((guard - revised).volume()),
            full_wire_TPU_overlap_mm3=float((A.D.outer_envelope(full, .30) ^ revised).volume()),
        )
        recuts.append(result)
        if (len(result["component_volumes_mm3"]) == 1 and result["full_wire_TPU_overlap_mm3"] <= .001
                and "mesh" not in recuts[0]):
            actual, meta = export_lattice(revised, folder / "FJ1468_proximal_transition_TPU.obj", start)
            full_wire = A.D.outer_envelope(full, .30)
            result["readback_wire_overlap_mm3"] = float((full_wire ^ actual).volume())
            if result["readback_wire_overlap_mm3"] > .001:
                raise ValueError("Channel clearance did not reproduce after mesh export")
            result["mesh"] = meta
            break
    save(folder / "FJ1468_proximal_transition_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_deltoid_proximal_transition.py"),
        source=provenance(source), original_min_radius_mm=float(min(baseline[2])),
        original_turn_max_deg=float(np.rad2deg(max(baseline[1]))), trials=trials, recuts=recuts,
        scope="Neutral FJ1468 centreline and matching proximal mount. The original bore axis "
              "and downstream points beyond the join are retained. Any exported TPU is a neutral "
              "passage candidate. Clipped wall, multi-muscle and loaded-motion checks remain open.",
    ))
    print("clear", len(clear), "best", max((r["transition_min_radius_mm"] for r, _, _ in clear), default=0.),
          flush=True)
    for row in recuts:
        print(row["join_index"], row["transition_min_radius_mm"], len(row["component_volumes_mm3"]),
              row["missing_nominal_wall_mm3"], row["full_wire_TPU_overlap_mm3"], flush=True)


if __name__ == "__main__":
    main()
