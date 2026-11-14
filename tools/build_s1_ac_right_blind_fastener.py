"""Build a short right AC screw on a screened accessible axis."""

import numpy as np

from audit_elbow_slider_short_guides import export_lattice
from audit_s1_infraspinatus_split_installation import insertion
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_AC_assembly_pose_screen.json"
    pose = read(source)
    access = load_reference(pose["source"])
    closure = load_reference(access["source"])
    entry = load_reference(closure["source"])
    scene = {key: load_component(meta) for key, meta in entry["sources"].items()}
    hardware = {key.removeprefix("AC:AC_"): load_component(meta)
                for key, meta in closure["hardware_sources"].items()}
    axis = np.array(closure["axis_source"])
    vertices = np.asarray(hardware["bolt_right"].to_mesh64().vert_properties)[:, :3]
    axial = vertices @ axis
    reference = vertices[axial >= axial.max() - .001].mean(axis=0) + .05 * axis
    cap = hardware["cap_right"]
    bone = scene["bone:scapula"]
    fixed = {key: part for key, part in scene.items() if key != "bone:scapula"}
    fixed.update({key: part for key, part in hardware.items() if key.endswith("_left")})
    fixed["bottom_washer_left"] = load_component(access["left_washer"]["mesh"])
    trials = []
    chosen = None
    for index, row in enumerate(pose["parallel_offset_diagnostics"]):
        if row["driver_hits_mm3"] or not any("cap_right" in s["stock_mm3"]
                                             for s in row.get("inboard_stock_slices", [])):
            continue
        point = reference + np.array(row["offset_source_mm"])

        def pin(start, end, radius):
            return rod(point - start * axis, point - end * axis, radius, 64)

        # Measure the outer cap surface with a narrow cylindrical witness.
        lo, hi = 2., 5.
        for _ in range(14):
            mid = (lo + hi) / 2.
            if (pin(0., mid, .2) ^ cap).volume() > 1e-7:
                hi = mid
            else:
                lo = mid
        seat = hi
        stock_depths = []
        for depth in np.arange(7.25, 16., .25):
            ring = pin(depth, depth + .25, 1.8) - pin(depth - .05, depth + .30, 1.25)
            stock_depths.append(float((ring ^ bone).volume() / ring.volume()))
        windows = [(min(stock_depths[i:i + 12]), i) for i in range(len(stock_depths) - 11)]
        stock_min, stock_index = max(windows, key=lambda pair: (round(pair[0], 6), -pair[1]))
        engagement_start = 7.25 + .25 * stock_index
        pilot_end = engagement_start + 3.25
        cap_hole = pin(-.5, 7.2, 1.65)
        pilot = pin(seat, pilot_end, 1.25)
        head_radius = 2.55
        head_seat = seat + .15
        head_counterbore = pin(-.5, head_seat, head_radius + .08)
        revised_cap = (cap - cap_hole - head_counterbore).as_original().simplify(1e-6)
        revised_bone = (bone - pilot).as_original().simplify(1e-6)
        screw = pin(head_seat - 1.8, head_seat, head_radius) + pin(head_seat - .15, pilot_end - .25, 1.19)
        witness = screw.translate(-.05 * axis)
        seat_volume = float((witness ^ revised_cap).volume())
        thread_zone = pin(engagement_start, engagement_start + 3., 1.8) - pilot
        thread_stock_fraction = float((thread_zone ^ revised_bone).volume() / thread_zone.volume())
        trial = dict(axis_index=index, point_source_mm=point.tolist(), head_seat_depth_mm=seat,
                     cap_removed_mm3=float((cap - revised_cap).volume()),
                     scapula_removed_mm3=float((bone - revised_bone).volume()),
                     head_bearing_probe_mm3=seat_volume, pilot_wall_stock_fraction=thread_stock_fraction,
                     engagement_start_mm=engagement_start, pilot_end_mm=pilot_end,
                     minimum_slice_stock_fraction=stock_min, stock_slice_fractions=stock_depths)
        trials.append(trial)
        if thread_stock_fraction < .98:
            # Keep the collar below the cap face and blend it into the
            # existing scapula along the full measured engagement zone.
            collar = pin(7.2, pilot_end + .5, 2.1)
            revised_bone = ((bone + collar) - pilot).as_original().simplify(1e-6)
        addition = revised_bone - bone
        trial["added_support_mm3"] = float(addition.volume())
        trial["added_support_hits_mm3"] = overlaps(addition, bounded({**fixed, "cap_right": revised_cap}))
        trial["revised_pilot_wall_stock_fraction"] = float((thread_zone ^ revised_bone).volume() / thread_zone.volume())
        trial["scapula_components"] = sum(p.volume() > 1e-6 for p in revised_bone.decompose())
        if (trial["revised_pilot_wall_stock_fraction"] < .98 or seat_volume < .05
                or trial["added_support_hits_mm3"] or trial["scapula_components"] != 1):
            continue
        world = bounded({**fixed, "cap_right": revised_cap, "bone:scapula": revised_bone})
        trial["entry"] = insertion(screw, world, axis, 16.)
        tool_start = point - (head_seat - 1.8 - .05) * axis
        tool = rod(tool_start, tool_start + 20. * axis, 1.5, 64)
        trial["driver_hits_mm3"] = overlaps(tool, world)
        if trial["entry"]["hits"] or trial["driver_hits_mm3"]:
            continue
        scapula_case = entry["cases"]["clavicle_seated_AC_hardware_deferred"]
        fixed_entry = bounded({**{key: scene[key] for key in scapula_case["fixed"]}, "cap_right": revised_cap})
        trial["support_on_scapula_entry"] = [
            insertion(addition, fixed_entry, np.asarray(trace["outward_direction_source"]), 24.)
            for trace in scapula_case["fine"] if not trace["hits"]
        ]
        if not any(not trace["hits"] for trace in trial["support_on_scapula_entry"]):
            continue
        parts = {}
        meshes = {}
        for key, solid in (("cap_right", revised_cap), ("bone:scapula", revised_bone), ("screw_right", screw)):
            writer = export_lattice if key == "bone:scapula" else export
            parts[key], meshes[key] = writer(solid, folder / f"AC_blind_{key.replace(':', '_')}.obj", point)
        world = bounded({**fixed, "cap_right": parts["cap_right"], "bone:scapula": parts["bone:scapula"]})
        trial["readback_entry"] = insertion(parts["screw_right"], world, axis, 16.)
        trial["readback_driver_hits_mm3"] = overlaps(tool, world)
        if trial["readback_entry"]["hits"] or trial["readback_driver_hits_mm3"]:
            raise ValueError("Right AC screw access changed after export")
        fixed_entry = bounded({**{key: scene[key] for key in scapula_case["fixed"]},
                               "cap_right": parts["cap_right"]})
        trial["readback_scapula_entry"] = []
        for prior in scapula_case["fine"]:
            if prior["hits"]:
                continue
            direction = np.asarray(prior["outward_direction_source"])
            traces = {"bone:scapula": insertion(parts["bone:scapula"], fixed_entry, direction, 24.)}
            for key in scapula_case["moving"]:
                if key != "bone:scapula":
                    traces[key] = insertion(scene[key], bounded({"cap_right": parts["cap_right"]}), direction, 24.)
            trial["readback_scapula_entry"].append(dict(traces=traces,
                sampled_clear=all(not trace["hits"] for trace in traces.values())))
        if not any(row["sampled_clear"] for row in trial["readback_scapula_entry"]):
            raise ValueError("No full scapula insertion path with the revised right cap held")
        chosen = dict(**trial, meshes=meshes)
        break
    save(folder / "FJ1504_AC_right_blind_fastener.json", dict(
        generator=provenance(ROOT / "tools/build_s1_ac_right_blind_fastener.py"), source=provenance(source),
        trials=trials, selected=chosen, scope="Short right screw seated through the existing right cap "
        "into a new scapular blind pilot. A smooth minor-diameter envelope is used for screw and pilot. "
        "The right cap is held at its registered pose. Thread mechanics, joint motion and tolerance "
        "remain to be checked.",
    ))
    print("selected", None if chosen is None else chosen["axis_index"],
          [(r["axis_index"], round(r["pilot_wall_stock_fraction"], 3), round(r["head_bearing_probe_mm3"], 3),
            len(r.get("entry", {}).get("hits", [])), r.get("driver_hits_mm3")) for r in trials], flush=True)


if __name__ == "__main__":
    main()
