"""Screen a separate spring lane and captured rail ends for the FJ1477 slider."""

from pathlib import Path

import numpy as np

from audit_elbow_cap_travel import load_geometry
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_guide import unit
from build_s1_teres_origin_mount import A, export, load_component, mesh_geom, rod, save
from study_elbow_series_package import ADAPTER, PREFIX, bounded, overlaps, union


def connected(part):
    return sum(piece.volume() > 1e-5 for piece in part.decompose()) == 1


def nearest_link(base, vertices, position):
    nearest = vertices[np.linalg.norm(vertices - position, axis=1).argmin()]
    vector = nearest - position
    distance = np.linalg.norm(vector)
    if distance > 18. or distance < 1e-8:
        return None
    link = rod(position, nearest + 2. * vector / distance, 2.5)
    return link if A.I.intersection_volume(base, link) > 1. else None


def eye(base, vertices, centre, pin_axis):
    # The hook occupies the free portion of the pin, clear of the side tab.
    inner, outer = centre + 1.3 * pin_axis, centre + 4. * pin_axis
    lug = rod(inner, outer, 2.4)
    link = nearest_link(base, vertices, (inner + outer) / 2.)
    if link is None:
        return None
    stock = lug + link
    bore = rod(centre - 1.4 * pin_axis, outer + .3 * pin_axis, .97)
    if A.I.intersection_volume(base, bore) > .001:
        return None
    part = (base + stock) - bore
    if not connected(part):
        return None
    pin = rod(centre - 1.3 * pin_axis, outer + .2 * pin_axis, .9)
    return part, pin


def spring_envelope(start, end, axis):
    body = rod(start + 4. * axis, end - 4. * axis, 3.25)
    hooks = union([rod(start + .9 * axis, start + 4.1 * axis, .8),
                   rod(end - 4.1 * axis, end - .9 * axis, .8)])
    return body + hooks


def main():
    folder = DEST.parent / "elbow_local_muscle_load"
    guide_file = folder / "elbow_series_guide.json"
    guide = read(guide_file)
    package_file = ROOT / guide["package"]["file"]
    assembly_file = ROOT / guide["stop_assembly"]["file"]
    for path, reference in ((package_file, guide["package"]), (assembly_file, guide["stop_assembly"])):
        if provenance(path) != reference:
            raise RuntimeError(f"Input changed: {path.name}")
    package = read(package_file)
    check = read(DEST / "elbow_passive_model/model_check.json")
    solids, _ = load_geometry(check, read(assembly_file))
    original_parts = guide["selected"]["sources"]
    parts = {name: load_component(source) for name, source in original_parts.items()}
    origin = np.asarray(package["bore_origin_mm"])
    axis = np.asarray(guide["axis_source"])
    travel = guide["travel_mm"]
    transverse = unit(np.cross(axis, [1., 0., 0.]))
    radial = np.cross(axis, transverse)
    bases = np.asarray(guide["selected"]["rail_origins_source_mm"])
    adapter, cap = parts["adapter"], parts["moving_cap"]
    adapter_points = mesh_geom.manifold_mesh_to_arrays(adapter.to_mesh64())[0]
    cap_points = mesh_geom.manifold_mesh_to_arrays(cap.to_mesh64())[0]
    moving_old = {name: solids[name][1] for name in package["proximal_moving_parts"]
                  if name != PREFIX + "proximal_rigid_cap"}
    old_fixed = {name: part for name, (_, part) in solids.items()
                 if name not in package["proximal_moving_parts"] and name != ADAPTER}
    fixed = bounded({**old_fixed, **{name: parts[name] for name in ("guide_a", "guide_b", "collar_a", "collar_b")}})
    full_fixed = {**fixed, **bounded({"adapter": adapter})}
    moving_fasteners = bounded(moving_old)
    bare_cap = union([cap, *moving_old.values()])
    travel_samples = np.linspace(0., travel, int(np.ceil(travel)) + 1)
    mf = A.I._manifold()
    bridge_clearance = mf.Manifold.cube([.25] * 3, True)
    rail_sweep_for_moving_lug = union(
        mf.Manifold.batch_hull([parts[name], parts[name].translate(-travel * axis)])
        for name in ("guide_a", "guide_b", "collar_a", "collar_b"))

    retain = {}
    for side, base in zip(("a", "b"), bases, strict=True):
        flange = rod(base + (travel - 10.) * axis, base + (travel - 6.) * axis, 4.)
        shaft = parts[f"guide_{side}"] + flange
        hole = rod(base - 27.1 * axis, base - 22. * axis, 1.23)
        shaft = shaft - hole
        head = rod(base - 29. * axis, base - 27. * axis, 3.1)
        shank = rod(base - 27.1 * axis, base - 22.2 * axis, 1.15)
        screw = head + shank
        check_against = {**{key: value for key, value in full_fixed.items()
                            if key not in (f"guide_{side}", f"collar_{side}")},
                         **bounded({"moving_cap": bare_cap, "flanged_shaft": shaft})}
        hit = overlaps(screw, check_against)
        if hit:
            retain[side] = dict(intersections_mm3=hit)
            continue
        retain[side] = dict(intersections_mm3={}, axial_seat_beyond_socket_mm=2.,
                            sources={})
        for name, solid in (("shaft", shaft), ("retainer_screw", screw)):
            _, meta = export(solid.as_original().simplify(1e-6),
                             folder / f"series_retention_{side}_{name}.obj", origin)
            retain[side]["sources"][name] = meta

    trials = []
    chosen = None
    positions = [(r, a, rear, front) for rear, front in ((-36., -10.), (-32., -10.), (-40., -12.),
                                                         (-32., -6.), (-30., -5.))
                 for r in (10., 14., 18., 22., 26.) for a in range(0, 360, 10)]
    retainers = {side: load_component(info["sources"]["retainer_screw"])
                 for side, info in retain.items() if "sources" in info}
    for radius, phase, rear, front in positions:
        angle = np.deg2rad(phase)
        radial_direction = np.cos(angle) * radial + np.sin(angle) * transverse
        base = origin + radius * radial_direction
        start, end = base + rear * axis, base + front * axis
        if min(np.linalg.norm(base - point) for point in bases) < 8.:
            continue
        spring_sweep = spring_envelope(start, end + travel * axis, axis)
        obstruction = overlaps(spring_sweep, full_fixed)
        row = dict(radial_offset_mm=radius, azimuth_deg=phase, rear_station_mm=rear, front_station_mm=front)
        trials.append(row)
        if obstruction:
            row["spring_sweep_hit"] = obstruction
            continue
        # Fit the lug on either side of the nominal spring line.
        for pin_sign in (-1., 1.):
            pin_axis = pin_sign * unit(np.cross(axis, radial_direction))
            fixed_eye = eye(adapter, adapter_points, start, pin_axis)
            moving_eye = eye(cap, cap_points, end, pin_axis)
            tag = f"pin_{int(pin_sign):+d}"
            if fixed_eye is None or moving_eye is None:
                row[tag + "_no_support"] = [name for name, value in
                                            (("adapter", fixed_eye), ("cap", moving_eye)) if value is None]
                continue
            new_adapter, fixed_pin = fixed_eye
            new_cap, moving_pin = moving_eye
            new_adapter -= spring_sweep.minkowski_sum(bridge_clearance) - adapter
            spring_in_moving_frame = spring_sweep.translate(-travel * axis)
            new_cap -= union([spring_in_moving_frame.minkowski_sum(bridge_clearance),
                              rail_sweep_for_moving_lug.minkowski_sum(bridge_clearance)]) - cap
            adapter_added, cap_added = new_adapter - adapter, new_cap - cap
            if not connected(new_adapter) or not connected(new_cap):
                row[tag + "_relief_breaks_support"] = True
                continue
            mount_hits = (overlaps(adapter_added + fixed_pin, {**fixed, **bounded(retainers)})
                          or overlaps(cap_added + moving_pin, moving_fasteners))
            pin_hits = [float(A.I.intersection_volume(fixed_pin, new_adapter)),
                        float(A.I.intersection_volume(moving_pin, new_cap))]
            if mount_hits or max(pin_hits) > .001:
                row[tag + "_mount_hits"] = dict(other_parts=mount_hits, pin_to_mount_mm3=pin_hits)
                continue
            fitted_fixed = {**fixed, **bounded({"adapter": new_adapter, "spring_pin": fixed_pin, **retainers})}
            moving = union([new_cap, moving_pin, *moving_old.values()])
            for position in travel_samples:
                translation = position * axis
                spring = spring_envelope(start, end + translation, axis)
                posed_moving = moving.translate(translation)
                spring_hits = overlaps(spring, {**fitted_fixed, **bounded({"moving_unit": posed_moving})})
                moving_hits = overlaps(posed_moving, fitted_fixed)
                if spring_hits or moving_hits:
                    row[tag + "_first_block"] = dict(
                        travel_mm=float(position), spring_hits=spring_hits, moving_hits=moving_hits)
                    break
            else:
                sources = {}
                for name, solid in (("adapter", new_adapter), ("cap", new_cap),
                                    ("fixed_pin", fixed_pin), ("moving_pin", moving_pin)):
                    for tolerance in (1e-6, 1e-5, 1e-4):
                        try:
                            _, meta = export(solid.as_original().simplify(tolerance),
                                             folder / f"series_spring_{name}.obj", origin)
                        except (AssertionError, ValueError):
                            continue
                        sources[name] = meta
                        break
                    else:
                        row[tag + "_export_failure"] = name
                        break
                else:
                    chosen = dict(**row, pin_axis_source=pin_axis.tolist(), sources=sources,
                                  fixed_eye_source_mm=start.tolist(), moving_eye_neutral_source_mm=end.tolist(),
                                  neutral_pin_span_mm=front - rear, peak_pin_span_mm=front - rear + travel)
                    break
        if chosen:
            break
    save(folder / "elbow_series_spring_mount.json", dict(
        generator=provenance(Path(__file__).resolve()), guide=provenance(guide_file),
        retention=retain, spring_mount=chosen, trials=trials, neutral_pose_travel_samples=len(travel_samples),
        travel_mm=travel, coil_envelope_diameter_mm=6.5, assumed_spring_stiffness_N_per_m=73.,
        assumed_preload_N=.25, full_extension_elastic_force_N=.25 + .073 * travel,
        scope="Spring is a geometric coil-body and hook-tail envelope, not a chosen hardware spring. "
        "Spring hooks, threaded engagement, pin retention, tooling, strength and shoulder motion "
        "are unresolved. Retention uses a nominal integral flanged shaft and proximal screw. "
        "Checks use the neutral elbow pose with prescribed slide travel."))
    print("retention", {k: v["intersections_mm3"] for k, v in retain.items()}, flush=True)
    print("spring candidate", chosen, "lanes tried", len(trials), flush=True)


if __name__ == "__main__":
    main()
