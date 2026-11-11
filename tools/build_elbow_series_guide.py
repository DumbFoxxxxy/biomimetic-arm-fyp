"""Build and screen a local FJ1477 series-slider guide on the relieved adapter."""

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np

from audit_elbow_cap_travel import load_geometry
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, mesh_geom, rod, save
from study_elbow_series_package import ADAPTER, PREFIX, bounded, maximum_rigid_travel, overlaps, union


@dataclass
class Lane:
    report: dict
    base: np.ndarray
    shaft: object
    collar: object
    socket_stock: object
    socket_bore: object
    cap_addition: object
    sleeve: object
    bridge: object
    protected_socket: object


def unit(vector):
    return vector / np.linalg.norm(vector)


def lane_parts(origin, axis, offset, travel):
    base = origin + offset

    def cylinder(start, stop, radius):
        return rod(base + start * axis, base + stop * axis, radius)

    guide = cylinder(-27., travel - 6., 2.)
    socket = cylinder(-27., -17., 4.)
    collar = cylinder(travel - 10., travel - 6., 4.) - cylinder(travel - 10., travel - 6., 2.)
    sleeve = cylinder(-17., -10., 4.) - cylinder(-17., -10., 2.18)
    return base, guide, socket, collar, sleeve


def attachment_link(base, axis, adapter_vertices):
    centre = base - 22. * axis
    nearest = adapter_vertices[np.linalg.norm(adapter_vertices - centre, axis=1).argmin()]
    distance = float(np.linalg.norm(nearest - centre))
    if distance > 13. or distance < 1e-8:
        return None, distance
    into_adapter = unit(nearest - centre)
    return rod(centre, nearest + 2. * into_adapter, 3.), distance


def moving_link(origin, base, axis, cap, existing_fasteners):
    transverse = unit(base - origin)
    root = origin - 13.5 * axis + 6. * transverse
    shoulder = base - 13.5 * axis
    for radius in (2.4, 2.8):
        bridge = rod(root, shoulder, radius)
        if (float(A.I.intersection_volume(bridge, cap)) > 1.
                and not overlaps(bridge - cap, existing_fasteners)):
            return bridge
    return None


def trimmed_root(root, moving_parts, travel, axis):
    mf = A.I._manifold()
    swept = union(mf.Manifold.batch_hull([part, part.translate(travel * axis)])
                  for part in moving_parts)
    return root - swept.minkowski_sum(mf.Manifold.cube([.3] * 3, True))


def main():
    folder = DEST.parent / "elbow_local_muscle_load"
    study_file = folder / "elbow_series_package.json"
    study = read(study_file)
    selected = next(row for row in study["relief_trials"]
                    if row.get("exported_adapter_cycle", {}).get("first_failure", True) is None)
    case = next(row for row in study["tilted_trials"]
                if (row["tilt_deg"], row["azimuth_deg"]) == (selected["tilt_deg"], selected["azimuth_deg"]))
    axis = unit(np.asarray(case["direction_source"]))
    origin = np.asarray(study["bore_origin_mm"])
    travel = selected["travel_mm"]
    side = unit(np.cross(axis, [1., 0., 0.]))
    radial = np.cross(axis, side)
    assembly_path = ROOT / study["stop_assembly"]["file"]
    if provenance(assembly_path) != study["stop_assembly"]:
        raise RuntimeError("Stop assembly fingerprint changed")
    check = read(DEST / "elbow_passive_model/model_check.json")
    solids, _ = load_geometry(check, read(assembly_path))
    adapter = load_component(selected["prototype_relief_source"])
    cap_name = PREFIX + "proximal_rigid_cap"
    cap = solids[cap_name][1]
    moving_names = study["proximal_moving_parts"]
    package = union(solids[name][1] for name in moving_names)
    fasteners = bounded({name: solids[name][1] for name in moving_names if name != cap_name})
    other = bounded({name: part for name, (_, part) in solids.items()
                     if name not in moving_names and name != ADAPTER})
    adapter_vertices, _ = mesh_geom.manifold_mesh_to_arrays(adapter.to_mesh64())
    mf = A.I._manifold()
    seat_keepout = [adapter ^ solids[PREFIX + f"bone_{index}_front_washer"][1].hull().minkowski_sum(
        mf.Manifold.cube([3.] * 3, True)) for index in (0, 1)]
    trials = []
    viable = []
    for offset in (13., 16., 19., 22., 25., 28.):
        for angle in range(0, 360, 10):
            direction = np.cos(np.deg2rad(angle)) * radial + np.sin(np.deg2rad(angle)) * side
            base, guide, socket, collar, sleeve = lane_parts(origin, axis, offset * direction, travel)
            info = dict(offset_mm=offset, azimuth_deg=angle)
            trials.append(info)
            shaft_other = overlaps(union([guide, socket, collar]), other)
            if shaft_other:
                info["fixed_obstacles"] = shaft_other
                continue
            link, distance = attachment_link(base, axis, adapter_vertices)
            if link is None:
                info["root_link_distance_mm"] = distance
                continue
            added_stock = union([socket, link])
            bore = rod(base - 27.5 * axis, base - 16.95 * axis, 2.04)
            root = (adapter + added_stock) - bore
            if (float(A.I.intersection_volume(added_stock, adapter)) < 2.
                    or len([part for part in root.decompose() if part.volume() > 1e-5]) != 1
                    or any((seat - root).volume() > .001 for seat in seat_keepout)):
                info["root_connection_or_seat"] = True
                continue
            addition_hits = overlaps(root - adapter, other)
            fixed_guide_hits = overlaps(guide, bounded({ADAPTER: root}))
            if addition_hits or fixed_guide_hits:
                info["root_overlap"] = addition_hits or fixed_guide_hits
                continue
            bridge = moving_link(origin, base, axis, cap, fasteners)
            if bridge is None:
                info["moving_cap_attachment"] = False
                continue
            moving_addition = (sleeve + bridge) - rod(base - 17.1 * axis, base - 9.9 * axis, 2.18)
            moving = cap + moving_addition
            if len([part for part in moving.decompose() if part.volume() > 1e-5]) != 1:
                info["moving_cap_connected"] = False
                continue
            seat_shaft = root ^ rod(base - 27. * axis, base - 17.4 * axis, 3.9)
            fitted_root = trimmed_root(root, [sleeve, bridge], travel, axis)
            if (len([part for part in fitted_root.decompose() if part.volume() > 1e-5]) != 1
                    or any((seat - fitted_root).volume() > .001 for seat in [*seat_keepout, seat_shaft])):
                info["added_clearance_breaks_support"] = True
                continue
            info["added_clearance_removed_mm3"] = float((root - fitted_root).volume())
            moving_package = package + moving_addition
            root = fitted_root
            fixed = {**other, **bounded({ADAPTER: root, "rail": guide, "end_collar": collar})}
            blocked = maximum_rigid_travel(moving_package, fixed, travel, axis)
            if blocked:
                info["travel_block"] = blocked
                continue
            info["single_lane_clear"] = True
            viable.append(Lane(info, base, guide, collar, added_stock, bore, moving_addition,
                               sleeve, bridge, seat_shaft))
    print("Single guide lanes clear", len(viable), "out of", len(trials), flush=True)
    pairs = []
    choice = None
    for left, right in combinations(viable, 2):
        centre_distance = float(np.linalg.norm(left.base - right.base))
        if centre_distance < 10.:
            continue
        root = union([adapter, left.socket_stock, right.socket_stock]) - union([left.socket_bore, right.socket_bore])
        root = trimmed_root(root, [left.sleeve, left.bridge, right.sleeve, right.bridge], travel, axis)
        cap_body = union([cap, left.cap_addition, right.cap_addition])
        guides = {"guide_a": left.shaft, "guide_b": right.shaft, "collar_a": left.collar, "collar_b": right.collar}
        row = dict(lanes=[left.report, right.report], centre_distance_mm=centre_distance)
        pairs.append(row)
        if (len([part for part in root.decompose() if part.volume() > 1e-5]) != 1
                or any((seat - root).volume() > .001
                       for seat in [*seat_keepout, left.protected_socket, right.protected_socket])):
            row["pair_clearance_breaks_support"] = True
            continue
        blocked = maximum_rigid_travel(union([package, cap_body]), {**other, **bounded({
            ADAPTER: root, **guides})}, travel, axis)
        if blocked:
            row["travel_block"] = blocked
            continue
        choice = row
        choice["sources"] = {}
        for name, solid in {"adapter": root, "moving_cap": cap_body, **guides}.items():
            _, choice["sources"][name] = export(solid.as_original().simplify(1e-6),
                                                 folder / f"series_guide_{name}.obj", origin)
        choice["rail_origins_source_mm"] = [left.base.tolist(), right.base.tolist()]
        break
    save(folder / "elbow_series_guide.json", dict(
        generator=provenance(Path(__file__).resolve()), package=provenance(study_file),
        stop_assembly=study["stop_assembly"], travel_mm=travel, axis_source=axis.tolist(),
        trials=trials, pair_trials=pairs, selected=choice, shaft_diameter_mm=4.,
        bushing_bore_diameter_mm=4.36, local_screen_step_max_mm=1.,
        scope="Neutral pose shaft packaging. Two parallel shafts intended to carry transverse reaction "
        "and prevent rotation. Added sockets connect to relieved adapter, sleeves to original cap. "
        "Nominal collars define end travel. No spring, wire or off-neutral collision acceptance."))
    if choice is not None:
        print("Paired guide exported", choice["rail_origins_source_mm"], flush=True)


if __name__ == "__main__":
    main()
