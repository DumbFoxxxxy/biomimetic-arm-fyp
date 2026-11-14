"""Screen right AC hardware access with the scapula temporarily rotated."""

from itertools import product

import numpy as np

from audit_s1_infraspinatus_split_installation import insertion
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_AC_fastener_access.json"
    access = read(source)
    closure = load_reference(access["source"])
    entry = load_reference(closure["source"])
    stage = entry["cases"]["clavicle_seated_AC_hardware_deferred"]
    pose_file = DEST / "elbow_passive_model/model_check.json"
    pivot = 1000. * np.asarray(read(pose_file)["neutral_origins_m"]["scapula"])
    assembly = {key: load_component(row) for key, row in entry["sources"].items()}
    hardware = {key.removeprefix("AC:AC_"): load_component(row)
                for key, row in closure["hardware_sources"].items()}
    axis = np.asarray(closure["axis_source"])
    bolt = hardware["bolt_right"]
    verts = np.asarray(bolt.to_mesh64().vert_properties)[:, :3]
    stations = verts @ axis
    driver_start = verts[stations >= stations.max() - .001].mean(axis=0) + .05 * axis
    driver = rod(driver_start, driver_start + 20. * axis, 1.5, 48)
    fixed = {key: assembly[key] for key in stage["fixed"]}
    mobile = {key: assembly[key] for key in stage["moving"]}
    mobile.update(hardware)
    fixed_bounds = bounded(fixed)
    radial = A.D.unit(np.cross(axis, np.eye(3)[np.argmin(abs(axis))]))
    tangent = np.cross(axis, radial)
    obstacles = bounded({**assembly, **{key: part for key, part in hardware.items()
                                        if not key.endswith("_right") or key == "cap_right"}})
    offsets = []
    for radius, azimuth in product((2., 4., 6., 8., 10.), np.linspace(0., 2. * np.pi, 24, endpoint=False)):
        delta = radius * (np.cos(azimuth) * radial + np.sin(azimuth) * tangent)
        point = driver_start + delta
        driver_hits = overlaps(rod(point, point + 20. * axis, 1.5, 48), obstacles)
        # The original bolt translated rigidly is a diagnostic for whether the
        # same through-bolt can simply move. It does not add a bearing or bore.
        fastener_hits = overlaps(bolt.translate(delta), obstacles)
        offsets.append(dict(offset_source_mm=delta.tolist(), radius_mm=radius,
                            driver_hits_mm3=driver_hits, fastener_hits_mm3=fastener_hits))
        if not driver_hits:
            # Locate possible near-side bearing stock independently of the
            # original through-bolt length. Each station is a short axial slice.
            parts = bounded({key: obstacles[key][0]
                             for key in ("bone:scapula", "bone:clavicle", "cap_right")})
            offsets[-1]["inboard_stock_slices"] = [
                dict(depth_mm=float(depth), stock_mm3=overlaps(
                    rod(point - (depth + .25) * axis, point - depth * axis, 1.5, 32), parts))
                for depth in np.arange(0., 16., .25)
            ]
    results = []
    for coord, angle in product(np.eye(3), (0., -2., 2., -5., 5., -10., 10., -15., 15.)):
        if angle == 0. and results:
            continue
        transform = A.I.rotation_about(pivot, coord, angle)
        posed = {key: part.transform(transform) for key, part in mobile.items()}
        rigid_hits = {key: found for key, part in posed.items() if key != "FJ1504/TPU"
                      and (found := overlaps(part, fixed_bounds))}
        soft_hits = overlaps(posed["FJ1504/TPU"], fixed_bounds)
        driver_hits = overlaps(driver.transform(transform), bounded({
            **fixed, **{key: part for key, part in posed.items() if key != "bolt_right"}}))
        row = dict(axis_source=coord.tolist(), angle_deg=angle, rigid_hits_mm3=rigid_hits,
                   soft_hits_mm3=soft_hits, driver_hits_mm3=driver_hits)
        if not rigid_hits and not driver_hits:
            direction = transform[:, :3] @ axis
            targets = bounded({**fixed, **{key: part for key, part in posed.items()
                                          if key != "bolt_right"}})
            row["bolt_entry"] = insertion(posed["bolt_right"], targets, direction, 22.)
        results.append(row)
    save(folder / "FJ1504_AC_assembly_pose_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_ac_assembly_pose.py"),
        source=provenance(source), pose_reference=provenance(pose_file), pivot_source_mm=pivot.tolist(),
        results=results, parallel_offset_diagnostics=offsets,
        scope="Single-axis candidate service poses. AC hardware moves with scapula "
        "as in the assembled reference. Sparse pose screen only; access and no rigid intersection "
        "at a pose do not establish a path there or a closure sequence.",
    ))
    for row in results:
        if not row["rigid_hits_mm3"] or not row["driver_hits_mm3"]:
            print(row["axis_source"], row["angle_deg"], list(row["rigid_hits_mm3"]),
                  list(row["driver_hits_mm3"]), len(row.get("bolt_entry", {}).get("hits", [])))
    print("poses", len(results), "with clear driver", sum(not row["driver_hits_mm3"] for row in results))
    print("parallel offsets", [(r, sum(row["radius_mm"] == r and not row["driver_hits_mm3"] for row in offsets),
                                 sum(row["radius_mm"] == r and not row["driver_hits_mm3"]
                                     and not row["fastener_hits_mm3"] for row in offsets))
                                for r in (2., 4., 6., 8., 10.)])


if __name__ == "__main__":
    main()
