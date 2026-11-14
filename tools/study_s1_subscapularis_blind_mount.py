"""Screen FJ1504 same-side blind mounts against FJ1500 belly insertion."""

from itertools import product

import numpy as np

from audit_elbow_slider_short_guides import export_lattice
from audit_s1_infraspinatus_split_installation import insertion
from audit_s1_series_elbow_integration import assert_source
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    audit_file = folder / "FJ1504_hardware_after_split_belly.json"
    neighbor = read(audit_file)
    split_audit = load_reference(neighbor["split_audit"])
    split = load_reference(split_audit["split"])
    side = load_reference(load_reference(split["parent"])["parent"])
    removable = load_reference(side["source"])
    feed = load_reference(removable["integrated_sleeve"])
    assembly = load_reference(feed["mount"])
    revision = read(assert_source(assembly["scapular_mount_revision"]))
    routes = read(assert_source(revision["route_reference"]))["routes"]
    fixation = read(assert_source(revision["fixation_geometry"]))
    route = next(row for row in routes if row["head"] == "FJ1504")
    holes = next(row["holes"] for row in fixation["rows"] if row["head"] == "FJ1504")
    points = np.array([hole["axis_point_source_mm"] for hole in holes])
    rear = np.array([hole["bone_crossings"][1]["position_along_cap_normal_mm"] for hole in holes])
    normal = np.array(route["ports"]["origin"]["normal_source"])
    normal /= np.linalg.norm(normal)
    origin = np.array(route["ports"]["origin"]["inner_face_source_mm"])
    prefix = "origin:FJ1504/"
    removed = {prefix + name for name in ("front_support", "backing_plate", "bolt_0", "bolt_1",
                                          "nut_0", "nut_1", "back_washer_0", "back_washer_1")}
    scene = {key: load_component(meta) for key, meta in assembly["sources"].items()}
    scene["cap:FJ1500/origin"] = load_component(removable["selected"]["cap"])
    bone = scene["bone:scapula"]
    front = scene[prefix + "front_support"]
    targets = {key: part for key, part in scene.items() if key not in removed and key != "bone:scapula"}
    wires = {row["head"]: A.D.outer_envelope(np.asarray(row["full_wire_path_source_mm"]), .3)
             for row in routes if row["head"] not in ("FJ1500", "FJ1504")}
    wires["FJ1500"] = A.D.outer_envelope(
        np.array(load_reference(feed["route"])["exported"]["trial"]["new_path_source_mm"]), .3)
    muscle_revision = load_reference(assembly["FJ1504_TPU_revision"])
    cap_revision = load_reference(assembly["FJ1504_cap_revision"])
    for label, path in (("FJ1504_TPU_wire", muscle_revision["selected"]["path_source_mm"]),
                        ("FJ1504_cap_wire", cap_revision["selected"]["local_wire_path_source_mm"])):
        wires[label] = A.D.outer_envelope(np.asarray(path), .3)
    # The cap and support belong to the same assembly used in the insertion audit.
    soft = load_component(assembly["flexible_geometry_sources"]["tpu:subscapularis"])
    fixed = bounded({**targets, **wires, "FJ1504/TPU": soft})
    bellies = [(row, load_component(row["belly_mesh"])) for row in split_audit["cases"]]
    trials = []
    selected = None
    for radius, inset in product((3.4, 3., 2.6, 2.2), (.3, .7, 1.1)):
        floor = rear + inset
        boss = union(rod(point + z * normal, point - 13.8 * normal, radius)
                     for point, z in zip(points, floor))
        bores = union(rod(point + (z + .25) * normal, point - 12.5 * normal, 1.25)
                      for point, z in zip(points, floor))
        revised = ((bone + front + boss) - bores).as_original().simplify(1e-6)
        added = revised - bone - front
        hits = overlaps(added, fixed)
        connected = sum(piece.volume() > 1e-6 for piece in revised.decompose())
        row = dict(radius_mm=radius, inset_from_centreline_rear_surface_mm=inset,
                   floor_from_cap_mm=floor.tolist(), bone_components=connected,
                   new_stock_overlap_mm3=hits, added_stock_mm3=float(added.volume()))
        trials.append(row)
        if hits or connected != 1:
            continue
        bolts = {}
        for index, (point, z) in enumerate(zip(points, floor)):
            head_start = -4.75
            shaft = rod(point + (z + .45) * normal, point + (head_start + .1) * normal, 1.19)
            head = rod(point + head_start * normal, point + (head_start + 3.) * normal, 2.75)
            bolts[prefix + f"bolt_{index}"] = shaft + head
        solids = {"bone:scapula": revised, **bolts}
        bolt_targets = bounded({**targets, **wires, "FJ1504/TPU": soft, "scapula": revised})
        row["bolt_hits_mm3"] = {key: overlaps(part, bolt_targets) for key, part in bolts.items()}
        if any(row["bolt_hits_mm3"].values()):
            continue
        def belly_checks(parts):
            world = bounded({**targets, **parts})
            return [dict(station_mm=case["split_station_mm"], traces=[
                insertion(belly, world, np.asarray(trace["offset_direction_source"]), 8.)
                for trace in case["results"]["FJ1504_origin_deferred"]["belly_insertion"]
                if not trace["hits"]
            ]) for case, belly in bellies]

        row["belly_insertion"] = belly_checks(solids)
        if not any(not trace["hits"] for case in row["belly_insertion"] for trace in case["traces"]):
            continue
        reloaded = {}
        meshes = {}
        for key, part in solids.items():
            writer = export_lattice if key == "bone:scapula" else export
            name = key.replace(":", "_").split("/")[-1]
            reloaded[key], meshes[key] = writer(part, folder / f"FJ1504_blind_{name}.obj", origin)
        row["readback_belly_insertion"] = belly_checks(reloaded)
        row["readback_new_stock_overlap_mm3"] = overlaps(reloaded["bone:scapula"] - bone - front, fixed)
        row["readback_bolt_hits_mm3"] = {key: overlaps(reloaded[key], bounded({
            **targets, **wires, "FJ1504/TPU": soft, "scapula": reloaded["bone:scapula"]}))
            for key in bolts}
        if (row["readback_new_stock_overlap_mm3"] or any(row["readback_bolt_hits_mm3"].values())
                or not any(not t["hits"] for c in row["readback_belly_insertion"] for t in c["traces"])):
            raise ValueError("Blind mount did not reproduce after export")
        selected = dict(**row, meshes=meshes, nominal_engagement_mm=(-13.8 - floor - .45).tolist())
        selected["bolt_entry"] = {}
        for label, world in (("on_thorax", targets), ("without_thorax_reference", {
                key: part for key, part in targets.items() if not key.startswith("thorax_reference:")})):
            assembled = {**world, "scapula": reloaded["bone:scapula"], "FJ1504/TPU": soft}
            traces = {}
            for key in bolts:
                traces[key] = insertion(reloaded[key], bounded(assembled), normal, 25.)
                assembled[key] = reloaded[key]
            selected["bolt_entry"][label] = traces
        break

    save(folder / "FJ1504_blind_mount_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_subscapularis_blind_mount.py"),
        parent=provenance(audit_file), retired_hardware=sorted(removed), trials=trials, selected=selected,
        scope="Nominal same-side hardware at neutral. The drilled substrate is an engineering geometry. "
              "Thread, insert retention, scapular fabrication, and loads are unvalidated.",
    ))
    print("trials", len(trials), "selected", bool(selected), flush=True)


if __name__ == "__main__":
    main()
