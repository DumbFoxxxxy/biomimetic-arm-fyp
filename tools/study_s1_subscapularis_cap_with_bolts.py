"""Test a preloaded FJ1504 cap and locate its remaining assembly limits."""

import numpy as np

from audit_s1_infraspinatus_split_installation import insertion
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_blind_mount_screen.json"
    blind = read(source)
    split_audit = load_reference(load_reference(blind["parent"])["split_audit"])
    split = load_reference(split_audit["split"])
    side = load_reference(load_reference(split["parent"])["parent"])
    removable = load_reference(side["source"])
    assembly = load_reference(load_reference(removable["integrated_sleeve"])["mount"])
    cap_key = "cap:FJ1504/origin"
    cap = load_component(assembly["sources"][cap_key])
    scene = {key: load_component(meta) for key, meta in assembly["sources"].items()
             if key not in {*blind["retired_hardware"], cap_key, "bone:scapula"}
             and not key.startswith("thorax_reference:")}
    scene["cap:FJ1500/origin"] = load_component(removable["selected"]["cap"])
    scene["bone:scapula"] = load_component(blind["selected"]["meshes"]["bone:scapula"])
    soft = load_component(assembly["flexible_geometry_sources"]["tpu:subscapularis"])
    scene["FJ1504/TPU"] = soft
    bolt_sources = {key: meta for key, meta in blind["selected"]["meshes"].items()
                    if key.startswith("origin:FJ1504/bolt_")}
    normal = np.asarray(next(iter(blind["selected"]["bolt_entry"]["on_thorax"].values()))
                        ["offset_direction_source"])
    bolts = {key: load_component(meta) for key, meta in bolt_sources.items()}
    washers = {key: scene.pop(key) for key in list(scene) if key.startswith("origin:FJ1504/front_washer_")}
    cap_and_washers = cap
    for part in washers.values():
        cap_and_washers += part
    travel = 12.
    cap_world = bounded({cap_key: cap, **washers})
    preloads = {}
    for key, bolt in bolts.items():
        preloads[key] = {
            "from_final_toward_head": insertion(bolt, cap_world, normal, 25.),
            "from_final_toward_tip": insertion(bolt, cap_world, -normal, 12.),
        }

    world = bounded(scene)
    assembly_cases = []
    for shift in (0., 1., 2., 3., 4.):
        carried = {key: part.translate(shift * normal) for key, part in bolts.items()}
        internal = {key: overlaps(part, cap_world) for key, part in carried.items()}
        if any(internal.values()):
            assembly_cases.append(dict(retraction_mm=shift, internal_hits_mm3=internal))
            continue
        group = cap_and_washers
        for part in carried.values():
            group += part
        placement = insertion(group, world, normal, travel)
        cap_seated = bounded({**scene, cap_key: cap, **washers})
        tightening = ({key: insertion(part, cap_seated, normal, shift) for key, part in bolts.items()}
                      if shift else {})
        assembly_cases.append(dict(retraction_mm=shift, internal_hits_mm3=internal,
                                   placement=placement, tightening=tightening))

    cap_trace = insertion(cap_and_washers, world, normal, travel)
    bolt_key = "origin:FJ1504/bolt_1"
    blocked = preloads[bolt_key]["from_final_toward_head"]["hits"]
    relief = None
    if blocked:
        # Use only the local intersection swept by the real bolt. Keep the
        # clearances around the guide lumen and final head seat as separate gates.
        offsets = np.arange(0., max(hit["offset_mm"] for hit in blocked) + .11, .05)
        tool_margin = A.I._manifold().Manifold.cube([.16] * 3, True)
        bolt = bolts[bolt_key]
        contact = []
        for offset in offsets:
            hit = bolt.translate(offset * normal) ^ cap
            if hit.volume() > 1e-8:
                contact.append(hit.minkowski_sum(tool_margin))
        cut = union(contact)
        revised = (cap - cut).as_original().simplify(1e-6)
        removed = cap - revised
        cap_revision = load_reference(assembly["FJ1504_cap_revision"])
        path = np.array(cap_revision["selected"]["local_wire_path_source_mm"])
        wire_guard = A.D.outer_envelope(path, cap_revision["bore_radius_mm"] + .55)
        washer = washers["origin:FJ1504/front_washer_1"]
        bearing = washer.translate(-.25 * normal) ^ cap
        relief = dict(removed_mm3=float(removed.volume()),
                      removed_guide_guard_mm3=float((removed ^ wire_guard).volume()),
                      washer_axial_bearing_probe_shift_mm=.25,
                      washer_to_cap_bearing_witness_mm3=float(bearing.volume()),
                      removed_washer_seat_witness_mm3=float((removed ^ bearing).volume()),
                      components=sum(part.volume() > 1e-6 for part in revised.decompose()))
        revised_targets = bounded({cap_key: revised, **washers})
        relief["bolt_entry"] = insertion(bolt, revised_targets, normal, 25.)
        if (relief["components"] == 1 and relief["removed_guide_guard_mm3"] <= .001
                and relief["washer_to_cap_bearing_witness_mm3"] > .001
                and relief["removed_washer_seat_witness_mm3"] <= .001 and not relief["bolt_entry"]["hits"]):
            origin = np.array(cap_revision["selected"]["cap_inner_source_mm"])
            actual, relief["mesh"] = export(revised, folder / "FJ1504_cap_bolt_entry_relief.obj", origin)
            relief["readback_bolt_entry"] = insertion(bolt, bounded({cap_key: actual, **washers}), normal, 25.)
            group = union([actual, *washers.values(), *bolts.values()])
            relief["readback_group_insertion"] = insertion(group, world, normal, travel)
            thorax = {key: load_component(meta) for key, meta in assembly["sources"].items()
                      if key.startswith("thorax_reference:")}
            relief["readback_group_insertion_on_thorax"] = insertion(
                group, bounded({**scene, **thorax}), normal, travel)
            if relief["readback_bolt_entry"]["hits"] or relief["readback_group_insertion"]["hits"]:
                raise ValueError("Relieved cap assembly did not reproduce after mesh export")
    output = dict(
        generator=provenance(ROOT / "tools/study_s1_subscapularis_cap_with_bolts.py"),
        source=provenance(source), free_cap_preload=preloads,
        cap_and_washers_reference=cap_trace, assembly_cases=assembly_cases, local_cap_relief=relief,
        scope="FJ1504 axial placement without thorax reference, 0.1 mm samples, "
              "including its stationary neutral TPU. Preloaded poses require a separate feasible "
              "preload route. Wire is threaded afterward; this probe does not certify threading.",
    )
    save(folder / "FJ1504_cap_with_bolts_screen.json", output)
    for key, directions in preloads.items():
        print(key, {label: len(trace["hits"]) for label, trace in directions.items()}, flush=True)
    print("cap only", len(cap_trace["hits"]), flush=True)
    for case in assembly_cases:
        print(case["retraction_mm"], case["internal_hits_mm3"],
              len(case.get("placement", {}).get("hits", [])), flush=True)
    if relief is not None:
        print({key: value for key, value in relief.items() if not isinstance(value, dict)}, flush=True)
        print("relief residual", len(relief["bolt_entry"]["hits"]), "mesh", "mesh" in relief, flush=True)


if __name__ == "__main__":
    main()
