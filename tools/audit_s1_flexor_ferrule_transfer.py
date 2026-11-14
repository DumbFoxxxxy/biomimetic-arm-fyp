"""Screen legacy ferrule seats and hardware on the extended-rail elbow."""

import argparse
import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_spring_integration import load_candidate
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, export, load_component, mesh_geom, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channel", choices=("FJ1486", "FJ1477"), default="FJ1486")
    parser.add_argument("--relief", action="store_true", help="Clear the FJ1477 adapter sweep from added stock")
    args = parser.parse_args()
    if args.relief and args.channel != "FJ1477":
        parser.error("--relief requires --channel FJ1477")
    folder = DEST.parent / "integration_review/s1_series_elbow"
    anchor_ref = provenance(folder / "elbow_physical_anchor_audit.json")
    anchor = load_reference(anchor_ref)
    _, extended, _, _, sources, owners, revisions, retention = load_candidate(folder)
    sources = {**sources, **revisions}
    head = args.channel
    channel = anchor["channels"][head]
    cap_key = f"elbow:{head}_proximal_rigid_cap"
    owner = owners[cap_key]
    if channel["inlet_cap"] != sources[cap_key] or anchor["retention"] != retention:
        raise ValueError("Anchor report refers to different elbow hardware")
    legacy_ports = load_reference(channel["old_ferrule_source"])
    legacy = load_reference(legacy_ports["assembly"])
    legacy_parent = legacy["parent"]
    legacy_base = load_reference(legacy_parent)
    old_cap = load_component(legacy_base["sources"][cap_key])
    legacy_cap = load_component(legacy["sources"][cap_key])
    current = {key: load_component(meta) for key, meta in sources.items()}
    added = legacy_cap - old_cap
    removed = old_cap - legacy_cap
    revised = ((current[cap_key] + added) - removed).as_original().simplify(1e-7)
    axis = np.asarray(channel["inlet_axis_source"])
    seat = np.asarray(channel["old_reaction_xyz_mm"])
    hardware = {key: load_component(meta) for key, meta in legacy["sources"].items()
                if key.startswith(f"supply:{head}_")}
    current_targets = {key: part for key, part in current.items() if key != cap_key}
    study_ref = provenance(folder / "nominal_spring_spool_seating.json")
    study = load_reference(study_ref)
    case = next(row for row in study["cases"] if row["case"] == "seating_k5")
    if study["drive_model"] != extended["model"]:
        raise ValueError("Trace and assembly model differ")
    trace = checked_trace(case["trace"])
    model_path = ROOT / extended["model"]["file"]
    if provenance(model_path) != extended["model"]:
        raise ValueError("Drive model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    relief = None
    suffix = "_relief" if args.relief else ""
    if args.relief:
        adapter = current["elbow:FJ1477_bone_adapter"]
        stock = revised - current[cap_key]
        path = model.joint("FJ1477_series_travel").axis
        limits = model.joint("FJ1477_series_travel").range * 1000.
        samples = []
        for shift in np.linspace(*limits, int(np.ceil(np.ptp(limits) / .1)) + 1):
            hit = stock ^ adapter.translate(-shift * path)
            if hit.volume() > 1e-8:
                vertices, _ = mesh_geom.manifold_mesh_to_arrays(hit.to_mesh64())
                samples.append(vertices)
        if not samples:
            raise ValueError("The current parts no longer show the expected adapter interference")
        mf = A.I._manifold()
        cutter = mf.Manifold.hull_points(np.ascontiguousarray(np.vstack(samples)))
        cutter = cutter.minkowski_sum(mf.Manifold.cube([.30] * 3, True))
        trimmed = stock ^ cutter
        guard = rod(seat - 8. * axis, seat + 8. * axis, 4.)
        guard_hit = float((trimmed ^ guard).volume())
        hardware_hit = overlaps(trimmed, bounded({
            key: part.minkowski_sum(mf.Manifold.cube([.50] * 3, True))
            for key, part in hardware.items()}))
        if guard_hit > .001 or hardware_hit:
            raise ValueError(f"Relief enters protected ferrule or fastener stock: {guard_hit:g}, {hardware_hit}")
        revised = (revised - trimmed).as_original().simplify(1e-7)
        relief = dict(removed_added_stock_mm3=float(trimmed.volume()), guarded_radius_mm=4.,
                      guard_hit_mm3=guard_hit, hardware_guard_hits_mm3=hardware_hit,
                      nominal_slide_limits_mm=limits.tolist(), maximum_sampling_step_mm=.1,
                      nominal_full_slide_hit_samples=len(samples))
    revised, cap_mesh = export(revised, folder / f"{head}_ferrule_transfer_cap{suffix}.obj", seat)
    parts = {cap_key: revised, **hardware}
    local = bounded(parts)
    own_hits = {key: overlaps(part, {other: val for other, val in local.items() if other < key})
                for key, part in parts.items()}
    if args.relief:
        stock = revised - current[cap_key]
        residual = [float((stock ^ adapter.translate(-shift * path)).volume())
                    for shift in np.linspace(*limits, int(np.ceil(np.ptp(limits) / .1)) + 1)]
        relief.update(full_slide_readback_hit_samples=sum(value > .001 for value in residual),
                      full_slide_peak_readback_overlap_mm3=max(residual))
    first = {}
    hit_counts = {}
    max_hits = {}
    peak_cap_contact = None
    cap_added_stock = revised - current[cap_key]
    for row in trace:
        data.qpos[:] = [float(row[f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        frames = transforms(model, data, neutral)
        posed = bounded({key: part.transform(frames[owners[key]]) if owners[key] in frames else part
                        for key, part in current_targets.items()})
        for key, part in parts.items():
            new_hits = overlaps(part.transform(frames[owner]), posed)
            if key == cap_key and new_hits:
                old_hits = overlaps(current[cap_key].transform(frames[owner]), posed)
                new_hits = {name: volume for name, volume in new_hits.items()
                            if volume > old_hits.get(name, 0.) + .001}
            for name, volume in new_hits.items():
                first.setdefault(key, {}).setdefault(name, dict(time_s=float(row["time_s"]), overlap_mm3=volume))
                pair = f"{key} / {name}"
                hit_counts[pair] = hit_counts.get(pair, 0) + 1
                max_hits[pair] = max(volume, max_hits.get(pair, 0.))
                if key == cap_key and (peak_cap_contact is None or volume > peak_cap_contact["overlap_mm3"]):
                    frame = frames[owner]
                    target = posed[name][0]
                    added_hit = cap_added_stock.transform(frame) ^ target
                    inverse = np.c_[frame[:, :3].T, -frame[:, :3].T @ frame[:, 3]]
                    peak_cap_contact = dict(target=name, time_s=float(row["time_s"]), slide_mm=float(row["slide_mm"]),
                                            overlap_mm3=volume, new_stock_overlap_mm3=float(added_hit.volume()),
                                            new_stock_contact_bounds_neutral_mm=list(added_hit.transform(inverse).bounding_box())
                                            if not added_hit.is_empty() else None)
    wire = rod(seat - 12. * axis, seat + 8. * axis, .30)
    wire_hits = overlaps(wire, local)
    save(folder / f"{head}_ferrule_transfer_audit{suffix}.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_flexor_ferrule_transfer.py"), anchors=anchor_ref,
        legacy_ports=channel["old_ferrule_source"], model=extended["model"], trace=case["trace"],
        sampled_frames=len(trace), tested_hardware=list(parts), transferred_cap_mesh=cap_mesh,
        hardware_sources={key: legacy["sources"][key] for key in hardware},
        sheath_reaction_seat_source_mm=seat.tolist(), owner=owner,
        transferred_cap_component_volumes_mm3=sorted((float(piece.volume()) for piece in revised.decompose()), reverse=True),
        cap_added_mm3=float((revised - current[cap_key]).volume()),
        cap_removed_mm3=float((current[cap_key] - revised).volume()),
        internal_hardware_hits_mm3={key: value for key, value in own_hits.items() if value},
        new_or_increased_external_pair_first_hits=first, steel_probe_hits_mm3=wire_hits,
        external_pair_hit_counts=hit_counts, external_pair_peak_overlap_mm3=max_hits,
        peak_cap_contact=peak_cap_contact, relief=relief,
        scope="Exact legacy cap material delta and matching rigid hardware against saved dynamics. "
              "Nominal steel probe is local. No external sheath motion, TPU or retention loads checked.",
    ))
    print(f"Checked {len(trace)} frames; external hits {first}; steel probe {wire_hits}", flush=True)


if __name__ == "__main__":
    main()
