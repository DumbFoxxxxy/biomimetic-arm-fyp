"""Compare historical ferrule locations with the current elbow port and crimp solids."""

import mujoco
import numpy as np

from audit_s1_series_elbow_spring_integration import load_candidate
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def axis_measurement(part, point, axis):
    vertices, _ = mesh_geom.manifold_mesh_to_arrays(part.to_mesh64())
    axial = (vertices - point) @ axis
    lo, hi = float(axial.min()), float(axial.max())
    # Rays on two transverse directions classify whether stock surrounds the
    # cable axis locally. This does not infer a ferrule shoulder from a bore.
    radial = np.cross(axis, np.eye(3)[np.argmin(abs(axis))])
    radial /= np.linalg.norm(radial)
    transverse = np.array([radial, np.cross(axis, radial)])
    samples = []
    for station in np.unique(np.r_[0., np.arange(np.ceil(lo), hi, 1.)]):
        center = point + station * axis
        entries = []
        for direction in np.vstack((transverse, -transverse)):
            hits = part.ray_cast(center, center + 15. * direction)
            entries.append(float((np.asarray(hits[0].position) - center) @ direction) if hits else None)
        samples.append(dict(station_mm=float(station), first_radial_surface_mm=entries))
    return dict(axial_bounds_from_reference_mm=[lo, hi], radial_probe_samples=samples)


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    register_ref = provenance(DEST / "s1_cable_anchor_register.json")
    register = load_reference(register_ref)
    prior, extended, _, _, sources, owners, revisions, retention = load_candidate(folder)
    assembly_ref = load_reference(prior["passive"])["assembly"]
    assembly = load_reference(assembly_ref)
    sources = {**sources, **revisions}
    coupling = load_reference(next(row for row in register["channels"]
                                   if row["id"] == "FJ1477")["route_reference"])
    if coupling["model"] != extended["model"]:
        raise ValueError("Anchor and slide motion reports refer to different models")
    model_path = ROOT / extended["model"]["file"]
    if provenance(model_path) != extended["model"]:
        raise ValueError("Drive model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    slide_axis = model.joint("FJ1477_series_travel").axis
    slide_mm = np.array([float(row["slide_mm"]) for row in checked_trace(coupling["trace"])])
    channels = {}
    for muscle in assembly["muscles"].values():
        head = muscle["anatomical_id"]
        entry_record = next(row for row in register["channels"] if row["id"] == head)
        ferrule_ref = entry_record["historical_ferrule_candidate"]["reference"]
        old = next(row for row in load_reference(ferrule_ref)["routes"] if row["id"] == head)
        path = np.array(muscle["wire_path_source_mm"])
        first, last = path[1], path[-2]
        np.testing.assert_allclose([first, last], [entry_record["entry_xyz_mm"], entry_record["terminal_xyz_mm"]],
                                   atol=1e-6, rtol=0.)
        cap_key = f"elbow:{head}_proximal_rigid_cap"
        barrel_key = f"elbow:{head}_active_wire_crimp_barrel"
        cap, barrel = (load_component(sources[key]) for key in (cap_key, barrel_key))
        inlet_axis = path[0] - first
        inlet_axis /= np.linalg.norm(inlet_axis)
        exit_axis = path[-1] - last
        exit_axis /= np.linalg.norm(exit_axis)
        old_seat = np.array(old["physical_sheath_reaction_seat_source_mm"])
        lead = old_seat - first
        shift = float(lead @ inlet_axis)
        old_bearing = np.array(old["physical_ferrule_bearing_source_mm"])
        seat_measure = axis_measurement(cap, old_seat, inlet_axis)
        barrel_measure = axis_measurement(barrel, last, exit_axis)
        inlet = A.D.outer_envelope(np.array([first, old_seat]), .30)
        termination = A.D.outer_envelope(path[-2:], .30)
        channels[head] = dict(
            old_ferrule_source=ferrule_ref, inlet_cap=sources[cap_key], inlet_body=owners[cap_key],
            crimp_body=owners[barrel_key], crimp_barrel=sources[barrel_key],
            inlet_axis_source=inlet_axis.tolist(), old_reaction_xyz_mm=old_seat.tolist(),
            old_reaction_axial_distance_from_model_point_mm=shift,
            old_reaction_offset_from_wire_axis_mm=float(np.linalg.norm(lead - shift * inlet_axis)),
            old_reaction_distance_to_cap_axial_end_mm=-seat_measure["axial_bounds_from_reference_mm"][1],
            old_bearing_distance_past_cap_axial_end_mm=float(
                (old_bearing - old_seat) @ inlet_axis - seat_measure["axial_bounds_from_reference_mm"][1]),
            inlet_cap_probe=seat_measure, terminal_crimp_probe=barrel_measure,
            inlet_lead_cap_hits_mm3=overlaps(inlet, bounded({cap_key: cap})),
            terminal_lead_crimp_overlap_mm3=float((termination ^ barrel).volume()),
            crimp_bore_radial_samples_are_geometric=True,
        )
        if owners[cap_key] == "FJ1477_carriage":
            channels[head]["carriage_relative_to_humerus"] = dict(
                axis_source=slide_axis.tolist(),
                travel_mm=[float(slide_mm.min()), float(slide_mm.max())],
                historical_seat_xyz_at_peak_if_carriage_owned_mm=(old_seat + slide_axis * slide_mm.max()).tolist(),
            )
    save(folder / "elbow_physical_anchor_audit.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_elbow_cable_anchors.py"), register=register_ref,
        retention=retention, model=extended["model"], wire_assembly=assembly_ref, channels=channels,
        scope="Neutral pose of the extended-rail candidate. Axial material bounds are projections, "
              "not a finding of a bearing face. The named crimp barrels are geometric solids, "
              "not evidence of cable retention. Historical ferrules remain reference only.",
    ))
    for name, row in channels.items():
        print(name, row["inlet_body"], "seat offset", row["old_reaction_offset_from_wire_axis_mm"],
              "cap end-to-seat", row["old_reaction_distance_to_cap_axial_end_mm"],
              "barrel span", row["terminal_crimp_probe"]["axial_bounds_from_reference_mm"], flush=True)


if __name__ == "__main__":
    main()
