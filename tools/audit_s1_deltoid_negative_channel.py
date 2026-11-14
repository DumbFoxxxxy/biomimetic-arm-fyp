"""Screen FJ1468 centreline deformation with both end tangents retained."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_coupled_scapula import ShoulderAssembly
from study_s1_infraspinatus_wire_travel import route_metrics
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_negative_path_local_stroke.json"
    stroke = read(source)
    integration = load_reference(stroke["path"])
    probes = integration["cases"]["FJ1504_local_relief"]["negative_half_degree_endpoint_probes"]
    path = min((p["extension"] for p in probes if p.get("extension", {}).get("outbound_and_inbound_clear")),
               key=lambda row: row["last_angle_deg"])
    route = next(row for row in load_reference(stroke["inventory"])["routes"] if row["id"] == "FJ1468")
    points = np.array(route["local_route"]["points_source_mm"])
    length, _, _ = route_metrics(points)
    station = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    origins = load_reference(stroke["origins"])["neutral_origins_m"]
    model = ShoulderAssembly({}, {}, origins)
    design_dir = DEST.parent / "shoulder_design"
    proximal_file = design_dir / "s1_anterior_mount/assembly_check.json"
    distal_file = design_dir / "s1_deltoid_distal_mount/assembly_check.json"
    proximal = read(proximal_file)["sources"]
    distal = read(distal_file)["sources"]
    rigid = {
        "clavicle": {key: load_component(proximal[key]) for key in (
            "bone:clavicle", "cap:FJ1468/origin", "origin:FJ1468/front_support")},
        "humerus": {key: load_component(distal[key]) for key in (
            "bone:humerus", "cap:FJ1468/distal", "distal:FJ1468/front_support")},
    }
    results = []
    for fixed_length in (10., 20., 30.):
        u = np.clip((station - fixed_length) / (length - 2. * fixed_length), 0., 1.)
        weight = (u * u * (3. - 2. * u))[:, None]
        records = []
        for phase in ("outbound", "inbound"):
            for index, segment in enumerate(path[phase]):
                for frame_index, frame in enumerate(segment):
                    angle = frame["scapular_coronal_deg"]
                    sc = np.array(frame["SC_rotvec_deg"])
                    clavicle, humerus, _ = model.kinematics(angle, sc)
                    pa = points @ clavicle[:, :3].T + clavicle[:, 3]
                    pb = points @ humerus[:, :3].T + humerus[:, 3]
                    centreline = (1. - weight) * pa + weight * pb
                    wire_length, turns, radii = route_metrics(centreline)
                    row = dict(phase=phase, angle_deg=angle, arc_mm=wire_length,
                               delta_mm=wire_length - length, min_sampled_radius_mm=float(radii.min()),
                               largest_turn_deg=float(np.rad2deg(turns.max())),
                               sharpest_radius_vertex=int(radii.argmin()) + 1)
                    if frame_index == 0 or frame_index == len(segment) - 1:
                        obstacles = {key: part.transform(transform)
                                     for owner, transform in (("clavicle", clavicle), ("humerus", humerus))
                                     for key, part in rigid[owner].items()}
                        wire = A.D.outer_envelope(centreline, .3)
                        row["wire_contact_at_node_mm3"] = overlaps(wire, bounded(obstacles))
                    records.append(row)
        results.append(dict(end_region_mm=fixed_length, records=records))
        print(fixed_length, min(row["delta_mm"] for row in records),
              min(row["min_sampled_radius_mm"] for row in records),
              records[0]["wire_contact_at_node_mm3"],
              sum(bool(row.get("wire_contact_at_node_mm3")) for row in records), flush=True)
    save(folder / "FJ1468_negative_channel_screen.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_deltoid_negative_channel.py"),
        path_stroke=provenance(source), proximal_mount=provenance(proximal_file),
        distal_mount=provenance(distal_file), neutral_length_mm=length, results=results,
        scope="Smooth coordinate blend of the inherited FJ1468 polyline, with rigid port "
              "sections attached to the clavicle and humerus. Separate end-mount candidates "
              "are used as local obstacles at path nodes. No loaded equilibrium, continuous "
              "swept collision proof, or TPU solid deformation model.",
    ))


if __name__ == "__main__":
    main()
