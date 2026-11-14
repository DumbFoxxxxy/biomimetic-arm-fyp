"""Measure local attachment travel on the verified negative SC/AC trajectory."""

import csv

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save
from study_s1_coupled_scapula import ShoulderAssembly
from study_s1_local_force_geometry import load_reference


def attachment_points(route):
    if "current_free_strand" in route:
        strand = route["current_free_strand"]
        return (strand["endpoint_A_source_mm"], strand["endpoint_B_source_mm"],
                "historical_neutral_free_terminal_chord")
    if route["region"] == "elbow":
        points = route["local_route"]["points_source_mm"]
        return points[0], points[-1], "reference_local_wire_endpoints"
    ports = route["endpoint_ports_reference"]
    return (ports["origin"]["outer_face_source_mm"], ports["distal"]["outer_face_source_mm"],
            "reference_outer_port_chord")


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1504_negative_relief_integration.json"
    case = read(path_file)
    candidates = case["cases"]["FJ1504_local_relief"]["negative_half_degree_endpoint_probes"]
    path = min((p["extension"] for p in candidates if "extension" in p
                and p["extension"]["outbound_and_inbound_clear"]), key=lambda row: row["last_angle_deg"])
    inventory_file = DEST / "route_inventory.json"
    inventory = read(inventory_file)
    relief = load_reference(case["source"])
    coupled = load_reference(relief["source"])
    passive_ref = coupled["pose_reference"]
    origins = load_reference(passive_ref)["neutral_origins_m"]
    model = ShoulderAssembly({}, {}, origins)
    summaries = []
    rows = []
    for route in inventory["routes"]:
        a, b, point_kind = attachment_points(route)
        a, b = np.asarray(a), np.asarray(b)
        origin_owner = route["anatomical_origin_body"]
        distal_owner = route["distal_body"]
        neutral = float(np.linalg.norm(b - a))
        neutral_direction = (b - a) / neutral
        phases = {}
        max_swing = 0.
        for phase in ("outbound", "inbound"):
            history = []
            for segment in path[phase]:
                for frame in segment:
                    if frame["hits"] or frame["limit"]:
                        raise ValueError("Path includes a blocked pose")
                    angle = frame["scapular_coronal_deg"]
                    sc = np.asarray(frame["SC_rotvec_deg"])
                    clavicle, scapula, _ = model.kinematics(angle, sc)
                    # GH and elbow are fixed, so their links share the same
                    # rigid world transform as the scapula on this trajectory.
                    transforms = {"clavicle": clavicle, "scapula": scapula,
                                  "humerus": scapula, "forearm": scapula}
                    first, second = transforms[origin_owner], transforms[distal_owner]
                    start = first[:, :3] @ a + first[:, 3]
                    end = second[:, :3] @ b + second[:, 3]
                    chord = float(np.linalg.norm(end - start))
                    direction = second[:, :3].T @ (end - start) / chord
                    swing = np.rad2deg(np.arctan2(np.linalg.norm(np.cross(neutral_direction, direction)),
                                                  neutral_direction @ direction))
                    max_swing = max(max_swing, float(swing))
                    history.append(chord - neutral)
                    rows.append(dict(muscle=route["id"], phase=phase, angle_deg=angle,
                                     chord_mm=chord, chord_change_mm=chord - neutral,
                                     nominal_chord_strain=(chord - neutral) / neutral,
                                     distal_frame_chord_swing_deg=swing,
                                     distal_frame_chord_x=direction[0], distal_frame_chord_y=direction[1],
                                     distal_frame_chord_z=direction[2]))
            phases[phase] = np.asarray(history)
        replay_delta = np.max(np.abs(phases["outbound"] - phases["inbound"][::-1]))
        excursion = float(np.ptp(phases["outbound"]))
        if (origin_owner != "clavicle" and excursion > 1e-8) or replay_delta > 1e-8:
            raise ValueError("Rigid-group invariance or reverse replay disagrees")
        summary = dict(muscle=route["id"], point_kind=point_kind, origin_body=origin_owner,
                       distal_body=distal_owner, endpoints_source_mm=[a.tolist(), b.tolist()],
                       neutral_chord_mm=neutral, minimum_delta_mm=float(min(phases["outbound"])),
                       maximum_delta_mm=float(max(phases["outbound"])),
                       endpoint_delta_mm=float(phases["outbound"][-1]), excursion_mm=excursion,
                       max_distal_frame_chord_swing_deg=max_swing,
                       replay_max_difference_mm=float(replay_delta))
        summaries.append(summary)
        print(route["id"], neutral, excursion, summary["endpoint_delta_mm"], flush=True)
    table = folder / "FJ1504_negative_path_local_stroke.csv"
    with table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save(folder / "FJ1504_negative_path_local_stroke.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_negative_path_local_stroke.py"),
        path=provenance(path_file), inventory=provenance(inventory_file), origins=passive_ref,
        table=provenance(table), endpoint_deg=path["last_angle_deg"],
        samples_each_way=sum(map(len, path["outbound"])), routes=summaries,
        scope="Prescribed motion and neutral port coordinates. Chords isolate local attachment "
              "differential travel. They are not routed wire payout, actual TPU strain or a force model. "
              "Free terminals are illustrative co-moving neutral references. Supply-side wire paths "
              "and contact/deformation are excluded.",
    ))


if __name__ == "__main__":
    main()
