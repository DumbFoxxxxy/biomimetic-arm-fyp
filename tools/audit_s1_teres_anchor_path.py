"""Compare shoulder port-to-port straight pull with the inventory wire route."""

import argparse
import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_ac_blind_coupled_motion import load_assembly
from audit_s1_infraspinatus_distributed_deformation import blend
from audit_s1_series_elbow_integration import MOBILE, compose
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_coupled_scapula import ShoulderAssembly, affine_at
from study_s1_local_force_geometry import load_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channel", choices=("FJ1507", "FJ1513", "FJ1468"), default="FJ1507")
    args = parser.parse_args()
    local = DEST.parent / "s1_local_muscle_load"
    register_ref = provenance(DEST / "s1_cable_anchor_register.json")
    register = load_reference(register_ref)
    inventory = load_reference(register["inventory"])
    item = next(row for row in inventory["routes"] if row["id"] == args.channel)
    origin_body = "clavicle" if args.channel == "FJ1468" else "scapula"
    if (item["anatomical_origin_body"], item["distal_body"]) != (origin_body, "humerus"):
        raise ValueError("Unexpected route owners")
    rigid_ref = provenance(local / "FJ1468_girdle_flexion_clearance.json")
    rigid = load_reference(rigid_ref)
    girdle = load_reference(rigid["girdle_path"])
    endpoints = load_reference(load_reference(rigid["span_study"])["local_endpoints"])
    model = ShoulderAssembly({}, {}, load_reference(girdle["pose_reference"])["neutral_origins_m"])
    assembly, _, original = load_assembly(local / "FJ1500_cap_relief")
    sources = register["neutral_rigid_sources"]
    owners = {key: assembly["owners"].get(key, "scapula") for key in sources}
    if "cap:FJ1507/origin" not in original or sources["cap:FJ1507/origin"] != original["cap:FJ1507/origin"]:
        raise ValueError("Anchor register uses another FJ1507 origin cap")
    if args.channel == "FJ1513":
        # Compare like for like with the prior neighbor screen.
        neighbor_ref = provenance(local / "FJ1468_girdle_flexion_channel_neighbor_wires.json")
        neighbor = load_reference(neighbor_ref)
        guide = load_reference(neighbor["guide_motion"])
        front = np.asarray(neighbor["neutral_path_source_mm"])
        axis = front[-1] - front[0]
        axis /= np.linalg.norm(axis)
    solids = {key: load_component(metadata) for key, metadata in sources.items()}
    wire = np.asarray(item["local_route"]["points_source_mm"])
    motion = dict(guided_weight_end_zones_mm=6.)
    if args.channel == "FJ1468":
        front_ref = provenance(local / "FJ1468_girdle_flexion_channel_liner_assembly.json")
        front_report = load_reference(front_ref)
        if front_report["rigid_audit"] != rigid_ref or front_report["adapter_mesh"] != sources["elbow:FJ1486_bone_adapter"]:
            raise ValueError("Front route does not match the anchor assembly")
        guide = load_reference(front_report["guide_motion"])
        wire = np.asarray(front_report["neutral_path_source_mm"])
        axis = wire[-1] - wire[0]
        axis /= np.linalg.norm(axis)
        motion = dict(revised_route=front_ref, blend_window_mm=guide["deformation_window_mm"])
    ports = item["endpoint_ports_reference"]
    origin, distal = (np.asarray(ports[name]["outer_face_source_mm"]) for name in ("origin", "distal"))
    path_axis = distal - origin
    length = np.linalg.norm(path_axis)
    path_axis /= length
    origin_axis = -np.asarray(ports["origin"]["normal_source"])
    lead_angle = float(np.rad2deg(np.arccos(np.clip(path_axis @ origin_axis, -1., 1.))))
    # Keep the cap regions with their respective bones. The intermediate
    # guide positions are a geometric hypothesis between the two end zones.
    weights, _ = blend((wire - origin) @ path_axis, 6., length - 6., .1)
    if args.channel == "FJ1513":
        weights, _ = blend((wire - front[0]) @ axis, *guide["deformation_window_mm"], .1)
        motion = dict(neighbor_replay=neighbor_ref, blend_window_mm=guide["deformation_window_mm"])
    elif args.channel == "FJ1468":
        weights, _ = blend((wire - wire[0]) @ axis, *guide["deformation_window_mm"], .1)
    neutral_stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(wire, axis=0), axis=1))]
    flex_axis = np.asarray(endpoints["axes_source"][0])
    rows = []
    poses = sorted({tuple(pose) for case in rigid["cases"].values() for pose in case["clear_path_deg"]})
    if not poses or any(case["reverse_replay_hits"] for case in rigid["cases"].values()):
        raise ValueError("Expected checked shoulder paths")
    for pose in poses:
        angle, *sc, flex = pose
        clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
        gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * flex_axis).as_matrix(), model.gh)
        humerus = compose(scapula, gh)
        frames = {"clavicle": clavicle, "scapula": scapula, **dict.fromkeys(MOBILE, humerus)}
        targets = bounded({key: solid.transform(frames[owners[key]]) if owners[key] in frames else solid
                           for key, solid in solids.items()})
        start = frames[origin_body]
        a, ta = start[:, :3], start[:, 3]
        b, tb = humerus[:, :3], humerus[:, 3]
        chord = np.array([a @ origin + ta, b @ distal + tb])
        guided = wire @ a.T + ta + weights[:, None] * (wire @ (b - a).T + tb - ta)
        segments = np.diff(guided, axis=0)
        lengths = np.linalg.norm(segments, axis=1)
        if np.any(lengths < 1e-8):
            raise ValueError(f"Collapsed cable segment at {pose}")
        tangents = segments / lengths[:, None]
        reaction = np.diff(tangents, axis=0)
        strongest = int(np.argmax(np.linalg.norm(reaction, axis=1)))
        node = strongest + 1
        rows.append(dict(pose_deg=pose, chord_mm=float(np.linalg.norm(chord[1] - chord[0])),
                         guided_arc_mm=float(lengths.sum()),
                         peak_turn_node=node, peak_turn_neutral_station_mm=float(neutral_stations[node]),
                         peak_turn_neutral_xyz_mm=wire[node].tolist(),
                         peak_turn_reaction_per_N=np.linalg.norm(reaction[strongest]).item(),
                         peak_turn_distal_blend_weight=float(weights[node]),
                         chord_hits_mm3=overlaps(A.D.outer_envelope(chord, .30), targets),
                         guided_hits_mm3=overlaps(A.D.outer_envelope(guided, .30), targets)))
    if args.channel == "FJ1468":
        recorded = {tuple(row["pose_deg"]): row for case in front_report["cases"].values() for row in case}
        if set(recorded) != set(poses):
            raise ValueError("Front route study uses a different pose grid")
        np.testing.assert_allclose([row["guided_arc_mm"] for row in rows],
                                   [recorded[pose]["arc_mm"] for pose in poses], atol=1e-8, rtol=0.)
    save(local / f"{args.channel}_anchor_path_audit.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_teres_anchor_path.py"),
        anchor_register=register_ref, rigid_path=rigid_ref, chord_to_inward_port_axis_deg=lead_angle,
        channel=args.channel, cases=rows, **motion,
        scope="1:1 rigid meshes with no excluded pierce pairs, neutral elbow, sampled shoulder paths. "
              "Port-to-port chord vs nominal cable path with axial guide interpolation, "
              "including its incoming lead. The change in adjacent unit tangents is a discrete "
              "support-load indicator under equal tension, not a contact result. No TPU mesh "
              "or tension equilibrium evaluated.",
    ))
    print(f"Checked {len(rows)} shoulder poses; neutral port approach angle {lead_angle:.2f} degrees", flush=True)


if __name__ == "__main__":
    main()
