"""Replay the revised FJ1468 guide and its actual TPU mesh under prescribed motion."""

import numpy as np

from audit_s1_infraspinatus_distributed_deformation import blend
from audit_elbow_slider_short_guides import export_lattice
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from study_elbow_series_package import bounded, overlaps
from study_s1_coupled_scapula import ShoulderAssembly
from study_s1_infraspinatus_wire_travel import route_metrics
from study_s1_local_force_geometry import load_reference


def shortcut_route(points, start, end, targets):
    """Shorten the near-port route with sampled collision-free chords."""
    chosen = [start]
    rejected = []
    while chosen[-1] < end:
        i = chosen[-1]
        for j in range(end, i, -1):
            swept = A.D.outer_envelope(points[[i, j]], .3)
            hit = overlaps(swept, targets)
            if not hit:
                chosen.append(j)
                break
            rejected.append(dict(from_index=i, to_index=j, hits_mm3=hit))
        else:
            return dict(valid=False, reached_index=i, rejected_chords=rejected)
    route = np.vstack([points[:start], points[chosen], points[end + 1:]])
    # A unioned sweep can expose small contacts spanning more than one chord.
    final_hits = overlaps(A.D.outer_envelope(route, .3), targets)
    length, turns, radii = route_metrics(route)
    # Include the junction angles at either end of the edited section.
    local_turn = float(turns[start - 1:start + len(chosen) - 1].sum())
    return dict(valid=not final_hits, accepted_indices=chosen, rejected_chords=rejected,
                full_route_overlap_mm3=final_hits, arc_mm=length,
                minimum_sampled_radius_mm=float(radii.min()), near_port_turn_rad=local_turn,
                centreline_source_mm=route.tolist(),
                ideal_equal_tension_reaction_N_at_10N=float(20. * np.sin(turns.max() / 2.)),
                hypothetical_full_sliding_capstan=[
                    dict(mu=mu, input_N_for_10N_exit=10. * float(np.exp(mu * local_turn)))
                    for mu in (0., .1, .2, .3)
                ])


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    audit_file = folder / "FJ1468_channel_wall_audit.json"
    audit = read(audit_file)
    transition = load_reference(audit["source"])
    selected = next(row for row in transition["recuts"] if "mesh" in row)
    prior = load_reference(transition["source"])
    stroke = load_reference(prior["path_stroke"])
    integration = load_reference(stroke["path"])
    probes = integration["cases"]["FJ1504_local_relief"]["negative_half_degree_endpoint_probes"]
    trajectory = min((p["extension"] for p in probes if p.get("extension", {}).get("outbound_and_inbound_clear")),
                     key=lambda row: row["last_angle_deg"])
    model = ShoulderAssembly({}, {}, load_reference(stroke["origins"])["neutral_origins_m"])
    points = np.asarray(selected["centreline_source_mm"])
    mf = A.I._manifold()
    proximal = load_reference(prior["proximal_mount"])["sources"]
    distal = load_reference(prior["distal_mount"])["sources"]
    original = load_component(proximal["tpu:deltoid"])
    revised = load_component(selected["mesh"])
    downstream = load_component(distal["tpu:deltoid"])
    added, removed = revised - original, original - revised
    tpu = ((downstream + added) - removed).as_original().simplify(1e-6)
    tpu, merged_mesh = export_lattice(tpu, folder / "FJ1468_transition_on_distal_mount_TPU.obj", points[0])
    inventory = load_reference(stroke["inventory"])
    original_point_count = next(len(route["local_route"]["points_source_mm"])
                                for route in inventory["routes"] if route["id"] == "FJ1468")
    neutral_wires = {}
    for route in inventory["routes"]:
        key = route["id"]
        if key in ("FJ1468", "FJ1467", "FJ1513"):
            wire_points = points if key == "FJ1468" else np.array(route["local_route"]["points_source_mm"])
            neutral_wires[key] = float((A.D.outer_envelope(wire_points, .3) ^ tpu).volume())
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(tpu.to_mesh64())
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    groups = [
        {key: load_component(proximal[key]) for key in ("bone:clavicle", "cap:FJ1468/origin",
                                                       "origin:FJ1468/front_support")},
        {key: load_component(distal[key]) for key in ("bone:humerus", "cap:FJ1468/distal",
                                                     "distal:FJ1468/front_support")},
    ]
    origin = points[0]
    axis = A.D.unit(points[-1] - origin)
    span = (points[-1] - origin) @ axis
    guards = []
    for group in groups:
        contact = []
        for key, part in group.items():
            if key.startswith("bone:"):
                continue
            padded = part.minkowski_sum(mf.Manifold.cube([1.] * 3, True))
            patch = padded ^ tpu
            if patch.volume() > .001:
                patch_vertices, _ = mesh_geom.manifold_mesh_to_arrays(patch.to_mesh64())
                station = (patch_vertices - origin) @ axis
                contact.append(dict(part=key, range_mm=[float(station.min()), float(station.max())]))
        guards.append(contact)
    lower = max([12.] + [row["range_mm"][1] for row in guards[0]])
    upper = min([span - 12.] + [row["range_mm"][0] for row in guards[1]])
    if upper - lower < 5.:
        raise ValueError("End guards leave too little space for the assumed blend")
    v_station = (vertices - origin) @ axis
    p_station = (points - origin) @ axis
    v_weight, derivative = blend(v_station, lower, upper, .1)
    p_weight, _ = blend(p_station, lower, upper, .1)
    neutral_length = route_metrics(points)[0]
    rows = []
    for segment in trajectory["outbound"]:
        for frame in segment:
            angle = frame["scapular_coronal_deg"]
            if abs(angle * 2 - round(angle * 2)) > 1e-9 or any(row["angle_deg"] == angle for row in rows):
                continue
            clavicle, humerus, _ = model.kinematics(angle, np.asarray(frame["SC_rotvec_deg"]))
            a, b = clavicle[:, :3], humerus[:, :3]
            ta, tb = clavicle[:, 3], humerus[:, 3]
            moved = vertices @ a.T + ta
            delta = vertices @ (b - a).T + (tb - ta)
            moved += v_weight[:, None] * delta
            jacobian = (a + v_weight[:, None, None] * (b - a)
                        + derivative[:, None, None] * delta[:, :, None] * axis)
            determinant = float(np.linalg.det(jacobian).min())
            shaped = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
            if shaped.status() != mf.Error.NoError or determinant <= 0.:
                raise ValueError(f"Invalid deformation at {angle} degrees")
            wire_path = points @ a.T + ta + p_weight[:, None] * (points @ (b - a).T + tb - ta)
            length, turns, radii = route_metrics(wire_path)
            obstacles = {key: part.transform(transform)
                         for group, transform in zip(groups, (clavicle, humerus))
                         for key, part in group.items()}
            rigid = bounded(obstacles)
            passage = {}
            for radius in (.3, .5):
                wire = A.D.outer_envelope(wire_path, radius)
                passage[f"{radius:.1f}"] = dict(rigid_mm3=overlaps(wire, rigid),
                                                  TPU_mm3=float((wire ^ shaped).volume()))
            row = dict(angle_deg=angle, arc_change_mm=length - neutral_length,
                       min_radius_mm=float(radii.min()), peak_turn_deg=float(np.rad2deg(turns.max())),
                       minimum_deformation_determinant=determinant, passage=passage,
                       TPU_rigid_overlap_mm3=overlaps(shaped, rigid))
            if angle in (0., -2., -4.):
                local_end = len(points) - original_point_count + selected["join_index"]
                cable_targets = bounded({**obstacles, "tpu:deltoid": shaped})
                result = shortcut_route(wire_path, 3, local_end, cable_targets)
                if "arc_mm" in result:
                    result["length_change_from_same_pose_centreline_mm"] = result["arc_mm"] - length
                row["near_port_shortcut_screen"] = result
                print("shortcut", angle, result.get("valid"), result.get("length_change_from_same_pose_centreline_mm"),
                      len(result["rejected_chords"]), flush=True)
            rows.append(row)
            print(angle, row["arc_change_mm"], passage, row["TPU_rigid_overlap_mm3"], flush=True)
    save(folder / "FJ1468_revised_channel_motion.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_deltoid_revised_channel_motion.py"),
        wall_audit=provenance(audit_file), protection_intersections=guards,
        merged_TPU=merged_mesh, neutral_three_wire_overlap_mm3=neutral_wires,
        deformation_window_mm=[lower, upper], neutral_arc_mm=neutral_length, cases=rows,
        scope="Prescribed half-degree waypoints from the outbound SC/AC path. Same axial blend for "
              "the lattice and wire. Return has the same kinematics in this model. The proximal "
              "Boolean change is applied to the distal-mount TPU candidate. The mesh and wire "
              "are not contact- or load-solved. Selected waypoints include an independent "
              "geometric shortcut screen and illustrative capstan sensitivity. Shortcut vertices "
              "are numerical samples, not physical guide points or contact equilibrium.",
    ))


if __name__ == "__main__":
    main()
