"""Record local cable anchors and screen the straight-pull option at neutral."""

import numpy as np

from audit_s1_ac_blind_coupled_motion import load_assembly
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    local = DEST.parent / "s1_local_muscle_load"
    inventory_ref = provenance(DEST / "route_inventory.json")
    inventory = load_reference(inventory_ref)
    posterior_ref = provenance(local / "FJ1467_distal_detour_route.json")
    posterior = load_reference(posterior_ref)
    mount_ref = posterior["mount"]
    mount = load_reference(mount_ref)
    proximal_ref = posterior["proximal_exit"]
    proximal = load_reference(proximal_ref)
    liner_ref = provenance(local / "FJ1468_girdle_flexion_channel_liner_assembly.json")
    liner = load_reference(liner_ref)
    combined_ref = provenance(DEST.parent / "integration_review/s1_series_elbow/nominal_spring_spool_seating_k5_coupling.json")
    combined = load_reference(combined_ref)
    ferrule_ref = provenance(DEST / "elbow_supply_interface/route_endpoints.json")
    ferrules = {row["id"]: row for row in load_reference(ferrule_ref)["routes"]}
    folder = local / "FJ1500_cap_relief"
    _, _, sources = load_assembly(folder)
    sources = {**sources, "elbow:FJ1486_bone_adapter": liner["adapter_mesh"], **mount["meshes"],
               "support": proximal["matching_port"]["mesh"]}
    targets = bounded({key: load_component(row) for key, row in sources.items()})

    channels = []
    for item in inventory["routes"]:
        name = item["id"]
        row = dict(id=name, region=item["region"], muscle=item["module"],
                   entry_body=item["anatomical_origin_body"], terminal_body=item["distal_body"],
                   motor_anchor=None, upstream_reaction_anchor=None,
                   upstream_supply_status="Deferred to arm integration; tabletop thorax is the ground.")
        if "current_free_strand" in item:
            strand = item["current_free_strand"]
            points = [strand["endpoint_A_source_mm"], strand["endpoint_B_source_mm"]]
            row.update(entry_body="floating TPU terminal", route_reference=strand["source"],
                       local_status="Conditional free span only. Proximal reaction seat to be resolved.",
                       local_reaction_anchor=None)
        elif item["region"] == "shoulder":
            ports = item["endpoint_ports_reference"]
            points = [ports[key]["outer_face_source_mm"] for key in ("origin", "distal")]
            row.update(route_reference=item["local_route_reference"], local_reaction_anchor=None,
                       local_status="Nominal port pair. Confirm clamp, stop and direction on the selected assembly.")
            if name == "FJ1467":
                row.update(route_reference=posterior_ref,
                           local_status="Candidate stepped support and nominal guided route. "
                                        "Localized turns require support on the scapular and humeral sides.")
                row["candidate_guide_nodes_source_mm"] = [
                    posterior["neutral_candidate_mm"][i] for i in posterior["edited_nodes"]]
            if name == "FJ1468":
                row.update(route_reference=liner_ref,
                           local_status="Candidate local route with rigid clearance on the saved shoulder paths. "
                                        "Verify reaction seat and assembled shoulder-and-elbow motion.")
        else:
            sites = combined["channels"][name]["equivalent_sites"]
            points = [site["reference_xyz_mm"] for site in sites]
            row.update(entry_body=sites[0]["body"], terminal_body=sites[1]["body"],
                       route_reference=combined_ref,
                       local_status="Equivalent sites for the series-rail replay; confirm actual feed and crimp.",
                       local_reaction_anchor=None,
                       historical_ferrule_candidate=dict(reference=ferrule_ref,
                           body=ferrules[name]["owner"],
                           seat_xyz_mm=ferrules[name]["physical_sheath_reaction_seat_source_mm"]))
        row.update(entry_xyz_mm=points[0], terminal_xyz_mm=points[1])
        row["neutral_chord_mm"] = float(np.linalg.norm(np.subtract(*points)))
        if item["region"] == "shoulder" and "current_free_strand" not in item:
            row["neutral_straight_pull_rigid_hits_mm3"] = overlaps(
                A.D.outer_envelope(np.asarray(points), .30), targets)
        channels.append(row)
    save(DEST / "s1_cable_anchor_register.json", dict(
        generator=provenance(ROOT / "tools/record_s1_cable_anchors.py"), inventory=inventory_ref,
        coordinate_frame="BP3D source mm; model reference pose", channels=channels,
        neutral_rigid_sources=sources,
        scope="Anchor candidates and diagnostic chords. Only four fixed-end shoulder chords are "
              "screened against the current neutral rigid assembly, with no excluded pair. "
              "Port numbers do not establish reaction-seat fit, dynamic clearance or tendon preload.",
    ))
    print(f"Recorded {len(channels)} local cable channels")


if __name__ == "__main__":
    main()
