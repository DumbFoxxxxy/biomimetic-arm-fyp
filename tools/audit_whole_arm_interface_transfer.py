"""Audit transfer of the existing shoulder and elbow packages to the shared chain."""

from itertools import combinations

import numpy as np

from audit_s1_ac_blind_coupled_motion import load_assembly
from audit_s1_series_elbow_integration import find_hits
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save


def main():
    shared_folder = DEST / "unified_SC_AC"
    assembly_file = shared_folder / "assembly.json"
    check_file = shared_folder / "geometry_mass_body_chain_check.json"
    shared = read(assembly_file)
    check = read(check_file)
    if check["reference"] != provenance(assembly_file):
        raise ValueError("The motion chain uses a different shared assembly")
    output = DEST.parent / "whole_arm_first_pass"
    package, _, _ = load_assembly(DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief")
    register_file = DEST / "s1_cable_anchor_register.json"
    register = read(register_file)
    shared_sources = {**shared["sources"], **shared["joint_context_sources"]}
    old_sources = register["neutral_rigid_sources"]
    shared_owners = {row["source_key"]: row["owner"] for row in check["geom_owners"]}

    # Use the shared skeleton and its matching joint hardware as a set.
    # Test only the muscle mounts and actuation hardware added to that set.
    added = {}
    for key, source in old_sources.items():
        if (key in shared_sources or key.startswith(("SC:", "AC:", "bone:", "thorax_reference:"))):
            continue
        added[key] = source
    seats_file = output / "girdle_seat_pair.json"
    seats = read(seats_file)
    if seats["assembly"] != provenance(assembly_file):
        raise ValueError("Reaction seat prototype uses a different shared assembly")
    owners = dict(package["owners"])
    owners.update(shared_owners)
    owners.update(dict.fromkeys(("fastener_0", "fastener_1", "support"), "scapula"))
    for name, info in seats["seats"].items():
        added[name] = info["source"]
        owners[name] = info["owner"]
    # Ownership is diagnostic at neutral; any unassigned parts require explicit follow-up.
    sources = {**shared_sources, **added}
    keys = list(sources)
    solids = [load_component(source) for source in sources.values()]
    pairs = np.array([(a, b) for a, b in combinations(range(len(keys)), 2)
                      if keys[a] in added or keys[b] in added])
    hits = find_hits(solids, keys, pairs)
    for hit in hits:
        hit["left_owner"] = owners.get(hit["left"])
        hit["right_owner"] = owners.get(hit["right"])
    shared_solids = dict(zip(keys, solids))
    original_humerus = shared_solids["bone:humerus"]
    donor_humerus = load_component(old_sources["bone:humerus"])
    # Exchange only the forearm-side stock. Confirm that both versions share
    # material in the splice band before retaining the distal interface.
    cube = A.I._manifold().Manifold.cube
    for height in (1180., 1170., 1190., 1160., 1200.):
        band = cube([600., 600., 2.]).translate([-400., -350., height - 1.])
        stock_a, stock_b = original_humerus ^ band, donor_humerus ^ band
        difference = float((stock_a - stock_b).volume() + (stock_b - stock_a).volume())
        if difference < .001 and (stock_a ^ stock_b).volume() > 10.:
            break
    else:
        height = None
    revisions = {}
    if height is not None:
        upper = cube([600., 600., 600.]).translate([-400., -350., height - .5])
        lower = cube([600., 600., 600.]).translate([-400., -350., height - 599.5])
        humerus = ((original_humerus ^ upper) + (donor_humerus ^ lower)).as_original().simplify(1e-6)
        mount_file = DEST / "elbow_mount_candidate/assembly.json"
        # The four muscle bracket bores lie above the elbow-stop splice.
        for seat in read(mount_file)["mount_seats"]:
            front = np.asarray(seat["front_seat_source_mm"])
            back = np.asarray(seat["back_seat_source_mm"])
            normal = np.asarray(seat["axis_source"])
            humerus -= rod(back - 20. * normal, front + 20. * normal, seat["bore_radius_mm"])
        rebuilt, source = export(humerus.as_original().simplify(1e-6), output / "humerus_elbow_transfer.obj",
                                 np.asarray(shared_sources["bone:humerus"]["translation_source_mm"]))
        revisions["bone:humerus"] = source
        shared_solids["bone:humerus"] = rebuilt
    stock_file = DEST.parent / "s1_local_muscle_load/FJ1467_current_scapula_blind_axes.json"
    stock = read(stock_file)
    normal = np.asarray(stock["normal_source"])
    mounts = stock["compact_pair_nearest_centre"]
    scapula = shared_solids["bone:scapula"]
    for row in (mounts["first"], mounts["second"]):
        axis = np.asarray(row["axis_source_mm"])
        entry = np.asarray(row["bone_entry_mm"])
        scapula -= rod(entry - 3. * normal, axis + 2.1 * normal, 1.15, 48)
    rebuilt, source = export(scapula.as_original().simplify(1e-7),
                             output / "scapula_posterior_bore_transfer.obj",
                             np.asarray(shared_sources["bone:scapula"]["translation_source_mm"]))
    revisions["bone:scapula"] = source
    shared_solids["bone:scapula"] = rebuilt
    revision_pairs = np.array([(a, b) for a, b in combinations(range(len(keys)), 2)
                                if keys[a] in added or keys[b] in added
                                or keys[a] in revisions or keys[b] in revisions])
    revised_hits = find_hits(list(shared_solids.values()), keys, revision_pairs)
    save(output / "interface_transfer_audit.json", dict(
        shared_assembly=provenance(assembly_file), shared_chain=provenance(check_file),
        s1_elbow_register=provenance(register_file), girdle_seats=provenance(seats_file),
        shared_sources=shared_sources, added_sources=added,
        owner_by_source={key: owners.get(key) for key in keys},
        unassigned=[key for key in added if key not in owners],
        changed_common_key_geometry=[key for key, source in shared_sources.items()
                                     if key in old_sources and old_sources[key] != source],
        reference_channels=register["channels"], checked_added_pairs=len(pairs),
        overlaps_mm3=hits, revisions=revisions, revised_overlaps_mm3=revised_hits,
        transferred_humerus_joint_source=old_sources["bone:humerus"],
        humerus_splice_z_mm=height, humerus_splice_band_difference_mm3=difference,
        humerus_mount_bore_reference=provenance(DEST / "elbow_mount_candidate/assembly.json"),
        scapular_bore_reference=provenance(stock_file), revised_checked_pairs=len(revision_pairs),
        generator=provenance(ROOT / "tools/audit_whole_arm_interface_transfer.py"),
        scope="Neutral rigid fit transfer, including same-body overlaps for diagnosis. Existing "
              "SC/AC/GH and elbow joint context uses the shared assembly. Potential intended "
              "contacts remain in the report. No flexibles, wire dynamics or package motion verified."))
    print(f"{len(shared_sources)} shared, {len(added)} added; {len(hits)} overlaps across {len(pairs)} checked pairs")
    for hit in hits[:10]:
        print(hit)
    print(f"After local transfer: {len(revised_hits)} overlaps; humerus splice at source Z {height}")
    for hit in revised_hits[:10]:
        print(hit)


if __name__ == "__main__":
    main()
