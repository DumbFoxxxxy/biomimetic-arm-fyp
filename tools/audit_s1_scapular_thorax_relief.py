"""Independently check the low-profile mounts in the shoulder-elbow assembly."""

from itertools import combinations
from pathlib import Path

import numpy as np

from audit_s1_series_elbow_integration import MOBILE, assert_source, find_hits
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from study_elbow_series_package import bounded, overlaps


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    revision_file = folder / "scapular_thorax_relief.json"
    revision = read(revision_file)
    parent_path = assert_source(revision["baseline"])
    parent = read(parent_path)
    replacements = revision["replacements"]
    removed = set(revision["replaced_nuts_and_washers"])
    if not all(row.get("selected") and row.get("exported_hits_mm3") == {} for row in revision["reports"]):
        raise RuntimeError("Both scapular mounts need verified exports for combined assessment")
    sources = {**{key: row for key, row in parent["sources"].items() if key not in removed}, **replacements}
    owners = {**{key: owner for key, owner in parent["owners"].items() if key in sources},
              **{key: "scapula" for key in replacements}}
    solids = {key: load_component(row) for key, row in sources.items()}
    keys = list(solids)
    neutral_pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                                if owners[keys[i]] != owners[keys[j]]
                                or keys[i] in replacements or keys[j] in replacements])
    neutral = find_hits(list(solids.values()), keys, neutral_pairs)
    new_pair_hits = [row for row in neutral if row["left"] in replacements or row["right"] in replacements]
    passive = read(DEST / "elbow_passive_model/model_check.json")
    centres = passive["neutral_origins_m"]
    chains = {"SC": (np.asarray(centres["clavicle"]) * 1000., MOBILE | {"clavicle", "scapula"}),
              "AC": (np.asarray(centres["scapula"]) * 1000., MOBILE | {"scapula"})}
    mount_geometry = {key: solids[key] for key in replacements}
    motions = []
    for label, (pivot, moving) in chains.items():
        fixed = {key: solid for key, solid in solids.items() if owners[key] not in moving}
        fixed_bounds = bounded(fixed)
        cross_pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                                   if (owners[keys[i]] in moving) != (owners[keys[j]] in moving)])
        for axis in np.eye(3):
            for sign in (-1., 1.):
                first_whole, first_mount = None, None
                last_whole, last_mount = 0., 0.

                def evaluate(angle, only_mount):
                    matrix = A.I.rotation_about(pivot, axis, sign * angle)
                    if only_mount:
                        return {key: hit for key, solid in mount_geometry.items()
                                if (hit := overlaps(solid.transform(matrix), fixed_bounds))}
                    placed = [solid.transform(matrix) if owners[key] in moving else solid
                              for key, solid in solids.items()]
                    return find_hits(placed, keys, cross_pairs)

                for angle in range(1, 11):
                    if first_whole is None:
                        hit = evaluate(angle, False)
                        if hit:
                            low, high = last_whole, float(angle)
                            while high - low > .1:
                                mid = (low + high) / 2.
                                smaller = evaluate(mid, False)
                                if smaller:
                                    high, hit = mid, smaller
                                else:
                                    low = mid
                            last_whole = low
                            first_whole = dict(angle_deg=sign * high, hits=hit)
                        else:
                            last_whole = float(angle)
                    if first_mount is None:
                        hit = evaluate(angle, True)
                        if hit:
                            low, high = last_mount, float(angle)
                            while high - low > .1:
                                mid = (low + high) / 2.
                                smaller = evaluate(mid, True)
                                if smaller:
                                    high, hit = mid, smaller
                                else:
                                    low = mid
                            last_mount = low
                            first_mount = dict(angle_deg=sign * high, hits=hit)
                        else:
                            last_mount = float(angle)
                    if first_whole and first_mount:
                        break
                motions.append(dict(joint=label, axis_source=axis.tolist(), direction=sign,
                                    global_last_sample_clear_deg=sign * last_whole, global_first_hit=first_whole,
                                    mount_last_sample_clear_deg=sign * last_mount, mount_first_hit=first_mount))
    routes = {row["head"]: row for row in read(assert_source(revision["route_reference"]))["routes"]}
    thorax = bounded({key: solid for key, solid in solids.items() if key.startswith("thorax_reference:")})
    inlet = []
    original_tpu_mount_overlaps = {}
    for row in revision["reports"]:
        head = row["head"]
        port = routes[head]["ports"]["origin"]
        normal = np.asarray(port["normal_source"])
        shift = row["selected"].get("shift_mm", 0.)
        old = np.asarray(port["outer_face_source_mm"])
        moved = old + shift * normal
        first = None
        for distance in np.arange(.1, 10.01, .1):
            hits = overlaps(rod(moved, moved + distance * normal, .3), thorax)
            if hits:
                first = dict(probe_length_mm=float(distance), hits_mm3=hits)
                break
        inlet.append(dict(head=head, previous_outer_face_mm=old.tolist(), revised_outer_face_mm=moved.tolist(),
                          normal_source=normal.tolist(), first_thorax_overlap=first))
        muscle_key = "tpu:" + routes[head]["module"]
        tpu = load_component(parent["flexible_geometry_sources"][muscle_key])
        mount_and_thorax = {**thorax, **bounded({
            key: solid for key, solid in solids.items()
            if key.startswith(f"origin:{head}/") or key == f"cap:{head}/origin"
        })}
        original_tpu_mount_overlaps[muscle_key] = overlaps(tpu, mount_and_thorax)
    target = folder / "assembly_scapular_relief.json"
    save(target, dict(**{key: value for key, value in parent.items() if key not in ("sources", "owners", "scope")},
                      sources=sources, owners=owners, scapular_mount_revision=provenance(revision_file),
                      scope="Combined rigid neutral candidate with low-profile scapular mounts. Original TPU "
                      "is retained as a geometry reference pending adaptation to revised FJ1504 cap."))
    save(folder / "scapular_thorax_relief_audit.json", dict(
        generator=provenance(Path(__file__).resolve()), assembly=provenance(target), previous=provenance(parent_path),
        rigid_count=len(solids), neutral_pairs_checked=len(neutral_pairs), neutral_overlaps=neutral,
        neutral_revised_mount_overlaps=new_pair_hits, shoulder_pose_probes=motions,
        local_inlet_axial_probe=inlet, original_tpu_neutral_mount_and_thorax_overlaps_mm3=original_tpu_mount_overlaps,
        sampling_deg=1., boundary_step_deg=.1, overlap_threshold_mm3=.001,
        scope="Neutral cross-body and revised same-body material checked. Each SC/AC probe keeps the "
        "other joints neutral; simulated pose angles do not validate continuous or powered motion. "
        "Inlet probes are straight physical clearance measurements from the cap face against thorax only. "
        "Actual redirected wire, flexible pad, fastener retention and installation remain pending."))
    print("neutral", len(neutral), "with revised mounts", len(new_pair_hits), "inlets", inlet,
          "original TPU", original_tpu_mount_overlaps, flush=True)
    for row in motions:
        print(row, flush=True)


if __name__ == "__main__":
    main()
