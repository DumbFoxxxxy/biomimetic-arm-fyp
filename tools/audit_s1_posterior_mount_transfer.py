"""Check transfer of the existing FJ1467 side-nut mount into the current shoulder assembly."""

from itertools import combinations

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_ac_blind_coupled_motion import load_assembly
from audit_s1_infraspinatus_distributed_deformation import blend
from audit_s1_series_elbow_integration import MOBILE, compose, find_hits
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_coupled_scapula import ShoulderAssembly, affine_at
from study_s1_local_force_geometry import load_reference


def main():
    local = DEST.parent / "s1_local_muscle_load"
    folder = local / "FJ1500_cap_relief"
    mount_ref = provenance(A.DESIGN / "s1_middle_side_nut_mount/assembly_check.json")
    mount = load_reference(mount_ref)
    parent = load_reference(dict(file=mount["parent_report"], sha256=mount["parent_sha256"]))
    locations_ref = provenance(local / "FJ1467_posterior_rigid_locations.json")
    locations = load_reference(locations_ref)
    neighbor = load_reference(locations["neighbor_replay"])
    rigid = load_reference(neighbor["rigid_audit"])
    path = load_reference(rigid["girdle_path"])
    endpoints = load_reference(load_reference(rigid["span_study"])["local_endpoints"])
    guide = load_reference(neighbor["guide_motion"])
    assembly, _, sources = load_assembly(folder)
    sources["elbow:FJ1486_bone_adapter"] = neighbor["adapter_mesh"]
    owners = {key: assembly["owners"].get(key, "scapula") for key in sources}
    base = {key: load_component(row) for key, row in sources.items()}
    scapula_key = "bone:scapula"
    source_before, source_after = (load_component(report["sources"][scapula_key]) for report in (parent, mount))
    parts = {key: load_component(row) for key, row in mount["sources"].items()
             if key.startswith(("origin:FJ1467/", "cap:FJ1467/"))}
    if not parts or parts.keys() & base.keys():
        raise ValueError("Expected absent posterior support hardware")
    # Transfer the removed stock only, retaining the current S1 engineering
    # scapula and its other interface revisions.
    pocket_stock = source_before - source_after
    current = base[scapula_key]
    installed = (current - pocket_stock).as_original().simplify(1e-7)
    donor_bounds = bounded({scapula_key: source_after})
    recipient_bounds = bounded({scapula_key: installed})
    fit = {key: dict(donor_scapula_hits_mm3=overlaps(part, donor_bounds),
                     recipient_scapula_hits_mm3=overlaps(part, recipient_bounds))
           for key, part in parts.items()}
    owners.update(dict.fromkeys(parts, "scapula"))
    migrated = {**base, **parts, scapula_key: installed}
    keys = list(migrated)
    pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                      if keys[i] in parts or keys[j] in parts])
    model = ShoulderAssembly({}, {}, load_reference(path["pose_reference"])["neutral_origins_m"])
    axis = np.asarray(endpoints["axes_source"][0])
    cable = np.asarray(locations["neutral_source_route_mm"])
    front = np.asarray(neighbor["neutral_path_source_mm"])
    axial = front[-1] - front[0]
    axial /= np.linalg.norm(axial)
    weights, _ = blend((cable - front[0]) @ axial, *guide["deformation_window_mm"], .1)
    rows = []
    for sample in locations["cases"]:
        angle, *sc, flex = sample["pose_deg"]
        clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
        gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
        humerus = compose(scapula, gh)
        frames = {"clavicle": clavicle, "scapula": scapula, **dict.fromkeys(MOBILE, humerus)}
        placed = {key: part.transform(frames[owners[key]]) if owners[key] in frames else part
                  for key, part in migrated.items()}
        wire = cable @ scapula[:, :3].T + scapula[:, 3]
        distal = cable @ humerus[:, :3].T + humerus[:, 3]
        envelope = A.D.outer_envelope(wire + weights[:, None] * (distal - wire), .30)
        hits = find_hits(list(placed.values()), keys, pairs)
        wire_hits = overlaps(envelope, bounded(placed))
        rows.append(dict(pose_deg=sample["pose_deg"], new_mount_hits=hits, wire_hits_mm3=wire_hits))
    out = local / "FJ1467_posterior_mount_transfer.json"
    save(out, dict(generator=provenance(ROOT / "tools/audit_s1_posterior_mount_transfer.py"),
                   historical_mount=mount_ref, posterior_locations=locations_ref,
                   scapula_removed_mm3=float((current - installed).volume()),
                   scapula_component_volumes_mm3=sorted((float(part.volume()) for part in installed.decompose()),
                                                      reverse=True),
                   transferred_part_keys=list(parts), donor_to_recipient_fit=fit,
                   checked_mount_pairs=len(pairs), cases=rows,
                   scope="Transfer screen in the current rigid short-bolt assembly. Transfers the historic "
                         "bone pocket cut and matching rigid parts without changing nominal cable. "
                         "TPU, assembly access, screw preload and strength remain outside this check."))
    print(f"Checked {len(rows)} poses and {len(pairs)} mount pairs", flush=True)


if __name__ == "__main__":
    main()
