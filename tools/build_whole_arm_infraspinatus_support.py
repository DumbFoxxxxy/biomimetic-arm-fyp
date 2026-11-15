"""Join the raised infraspinatus bearing to the existing elbow adapter."""

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import transforms
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from project_runtime import load_mujoco_model
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    audit_ref = provenance(folder / "whole_arm_offset_infraspinatus_audit.json")
    audit = load_reference(audit_ref)
    motion = load_reference(audit["motion"])
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    selected = load_reference(audit["humeral_revision"])["selected"]
    kinematics = load_reference(load_reference(audit["original_kinematics"])["kinematics"])
    spec = next(row for row in kinematics["models"] if row["head"] == "FJ1500")
    pin = np.array(selected["parameters"]["new_pin_source_mm"])
    axis = np.array(spec["bone_pin_axis_source"])
    bearing = rod(pin - 1.5 * axis, pin + 1.5 * axis, 3.2) - rod(pin - 5. * axis, pin + 5. * axis, 1.1)
    adapter_key = "elbow:FJ1477_bone_adapter"
    jaw_key = "collar:jaw_upper"
    fixed = {key: meta for key, meta in selected["meshes"].items()
             if not key.endswith(("bone_groove_ring", "bone_uncrimped_sleeve"))}
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"], **fixed}
    solids = {key: load_component(meta) for key, meta in sources.items()}
    adapter, jaw = solids[adapter_key], solids.pop(jaw_key)
    mf = A.I._manifold()
    face_pad = mf.Manifold.cube([.4] * 3, True)
    # Remove only jaw material at the shared adapter fasteners. Retain the
    # original adapter and its fastener interfaces, and fuse the two solids.
    fastener_keys = [key for key in audit["fixed_revision_overlap_mm3"][jaw_key] if key != adapter_key]
    for key in fastener_keys:
        jaw -= solids[key].minkowski_sum(face_pad)
    interface_volume = float((jaw ^ adapter).volume())
    candidate = (adapter + jaw).as_original().simplify(1e-6)
    bearing_loss = max(0., float((bearing - candidate).volume()))
    if interface_volume <= 1. or bearing_loss > .001:
        raise ValueError("Shared support has insufficient joint material or a damaged bearing")
    actual, mesh = export(candidate, folder / "FJ1500_FJ1477_shared_support.obj", pin)
    added = actual - adapter
    baseline = {key: part for key, part in solids.items() if key != adapter_key}
    neutral = overlaps(actual, bounded(baseline))
    if neutral:
        print("Neutral support overlaps:", neutral, flush=True)
    owners = transfer["owner_by_source"]
    origins = {key: np.asarray(scene["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    xml = ROOT / scene["model"]["file"]
    if provenance(xml) != scene["model"]:
        raise ValueError("Shared model was modified after motion screening")
    model = load_mujoco_model(xml)
    data = mujoco.MjData(model)
    checked = []
    for row in motion["paths"]["K7_SC1_pos"]["accepted"]:
        data.qpos[:] = row["qpos"]
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        neighbors = bounded({k: part.transform(frames[owners[k]]) for k, part in baseline.items()})
        checked.append(dict(input_deg=row["input_deg"],
                            incremental_hits_mm3=overlaps(added.transform(frames["humerus"]), neighbors)))
    save(folder / "whole_arm_infraspinatus_support.json", dict(
        audit=audit_ref, generator=provenance(ROOT / "tools/build_whole_arm_infraspinatus_support.py"),
        source=mesh, replaced_adapter=sources[adapter_key], retained_pin_source_mm=pin.tolist(),
        bearing_material_loss_mm3=bearing_loss, joint_intersection_mm3=interface_volume,
        local_fastener_clearance_mm=.2, cut_against=fastener_keys,
        existing_adapter_material_lost_mm3=max(0., float((adapter - actual).volume())),
        neutral_hits_mm3=neutral, K7_added_material_samples=checked,
        scope="Connected replacement support retains the original adapter; the existing jaw is joined "
              "after local clearance around shared fasteners. This checks fixed geometry only. "
              "Assembly sequence and the moving eye and sleeve must be rechecked with this support."))
    print(len(checked), "K7 poses;", sum(bool(row["incremental_hits_mm3"]) for row in checked), "blocked", flush=True)


if __name__ == "__main__":
    main()
