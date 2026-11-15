"""Generate a conformal T12 mounting socket and tabletop foot for the shared arm."""

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, save
from project_runtime import load_mujoco_model


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    check_file = folder / "whole_arm_shared_scene_check.json"
    check = read(check_file)
    motion_file = folder / "whole_arm_shared_motion.json"
    motion = read(motion_file)
    model_file = ROOT / check["model"]["file"]
    if provenance(model_file) != check["model"] or motion["model"] != check["model"]:
        raise ValueError("The base must be checked on the current shared scene and motion")
    transfer = read(ROOT / check["transfer"]["file"])
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"]}
    solids = {key: load_component(source) for key, source in sources.items()}
    owners = transfer["owner_by_source"]
    origins = {key: np.asarray(check["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    receiver = "thorax_reference:FJ3156"
    t12 = solids[receiver]
    mf = A.I._manifold()
    plate_size = (310., 262., 10.)
    plate_min = np.array([-155., -232., 929.])
    post_size = (32., 24., 109.)
    post_min = np.array([-16., -56., 938.])
    seat_size = (26., 18., 35.)
    seat_min = np.array([-13., -51., 1035.])
    plate = mf.Manifold.cube(plate_size).translate(plate_min)
    post = mf.Manifold.cube(post_size).translate(post_min)
    seat_stock = mf.Manifold.cube(seat_size).translate(seat_min)
    # Offset the conformal pocket slightly to leave room for a bond line and
    # avoid numerically coincident bone/socket triangles in the overlap check.
    pocket = t12.minkowski_sum(mf.Manifold.sphere(.2, 12))
    seat = seat_stock - pocket
    base = ((plate + post + seat_stock) - pocket).as_original().simplify(1e-6)
    rebuilt, meta = export(base, folder / "tabletop_T12_socket.obj", plate_min)
    t12_contact_volume = float((t12 ^ seat_stock).volume())
    if t12_contact_volume < 100. or seat.volume() <= 0.:
        raise ValueError("The conformal seat does not engage the T12 surface")
    fixed_hits = {key: volume for key, part in solids.items() if owners[key] == "thorax"
                  if (volume := A.I.intersection_volume(rebuilt, part)) > .001}
    if fixed_hits:
        raise ValueError(f"Tabletop base overlaps the fixed thorax: {fixed_hits}")

    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    moving = [(key, part) for key, part in solids.items() if owners[key] != "thorax"]
    paths = {}
    for name, branch in motion["paths"].items():
        min_gap = 15.
        first_hit = None
        for row in branch["accepted"]:
            data.qpos[:] = row["qpos"]
            mujoco.mj_forward(model, data)
            pose = transforms(model, data, origins)
            for key, solid in moving:
                placed = solid.transform(pose[owners[key]])
                volume = A.I.intersection_volume(rebuilt, placed)
                if volume > .001 and first_hit is None:
                    first_hit = dict(input_deg=row["input_deg"], part=key, overlap_mm3=volume)
                if min_gap > 0.:
                    min_gap = min(min_gap, A.I.gap(rebuilt, placed, min_gap))
        paths[name] = dict(checked_outward_samples=len(branch["accepted"]),
                           same_pose_return_samples=branch["checked_return"],
                           minimum_clearance_capped_15mm=min_gap, first_hit=first_hit)

    save(folder / "tabletop_base_check.json", dict(
        generator=provenance(ROOT / "tools/build_whole_arm_table_base.py"),
        scene=provenance(check_file), motion=provenance(motion_file), source=meta,
        bone_receiver=receiver, bone_receiver_source=sources[receiver],
        fixed_owner="thorax", table_plane_z_source_mm=plate_min[2],
        plate=dict(size_mm=plate_size, lower_corner_source_mm=plate_min.tolist()),
        post=dict(size_mm=post_size, lower_corner_source_mm=post_min.tolist()),
        seat=dict(size_mm=seat_size, lower_corner_source_mm=seat_min.tolist(),
                  T12_stock_intersection_mm3=t12_contact_volume, nominal_pocket_offset_mm=.2),
        neutral_fixed_overlaps_mm3=fixed_hits, motion_checks=paths,
        scope="Nominal bonded conformal T12 socket. The visual thorax is one fixed body. "
              "Thoracic inter-bone connections, adhesive bond strength, installation sequence, "
              "tip resistance and load limits are pending. Return reuses audited outward poses. "
              "Clearance is to moving rigid geometry; flexible mesh deformation and cable loading are excluded."))
    print(f"Table base exported; {sum(v['checked_outward_samples'] for v in paths.values())} poses; "
          f"{sum(v['first_hit'] is not None for v in paths.values())} paths with base interference", flush=True)


if __name__ == "__main__":
    main()
