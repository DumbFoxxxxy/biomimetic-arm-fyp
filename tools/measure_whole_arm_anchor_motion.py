"""Measure signed cable work directions on verified rigid whole-arm paths."""

import sys
import json

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read

sys.path.insert(0, str(ROOT))
from project_runtime import load_mujoco_model


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    motion_file = folder / "whole_arm_shared_motion.json"
    motion = read(motion_file)
    scene_file = folder / "whole_arm_shared_scene_check.json"
    scene = read(scene_file)
    model_file = ROOT / scene["model"]["file"]
    if motion["model"] != scene["model"] or provenance(model_file) != scene["model"]:
        raise ValueError("Motion does not correspond to the shared scene")
    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    transfer = read(ROOT / scene["transfer"]["file"])
    seat_file = folder / "girdle_seat_pair.json"
    seats = read(seat_file)["seats"]
    definitions = {}
    unresolved = {}
    for row in transfer["reference_channels"]:
        owners = (row["entry_body"], row["terminal_body"])
        if not all(owner in scene["neutral_body_origins_mm"] for owner in owners):
            unresolved[row["id"]] = "Floating TPU terminal requires the actual local guide and tensioned pose."
            continue
        definitions[row["id"]] = (owners, (row["entry_xyz_mm"], row["terminal_xyz_mm"]), "equivalent_port_chord")
    for key in ("FJ1459_g3", "FJ1536_g2"):
        endpoints = [seats[f"{key}_{role}"] for role in ("origin", "insertion")]
        definitions[key] = ([row["owner"] for row in endpoints], [row["port_source_mm"] for row in endpoints],
                            "straight_cable_through_prototype_eyes")
    routes = {}
    for key, (owners, points, basis) in definitions.items():
        local = [(np.asarray(point) - scene["neutral_body_origins_mm"][owner]) * .001
                 for owner, point in zip(owners, points)]
        ids = [model.body(owner).id for owner in owners]
        routes[key] = dict(owners=owners, neutral_source_mm=points, basis=basis, motion={})
        for name, path in motion["paths"].items():
            rows = []
            for pose in path["accepted"]:
                data.qpos[:] = pose["qpos"]
                mujoco.mj_forward(model, data)
                positions, jacobians = [], []
                for body_id, anchor in zip(ids, local):
                    xyz = data.xpos[body_id] + data.xmat[body_id].reshape(3, 3) @ anchor
                    jacobian = np.empty((3, model.nv))
                    mujoco.mj_jac(model, data, jacobian, None, xyz, body_id)
                    positions.append(xyz)
                    jacobians.append(jacobian)
                span = positions[1] - positions[0]
                length = np.linalg.norm(span)
                moment = (span / length) @ (jacobians[0] - jacobians[1])
                rows.append(dict(input_deg=pose["input_deg"], length_mm=length.item() * 1000.,
                                 moment_arm_mm=(moment * 1000.).tolist()))
            derivatives = -np.diff([row["length_mm"] for row in rows]) / np.diff(
                np.deg2rad([row["input_deg"] for row in rows]))
            routes[key]["motion"][name] = dict(
                samples=rows, length_excursion_mm=float(np.ptp([row["length_mm"] for row in rows])),
                min_leverage_along_input_mm=float(derivatives.min()) if derivatives.size else None,
                max_leverage_along_input_mm=float(derivatives.max()) if derivatives.size else None)
    dofs = [model.joint(i).name for i in model.dof_jntid]
    target = folder / "whole_arm_anchor_motion.json"
    target.write_text(json.dumps(dict(motion=provenance(motion_file), seats=provenance(seat_file),
                                      dof_coordinates=dofs, routes=routes, unresolved=unresolved,
                                      scope="Signed derivative of equivalent port chord lengths only. "
                                            "Shoulder and elbow wires can require guided routing; "
                                            "actual routed excursion, friction and tension remain open."),
                                 ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(routes)} route proxies; {len(unresolved)} floating endpoints recorded separately")


if __name__ == "__main__":
    main()
