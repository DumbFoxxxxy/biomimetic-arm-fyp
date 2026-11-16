"""Export audited MuJoCo poses and rigid geometry for the 360-degree replay."""

import json
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_web_model import write_scene
from render_current_arm_actions import DEST, FOLDER, FPS, ROOT, Replay, provenance, read


VIEWER = ROOT / "docs/viewer"
STEM = "whole-arm-actions"


def main():
    replay = Replay()
    record_path = DEST / "current_actions.json"
    record = read(record_path)
    model, data = replay.model, replay.data
    xml = ET.parse(ROOT / record["assembly"]["file"])
    transfer = read(FOLDER / "interface_transfer_audit.json")
    source_keys = list({**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"]})
    parts, owners = [], {}
    for index in range(model.ngeom):
        if not replay.options.geomgroup[model.geom_group[index]]:
            continue
        geom = model.geom(index)
        mesh = model.mesh(geom.dataid.item())
        start, end = mesh.vertadr.item(), mesh.vertadr.item() + mesh.vertnum.item()
        first_face = mesh.faceadr.item()
        name = geom.name
        if name.startswith("whole_arm_mesh_"):
            name = source_keys[int(name.removeprefix("whole_arm_mesh_"))]
        points = model.mesh_vert[start:end] @ data.geom_xmat[index].reshape(3, 3).T + data.geom_xpos[index]
        parts.append((name, points * 1000., model.mesh_face[first_face:first_face + mesh.facenum.item()]))
        owners[name] = geom.bodyid.item()
    if len(owners) != len(parts):
        raise ValueError("Multiple rigid parts have the same name")
    output = VIEWER / f"{STEM}.json"
    write_scene(parts, STEM, "Audited whole-arm motion", ROOT / record["assembly"]["file"],
                scope="Rigid geometry from the same MuJoCo assembly and relieved cap as the nine action replays.",
                evidence={"replay": record_path, "cap": ROOT / record["cap"]["file"], "routing": ROOT / record["routing"]["file"]})
    manifest = read(output)
    body_ids = sorted(set(owners.values()))
    pivot = np.asarray(manifest["source_origin_mm"]) * .001
    initial = data.xpos[body_ids].copy()
    initial_axes = data.xmat[body_ids].reshape(-1, 3, 3).copy()
    export = {}
    for key, action in replay.actions.items():
        samples = []
        for frame in replay.frames(key):
            replay.set_frame(frame)
            rotations = data.xmat[body_ids].reshape(-1, 3, 3) @ initial_axes.transpose(0, 2, 1)
            shifts = data.xpos[body_ids] - pivot - (rotations @ (initial - pivot)[..., None])[..., 0]
            poses = []
            for rotation, shift in zip(rotations, shifts, strict=True):
                quaternion = np.empty(4)
                mujoco.mju_mat2Quat(quaternion, rotation.ravel())
                poses.append([*shift, *quaternion[[1, 2, 3, 0]]])
            points, support = replay.points()
            lines = [[*points], [support, points[2]]]
            lines.extend(path for path, _, _, _ in replay.pulling_lines(key))
            if key in ("K4", "K5"):
                lines.append(replay.rotation_marker())
            samples.append([round(frame["input_value"], 5),
                            np.round(poses, 7).tolist(),
                            [np.round(np.asarray(line) - pivot, 7).tolist() for line in lines]])
        active = action["status"] == "active_return"
        sample_interval = frame["time_s"] / (len(samples) - 1) if active else 1. / FPS
        export[key] = dict(label=action["label"], units=action["units"], active=active,
                           extent=round(action["screened_extent"], 2), sample_interval=sample_interval,
                           samples=samples)
    manifest.update(part_body=[body_ids.index(owners[p["key"]]) for p in manifest["parts"]],
                    motion_file="whole-arm-motion.json")
    output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    motion_path = VIEWER / manifest["motion_file"]
    motion_path.write_text(json.dumps(dict(mesh_sha256=manifest["mesh_sha256"], actions=export),
                                      ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"Exported {len(export)} actions across {len(body_ids)} rigid bodies")


if __name__ == "__main__":
    main()
