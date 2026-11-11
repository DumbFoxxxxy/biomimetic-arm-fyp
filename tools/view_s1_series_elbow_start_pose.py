"""Replay the checked gravity-balanced shoulder–elbow motion with bone meshes."""

import argparse
import csv
from time import sleep
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from study_s1_local_force_geometry import load_reference


def load_scene(report, case):
    baseline = load_reference(report["baseline"])
    model_path = ROOT / baseline["model"]["file"]
    if provenance(model_path) != baseline["model"]:
        raise ValueError("Source MuJoCo model changed")
    xml = model_path.read_text(encoding="utf-8")
    physics = mujoco.MjModel.from_xml_string(xml)
    neutral = mujoco.MjData(physics)
    mujoco.mj_kinematics(physics, neutral)
    source = load_reference(baseline["geometry_assembly"])
    root = ET.fromstring(xml)
    world = root.find("worldbody")
    asset = ET.SubElement(root, "asset")
    bodies = {row.get("name"): row for row in world.findall(".//body")}
    assets = {}
    for index, (key, row) in enumerate(source["sources"].items()):
        if not (key.startswith(("bone:", "thorax_reference:")) or key == "elbow:elbow_shoulder_pin"):
            continue
        mesh_file = ROOT / row["file"]
        if provenance(mesh_file)["sha256"] != row["sha256"]:
            raise ValueError(f"Source mesh changed: {key}")
        alias = f"replay_{index}{mesh_file.suffix.lower()}"
        assets[alias] = mesh_file.read_bytes()
        mesh_name = f"replay_{index}"
        ET.SubElement(asset, "mesh", name=mesh_name, file=alias, scale=".001 .001 .001")
        owner = source["owners"][key]
        body = bodies.get(owner, world)
        centre = neutral.body(owner).xpos if owner in bodies else np.zeros(3)
        offset = np.array(row["translation_source_mm"]) / 1000. - centre
        ET.SubElement(body, "geom", name=mesh_name, type="mesh", mesh=mesh_name,
                      pos=" ".join(map(str, offset)), mass="0", contype="0", conaffinity="0",
                      rgba=(".58 .68 .74 .40" if owner == "thorax" else ".89 .81 .64 1"))
    trace = report["cases"][case]
    trace_path = ROOT / trace["trace"]["file"]
    if provenance(trace_path) != trace["trace"]:
        raise ValueError("Source trajectory changed")
    with trace_path.open(newline="", encoding="utf-8") as file:
        samples = list(csv.DictReader(file))
    first = trace["first_collision"]
    rows = [row for row in samples if first is None or float(row["time_s"]) < first["time_s"]]
    frames = np.array([[float(row[f"qpos_{i}"]) for i in range(physics.nq)] for row in rows])
    model = mujoco.MjModel.from_xml_string(ET.tostring(root, encoding="unicode"), assets)
    if model.nq != physics.nq or model.nv != physics.nv or not len(frames) or not np.isfinite(frames).all():
        raise ValueError("Replay geometry changed the recorded state space")
    return model, frames, len(assets)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default="elbow_load_release",
                        choices=("hold", "elbow_load_release", "shoulder_low", "shoulder_high",
                                 "combined_low", "combined_high", "combined_low_zero_rail_friction", "timed"))
    parser.add_argument("--check", action="store_true", help="Load all data without opening a window")
    args = parser.parse_args()
    filename = ("gravity_balanced_start.json" if args.case in ("hold", "elbow_load_release")
                else "gravity_balanced_release_timing.json" if args.case == "timed"
                else "gravity_balanced_active_cycles.json")
    report = read(DEST.parent / "integration_review/s1_series_elbow" / filename)
    case = report["selected"] if args.case == "timed" else args.case
    if case is None:
        raise ValueError("The timing study has no candidate that passed the rigid audit")
    model, frames, assets = load_scene(report, case)
    data = mujoco.MjData(model)
    for frame in (frames[0], frames[-1]):
        data.qpos[:] = frame
        mujoco.mj_kinematics(model, data)
        if not np.isfinite(data.geom_xpos).all():
            raise ValueError("Invalid replay mesh pose")
    print(f"Loaded {len(frames)} recorded frames and {assets} visual meshes", flush=True)
    if args.check:
        return
    from mujoco import viewer

    with viewer.launch_passive(model, data) as window:
        with window.lock():
            window.cam.lookat[:] = [-.12, -.04, 1.10]
            window.cam.distance = .90
            window.cam.azimuth = 145.
            window.cam.elevation = -12.
        while window.is_running():
            for frame in frames:
                if not window.is_running():
                    break
                with window.lock():
                    data.qpos[:] = frame
                    mujoco.mj_forward(model, data)
                window.sync()
                sleep(report.get("period_s", .01))


if __name__ == "__main__":
    main()
