"""Replay the checked gravity-balanced shoulder–elbow motion with bone meshes."""

import argparse
import csv
from time import sleep
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from study_s1_local_force_geometry import load_reference


INTEGRATION = DEST.parent / "integration_review/s1_series_elbow"
SPOOL_CASES = {
    "spool_adjusted": ("friction", "mu_30_adjusted"),
    "seating_k5": ("seating", "seating_k5"),
}
SPRING_CASES = ("extended", *SPOOL_CASES)


def load_scene(report, case):
    hardware_keys = set()
    if case in SPRING_CASES:
        from audit_s1_series_elbow_spring_integration import load_candidate

        _, extended, _, _, sources, owners, revisions, retention = load_candidate(INTEGRATION)
        sources = {**sources, **revisions}
        hardware_keys = set(revisions) | {key for key in sources if key.startswith("series:")}
        model_source = extended["model"]
        trace_source = extended["trajectory"]
        first = extended["first_interference"]
        if case in SPOOL_CASES:
            study, name = SPOOL_CASES[case]
            spool_ref = provenance(INTEGRATION / f"nominal_spring_spool_{study}.json")
            spool = load_reference(spool_ref)
            selected = next(row for row in spool["cases"] if row["case"] == name)
            audit = read(INTEGRATION / f"nominal_spring_spool_{name}_rigid_audit.json")
            trace_source = selected["trace"]
            if (spool["drive_model"] != model_source or audit["trace"] != trace_source
                    or audit["retention"] != retention or audit["spool_study"] != spool_ref
                    or audit["findings"]):
                raise ValueError("Missing matching clear rigid audit for the spool trajectory")
            first = None
    else:
        baseline = load_reference(report["baseline"])
        source = load_reference(baseline["geometry_assembly"])
        sources, owners = source["sources"], source["owners"]
        model_source = baseline["model"]
        trace_source = report["cases"][case]["trace"]
        first = report["cases"][case]["first_collision"]
    model_path = ROOT / model_source["file"]
    if provenance(model_path) != model_source:
        raise ValueError("Source MuJoCo model changed")
    xml = model_path.read_text(encoding="utf-8")
    physics = mujoco.MjModel.from_xml_string(xml)
    neutral = mujoco.MjData(physics)
    mujoco.mj_kinematics(physics, neutral)
    root = ET.fromstring(xml)
    world = root.find("worldbody")
    asset = ET.SubElement(root, "asset")
    bodies = {row.get("name"): row for row in world.findall(".//body")}
    assets = {}
    for index, (key, row) in enumerate(sources.items()):
        bone = key.startswith(("bone:", "thorax_reference:")) or key == "elbow:elbow_shoulder_pin"
        hardware = key in hardware_keys
        if not bone and not hardware:
            continue
        mesh_file = ROOT / row["file"]
        if provenance(mesh_file)["sha256"] != row["sha256"]:
            raise ValueError(f"Source mesh changed: {key}")
        alias = f"replay_{index}{mesh_file.suffix.lower()}"
        assets[alias] = mesh_file.read_bytes()
        mesh_name = f"replay_{index}"
        ET.SubElement(asset, "mesh", name=mesh_name, file=alias, scale=".001 .001 .001")
        owner = owners[key]
        body = bodies.get(owner, world)
        centre = neutral.body(owner).xpos if owner in bodies else np.zeros(3)
        offset = np.array(row["translation_source_mm"]) / 1000. - centre
        color = ".19 .37 .62 1" if hardware else (".58 .68 .74 .40" if owner == "thorax" else ".89 .81 .64 1")
        ET.SubElement(body, "geom", name=mesh_name, type="mesh", mesh=mesh_name,
                      pos=" ".join(map(str, offset)), mass="0", contype="0", conaffinity="0",
                      rgba=color)
    trace_path = ROOT / trace_source["file"]
    if provenance(trace_path) != trace_source:
        raise ValueError("Source trajectory changed")
    with trace_path.open(newline="", encoding="utf-8") as file:
        samples = list(csv.DictReader(file))
    if case in SPOOL_CASES and audit["sampled_frames"] != len(samples):
        raise ValueError("Spool trajectory and rigid audit cover different numbers of frames")
    rows = [row for row in samples if first is None or float(row["time_s"]) < first["time_s"]]
    frames = np.array([[float(row[f"qpos_{i}"]) for i in range(physics.nq)] for row in rows])
    model = mujoco.MjModel.from_xml_string(ET.tostring(root, encoding="unicode"), assets)
    if model.nq != physics.nq or model.nv != physics.nv or not len(frames) or not np.isfinite(frames).all():
        raise ValueError("Replay geometry changed the recorded state space")
    return model, frames, len(assets)


def draw_return_spring(scene, data, specification):
    from build_s1_series_elbow_return_spring import spring_path

    mount, design, neutral = specification
    start, end = (np.asarray(mount[key]) for key in ("fixed_eye_source_mm", "moving_eye_neutral_source_mm"))
    axis = (end-start) / np.linalg.norm(end-start)
    pin = np.asarray(mount["pin_axis_source"])
    travel = max(0., data.joint("FJ1477_series_travel").qpos.item() * 1000.)
    radius = .5 * design["wire_diameter_mm"]
    mean = .5 * (design["coil_OD_mm"] - design["wire_diameter_mm"])
    points = spring_path(np.linalg.norm(end-start)+travel, design["active_turns"]+2, radius, mean)
    points = points[np.r_[np.arange(0, len(points)-1, 4), len(points)-1]]
    h = data.body("humerus")
    points = ((points @ np.array([pin, np.cross(axis, pin), axis]) + start) / 1000. - neutral) @ h.xmat.reshape(3, 3).T + h.xpos
    if len(points)-1 > scene.maxgeom:
        raise ValueError("Spring display exceeds the MuJoCo scene capacity")
    connector = mujoco.mjtGeom.mjGEOM_CAPSULE
    scene.ngeom = len(points)-1
    for i, (a, b) in enumerate(zip(points[:-1], points[1:], strict=True)):
        geom = scene.geoms[i]
        mujoco.mjv_initGeom(geom, connector, np.zeros(3), np.zeros(3), np.eye(3).ravel(), np.array([.88, .46, .16, 1.]))
        mujoco.mjv_connector(geom, connector, radius/1000., a, b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default="elbow_load_release",
                        choices=("hold", "elbow_load_release", "shoulder_low", "shoulder_high",
                                 "combined_low", "combined_high", "combined_low_zero_rail_friction",
                                 "timed", *SPRING_CASES))
    parser.add_argument("--check", action="store_true", help="Load all data without opening a window")
    args = parser.parse_args()
    filename = ("gravity_balanced_start.json" if args.case in ("hold", "elbow_load_release")
                else "gravity_balanced_release_timing.json" if args.case == "timed"
                else "gravity_balanced_active_cycles.json")
    report = {} if args.case in SPRING_CASES else read(INTEGRATION / filename)
    case = report["selected"] if args.case == "timed" else args.case
    if case is None:
        raise ValueError("The timing study has no candidate that passed the rigid audit")
    model, frames, assets = load_scene(report, case)
    data = mujoco.MjData(model)
    spring = None
    if case in SPRING_CASES:
        from audit_s1_series_elbow_spring_integration import load_candidate

        _, _, mount, geometry, *_ = load_candidate(INTEGRATION)
        mujoco.mj_kinematics(model, data)
        spring = (mount, geometry["design"], data.body("humerus").xpos.copy())
    for frame in (frames[0], frames[-1]):
        data.qpos[:] = frame
        mujoco.mj_kinematics(model, data)
        if not np.isfinite(data.geom_xpos).all():
            raise ValueError("Invalid replay mesh pose")
    print(f"Loaded {len(frames)} recorded frames and {assets} visual meshes", flush=True)
    if spring is not None:
        address = model.joint("FJ1477_series_travel").qposadr.item()
        data.qpos[:] = frames[int(np.argmax(frames[:, address]))]
        mujoco.mj_kinematics(model, data)
        scene = mujoco.MjvScene(model, maxgeom=1000)
        draw_return_spring(scene, data, spring)
        if not all(np.isfinite(geom.pos).all() for geom in scene.geoms[:scene.ngeom]):
            raise ValueError("Invalid spring display at peak travel")
        print(f"Checked {scene.ngeom} spring display segments", flush=True)
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
                    if spring is not None:
                        draw_return_spring(window.user_scn, data, spring)
                window.sync()
                sleep(.005 if case in SPRING_CASES else report.get("period_s", .01))


if __name__ == "__main__":
    main()
