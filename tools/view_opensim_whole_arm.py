"""Show anatomical references, surface candidates and cables on the existing arm."""

import argparse
import time

import mujoco
import numpy as np

from compare_opensim_whole_arm import FOLDER, REPORT, ROOT, load_k1, pose_cable, world_points
from build_channel_first_inventory import provenance, read
from project_runtime import load_mujoco_model
from study_s1_local_force_geometry import load_reference


COLORS = {
    "reference": (.1, .75, 1., 1.),
    "engineering": (1., .55, .12, 1.),
    "anchors": (.9, .25, 1., 1.),
}


def overlay(scene, layer, points, *, markers=False):
    color = np.array(COLORS[layer], dtype=np.float32)
    if markers:
        for point in points:
            geom = scene.geoms[scene.ngeom]
            mujoco.mjv_initGeom(geom, mujoco.mjtGeom.mjGEOM_SPHERE,
                                np.full(3, .003), point * .001, np.eye(3).ravel(), color)
            scene.ngeom += 1
        return
    for start, end in zip(points, points[1:]):
        if np.linalg.norm(end - start) <= 1e-6:
            continue
        geom = scene.geoms[scene.ngeom]
        mujoco.mjv_initGeom(geom, mujoco.mjtGeom.mjGEOM_CAPSULE,
                            np.zeros(3), np.zeros(3), np.eye(3).ravel(), color)
        mujoco.mjv_connector(geom, mujoco.mjtGeom.mjGEOM_CAPSULE, .0007, start * .001, end * .001)
        scene.ngeom += 1


def active_reference_points(path, elbow_deg, shoulder_neutral):
    groups, points = [], []
    for point in path["points"]:
        if point["type"] == "PathPoint":
            points.append(point)
            continue
        if point["type"] == "ConditionalPathPoint" and point["active_range"]:
            coordinate = point["coordinate"].rsplit("/", 1)[-1]
            angle = elbow_deg if coordinate == "r_elbow_flex" else 0.
            lower, upper = np.fromstring(point["active_range"], sep=" ")
            if (coordinate == "r_elbow_flex" or shoulder_neutral) and lower <= np.deg2rad(angle) <= upper:
                points.append(point)
            continue
        if points:
            groups.append(points)
            points = []
    if points:
        groups.append(points)
    return groups


def main():
    report = read(REPORT)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route", choices=report["routes"], action="append",
                        help="Show a cable and its human reference; repeat to add routes")
    parser.add_argument("--action", default="neutral", choices=["neutral", *report["playback"]])
    parser.add_argument("--check", action="store_true", help="Check the complete replay without opening a window")
    args = parser.parse_args()
    chosen = args.route or ["FJ1486", "FJ1477"]
    guides = load_reference(report["guide_points"])["data"]
    assembly = load_reference(report["engineering"])
    if args.action == "K1":
        if provenance(FOLDER / "whole_arm_k1_cap_relief.json") != report["K1"]["source"]:
            raise ValueError("Rebuild the comparison because the K1 record changed")
        _, model = load_k1()
    else:
        ref = assembly["installed"] if args.action == "neutral" else assembly["neutral_reference"]
        path = ROOT / ref["file"]
        if provenance(path) != ref:
            raise ValueError("Comparison assembly changed")
        model = load_mujoco_model(path)
    data = mujoco.MjData(model)
    outbound = report["playback"].get(args.action, [model.qpos0.copy()])
    if not outbound:
        parser.error("No accepted replay frames for this action")
    timeline = list(range(len(outbound))) + list(range(len(outbound) - 2, -1, -1))
    flags = dict.fromkeys(COLORS, True)
    paused = False

    def frame(index, scene):
        source_index = timeline[index]
        data.qpos[:] = outbound[source_index]
        mujoco.mj_forward(model, data)
        scene.ngeom = 0
        for key in chosen:
            route = report["routes"][key]
            source = guides[key]
            if flags["engineering"]:
                if "owners" in source:
                    engineering = pose_cable(data, source)
                    overlay(scene, "engineering", engineering)
                    overlay(scene, "engineering", engineering[[0, -1]], markers=True)
                else:
                    # The cuff model defines a load-dependent free cap. The
                    # known rigid ports may be shown independently of the cap.
                    if args.action == "K7":
                        chords = report["shoulder_motion"]["K7"]["loaded_chords_at_10N"]
                        overlay(scene, "engineering", np.asarray(chords[source_index][key]))
                    elif args.action == "neutral":
                        overlay(scene, "engineering", np.asarray(source["loaded_neutral_tendon_tips_world_mm"]))
                    ports = [world_points(data, owner, points) for owner, points in
                             zip(route["engineering_owners"], route["engineering_ports_local_mm"], strict=True)]
                    overlay(scene, "engineering", ports, markers=True)
            for ref in route["anatomical_paths"]:
                if flags["reference"]:
                    elbow_deg = float(np.rad2deg(data.joint("elbow").qpos.item()))
                    for points in active_reference_points(ref, elbow_deg, args.action in ("neutral", "K1")):
                        overlay(scene, "reference",
                                [world_points(data, p["owner"], p["mapped_local_mm"]) for p in points])
                if flags["anchors"]:
                    overlay(scene, "anchors",
                            [world_points(data, p["owner"], p["surface_local_mm"]) for p in ref["points"]
                             if "surface_local_mm" in p], markers=True)

    print("R = OpenSim points | E = engineering cables and ports | A = surface candidates | Space = pause")
    if args.check:
        scene = mujoco.MjvScene(model, maxgeom=10000)
        peak = 0
        for index in range(len(timeline)):
            frame(index, scene)
            peak = max(peak, scene.ngeom)
        if peak == 0:
            raise ValueError("No overlay geometry was drawn")
        print(f"Checked {len(timeline)} frames, up to {peak} overlay geometries")
        return

    def key_callback(key):
        nonlocal paused
        layer = {ord("R"): "reference", ord("E"): "engineering", ord("A"): "anchors"}.get(key)
        if layer:
            flags[layer] = not flags[layer]
        if key == ord(" "):
            paused = not paused

    from mujoco import viewer
    with viewer.launch_passive(model, data, key_callback=key_callback) as view:
        view.cam.lookat[:] = [-.12, -.08, 1.15]
        view.cam.distance = .85
        view.cam.azimuth = 135
        view.cam.elevation = 10
        view.opt.geomgroup[3] = int(args.action == "neutral")
        time_position = 0.
        previous = time.perf_counter()
        while view.is_running():
            now = time.perf_counter()
            if not paused:
                time_position += now - previous
            previous = now
            with view.lock():
                frame(int(time_position * len(timeline) / 8.) % len(timeline), view.user_scn)
            view.sync()
            time.sleep(.016)


if __name__ == "__main__":
    main()
