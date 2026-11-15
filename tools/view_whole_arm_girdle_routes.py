"""Replay one sampled shoulder-girdle cable path on the whole-arm reference."""

import argparse
import sys
import time

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read

sys.path.insert(0, str(ROOT))
from project_runtime import load_mujoco_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route", help="Shoulder-girdle route ID, e.g. FJ1459_g3")
    parser.add_argument("--check", action="store_true", help="Verify the replay without opening a window")
    args = parser.parse_args()
    report = read(DEST.parent / "whole_arm_first_pass/girdle_route_screen.json")
    model_path = ROOT / report["model"]["file"]
    if provenance(model_path) != report["model"]:
        raise ValueError("The shared model changed since the route screen")
    paths = [(route, branch) for route in report["routes"] if not args.route or route["id"] == args.route
             for branch in route.get("motion", []) if len(branch["outgoing"]) > 1]
    if not paths:
        parser.error("No nonzero screened path for the requested route")
    route, branch = max(paths, key=lambda pair: abs(pair[1]["outgoing"][-1]["angle_deg"]))
    outbound = branch["outgoing"]
    timeline = outbound + outbound[-2::-1]
    model = load_mujoco_model(model_path)
    data = mujoco.MjData(model)
    sc = model.joint("SC").dofadr.item()
    ac = model.joint("AC").dofadr.item()
    local_points = np.asarray(route["neutral_selection"]["body_local_points_mm"]) * .001
    delta = np.zeros(model.nv)

    def frame(row):
        mujoco.mj_resetData(model, data)
        delta[sc + branch["axis"]] = np.deg2rad(row["angle_deg"])
        delta[ac:ac + 3] = np.deg2rad(row["ac_rotvec_deg"])
        mujoco.mj_integratePos(model, data.qpos, delta, 1.)
        mujoco.mj_forward(model, data)
        points = []
        for owner, point in zip(route["owners"], local_points):
            body = data.body(owner)
            points.append(body.xmat.reshape(3, 3) @ point + body.xpos)
        actual = np.linalg.norm(points[1] - points[0]) * 1000
        if not np.isclose(actual, row["length_mm"], atol=1e-7, rtol=0):
            raise ValueError("Replayed cable length differs from the screen")
        return points

    print(f'{route["id"]}: prescribed geometry, SC axis {branch["axis"]}, '
          f'{outbound[-1]["angle_deg"]:.3f} deg endpoint, {len(timeline)} outbound/return frames')
    if args.check:
        for row in timeline:
            frame(row)
        return

    from mujoco import viewer
    with viewer.launch_passive(model, data) as view:
        view.cam.lookat[:] = [-.08, -.08, 1.23]
        view.cam.distance = .8
        view.cam.azimuth = 125
        view.cam.elevation = 8
        start = time.perf_counter()
        while view.is_running():
            index = int((time.perf_counter() - start) * len(timeline) / 6.) % len(timeline)
            with view.lock():
                a, b = frame(timeline[index])
                cable = view.user_scn.geoms[0]
                mujoco.mjv_initGeom(cable, mujoco.mjtGeom.mjGEOM_CAPSULE,
                                    np.zeros(3), np.zeros(3), np.eye(3).ravel(),
                                    np.array([1., .45, .08, 1.], dtype=np.float32))
                mujoco.mjv_connector(cable, mujoco.mjtGeom.mjGEOM_CAPSULE, .0003, a, b)
                view.user_scn.ngeom = 1
            view.sync()
            time.sleep(.015)


if __name__ == "__main__":
    main()
