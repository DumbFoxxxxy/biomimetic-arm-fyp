"""Compare local stop-response poses with the existing rigid assembly."""

import argparse
import csv
from itertools import combinations
from pathlib import Path

import mujoco
import numpy as np

from audit_elbow_cap_travel import load_geometry, rigid_hits
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save
from simulate_elbow_antagonistic_load import JOINTS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attachment-options", action="store_true")
    args = parser.parse_args()
    output = DEST.parent / "elbow_local_muscle_load"
    prefix = "elbow_attachment" if args.attachment_options else "elbow_local_stop"
    response_path = output / f"{prefix}_response.json"
    response = read(response_path)
    check_path = ROOT / response["stop_reference"]["file"]
    check = read(check_path)
    passive_path = ROOT / check["parent"]["file"]
    assembly_path = ROOT / check["assembly"]["file"]
    xml_path = ROOT / response["local_model"]["file"]
    for path, expected in ((check_path, response["stop_reference"]), (passive_path, check["parent"]),
                           (assembly_path, check["assembly"]), (xml_path, response["local_model"])):
        if provenance(path) != expected:
            raise RuntimeError(f"Source fingerprint changed: {path.name}")
    model = mujoco.MjModel.from_xml_string(xml_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    qpos = [model.joint(name).qposadr.item() for name in JOINTS]
    solids, origins = load_geometry(read(passive_path), read(assembly_path))
    owners = [name for name, _ in solids.values()]
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(owners)), 2)
                        if owners[i] != owners[j] and (owners[i] in origins or owners[j] in origins)])
    cases = []
    for case in response["cases"]:
        trace_path = ROOT / case["trace"]["file"]
        if provenance(trace_path) != case["trace"]:
            raise RuntimeError(f"Trace fingerprint changed: {trace_path.name}")
        with trace_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        angles = np.asarray([[float(row[f"{name}_deg"]) for name in JOINTS] for row in rows])
        times = np.asarray([float(row["time_s"]) for row in rows])
        max_step = float(np.max(abs(np.diff(angles, axis=0))))
        if times[-1] < 4. or max_step > 1. or np.any(np.diff(times) <= 0.):
            raise RuntimeError(f"Invalid sampled trace: {trace_path.name}, angular step {max_step:g}")
        first_stop = None
        first_other = None
        last_clear = None
        for index, (time, pose) in enumerate(zip(times, angles, strict=True)):
            data.qpos[qpos] = np.deg2rad(pose)
            mujoco.mj_kinematics(model, data)
            hits = rigid_hits(data, solids, origins, pairs)
            stop_hits = [hit for hit in hits if set(hit["parts"]) == {"stop:arm", "bone:forearm"}]
            other_hits = [hit for hit in hits if hit not in stop_hits]
            if stop_hits and first_stop is None:
                first_stop = dict(time_s=float(time), pose_deg=pose.tolist(), intersections=stop_hits)
            if other_hits:
                first_other = dict(time_s=float(time), pose_deg=pose.tolist(), intersections=other_hits)
                break
            if first_stop is None:
                last_clear = dict(time_s=float(time), pose_deg=pose.tolist())
        result = dict(trace=case["trace"], poses_checked=index + 1, max_saved_joint_step_deg=max_step,
                      last_pose_before_first_overlap=last_clear, first_stop_overlap=first_stop,
                      first_other_overlap=first_other)
        cases.append(result)
        print(trace_path.name, result, flush=True)
    save(output / f"{prefix}_rigid_audit.json", dict(
        generator=provenance(Path(__file__).resolve()), response=provenance(response_path),
        assembly=provenance(assembly_path), rigid_parts=len(solids), varying_pairs=len(pairs),
        threshold_mm3=.001, cases=cases,
        scope="Saved dynamics poses through the first non-stop rigid overlap. The exact stop arm and "
        "forearm material pair is reported separately; its volume is still checked. No TPU, tool, "
        "cable, candidate connector geometry, or between-sample proof."))


if __name__ == "__main__":
    main()
