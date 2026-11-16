"""Recheck the selected scapular routing on the existing whole-arm action traces."""

from dataclasses import dataclass
from pathlib import Path

import mujoco
import numpy as np

from audit_whole_arm_k1_shared import CAP
from compare_opensim_whole_arm import REPORT, load_k1
from refine_scapular_routing import RESULT, restore_stock
from scapular_cable_control import (
    ASM, FOLDER, OUTPUT, ROOT, Geometry, active_model, load_component, provenance, read, save, simulate,
)
from project_runtime import load_mujoco_model


DESTINATION = FOLDER / "action_showcase"
RECORD = DESTINATION / "current_actions.json"


@dataclass
class Motion:
    label: str
    source: Path
    frames: list
    evidence: str
    units: str = "°"


def lift_reduced_pose(full, reduced, qpos):
    """Put the reduced SC/AC angles back into the same full-arm bone frames."""
    tangent = np.zeros(full.nv)
    for index in range(reduced.njnt):
        joint = reduced.joint(index)
        if joint.name not in ("SC", "AC") or joint.type.item() != mujoco.mjtJoint.mjJNT_HINGE:
            raise ValueError("Unexpected reduced shoulder joint")
        start = full.joint(joint.name).dofadr.item()
        tangent[start:start + 3] = joint.axis * qpos[joint.qposadr.item()]
    pose = full.qpos0.copy()
    mujoco.mj_integratePos(full, pose, tangent, 1.)
    return pose.tolist()


def accepted_trace(rows, metric):
    return [dict(qpos=row["qpos"], input_value=row[metric]) for row in rows]


def axial_trace(model, poses):
    start = model.joint("GH").dofadr.item()
    angular = np.zeros(model.nv)
    result = []
    for pose in poses:
        mujoco.mj_differentiatePos(model, angular, 1., model.qpos0, np.asarray(pose))
        result.append(dict(qpos=pose, input_value=float(np.rad2deg(np.linalg.norm(angular[start:start + 3])))))
    return result


def motions(geometry, routes):
    reference = read(REPORT)
    shoulder_file = FOLDER / "whole_arm_shoulder_elevation.json"
    shoulder = read(shoulder_file)["paths"]
    anatomy_file = FOLDER / "whole_arm_anatomical_actions.json"
    anatomy = read(anatomy_file)
    candidate = read(RESULT)
    path = ROOT / candidate["model"]["file"]
    if provenance(path) != candidate["model"]:
        raise ValueError("K6 replay model has changed")
    reduced_k6 = load_mujoco_model(path)
    k6 = [dict(qpos=lift_reduced_pose(geometry.model, reduced_k6, row["qpos"]),
               input_value=row["sc_deg"], time_s=row["time_s"], tension_N=row["tension_N"])
          for row in candidate["best"]["frames"]]
    k7_file = FOLDER / "whole_arm_shared_motion.json"
    k7 = read(k7_file)["paths"]["K7_SC1_pos"]["accepted"][:len(reference["playback"]["K7"])]
    elbow = geometry.model.joint("elbow").qposadr.item()
    result = {
        "K1": Motion("肘屈曲", REPORT, [dict(qpos=q, input_value=np.rad2deg(q[elbow]).item())
                                         for q in reference["playback"]["K1"]], "prescribed"),
        "K2": Motion("肩关节前屈", shoulder_file,
                     accepted_trace(shoulder["K2_GH_flexion"]["accepted"], "input_deg"), "prescribed"),
        "K3": Motion("肩关节外展", shoulder_file,
                     accepted_trace(shoulder["K3_GH_abduction"]["accepted"], "input_deg"), "prescribed"),
        "K4": Motion("肩关节内旋", REPORT, axial_trace(geometry.model, reference["playback"]["K4_internal"]),
                     "prescribed_with_separate_cuff_probes"),
        "K5": Motion("肩关节外旋", REPORT, axial_trace(geometry.model, reference["playback"]["K5_external"]),
                     "prescribed_with_separate_cuff_probes"),
        "K6": Motion("肩胛上回旋", RESULT, k6, "tension_driven"),
        "K7": Motion("肩肘联合抬升", k7_file, accepted_trace(k7, "humerus_direction_change_deg"), "prescribed"),
        "protraction": Motion("肩胛前伸", anatomy_file,
                             accepted_trace(anatomy["scapular"]["protraction"]["accepted"],
                                            "anterior_displacement_mm"), "prescribed", "mm"),
    }
    selected = read(OUTPUT / "scapular_cable_control.json")
    for key in ("protraction", "retraction_alternative"):
        action = anatomy["scapular"][key]
        before = selected["trials"][key]
        _, reduced = active_model(geometry, action, routes)
        driven = simulate(reduced, geometry, action, routes, before["driven"], before["target_deg"],
                          neutral_preload=True)
        direction = 1. if key == "protraction" else -1.
        frames = [dict(qpos=lift_reduced_pose(geometry.model, reduced, row["qpos"]),
                       input_value=direction * row["anterior_mm"], time_s=row["time_s"], tension_N=row["tension_N"])
                  for row in driven["frames"]]
        name = "protraction_pull" if key == "protraction" else "retraction"
        label = "肩胛前伸牵引试验" if key == "protraction" else "肩胛后缩"
        result[name] = Motion(label, OUTPUT / "scapular_cable_control.json", frames, "tension_driven", "mm")
    return result


def interpolate(model, left, right, fraction):
    velocity = np.zeros(model.nv)
    a = np.asarray(left["qpos"])
    b = np.asarray(right["qpos"])
    mujoco.mj_differentiatePos(model, velocity, 1., a, b)
    qpos = a.copy()
    mujoco.mj_integratePos(model, qpos, velocity, fraction)
    return dict(qpos=qpos.tolist(), input_value=(1 - fraction) * left["input_value"] + fraction * right["input_value"])


def check_motion(geometry, routes, motion):
    data = geometry.data

    def test(frame):
        data.qpos[:] = frame["qpos"]
        mujoco.mj_forward(geometry.model, data)
        if not np.isfinite(data.xpos).all():
            raise ValueError(f"Nonfinite pose in {motion.label}")
        hit = geometry.screen(data, routes)["hits"]
        geometry.cached.clear()
        return hit

    active = motion.evidence == "tension_driven"
    if active:
        # Keep the full simulated time series for replay; screen independent
        # 0.04 s samples and the exact recorded extrema.
        extrema = [max(range(len(motion.frames)), key=lambda index: motion.frames[index]["input_value"]),
                   min(range(len(motion.frames)), key=lambda index: motion.frames[index]["input_value"])] if motion.frames else []
        indices = sorted({*range(0, len(motion.frames), 2), *extrema, len(motion.frames) - 1}) if motion.frames else []
        frames = [motion.frames[index] for index in indices]
    else:
        frames = motion.frames
    accepted, first = [], None
    for frame in frames:
        hits = test(frame)
        if not hits:
            accepted.append(frame)
            continue
        first = dict(**frame, hits=hits)
        if accepted and not active:
            left, right = accepted[-1], frame
            while abs(right["input_value"] - left["input_value"]) > .1:
                mid = interpolate(geometry.model, left, right, .5)
                hits = test(mid)
                if hits:
                    right = mid
                    first = dict(**mid, hits=hits)
                else:
                    left = mid
                    accepted.append(mid)
        break
    if active:
        complete = (first is None and bool(frames) and abs(motion.frames[-1]["input_value"]) < .5 and
                    max(row["input_value"] for row in frames) > .05)
        timeline = motion.frames if complete else []
        status = "active_return" if complete else "blocked_active"
    else:
        for frame in reversed(accepted):
            if test(frame):
                raise ValueError(f"Return mismatch in {motion.label}")
        complete = len(accepted) > 1 and accepted[-1]["input_value"] > 0.
        timeline = accepted + accepted[-2::-1] if complete else []
        status = "geometric" if complete else "blocked"
    extent = max((row["input_value"] for row in accepted), default=0.)
    if timeline and active:
        extent = max(row["input_value"] for row in timeline)
    return dict(label=motion.label, status=status, source=provenance(motion.source), evidence=motion.evidence,
                units=motion.units, screened_extent=extent, original_extent=max((f["input_value"] for f in motion.frames), default=0.),
                checked_samples=len(accepted), first_limit=first, playback=timeline)


def main():
    geometry = Geometry()
    cap, common_model = load_k1()
    geometry.parts[CAP] = load_component(cap["cap_relief"]["selected_source"])
    geometry.model = common_model
    geometry.data = mujoco.MjData(common_model)
    source = read(OUTPUT / "scapular_cable_control.json")
    route = read(RESULT)["candidate"]
    selected = [route, source["routes"]["trapezius_superior"], source["routes"]["FJ1536_g2"]]
    restore_stock(geometry, selected)
    actions = {}
    for name, action in motions(geometry, selected).items():
        actions[name] = check_motion(geometry, selected, action)
        checked = actions[name]
        print(name, checked["status"], round(checked["screened_extent"], 4),
              (checked["first_limit"] or {}).get("hits", [])[:1], flush=True)
    DESTINATION.mkdir(parents=True, exist_ok=True)
    save(RECORD, dict(assembly=provenance(ASM), cap=cap["cap_relief"]["selected_source"], routing=provenance(RESULT),
                      actions=actions, scope="Full rigid assembly with K1 cap relief and the new serratus guide. "
                      "Only changed routes/supports were checked against all other rigid parts and the pre-existing local wires. "
                      "Archived action traces provide the prior geometric evidence. TPU and load-dependent cuff motion are outside this recheck."))


if __name__ == "__main__":
    main()
