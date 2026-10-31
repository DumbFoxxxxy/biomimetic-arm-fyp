"""Structural, geometric and kinematic checks for the BodyParts3D arm model.

The important measurement here is the moment-arm matrix: for every tendon and
every named joint axis, ``dL/dq`` at the neutral pose.  It says which wire can
drive which motion and in which direction, which is the only way to tell a
mis-routed muscle from a weak one.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mujoco
import numpy as np

import bp3d_common as B

AXIS_LIST = (
    ("sternoclavicular", "elevation"),
    ("sternoclavicular", "protraction"),
    ("sternoclavicular", "axial_rotation"),
    ("acromioclavicular", "upward_rotation"),
    ("acromioclavicular", "posterior_tilt"),
    ("acromioclavicular", "external_rotation"),
    ("glenohumeral", "abduction"),
    ("glenohumeral", "flexion"),
    ("glenohumeral", "internal_rotation"),
    ("elbow", "flexion"),
    ("proximal_radioulnar", "valgus"),
)


def quaternion_from_axis_angle(axis: np.ndarray, degrees: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    half = np.radians(degrees) * 0.5
    return np.concatenate([[np.cos(half)], axis * np.sin(half)])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=B.MODEL_PATH)
    parser.add_argument("--contract", type=Path, default=B.CONTRACT_PATH)
    parser.add_argument("--muscles", type=Path, default=B.MUSCLE_MANIFEST_PATH)
    parser.add_argument("--delta-deg", type=float, default=0.5)
    parser.add_argument(
        "--report", type=Path, default=B.REPORT_DIR / "check.json"
    )
    args = parser.parse_args()

    contract = B.load_json(args.contract)
    muscles = B.load_json(args.muscles)
    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)
    report: dict[str, object] = {}
    failures: list[str] = []

    # --- structure --------------------------------------------------------
    bodies = [
        mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i)
        for i in range(model.nbody)
    ]
    joints = [
        mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, i)
        for i in range(model.njnt)
    ]
    tendons = [
        mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_TENDON, i)
        for i in range(model.ntendon)
    ]
    expected_wires = sum(
        len(channel["paths"]) for channel in muscles["channels"]
    )
    structure = {
        "bones": [name for name in bodies if name in B.BONE_NAMES],
        "skeleton_joints": [name for name in joints if not name.startswith("slider_")],
        "wires": sum(1 for name in tendons if name.startswith("wire_")),
        "ligaments": sum(1 for name in tendons if name.startswith("lig_")),
        "sliders": sum(1 for name in joints if name.startswith("slider_")),
        "expected_wires": expected_wires,
    }
    report["structure"] = structure
    print("== structure ==")
    print(f"  bones           {len(structure['bones'])} {structure['bones']}")
    print(f"  skeleton joints {len(structure['skeleton_joints'])}")
    print(
        f"  wires           {structure['wires']} (manifest {expected_wires})"
    )
    print(f"  ligaments       {structure['ligaments']}")
    print(f"  sliders         {structure['sliders']}")
    if structure["wires"] != expected_wires:
        failures.append("wire count does not match the manifest")
    if len(structure["bones"]) != len(B.BONE_NAMES):
        failures.append("not every bone was built")

    # --- meshes -----------------------------------------------------------
    missing = [
        str(path)
        for path in (B.ASSET_DIR / f"{name}.obj" for name in B.BONE_NAMES)
        if not path.exists()
    ]
    report["missing_assets"] = missing
    if missing:
        failures.append(f"missing assets: {missing}")
    print(f"  meshes on disk  {len(B.BONE_NAMES) - len(missing)}/{len(B.BONE_NAMES)}")

    # --- joint geometry ---------------------------------------------------
    print("\n== joint geometry ==")
    geometry = {}
    for name, entry in contract["joints"].items():
        gap = float(entry["gap_mm"])
        geometry[name] = {
            "gap_mm": gap,
            "design_radius_mm": round(float(entry["design_radius_m"]) * 1000.0, 1),
            "fit_rms_mm": entry.get("fit_rms_mm"),
        }
        print(
            f"  {name:20s} gap={gap:5.2f} mm  "
            f"collar r={geometry[name]['design_radius_mm']:5.1f} mm"
            + (
                f"  fit rms={entry['fit_rms_mm']} mm"
                if entry.get("fit_rms_mm") is not None
                else ""
            )
        )
        if name in ("sternoclavicular", "acromioclavicular", "elbow") and gap > 3.0:
            failures.append(f"{name} bones are {gap} mm apart")
    report["geometry"] = geometry

    # --- neutral pose -----------------------------------------------------
    data.qpos[:] = 0.0
    data.ctrl[:] = 0.0
    for _ in range(3000):
        B.apply_pull_only_tension(model, data)
        mujoco.mj_step(model, data)
    finite = bool(np.isfinite(data.qpos).all())
    over = [
        mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_TENDON, index)
        for index in range(model.ntendon)
        if data.ten_length[index] > model.tendon_lengthspring[index, 0] + 1e-9
    ]
    report["neutral"] = {
        "finite": finite,
        "taut_tendons": over,
        "energy_j": float(data.energy[0] + data.energy[1]),
    }
    print("\n== neutral pose ==")
    print(f"  finite          {finite}")
    print(f"  taut tendons    {len(over)}")
    if not finite:
        failures.append("the model diverged at the neutral pose")

    # --- moment arms ------------------------------------------------------
    print("\n== moment arms  dL/dq  (millimetres per radian, sign included) ==")
    index_of = {}
    for index in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, index)
        index_of[name] = index
    tendon_names = [
        mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_TENDON, i)
        for i in range(model.ntendon)
    ]
    addressed = {
        name: int(model.jnt_qposadr[index_of[name]])
        for name in contract["joints"]
    }
    arm_table: dict[str, dict[str, float]] = {}
    for joint, axis_name in AXIS_LIST:
        axis = np.asarray(contract["joints"][joint]["axes"][axis_name], dtype=float)
        address = addressed[joint]
        lengths = {}
        for sign in (+1.0, -1.0):
            data.qpos[:] = 0.0
            if model.jnt_type[index_of[joint]] == mujoco.mjtJoint.mjJNT_BALL:
                data.qpos[address : address + 4] = quaternion_from_axis_angle(
                    axis, sign * args.delta_deg
                )
            else:
                data.qpos[address] = np.radians(sign * args.delta_deg)
            mujoco.mj_forward(model, data)
            lengths[sign] = np.array(data.ten_length, dtype=float)
        derivative = (lengths[1.0] - lengths[-1.0]) / np.radians(2.0 * args.delta_deg)
        for name, value in zip(tendon_names, derivative):
            arm_table.setdefault(name, {})[f"{joint}.{axis_name}"] = float(value)

    wire_rows = []
    for channel in muscles["channels"]:
        for path in channel["paths"]:
            name = f"wire_{path['id']}"
            row = arm_table.get(name, {})
            if not row:
                continue
            best = max(row.items(), key=lambda item: abs(item[1]))
            second = sorted(row.items(), key=lambda item: -abs(item[1]))[1]
            wire_rows.append(
                {
                    "wire": path["id"],
                    "channel": channel["id"],
                    "expected": channel["target_joint"],
                    "best_axis": best[0],
                    "best_mm_per_rad": round(best[1] * 1000.0, 2),
                    "second_axis": second[0],
                    "second_mm_per_rad": round(second[1] * 1000.0, 2),
                    "expected_mm_per_rad": round(
                        row.get(_expected_axis(channel["target_joint"]), 0.0) * 1000.0,
                        2,
                    ),
                }
            )
    report["wire_moment_arms"] = wire_rows
    print(f"  {'wire':30s} {'expected axis':32s} {'dL/dq':>9s}  {'strongest axis':32s} {'dL/dq':>9s}")
    for row in wire_rows:
        print(
            f"  {row['wire']:30s} {_expected_axis(row['expected']):32s} "
            f"{row['expected_mm_per_rad']:9.2f}  {row['best_axis']:32s} "
            f"{row['best_mm_per_rad']:9.2f}"
        )

    weak = [
        row["wire"]
        for row in wire_rows
        if abs(row["expected_mm_per_rad"]) < 3.0
    ]
    report["weak_wires"] = weak
    if weak:
        print(f"\n  {len(weak)} wires have less than 3 mm/rad on their target axis")

    # --- ligament slackness ----------------------------------------------
    print("\n== ligaments ==")
    data.qpos[:] = 0.0
    mujoco.mj_forward(model, data)
    taut_at_neutral = 0
    for index in range(model.ntendon):
        name = tendon_names[index]
        if not name.startswith("lig_"):
            continue
        if data.ten_length[index] > model.tendon_lengthspring[index, 0] + 1e-9:
            taut_at_neutral += 1
    print(f"  taut at the neutral pose: {taut_at_neutral}")
    report["ligaments_taut_at_neutral"] = taut_at_neutral
    if taut_at_neutral:
        failures.append("some ligaments are already taut at the neutral pose")

    report["failures"] = failures
    B.save_json(args.report, report)
    print(f"\n{'FAIL' if failures else 'PASS'} - {len(failures)} problem(s)")
    for item in failures:
        print(f"  ! {item}")
    print(f"report: {args.report}")


def _expected_axis(target_joint: str) -> str:
    """Manifest target name -> the joint axis it is supposed to drive."""

    table = {
        "glenohumeral_abduction": "glenohumeral.abduction",
        "glenohumeral_flexion": "glenohumeral.flexion",
        "glenohumeral_extension": "glenohumeral.flexion",
        "glenohumeral_internal_rotation": "glenohumeral.internal_rotation",
        "glenohumeral_external_rotation": "glenohumeral.internal_rotation",
        "glenohumeral_adduction": "glenohumeral.abduction",
        "scapulothoracic_upward_rotation": "acromioclavicular.upward_rotation",
        "scapulothoracic_downward_rotation": "acromioclavicular.upward_rotation",
        "sternoclavicular_depression": "sternoclavicular.elevation",
        "elbow_flexion": "elbow.flexion",
        "elbow_extension": "elbow.flexion",
    }
    return table.get(target_joint, target_joint)


if __name__ == "__main__":
    main()
