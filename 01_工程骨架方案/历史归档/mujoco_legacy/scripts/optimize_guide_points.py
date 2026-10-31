"""Choose one rigid wire-guide point per cable channel.

The guide stands in for the muscle's origin region: a fixture bolted to a bone
that the wire passes through.  Because it decides the direction of the wire
where it leaves the moving bone, it sets the moment arm the cable develops
about the target joint.

Candidate guides are only ever taken from the *surface* of an allowed bone,
offset outwards along the triangle normal by a clearance, so the fixture is a
point on a bone that could actually be machined.  The earlier version searched
a free box around the joint and picked points floating tens of millimetres off
the arm, which no fixture can reach and which made wires cross the thorax.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mujoco
import numpy as np

from cable_drive_common import (
    DEFAULT_EXPORT_DIR,
    DEFAULT_MANIFEST,
    DEFAULT_OUTPUT_XML,
    JOINT_SPECS,
    SKELETON_JOINT_NAMES,
    load_cable_manifest,
    read_binary_stl,
    vec3,
)


# Model body -> bone mesh.  The mesh rides on the last body of each chain, so
# e.g. the clavicle mesh is attached to clavicle_axial_rotation.
BODY_BONE = {"torso": "thorax"}
for _spec in JOINT_SPECS:
    if _spec["mesh"]:
        BODY_BONE[_spec["body"]] = _spec["mesh"]

# Each entry lists the bones a guide may sit on.  Shoulder muscles get their
# anatomical origin bones first and the thorax as a fallback; the elbow muscles
# may only run along the humerus.
#
# `origin_bones` is the anatomical origin for each shoulder muscle; `torso`
# stands in for a bracket bolted to the trunk frame.  Which set is used is
# chosen by --shoulder-origin.
TARGETS = (
    {
        "muscle": "supraspinatus",
        "joint": "shoulder_abduction",
        "sign": 1.0,
        "q_min_deg": -20.0,
        "q_max_deg": 100.0,
        "origin_bones": ("scapula_external_rotation",),
    },
    {
        "muscle": "infraspinatus",
        "joint": "shoulder_rotation",
        "sign": 1.0,
        "q_min_deg": -60.0,
        "q_max_deg": 60.0,
        "origin_bones": ("scapula_external_rotation",),
    },
    {
        "muscle": "deltoid_anterior",
        "joint": "shoulder_flexion",
        "sign": 1.0,
        "q_min_deg": -20.0,
        "q_max_deg": 100.0,
        "origin_bones": (
            "clavicle_axial_rotation",
            "scapula_external_rotation",
        ),
    },
    {
        "muscle": "deltoid_middle",
        "joint": "shoulder_abduction",
        "sign": 1.0,
        "q_min_deg": -20.0,
        "q_max_deg": 100.0,
        "origin_bones": ("scapula_external_rotation",),
    },
    {
        "muscle": "deltoid_posterior",
        "joint": "shoulder_flexion",
        "sign": -1.0,
        "q_min_deg": -20.0,
        "q_max_deg": 100.0,
        "origin_bones": ("scapula_external_rotation",),
    },
    {
        "muscle": "elbow_flexors",
        "joint": "elbow_flexion",
        "sign": 1.0,
        "q_min_deg": 0.0,
        "q_max_deg": 145.0,
        "guide_bodies": ("shoulder_abduction",),
    },
    {
        "muscle": "triceps",
        "joint": "elbow_flexion",
        "sign": -1.0,
        "q_min_deg": 0.0,
        "q_max_deg": 145.0,
        "guide_bodies": ("shoulder_abduction",),
    },
)


def name_id(model: mujoco.MjModel, obj_type, name: str) -> int:
    value = mujoco.mj_name2id(model, obj_type, name)
    if value < 0:
        raise ValueError(f"Object not found: {name}")
    return int(value)


def load_bone_triangles(bone: str) -> np.ndarray:
    path = DEFAULT_EXPORT_DIR / "bones" / f"{bone}.stl"
    if not path.exists():
        raise FileNotFoundError(f"Missing bone mesh export: {path}")
    triangles = read_binary_stl(path)
    if float(np.max(np.ptp(triangles.reshape(-1, 3), axis=0))) > 5.0:
        triangles = triangles * 0.001
    return triangles


def surface_candidates(
    bone: str, clearance: float, limit: int, seed: int
) -> np.ndarray:
    """Every triangle centroid pushed out along its own surface normal."""

    triangles = load_bone_triangles(bone)
    corner0, corner1, corner2 = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    centroids = (corner0 + corner1 + corner2) / 3.0
    normals = np.cross(corner1 - corner0, corner2 - corner0)
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    lengths[lengths < 1e-12] = 1.0
    points = centroids + (normals / lengths) * clearance
    # A hollow bone exports both an outer wall and a cavity wall, and the
    # cavity triangles push their candidates *into* the bone.  Keep only the
    # points that are outside the mesh volume.
    points = points[outside_mesh(points, triangles)]
    if len(points) > limit:
        rng = np.random.default_rng(seed)
        keep = rng.choice(len(points), size=limit, replace=False)
        points = points[np.sort(keep)]
    return points


def outside_mesh(
    points: np.ndarray, triangles: np.ndarray, chunk: int = 400
) -> np.ndarray:
    """Ray parity along +X: an odd crossing count means inside the mesh."""

    direction = np.array([1.0, 0.0, 0.0])
    corner0, corner1, corner2 = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    edge1 = corner1 - corner0
    edge2 = corner2 - corner0
    pvec = np.cross(np.broadcast_to(direction, edge2.shape), edge2)
    det = np.einsum("ij,ij->i", edge1, pvec)
    usable = np.abs(det) > 1e-12
    inv_det = np.zeros_like(det)
    inv_det[usable] = 1.0 / det[usable]
    inside = np.zeros(len(points), dtype=bool)
    for start in range(0, len(points), chunk):
        block = points[start : start + chunk]
        tvec = block[:, None, :] - corner0[None, :, :]
        u = np.einsum("cij,ij->ci", tvec, pvec) * inv_det[None, :]
        qvec = np.cross(tvec, np.broadcast_to(edge1, tvec.shape))
        v = (
            np.einsum("cij,j->ci", qvec, direction) * inv_det[None, :]
        )
        t = np.einsum("cij,ij->ci", qvec, edge2) * inv_det[None, :]
        hit = (
            usable[None, :]
            & (u >= 0.0)
            & (u <= 1.0)
            & (v >= 0.0)
            & (u + v <= 1.0)
            & (t > 1e-9)
        )
        inside[start : start + chunk] = np.count_nonzero(hit, axis=1) % 2 == 1
    return ~inside


def insertion_pairs(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    body_id: int,
    local: np.ndarray,
    qpos_address: int,
    q_values: np.ndarray,
    epsilon: float,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Insertion point at each sample angle, and at that angle plus epsilon."""

    pairs = []
    for q_value in q_values:
        positions = []
        for offset in (0.0, epsilon):
            data.qpos[:] = 0.0
            data.qpos[qpos_address] = q_value + offset
            mujoco.mj_forward(model, data)
            positions.append(
                data.xpos[body_id]
                + data.xmat[body_id].reshape(3, 3) @ local
            )
        pairs.append((positions[0], positions[1]))
    return pairs


def insertion_offset(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    body_id: int,
    local: np.ndarray,
    address: int,
    qpos_address: int,
    q_value: float,
    epsilon: float,
) -> np.ndarray:
    """Insertion point with one *other* joint knocked off its neutral value."""

    data.qpos[:] = 0.0
    data.qpos[qpos_address] = q_value
    data.qpos[address] = epsilon
    mujoco.mj_forward(model, data)
    return (
        data.xpos[body_id] + data.xmat[body_id].reshape(3, 3) @ local
    )


def optimize_target(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    muscle: dict,
    target: dict,
    clearance: float,
    candidate_limit: int,
    seed: int,
    coupling_weight: float = 0.5,
) -> dict:
    guide_bodies = target.get("guide_bodies") or tuple(
        target["origin_bones"]
    ) + ("torso",)
    insertion_body = muscle["insertion_body"]
    insertion_local = vec3(muscle["insertion_local_m"], "insertion_local_m")
    insertion_body_id = name_id(model, mujoco.mjtObj.mjOBJ_BODY, insertion_body)
    joint_id = name_id(model, mujoco.mjtObj.mjOBJ_JOINT, target["joint"])
    qpos_address = int(model.jnt_qposadr[joint_id])

    q_values = np.radians(np.linspace(target["q_min_deg"], target["q_max_deg"], 13))
    epsilon = 1e-4
    pairs = insertion_pairs(
        model, data, insertion_body_id, insertion_local, qpos_address, q_values, epsilon
    )
    # A wire that also torques another joint is a wire whose slider does not do
    # what its label says.  The glenohumeral axes share one centre and the
    # humeral twist is unopposed, so without this term the deltoid cables spend
    # their force twisting the arm to its limit instead of lifting it.
    other_joints = [
        name
        for name in SKELETON_JOINT_NAMES
        if name != target["joint"]
    ]
    other_addresses = {
        name: int(
            model.jnt_qposadr[
                name_id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
            ]
        )
        for name in other_joints
    }

    def worst_other_moments(
        candidates: np.ndarray, q_value: float
    ) -> np.ndarray:
        data.qpos[:] = 0.0
        data.qpos[qpos_address] = q_value
        mujoco.mj_forward(model, data)
        base = (
            data.xpos[insertion_body_id]
            + data.xmat[insertion_body_id].reshape(3, 3) @ insertion_local
        )
        base_length = np.linalg.norm(candidates - base, axis=1)
        worst = np.zeros(len(candidates))
        for address in other_addresses.values():
            moved = insertion_offset(
                model,
                data,
                insertion_body_id,
                insertion_local,
                address,
                qpos_address,
                q_value,
                epsilon,
            )
            moment = np.abs(
                (np.linalg.norm(candidates - moved, axis=1) - base_length)
                / epsilon
            )
            worst = np.maximum(worst, moment)
        return worst

    best: dict | None = None
    for guide_body in guide_bodies:
        bone = BODY_BONE[guide_body]
        candidates = surface_candidates(bone, clearance, candidate_limit, seed)
        moments = np.empty((len(q_values), len(candidates)))
        penalties = np.empty((len(q_values), len(candidates)))
        for index, (base, shifted) in enumerate(pairs):
            base_length = np.linalg.norm(candidates - base, axis=1)
            shifted_length = np.linalg.norm(candidates - shifted, axis=1)
            moments[index] = (
                -target["sign"] * (shifted_length - base_length) / epsilon
            )
            penalties[index] = (
                worst_other_moments(candidates, q_values[index])
                if coupling_weight > 0.0
                else 0.0
            )
        # The worst moment arm over the working range decides the score, so a
        # guide only wins if it pulls the right way at every sampled angle.
        worst = moments.min(axis=0)
        coupling = penalties.max(axis=0)
        best_target = float(worst.max())
        if best_target > 0.0 and coupling_weight > 0.0:
            # Two stages: keep only the points whose worst-case target moment
            # arm is close to the best available, then among those take the one
            # that torques the other joints least.  Penalising the coupling
            # directly would trade the target moment away, because the target
            # moment of a bone-hugging wire is small to begin with.
            threshold = coupling_weight * best_target
            keep = worst >= threshold
            score = np.where(keep, -coupling, -np.inf)
            if not keep.any():
                score = worst
        else:
            # No surface point drives the joint in the right direction; fall
            # back to the least bad one so the channel still has a guide.
            score = worst
        index = int(np.argmax(score))
        candidate = {
            "guide_body": guide_body,
            "guide_bone": bone,
            "guide_world": candidates[index],
            "minimum_signed_moment_arm_m": float(worst[index]),
            "worst_coupling_m": float(coupling[index]),
            "score": float(score[index]),
        }
        if (
            best is None
            or candidate["score"] > best["score"]
        ):
            best = candidate

    assert best is not None
    guide_body_id = name_id(model, mujoco.mjtObj.mjOBJ_BODY, best["guide_body"])
    data.qpos[:] = 0.0
    mujoco.mj_forward(model, data)
    guide_world = best["guide_world"]
    guide_local = data.xmat[guide_body_id].reshape(3, 3).T @ (
        guide_world - data.xpos[guide_body_id]
    )
    return {
        "muscle": muscle["id"],
        "joint": target["joint"],
        "guide_body": best["guide_body"],
        "guide_bone": best["guide_bone"],
        "guide_local_m": guide_local.tolist(),
        "guide_world_m": guide_world.tolist(),
        "minimum_signed_moment_arm_m": best["minimum_signed_moment_arm_m"],
        "worst_coupling_m": best["worst_coupling_m"],
        "coupling_weight": coupling_weight,
        "clearance_m": clearance,
        "penalty_weight": coupling_weight,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_OUTPUT_XML)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_MANIFEST.with_name("guide_optimization.json"),
    )
    parser.add_argument(
        "--only",
        type=str,
        default="",
        help="Comma-separated muscle ids to optimize; empty means all.",
    )
    parser.add_argument(
        "--clearance-mm",
        type=float,
        default=6.0,
        help="How far the fixture sits off the bone surface.",
    )
    parser.add_argument(
        "--candidates",
        type=int,
        default=4000,
        help="Maximum surface candidates sampled per allowed bone.",
    )
    parser.add_argument("--seed", type=int, default=3)
    parser.add_argument(
        "--coupling-weight",
        type=float,
        default=0.5,
        help=(
            "Penalty on the moment arm a guide creates about the other joints. "
            "A wire should mostly drive the joint its slider is labelled with; "
            "set 0 to score on the target moment arm alone."
        ),
    )
    parser.add_argument(
        "--shoulder-origin",
        choices=("torso", "girdle", "both"),
        default="torso",
        help=(
            "Where shoulder-muscle guides may sit.  'torso' bolts the fixture "
            "to the trunk frame, which acts only on the humerus; 'girdle' uses "
            "the anatomical origin bone (clavicle / scapula) and needs the "
            "shoulder girdle to have real muscles or stiff actuators, because "
            "the wire then loads the girdle several times harder than the arm."
        ),
    )
    args = parser.parse_args()

    manifest = load_cable_manifest(args.manifest)
    by_id = {item["id"]: item for item in manifest["cables"]}
    model = mujoco.MjModel.from_xml_path(str(args.model))
    data = mujoco.MjData(model)

    reports = []
    selected = {value.strip() for value in args.only.split(",") if value.strip()}
    clearance = args.clearance_mm / 1000.0
    for target in TARGETS:
        if selected and target["muscle"] not in selected:
            continue
        if "origin_bones" in target:
            origins = tuple(target["origin_bones"])
            if args.shoulder_origin == "torso":
                target = {**target, "guide_bodies": ("torso",)}
            elif args.shoulder_origin == "girdle":
                target = {**target, "guide_bodies": origins}
            else:
                target = {**target, "guide_bodies": origins + ("torso",)}
        muscle = by_id[target["muscle"]]
        report = optimize_target(
            model,
            data,
            muscle,
            target,
            clearance,
            args.candidates,
            args.seed,
            args.coupling_weight,
        )
        muscle["guide_mode"] = "explicit"
        muscle["guides"] = [
            {
                "body": report["guide_body"],
                "local_m": report["guide_local_m"],
                "radius_m": 0.004,
            }
        ]
        reports.append(report)
        print(
            f"{report['muscle']:20s} {report['joint']:20s} "
            f"guide on {report['guide_bone']:9s} "
            f"min_moment={report['minimum_signed_moment_arm_m']:.5f} m  "
            f"worst_coupling={report['worst_coupling_m']:.5f} m"
        )

    args.manifest.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    args.report.write_text(
        json.dumps(reports, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(f"updated manifest: {args.manifest}")
    print(f"report: {args.report}")


if __name__ == "__main__":
    main()
