"""Generate the ligament strips that define each joint's range of motion.

The source framework builds a joint by unfolding the joint cylinder into a
plane, laying a 2D branch-network pattern on it and projecting it back.  A
branch of that pattern behaves as a folded member: it unfolds freely and only
becomes a stiff spring once it is straight.  This script reproduces that
behaviour one dimension lower:

* candidate strips are laid at angles around the joint's primary axis with the
  two ends at different angles, so the strip is oblique to the axis the way a
  real collateral ligament is,
* the two ends are real bone surface vertices inside a band around the joint
  cylinder, so an attachment always lands on the correct bone,
* a wrap point on the joint cylinder keeps the strip outside the bones,
* the strip's free length makes it go taut at a fraction of the joint's design
  limit, so the last degrees of every motion stiffen gradually instead of
  hitting a hard stop,
* only strips that actually lengthen in a given direction are kept, because
  the ones on the far side could only ever go slack.

Outputs ``config/bp3d_ligament_manifest.json`` and ``config/bp3d_joint_rom.json``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import bp3d_common as B

# Stiffness of the printed 2D pattern once it has unfolded.  The folding is
# represented by the free length instead of by the modulus, so the modulus
# here is the printed TPU95A value: the pattern is soft while it unfolds and
# solid TPU once its branches are straight.
STRIP_MODULUS_PA = 26.0e6
STRIP_WIDTH_M = 0.008
STRIP_THICKNESS_M = 0.002
STRIP_MAX_FORCE_N = 150.0

# Candidate angles around the joint axis, and the angular offset between the
# two ends of a strip.  Without that offset the strip lies parallel to the
# axis and carries no moment at the neutral pose.
BASE_ANGLES = 6
SPIRAL_DEG = 75.0

MAX_ATTACH_PER_DIRECTION = 3
TAUT_FRACTIONS = (0.86, 1.0)
SHORT_LIMIT_DEG = 15.0

# Safety margin between the ligament-defined limit and the hard joint range.
# Keeping the hard stop slightly outside means widening the range further
# changes nothing, which is what makes the ligaments the real constraint.
HARD_RANGE_MARGIN = 1.05

DESIGN_LIMITS = {
    "sternoclavicular": {"elevation": (-25.0, 30.0), "protraction": (-25.0, 25.0)},
    "acromioclavicular": {
        "upward_rotation": (-8.0, 55.0),
        "posterior_tilt": (-20.0, 20.0),
    },
    "glenohumeral": {"abduction": (-20.0, 115.0), "flexion": (-30.0, 115.0)},
    "elbow": {"flexion": (-5.0, 150.0)},
    "proximal_radioulnar": {"valgus": (-8.0, 8.0)},
}


def load_bone_vertices(name: str) -> np.ndarray:
    vertices, _ = B.read_obj(B.ASSET_DIR / f"{name}.obj")
    return vertices


def direction_at(u: np.ndarray, v: np.ndarray, angle_deg: float) -> np.ndarray:
    theta = np.radians(angle_deg)
    return np.cos(theta) * u + np.sin(theta) * v


def band_attachment(
    vertices: np.ndarray,
    centre: np.ndarray,
    axis: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    angle_deg: float,
    r_design: float,
) -> np.ndarray | None:
    """Bone surface vertex that best matches an angle around the joint axis.

    Rays cast from the joint axis are unreliable here: near the shoulder the
    ray often starts inside the parent bone and then leaves it in a direction
    where nothing is hit.  Choosing the closest real surface vertex inside a
    band around the joint cylinder always lands on the bone that the ligament
    attaches to.
    """

    local = vertices - centre
    axial = local @ axis
    radial = local - np.outer(axial, axis)
    radius = np.linalg.norm(radial, axis=1)
    mask = (
        (radius >= 0.15 * r_design)
        & (radius <= 2.5 * r_design)
        & (np.abs(axial) <= 1.8 * r_design)
    )
    if int(mask.sum()) < 4:
        mask = radius <= 2.5 * r_design
    if int(mask.sum()) < 4:
        return None
    selected = radial[mask]
    angle = np.arctan2(selected @ v, selected @ u)
    delta = (angle - np.radians(angle_deg) + np.pi) % (2.0 * np.pi) - np.pi
    score = np.abs(delta) + 0.8 * np.abs(radius[mask] - 0.9 * r_design) / r_design
    return vertices[mask][int(np.argmin(score))]


def rotate(
    point: np.ndarray, centre: np.ndarray, axis: np.ndarray, degrees: float
) -> np.ndarray:
    angle = np.radians(degrees)
    local = point - centre
    cos, sin = np.cos(angle), np.sin(angle)
    return (
        centre
        + local * cos
        + np.cross(axis, local) * sin
        + axis * float(axis @ local) * (1.0 - cos)
    )


def path_length(
    parent_point: np.ndarray, wrap: np.ndarray, child_point: np.ndarray
) -> float:
    return float(
        np.linalg.norm(parent_point - wrap) + np.linalg.norm(wrap - child_point)
    )


def moment_arm(
    parent_point: np.ndarray,
    wrap: np.ndarray,
    child_point: np.ndarray,
    centre: np.ndarray,
    axis: np.ndarray,
    sample_deg: float,
) -> float:
    """dL/dphi about the joint axis, sampled over a small rotation."""

    forward = rotate(child_point, centre, axis, sample_deg)
    backward = rotate(child_point, centre, axis, -sample_deg)
    delta = path_length(parent_point, wrap, forward) - path_length(
        parent_point, wrap, backward
    )
    return delta / np.radians(2.0 * sample_deg)


def candidate_strip(
    centre: np.ndarray,
    axis: np.ndarray,
    radius: float,
    parent_vertices: np.ndarray,
    child_vertices: np.ndarray,
    parent_angle_deg: float,
    child_angle_deg: float,
    limits: tuple[float, float],
    sample_deg: float,
) -> dict[str, object] | None:
    """Geometry and length behaviour of one candidate ligament branch."""

    u, v = B._plane_basis(axis)
    parent_point = band_attachment(
        parent_vertices, centre, axis, u, v, parent_angle_deg, radius
    )
    child_point = band_attachment(
        child_vertices, centre, axis, u, v, child_angle_deg, radius
    )
    if parent_point is None or child_point is None:
        return None
    wrap = centre + radius * direction_at(
        u, v, 0.5 * (parent_angle_deg + child_angle_deg)
    )
    neutral = path_length(parent_point, wrap, child_point)
    if neutral < 0.02:
        return None
    low, high = limits
    stretch_high = (
        path_length(parent_point, wrap, rotate(child_point, centre, axis, high))
        - neutral
    )
    stretch_low = (
        path_length(parent_point, wrap, rotate(child_point, centre, axis, low))
        - neutral
    )
    return {
        "parent_angle_deg": float(parent_angle_deg),
        "child_angle_deg": float(child_angle_deg),
        "parent_point_m": parent_point,
        "wrap_m": wrap,
        "child_point_m": child_point,
        "neutral_length_m": neutral,
        "stretch_high_m": stretch_high,
        "stretch_low_m": stretch_low,
        "moment_arm_m": moment_arm(
            parent_point, wrap, child_point, centre, axis, sample_deg
        ),
    }


def with_taut_angle(
    strip: dict[str, object],
    centre: np.ndarray,
    axis: np.ndarray,
    taut_deg: float,
) -> dict[str, object]:
    """Copy of a strip whose free length makes it taut exactly at ``taut_deg``."""

    parent_point = np.asarray(strip["parent_point_m"], dtype=float)
    wrap = np.asarray(strip["wrap_m"], dtype=float)
    child_point = np.asarray(strip["child_point_m"], dtype=float)
    neutral = float(strip["neutral_length_m"])
    free = max(
        path_length(parent_point, wrap, rotate(child_point, centre, axis, taut_deg)),
        neutral,
    )
    stiffness = STRIP_MODULUS_PA * STRIP_WIDTH_M * STRIP_THICKNESS_M / max(
        free, 1e-6
    )
    updated = dict(strip)
    updated["taut_deg"] = float(taut_deg)
    updated["free_length_m"] = free
    updated["reserve_m"] = free - neutral
    updated["stiffness_n_per_m"] = stiffness
    return updated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=B.CONTRACT_PATH)
    parser.add_argument("--manifest", type=Path, default=B.LIGAMENT_MANIFEST_PATH)
    parser.add_argument("--rom", type=Path, default=B.ROM_PATH)
    parser.add_argument("--max-attachments", type=int, default=MAX_ATTACH_PER_DIRECTION)
    args = parser.parse_args()

    contract = B.load_json(args.contract)
    joints = contract["joints"]
    meshes: dict[str, np.ndarray] = {}

    strips: list[dict[str, object]] = []
    rom_report: dict[str, object] = {}
    print("== ligament strips ==")
    for joint_name, axis_limits in DESIGN_LIMITS.items():
        entry = joints[joint_name]
        centre = np.asarray(entry["origin_m"], dtype=float)
        radius = float(entry["design_radius_m"])
        parent_bone = entry["parent_bone"]
        child_bone = entry["child_bone"]
        for bone in (parent_bone, child_bone):
            if bone not in meshes:
                meshes[bone] = load_bone_vertices(bone)

        joint_report: dict[str, object] = {
            "hard_range_deg": {},
            "axes": {},
            "parent_bone": parent_bone,
            "child_bone": child_bone,
        }
        for axis_name, limits in axis_limits.items():
            axis = B.unit(
                np.asarray(entry["axes"][axis_name], dtype=float),
                f"{joint_name}.{axis_name}",
            )
            joint_report["hard_range_deg"][axis_name] = [
                round(limits[0] * HARD_RANGE_MARGIN, 2),
                round(limits[1] * HARD_RANGE_MARGIN, 2),
            ]
            sample_deg = float(
                np.clip(0.25 * max(abs(limits[0]), abs(limits[1])), 1.0, 10.0)
            )
            candidates: list[dict[str, object]] = []
            for index in range(BASE_ANGLES):
                base = 360.0 * index / BASE_ANGLES
                for spiral in (SPIRAL_DEG, -SPIRAL_DEG):
                    built = candidate_strip(
                        centre,
                        axis,
                        radius,
                        meshes[parent_bone],
                        meshes[child_bone],
                        base,
                        base + spiral,
                        limits,
                        sample_deg,
                    )
                    if built is not None:
                        candidates.append(built)

            axis_report: dict[str, object] = {"directions": {}, "k_phi_nm_per_rad": {}}
            for label, key, limit in (
                ("positive", "stretch_high_m", limits[1]),
                ("negative", "stretch_low_m", limits[0]),
            ):
                if abs(limit) < 1.0:
                    continue
                family = [
                    (float(strip[key]), strip)
                    for strip in candidates
                    if float(strip[key]) > 0.002
                ]
                family.sort(key=lambda item: item[0], reverse=True)
                chosen: list[tuple[float, dict[str, object]]] = []
                seen: set[tuple[int, int]] = set()
                for stretch, strip in family:
                    key_pair = (
                        int(
                            np.argmin(
                                np.linalg.norm(
                                    meshes[parent_bone]
                                    - np.asarray(strip["parent_point_m"], dtype=float),
                                    axis=1,
                                )
                            )
                        ),
                        int(
                            np.argmin(
                                np.linalg.norm(
                                    meshes[child_bone]
                                    - np.asarray(strip["child_point_m"], dtype=float),
                                    axis=1,
                                )
                            )
                        ),
                    )
                    if key_pair in seen:
                        continue
                    seen.add(key_pair)
                    chosen.append((stretch, strip))
                    if len(chosen) >= args.max_attachments:
                        break

                fractions = (
                    TAUT_FRACTIONS
                    if abs(limit) >= SHORT_LIMIT_DEG
                    else (TAUT_FRACTIONS[-1],)
                )
                entries = []
                for order, (stretch, strip) in enumerate(chosen):
                    for fraction in fractions:
                        strip_id = (
                            f"{joint_name}_{axis_name}_{label[0]}"
                            f"{order}{int(round(fraction * 100))}"
                        )
                        final = with_taut_angle(strip, centre, axis, limit * fraction)
                        final.update(
                            {
                                "id": strip_id,
                                "joint": joint_name,
                                "axis": axis_name,
                                "direction": label,
                                "axis_direction": axis,
                                "limit_deg": limit,
                                "parent_bone": parent_bone,
                                "child_bone": child_bone,
                                "width_m": STRIP_WIDTH_M,
                                "thickness_m": STRIP_THICKNESS_M,
                                "max_force_n": STRIP_MAX_FORCE_N,
                                "springlength_m": final["free_length_m"],
                                "sites": [
                                    {
                                        "body": parent_bone,
                                        "pos_m": final["parent_point_m"],
                                    },
                                    {"body": parent_bone, "pos_m": final["wrap_m"]},
                                    {
                                        "body": child_bone,
                                        "pos_m": final["child_point_m"],
                                    },
                                ],
                            }
                        )
                        final["damping"] = float(
                            2.0
                            * 0.6
                            * np.sqrt(final["stiffness_n_per_m"] * 0.05)
                        )
                        strips.append(final)
                        entries.append(
                            {
                                "id": strip_id,
                                "taut_deg": round(float(final["taut_deg"]), 2),
                                "free_mm": round(final["free_length_m"] * 1000.0, 2),
                                "stretch_at_limit_mm": round(
                                    float(stretch) * 1000.0, 2
                                ),
                                "stiffness_n_per_m": round(
                                    float(final["stiffness_n_per_m"]), 1
                                ),
                                "moment_arm_mm": round(
                                    float(final["moment_arm_m"]) * 1000.0, 2
                                ),
                            }
                        )
                axis_report["directions"][label] = {
                    "limit_deg": limit,
                    "strips": entries,
                    "k_phi_nm_per_rad": round(
                        float(
                            sum(
                                entry["stiffness_n_per_m"]
                                * (entry["moment_arm_mm"] / 1000.0) ** 2
                                for entry in entries
                            )
                        ),
                        2,
                    ),
                }
                axis_report["k_phi_nm_per_rad"][label] = axis_report["directions"][
                    label
                ]["k_phi_nm_per_rad"]
            joint_report["axes"][axis_name] = axis_report
        rom_report[joint_name] = joint_report
        for axis_name, axis_report in joint_report["axes"].items():
            counts = {
                label: len(entry["strips"])
                for label, entry in axis_report["directions"].items()
            }
            print(
                f"  {joint_name:20s} {axis_name:18s} strips={counts}  "
                f"k_phi(Nm/rad)={axis_report['k_phi_nm_per_rad']}"
            )

    manifest = {
        "schema_version": 1,
        "units": "meters",
        "notes": (
            "Each strip is one ligament branch of a 2D pattern projected onto "
            "the joint cylinder.  springlength is where the pattern has fully "
            "unfolded; beyond it the strip behaves as a linear tension spring."
        ),
        "strip_defaults": {
            "modulus_pa": STRIP_MODULUS_PA,
            "width_m": STRIP_WIDTH_M,
            "thickness_m": STRIP_THICKNESS_M,
            "max_force_n": STRIP_MAX_FORCE_N,
        },
        "hard_range_margin": HARD_RANGE_MARGIN,
        "strips": strips,
    }
    B.save_json(args.manifest, manifest)
    B.save_json(args.rom, rom_report)
    print(f"\n{len(strips)} ligament strips")
    print(f"manifest: {args.manifest}")
    print(f"rom:      {args.rom}")


if __name__ == "__main__":
    main()
