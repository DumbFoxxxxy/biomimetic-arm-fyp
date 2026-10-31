"""Shared kinematics, geometry and model-building interfaces for V8."""

from __future__ import annotations

import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation


ROOT = Path(__file__).resolve().parent
MUJOCO_ROOT = ROOT.parent
PROJECT_ROOT = MUJOCO_ROOT.parent
CONFIG_PATH = ROOT / "shoulder_chain_v8.json"
TARGETS_PATH = ROOT / "shoulder_targets_v8.json"
OPEN_CHAIN_PATH = ROOT / "v8_shoulder_open_chain.xml"
CONTACT_PATH = ROOT / "v8_shoulder_contact.xml"
EXPORT_ROOT = MUJOCO_ROOT / "exports" / "cable_design_v8"

ASSET_DIR = MUJOCO_ROOT / "assets_tendon_demo"
THORAX_MESH = ASSET_DIR / "thorax.stl"
CLAVICLE_MESH = ASSET_DIR / "clavicle.stl"
SCAPULA_MESH = ASSET_DIR / "scapula.stl"
CONTACT_BAND_MESH = MUJOCO_ROOT / "engineering" / "assets" / "thorax_contact_band.stl"
CONTACT_PAD_MESH = MUJOCO_ROOT / "engineering" / "assets" / "scapula_contact_pad.stl"

DESIGN_TO_MJ = np.array(
    [
        [1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0],
    ],
    dtype=float,
)

_MESH_VERTEX_CACHE: dict[str, np.ndarray] = {}


def load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def normalize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    length = float(np.linalg.norm(values))
    if length < 1e-12:
        raise ValueError("Cannot normalize a zero vector")
    return values / length


def vec(values: np.ndarray) -> str:
    return " ".join(f"{float(value):.9f}" for value in values)


def design_to_mj(point_mm: np.ndarray) -> np.ndarray:
    return (DESIGN_TO_MJ @ np.asarray(point_mm, dtype=float)) * 0.001


def design_vector_to_mj(vector: np.ndarray) -> np.ndarray:
    return normalize(DESIGN_TO_MJ @ np.asarray(vector, dtype=float))


def mj_to_design(point_m: np.ndarray) -> np.ndarray:
    return (DESIGN_TO_MJ.T @ np.asarray(point_m, dtype=float)) * 1000.0


def rotation_design_to_mj(rotation: np.ndarray) -> np.ndarray:
    rotation = np.asarray(rotation, dtype=float)
    return DESIGN_TO_MJ @ rotation @ DESIGN_TO_MJ.T


def rotation_matrix(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    axis = normalize(axis)
    x, y, z = axis
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ],
        dtype=float,
    )


def rotation_vector(rotation: np.ndarray) -> np.ndarray:
    rotation = np.asarray(rotation, dtype=float)
    angle = math.acos(
        float(np.clip((np.trace(rotation) - 1.0) * 0.5, -1.0, 1.0))
    )
    if angle < 1e-10:
        return np.array(
            [
                rotation[2, 1] - rotation[1, 2],
                rotation[0, 2] - rotation[2, 0],
                rotation[1, 0] - rotation[0, 1],
            ],
            dtype=float,
        ) * 0.5
    axis = np.array(
        [
            rotation[2, 1] - rotation[1, 2],
            rotation[0, 2] - rotation[2, 0],
            rotation[1, 0] - rotation[0, 1],
        ],
        dtype=float,
    ) / (2.0 * math.sin(angle))
    return axis * angle


def quat_wxyz_from_design_rotation(rotation: np.ndarray) -> np.ndarray:
    xyzw = Rotation.from_matrix(
        rotation_design_to_mj(rotation)
    ).as_quat()
    return np.array([xyzw[3], xyzw[0], xyzw[1], xyzw[2]], dtype=float)


def quat_from_z_axis(axis: np.ndarray) -> np.ndarray:
    axis = normalize(axis)
    local_z = np.array([0.0, 0.0, 1.0])
    cross = np.cross(local_z, axis)
    dot = float(np.clip(np.dot(local_z, axis), -1.0, 1.0))
    if np.linalg.norm(cross) < 1e-10:
        if dot > 0.0:
            return np.array([1.0, 0.0, 0.0, 0.0])
        return np.array([0.0, 1.0, 0.0, 0.0])
    angle = math.acos(dot)
    axis_rotation = normalize(cross)
    return np.array(
        [
            math.cos(angle * 0.5),
            axis_rotation[0] * math.sin(angle * 0.5),
            axis_rotation[1] * math.sin(angle * 0.5),
            axis_rotation[2] * math.sin(angle * 0.5),
        ],
        dtype=float,
    )


def points(config: dict) -> dict[str, np.ndarray]:
    return {
        name: np.asarray(config["geometry_mm"][name], dtype=float)
        for name in ("SC", "AC", "P", "SUP", "INF")
    }


def source_basis(config: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    point = points(config)
    centre = np.asarray(config["surface"]["center_xz_mm"], dtype=float)
    radial = np.array(
        [
            point["P"][0] - centre[0],
            0.0,
            point["P"][2] - centre[1],
        ],
        dtype=float,
    )
    w_axis = normalize(radial)
    raw_u = point["SUP"] - point["INF"]
    u_axis = normalize(raw_u - np.dot(raw_u, w_axis) * w_axis)
    v_axis = normalize(np.cross(w_axis, u_axis))
    return u_axis, v_axis, w_axis


def neutral_q(config: dict) -> np.ndarray:
    zero = config["zero_pose"]
    return np.radians(
        [
            float(zero["sc_elevation_deg"]),
            float(zero["sc_protraction_deg"]),
            float(zero["ac_u_deg"]),
            float(zero["ac_v_deg"]),
            float(zero["ac_w_deg"]),
        ]
    )


def rotation_from_q(config: dict, q: np.ndarray) -> np.ndarray:
    u0, v0, w0 = source_basis(config)
    sc_rotation = rotation_matrix(np.array([0.0, 0.0, 1.0]), q[0])
    sc_rotation = sc_rotation @ rotation_matrix(
        np.array([0.0, 1.0, 0.0]), q[1]
    )
    ball_rotation = (
        rotation_matrix(u0, q[2])
        @ rotation_matrix(v0, q[3])
        @ rotation_matrix(w0, q[4])
    )
    return sc_rotation @ ball_rotation


def surface_target(config: dict, s_mm: float, t_mm: float, phi_deg: float) -> dict:
    point = points(config)
    surface = config["surface"]
    centre = np.asarray(surface["center_xz_mm"], dtype=float)
    radius = float(surface["radius_mm"]) + float(surface["nominal_offset_mm"])
    p0 = point["P"]
    theta0 = math.atan2(p0[2] - centre[1], p0[0] - centre[0])
    delta = -float(s_mm) / float(surface["radius_mm"])
    theta = theta0 + delta
    u0, v0, w0 = source_basis(config)
    neutral_rotation = rotation_from_q(config, neutral_q(config))
    neutral_w = neutral_rotation @ w0
    arc_rotation = rotation_matrix(np.array([0.0, 1.0, 0.0]), -delta)
    w_axis = normalize(arc_rotation @ neutral_w)
    phi_rotation = rotation_matrix(w_axis, math.radians(float(phi_deg)))
    rotation = phi_rotation @ arc_rotation @ neutral_rotation
    contact = np.array(
        [
            centre[0] + radius * w_axis[0],
            p0[1] + float(t_mm),
            centre[1] + radius * w_axis[2],
        ],
        dtype=float,
    )
    ac = contact - rotation @ (p0 - point["AC"])
    return {
        "s_mm": float(s_mm),
        "t_mm": float(t_mm),
        "phi_deg": float(phi_deg),
        "contact_mm": contact.tolist(),
        "ac_mm": ac.tolist(),
        "rotation": rotation.tolist(),
        "neutral_rotation": neutral_rotation.tolist(),
        "u_axis": (rotation @ u0).tolist(),
        "v_axis": (rotation @ v0).tolist(),
        "w_axis": w_axis.tolist(),
        "theta_rad": theta,
        "theta0_rad": theta0,
    }


def generate_targets(config: dict) -> dict:
    fractions = [float(value) for value in config["targets"]["fractions"]]
    targets = [
        {
            "id": "ZERO",
            "action": "ZERO",
            "fraction": 0.0,
            **surface_target(config, 0.0, 0.0, 0.0),
        }
    ]
    for action in (
        "K6",
        "SCAP_EXT",
        "SCAP_INT",
        "SLIDE_FWD",
        "SLIDE_BACK",
        "AXIAL",
        "AXIAL_NEGATIVE",
    ):
        endpoint = config["targets"][action]
        for fraction in fractions:
            targets.append(
                {
                    "id": f"{action}_{fraction:.2f}",
                    "action": action,
                    "fraction": fraction,
                    **surface_target(
                        config,
                        float(endpoint["s_mm"]) * fraction,
                        float(endpoint["t_mm"]) * fraction,
                        float(endpoint["phi_deg"]) * fraction,
                    ),
                }
            )
    return {
        "schema_version": 8,
        "surface_definition": config["surface"],
        "targets": targets,
    }


def direct_kinematics(config: dict, q: np.ndarray) -> dict:
    point = points(config)
    u0, v0, w0 = source_basis(config)
    q = np.asarray(q, dtype=float)
    sc_rotation = rotation_matrix(np.array([0.0, 0.0, 1.0]), q[0])
    sc_rotation = sc_rotation @ rotation_matrix(
        np.array([0.0, 1.0, 0.0]), q[1]
    )
    ball_rotation = (
        rotation_matrix(u0, q[2])
        @ rotation_matrix(v0, q[3])
        @ rotation_matrix(w0, q[4])
    )
    body_rotation = sc_rotation @ ball_rotation
    ac = point["SC"] + sc_rotation @ (point["AC"] - point["SC"])
    contact = ac + body_rotation @ (point["P"] - point["AC"])
    return {
        "sc_rotation": sc_rotation,
        "ball_rotation": ball_rotation,
        "body_rotation": body_rotation,
        "ac_mm": ac,
        "contact_mm": contact,
        "u_axis": body_rotation @ u0,
        "v_axis": body_rotation @ v0,
        "w_axis": body_rotation @ w0,
    }


def target_rotation(target: dict) -> np.ndarray:
    return np.asarray(target["rotation"], dtype=float)


def target_neutral_rotation(target: dict) -> np.ndarray:
    return np.asarray(target["neutral_rotation"], dtype=float)


def orientation_projection(
    config: dict,
    body_rotation: np.ndarray,
    target: dict,
) -> tuple[float, float, np.ndarray]:
    neutral = target_neutral_rotation(target)
    relative = body_rotation @ neutral.T
    relative_vector = rotation_vector(relative)
    w_axis = np.asarray(target["w_axis"], dtype=float)
    projected_angle = float(np.dot(relative_vector, w_axis))
    target_angle = math.radians(float(target["phi_deg"]))
    error_angle = projected_angle - target_angle
    perpendicular = relative_vector - projected_angle * w_axis
    return error_angle, projected_angle, perpendicular


def pose_error(config: dict, q: np.ndarray, target: dict) -> tuple[np.ndarray, dict]:
    actual = direct_kinematics(config, q)
    position_error = actual["contact_mm"] - np.asarray(
        target["contact_mm"], dtype=float
    )
    phi_error, phi_actual, perpendicular = orientation_projection(
        config,
        actual["body_rotation"],
        target,
    )
    weight = float(config["ik"]["orientation_weight_mm_per_rad"])
    residual = np.concatenate(
        (
            position_error,
            np.array([phi_error * weight]),
            perpendicular * 0.05,
        )
    )
    metrics = {
        "position_error_mm": float(np.linalg.norm(position_error)),
        "orientation_error_deg": float(math.degrees(abs(phi_error))),
        "contact_error_vector_mm": position_error.tolist(),
        "phi_actual_deg": float(math.degrees(phi_actual)),
        "perpendicular_tilt_deg": float(
            math.degrees(np.linalg.norm(perpendicular))
        ),
    }
    return residual, metrics


def solve_ik(
    config: dict,
    target: dict,
    initial_q: np.ndarray | None = None,
) -> dict:
    sc_elevation = config["joints"]["sc_elevation"]["range_deg"]
    sc_protraction = config["joints"]["sc_protraction"]["range_deg"]
    ball_limit = float(
        config["joints"]["ac_ball"]["component_limit_deg"]
    )
    lower = np.radians(
        [
            float(sc_elevation[0]),
            float(sc_protraction[0]),
            -ball_limit,
            -ball_limit,
            -ball_limit,
        ]
    )
    upper = np.radians(
        [
            float(sc_elevation[1]),
            float(sc_protraction[1]),
            ball_limit,
            ball_limit,
            ball_limit,
        ]
    )

    def objective(values: np.ndarray) -> np.ndarray:
        residual, _ = pose_error(config, values, target)
        return residual

    starts = []
    if initial_q is not None:
        starts.append(np.clip(np.asarray(initial_q, dtype=float), lower, upper))
    starts.append(np.clip(neutral_q(config), lower, upper))
    starts.append(np.zeros(5, dtype=float))
    best = None
    for start in starts:
        result = least_squares(
            objective,
            start,
            bounds=(lower, upper),
            ftol=1e-13,
            xtol=1e-13,
            gtol=1e-13,
            max_nfev=int(config["ik"]["max_nfev"]),
        )
        _, metrics = pose_error(config, result.x, target)
        score = (
            metrics["position_error_mm"]
            + math.radians(metrics["orientation_error_deg"]) * 100.0
        )
        if best is None or score < best[0]:
            best = (
                score,
                result,
                metrics,
            )
        if (
            metrics["position_error_mm"]
            <= float(config["ik"]["position_tolerance_mm"])
            and metrics["orientation_error_deg"]
            <= float(config["ik"]["orientation_tolerance_deg"])
        ):
            break
    assert best is not None
    _, result, metrics = best
    return {
        "q_rad": result.x.tolist(),
        "q_deg": np.degrees(result.x).tolist(),
        "ac_u_deg": float(np.degrees(result.x[2])),
        "ac_v_deg": float(np.degrees(result.x[3])),
        "ac_w_deg": float(np.degrees(result.x[4])),
        "success": bool(result.success),
        "nfev": int(result.nfev),
        "cost": float(result.cost),
        **metrics,
    }


def body_local_from_design(design_vector_mm: np.ndarray) -> np.ndarray:
    return design_to_mj(design_vector_mm)


def build_shoulder_tree(config: dict, roller_count: int = 0) -> ET.ElementTree:
    if roller_count not in (0, 1, 2):
        raise ValueError("roller_count must be 0, 1 or 2")
    point = points(config)
    joint = config["joints"]
    contact = config["contact"]
    normal_slide_half_range = (
        float(contact["roller"]["normal_travel_mm"]) * 0.0005
    )
    stiffness = float(contact["roller"]["stiffness_N_per_m"])
    preload = float(contact["roller"]["preload_N"])
    spring_ref = -preload / stiffness
    roller_radius = float(contact["roller"]["radius_mm"])
    roller_half_width = float(contact["roller"]["width_mm"]) * 0.0005
    zero_target = surface_target(config, 0.0, 0.0, 0.0)
    zero_rotation = target_rotation(zero_target)
    zero_ac = np.asarray(zero_target["ac_mm"], dtype=float)
    zero_w = np.asarray(zero_target["w_axis"], dtype=float)
    zero_u = np.asarray(zero_target["u_axis"], dtype=float)
    zero_contact = np.asarray(zero_target["contact_mm"], dtype=float)
    roller_centres = []
    if roller_count == 2:
        spacing = float(contact["roller"]["dual_spacing_mm"]) * 0.5
        roller_centres = [
            zero_contact + zero_w * roller_radius + np.array([0.0, -spacing, 0.0]),
            zero_contact + zero_w * roller_radius + np.array([0.0, spacing, 0.0]),
        ]
    elif roller_count == 1:
        roller_centres = [zero_contact + zero_w * roller_radius]

    root = ET.Element("mujoco", {"model": f"v8_shoulder_five_dof_r{roller_count}"})
    ET.SubElement(
        root,
        "compiler",
        {
            "angle": "degree",
            "autolimits": "true",
            "meshdir": str(ASSET_DIR),
        },
    )
    ET.SubElement(
        root,
        "option",
        {
            "timestep": "0.0005",
            "gravity": "0 0 -9.81",
            "integrator": "implicitfast",
        },
    )
    visual = ET.SubElement(root, "visual")
    ET.SubElement(visual, "global", {"azimuth": "45", "elevation": "-16"})
    ET.SubElement(
        visual,
        "headlight",
        {"diffuse": "0.8 0.8 0.8", "ambient": "0.3 0.3 0.3"},
    )
    asset = ET.SubElement(root, "asset")
    ET.SubElement(asset, "mesh", {"name": "thorax", "file": str(THORAX_MESH)})
    ET.SubElement(asset, "mesh", {"name": "clavicle", "file": str(CLAVICLE_MESH)})
    ET.SubElement(asset, "mesh", {"name": "scapula", "file": str(SCAPULA_MESH)})
    ET.SubElement(
        asset,
        "mesh",
        {"name": "v8_contact_band", "file": str(CONTACT_BAND_MESH)},
    )
    ET.SubElement(
        asset,
        "mesh",
        {"name": "v8_contact_pad", "file": str(CONTACT_PAD_MESH)},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "thorax_mat", "rgba": "0.45 0.55 0.65 0.55"},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "clavicle_mat", "rgba": "0.72 0.55 0.34 1"},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "scapula_mat", "rgba": "0.30 0.72 0.52 1"},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "band_mat", "rgba": "0.10 0.75 0.78 0.42"},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "pad_mat", "rgba": "0.95 0.48 0.12 0.75"},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "roller_mat", "rgba": "0.18 0.65 0.82 1"},
    )

    world = ET.SubElement(root, "worldbody")
    ET.SubElement(
        world,
        "light",
        {
            "name": "key",
            "pos": "-0.4 -0.5 1.5",
            "dir": "0.4 0.3 -1",
            "directional": "true",
        },
    )
    ET.SubElement(
        world,
        "geom",
        {
            "name": "floor",
            "type": "plane",
            "size": "1.5 1.5 0.1",
            "pos": "0 0 -0.45",
            "rgba": "0.12 0.14 0.17 1",
            "contype": "0",
            "conaffinity": "0",
        },
    )
    thorax = ET.SubElement(world, "body", {"name": "thorax", "pos": "0 0 0"})
    ET.SubElement(
        thorax,
        "geom",
        {
            "name": "thorax_geom",
            "type": "mesh",
            "mesh": "thorax",
            "material": "thorax_mat",
            "contype": "0",
            "conaffinity": "0",
        },
    )
    band_geom = {
        "name": "v8_thorax_contact_band",
        "type": "mesh",
        "mesh": "v8_contact_band",
        "material": "band_mat",
        "density": "0",
        "group": "1",
    }
    if roller_count:
        band_geom.update(
            {
                "contype": "1",
                "conaffinity": "2",
                "condim": str(contact["band"]["condim"]),
                "friction": vec(contact["band"]["friction"]),
                "solref": vec(contact["band"]["solref"]),
                "solimp": vec(contact["band"]["solimp"]),
            }
        )
    else:
        band_geom.update({"contype": "0", "conaffinity": "0"})
    ET.SubElement(thorax, "geom", band_geom)

    sc1 = ET.SubElement(
        world,
        "body",
        {"name": "sc_elevation_body", "pos": vec(design_to_mj(point["SC"]))},
    )
    ET.SubElement(
        sc1,
        "joint",
        {
            "name": "sc_elevation",
            "type": "hinge",
            "axis": vec(design_vector_to_mj(joint["sc_elevation"]["axis_design"])),
            "range": vec(joint["sc_elevation"]["range_deg"]),
            "damping": str(joint["sc_elevation"]["damping"]),
            "armature": str(joint["sc_elevation"]["armature"]),
        },
    )
    ET.SubElement(
        sc1,
        "inertial",
        {"pos": "0 0 0", "mass": "0.005", "diaginertia": "1e-6 1e-6 1e-6"},
    )
    sc2 = ET.SubElement(sc1, "body", {"name": "sc_protraction_body", "pos": "0 0 0"})
    ET.SubElement(
        sc2,
        "joint",
        {
            "name": "sc_protraction",
            "type": "hinge",
            "axis": vec(
                design_vector_to_mj(joint["sc_protraction"]["axis_design"])
            ),
            "range": vec(joint["sc_protraction"]["range_deg"]),
            "damping": str(joint["sc_protraction"]["damping"]),
            "armature": str(joint["sc_protraction"]["armature"]),
        },
    )
    ET.SubElement(
        sc2,
        "inertial",
        {"pos": "0 0 0", "mass": "0.005", "diaginertia": "1e-6 1e-6 1e-6"},
    )
    clavicle = ET.SubElement(
        sc2,
        "body",
        {"name": "clavicle_body", "pos": "0 0 0"},
    )
    ET.SubElement(
        clavicle,
        "inertial",
        {"pos": "0 0 0", "mass": "0.02", "diaginertia": "2e-6 2e-6 2e-6"},
    )
    ET.SubElement(
        clavicle,
        "geom",
        {
            "name": "clavicle_mesh",
            "type": "mesh",
            "mesh": "clavicle",
            "pos": vec(
                body_local_from_design(
                    np.asarray(
                        config["geometry_mm"]["clavicle_mesh_pos_from_sc_mm"],
                        dtype=float,
                    )
                )
            ),
            "material": "clavicle_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "1100",
        },
    )
    scapula = ET.SubElement(
        clavicle,
        "body",
        {
            "name": "scapula",
            "pos": vec(design_to_mj(point["AC"] - point["SC"])),
        },
    )
    ET.SubElement(
        scapula,
        "joint",
        {
            "name": "ac_ball",
            "type": "ball",
            "damping": str(joint["ac_ball"]["damping"]),
            "armature": str(joint["ac_ball"]["armature"]),
        },
    )
    ET.SubElement(
        scapula,
        "inertial",
        {"pos": "0 0 0", "mass": "0.05", "diaginertia": "5e-6 5e-6 5e-6"},
    )
    ET.SubElement(
        scapula,
        "geom",
        {
            "name": "scapula_mesh",
            "type": "mesh",
            "mesh": "scapula",
            "pos": vec(
                body_local_from_design(
                    np.asarray(
                        config["geometry_mm"]["scapula_mesh_pos_from_ac_mm"],
                        dtype=float,
                    )
                )
            ),
            "material": "scapula_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "1100",
        },
    )
    ET.SubElement(
        scapula,
        "geom",
        {
            "name": "v8_scapula_contact_pad",
            "type": "mesh",
            "mesh": "v8_contact_pad",
            "material": "pad_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "0",
            "group": "1",
        },
    )
    ET.SubElement(
        scapula,
        "site",
        {
            "name": "v8_contact_marker",
            "type": "sphere",
            "pos": vec(design_to_mj(point["P"] - point["AC"])),
            "size": "0.003",
            "rgba": "0.95 0.55 0.15 0.9",
            "group": "3",
        },
    )
    for index, centre_world in enumerate(roller_centres, start=1):
        suffix = f"_{index}" if roller_count == 2 else ""
        local_centre = zero_rotation.T @ (centre_world - zero_ac)
        local_outward = zero_rotation.T @ zero_w
        local_roller_axis = zero_rotation.T @ np.array([0.0, 1.0, 0.0])
        axis_mj = design_vector_to_mj(local_roller_axis)
        quat = quat_from_z_axis(axis_mj)
        carriage = ET.SubElement(
            scapula,
            "body",
            {"name": f"v8_roller_carriage{suffix}", "pos": vec(design_to_mj(local_centre))},
        )
        ET.SubElement(
            carriage,
            "joint",
            {
                "name": f"roller_normal{suffix}",
                "type": "slide",
                "axis": vec(design_vector_to_mj(local_outward)),
                "range": f"{-normal_slide_half_range:.9f} {normal_slide_half_range:.9f}",
                "stiffness": str(stiffness),
                "springref": f"{spring_ref:.9f}",
                "damping": "2.0",
            },
        )
        ET.SubElement(
            carriage,
            "geom",
            {
                "name": f"v8_shoulder_roller{suffix}",
                "type": "cylinder",
                "size": f"{roller_radius * 0.001:.9f} {roller_half_width:.9f}",
                "quat": vec(quat),
                "material": "roller_mat",
                "contype": "2",
                "conaffinity": "1",
                "condim": "3",
                "friction": vec(contact["roller"]["friction"]),
                "solref": vec(contact["band"]["solref"]),
                "solimp": vec(contact["band"]["solimp"]),
                "density": "1400",
                "group": "1",
            },
        )
        ET.SubElement(
            carriage,
            "site",
            {
                "name": f"v8_roller_marker{suffix}",
                "type": "sphere",
                "pos": "0 0 0",
                "size": "0.003",
                "rgba": "0.1 0.85 1 0.9",
                "group": "3",
            },
        )
    ET.SubElement(root, "tendon")
    ET.SubElement(root, "actuator")
    ET.SubElement(root, "sensor")
    ET.indent(root, space="  ")
    return ET.ElementTree(root)


def write_model(tree: ET.ElementTree, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(path, encoding="utf-8", xml_declaration=True)


def compile_tree(tree: ET.ElementTree) -> mujoco.MjModel:
    xml = ET.tostring(tree.getroot(), encoding="unicode")
    return mujoco.MjModel.from_xml_string(xml)


def object_id(model: mujoco.MjModel, kind: mujoco.mjtObj, name: str) -> int:
    value = int(mujoco.mj_name2id(model, kind, name))
    if value < 0:
        raise KeyError(name)
    return value


def geom_mesh_vertices(model: mujoco.MjModel, name: str) -> np.ndarray:
    geom = object_id(model, mujoco.mjtObj.mjOBJ_GEOM, name)
    mesh = int(model.geom_dataid[geom])
    if mesh < 0:
        raise ValueError(f"{name} is not a mesh geom")
    start = int(model.mesh_vertadr[mesh])
    count = int(model.mesh_vertnum[mesh])
    return model.mesh_vert[start : start + count].copy()


def scapula_mesh_vertices_design(config: dict) -> np.ndarray:
    key = f"scapula:{CONFIG_PATH}"
    if key not in _MESH_VERTEX_CACHE:
        model = compile_tree(build_shoulder_tree(config, roller_count=0))
        data = mujoco.MjData(model)
        mujoco.mj_forward(model, data)
        geom = object_id(model, mujoco.mjtObj.mjOBJ_GEOM, "scapula_mesh")
        body = object_id(model, mujoco.mjtObj.mjOBJ_BODY, "scapula")
        raw = geom_mesh_vertices(model, "scapula_mesh")
        world = (
            data.geom_xmat[geom].reshape(3, 3) @ raw.T
        ).T + data.geom_xpos[geom]
        relative_mj = (
            data.xmat[body].reshape(3, 3).T
            @ (
                world - data.xpos[body]
            ).T
        ).T
        _MESH_VERTEX_CACHE[key] = (
            DESIGN_TO_MJ.T @ relative_mj.T
        ).T * 1000.0
    return _MESH_VERTEX_CACHE[key].copy()


def direct_mesh_clearance(config: dict, q: np.ndarray) -> float:
    kinematics = direct_kinematics(config, q)
    vertices = scapula_mesh_vertices_design(config)
    world = (
        kinematics["body_rotation"] @ vertices.T
    ).T + kinematics["ac_mm"]
    return float(np.min(thorax_clearance(world, config)))


def geom_world_vertices(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    geom_name: str,
    vertices_local: np.ndarray,
) -> np.ndarray:
    geom = object_id(model, mujoco.mjtObj.mjOBJ_GEOM, geom_name)
    rotation = data.geom_xmat[geom].reshape(3, 3)
    position = data.geom_xpos[geom]
    return (rotation @ vertices_local.T).T + position


def thorax_clearance(points_design: np.ndarray, config: dict) -> np.ndarray:
    points_design = np.asarray(points_design, dtype=float)
    x, y, z = points_design.T
    centre = np.asarray(config["surface"]["center_xz_mm"], dtype=float)
    radius = float(config["surface"]["radius_mm"])
    radius_gap = (
        np.linalg.norm(points_design[:, [0, 2]] - centre, axis=1) - radius
    )
    x_range = config["surface"]["x_range_mm"]
    y_range = config["surface"]["y_range_mm"]
    in_band = (
        (x >= float(x_range[0]))
        & (x <= float(x_range[1]))
        & (z >= 0.0)
    )
    ellipse_radius = np.sqrt((x / 140.0) ** 2 + (z / 106.0) ** 2)
    clearance = np.where(
        in_band,
        radius_gap,
        (ellipse_radius - 1.0) * 106.0,
    )
    return np.where(
        (y >= float(y_range[0])) & (y <= float(y_range[1])),
        clearance,
        20.0,
    )


def mesh_clearance(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    geom_name: str,
    vertices_local: np.ndarray,
    config: dict,
) -> float:
    world = geom_world_vertices(model, data, geom_name, vertices_local)
    design = (DESIGN_TO_MJ.T @ world.T).T * 1000.0
    return float(np.min(thorax_clearance(design, config)))


def set_joint_state(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    q: np.ndarray,
) -> None:
    q = np.asarray(q, dtype=float)
    for index, name in enumerate(("sc_elevation", "sc_protraction")):
        joint = object_id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
        qpos = int(model.jnt_qposadr[joint])
        dof = int(model.jnt_dofadr[joint])
        data.qpos[qpos] = q[index]
        data.qvel[dof] = 0.0
    ball = object_id(model, mujoco.mjtObj.mjOBJ_JOINT, "ac_ball")
    qpos = int(model.jnt_qposadr[ball])
    dof = int(model.jnt_dofadr[ball])
    u0, v0, w0 = source_basis(load_config())
    ball_rotation = (
        rotation_matrix(u0, q[2])
        @ rotation_matrix(v0, q[3])
        @ rotation_matrix(w0, q[4])
    )
    data.qpos[qpos : qpos + 4] = quat_wxyz_from_design_rotation(ball_rotation)
    data.qvel[dof : dof + 3] = 0.0


def actual_body_geometry(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    config: dict,
) -> dict:
    point = points(config)
    scapula = object_id(model, mujoco.mjtObj.mjOBJ_BODY, "scapula")
    rotation = data.xmat[scapula].reshape(3, 3)
    position_mj = data.xpos[scapula]
    rotation_design = DESIGN_TO_MJ.T @ rotation @ DESIGN_TO_MJ
    position_design = mj_to_design(position_mj)
    u0, v0, w0 = source_basis(config)
    contact = position_design + rotation_design @ (point["P"] - point["AC"])
    return {
        "rotation": rotation_design,
        "ac_mm": position_design,
        "contact_mm": contact,
        "u_axis": rotation_design @ u0,
        "v_axis": rotation_design @ v0,
        "w_axis": rotation_design @ w0,
    }


def contact_force(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    roller_names: list[str],
) -> tuple[float, int]:
    roller_geoms = {
        object_id(model, mujoco.mjtObj.mjOBJ_GEOM, name)
        for name in roller_names
    }
    band = object_id(
        model,
        mujoco.mjtObj.mjOBJ_GEOM,
        "v8_thorax_contact_band",
    )
    total = 0.0
    count = 0
    force = np.zeros(6, dtype=float)
    for index in range(int(data.ncon)):
        contact = data.contact[index]
        pair = {int(contact.geom1), int(contact.geom2)}
        if band in pair and roller_geoms.intersection(pair):
            mujoco.mj_contactForce(model, data, index, force)
            total += abs(float(force[0]))
            count += 1
    return total, count


def roller_gap_mm(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    marker_name: str,
    config: dict,
) -> float:
    site = object_id(model, mujoco.mjtObj.mjOBJ_SITE, marker_name)
    centre = mj_to_design(data.site_xpos[site])
    surface_center = np.asarray(config["surface"]["center_xz_mm"], dtype=float)
    radius = float(config["surface"]["radius_mm"])
    roller_radius = float(config["contact"]["roller"]["radius_mm"])
    return float(
        np.linalg.norm(centre[[0, 2]] - surface_center)
        - radius
        - roller_radius
    )

