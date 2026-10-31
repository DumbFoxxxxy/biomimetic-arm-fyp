"""Build and audit the V7 isolated scapulothoracic roller model."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import copy
import json
from project_runtime import WorkspacePath as Path
import sys
import xml.etree.ElementTree as ET

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
V6_ROOT = ROOT.parent / "cable_design_v6"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from evaluate_v2 import (  # noqa: E402
    geom_mesh_vertices,
    design_vector_to_mj,
    r400_center,
)
from v2_common import POINTS_MM, design_to_mj, mj_to_design, normalize, vec  # noqa: E402


V6_MODEL = V6_ROOT / "arm_cable_design_v6.xml"
V6_CONFIG = V6_ROOT / "cable_design_v6.json"
OUTWARD_DESIGN = normalize(np.array([0.6556, 0.0, 0.7552]))
U_AXIS = normalize(POINTS_MM["SUP"] - POINTS_MM["INF"])
ROLLER_RADIUS_MM = 8.0
ROLLER_WIDTH_MM = 16.0
SPRING_STIFFNESS_N_PER_M = 2000.0
TRAVEL_M = 0.002


class ShoulderGeometry:
    def __init__(self, model: mujoco.MjModel, config: dict):
        self.model = model
        self.scapula_visual_geom = mujoco.mj_name2id(
            model,
            mujoco.mjtObj.mjOBJ_GEOM,
            "scapula_mesh",
        )
        self.contact_pad_geom = mujoco.mj_name2id(
            model,
            mujoco.mjtObj.mjOBJ_GEOM,
            "v2_scapula_contact_pad",
        )
        self.contact_band_geom = mujoco.mj_name2id(
            model,
            mujoco.mjtObj.mjOBJ_GEOM,
            "v2_thorax_contact_band",
        )
        self.scapula_vertices = geom_mesh_vertices(model, "scapula_mesh")
        self.pad_vertices = geom_mesh_vertices(
            model,
            "v2_scapula_contact_pad",
        )
        self.band_centre = r400_center(config)
        band = config["contact"]["thorax_band"]
        self.band_radius = float(band["radius_mm"])
        self.band_x = (
            float(band["x_min_mm"]),
            float(band["x_max_mm"]),
        )
        self.band_y = (
            float(band["y_min_mm"]),
            float(band["y_max_mm"]),
        )

    def thorax_clearance(self, points_design: np.ndarray) -> np.ndarray:
        points = np.asarray(points_design, dtype=float)
        x, y, z = points.T
        radius_pad = (
            np.linalg.norm(
                points[:, [0, 2]] - self.band_centre,
                axis=1,
            )
            - self.band_radius
        )
        in_band = (
            (x >= self.band_x[0])
            & (x <= self.band_x[1])
            & (z >= 0.0)
        )
        ellipse_radius = np.sqrt((x / 140.0) ** 2 + (z / 106.0) ** 2)
        clearance = np.where(
            in_band,
            radius_pad,
            (ellipse_radius - 1.0) * 106.0,
        )
        return np.where(
            (y >= self.band_y[0]) & (y <= self.band_y[1]),
            clearance,
            20.0,
        )

    def visual_clearance(self, data: mujoco.MjData) -> float:
        rotation = data.geom_xmat[self.scapula_visual_geom].reshape(3, 3)
        position = data.geom_xpos[self.scapula_visual_geom]
        world = (rotation @ self.scapula_vertices.T).T + position
        design = np.asarray([mj_to_design(point) for point in world])
        return float(np.min(self.thorax_clearance(design)))

    def pad_gap(self, data: mujoco.MjData) -> float:
        rotation = data.geom_xmat[self.contact_pad_geom].reshape(3, 3)
        position = data.geom_xpos[self.contact_pad_geom]
        world = (rotation @ self.pad_vertices.T).T + position
        design = np.asarray([mj_to_design(point) for point in world])
        valid = (
            (design[:, 0] >= self.band_x[0])
            & (design[:, 0] <= self.band_x[1])
            & (design[:, 1] >= self.band_y[0])
            & (design[:, 1] <= self.band_y[1])
        )
        if not np.any(valid):
            return float("inf")
        radial = (
            np.linalg.norm(
                design[valid][:, [0, 2]] - self.band_centre,
                axis=1,
            )
            - self.band_radius
        )
        return float(np.min(radial))


def remove_body(root: ET.Element, name: str) -> None:
    for parent in root.iter():
        for child in list(parent):
            if child.tag == "body" and child.get("name") == name:
                parent.remove(child)
                return
    raise KeyError(name)


def clear_children(root: ET.Element, path: str) -> None:
    element = root.find(path)
    if element is None:
        element = ET.SubElement(root, path)
    for child in list(element):
        element.remove(child)


def roller_centres(count: int) -> list[np.ndarray]:
    spacing = 40.0 if count == 2 else 0.0
    offsets = [-spacing * 0.5, spacing * 0.5] if count == 2 else [0.0]
    return [
        POINTS_MM["P"]
        + U_AXIS * offset
        + OUTWARD_DESIGN * (ROLLER_RADIUS_MM + 0.35)
        for offset in offsets
    ]


def neutral_scapula_transform(config: dict) -> tuple[np.ndarray, np.ndarray]:
    model = _load_project_model(str(V6_MODEL))
    data = mujoco.MjData(model)
    for joint_name, offset in config.get(
        "zero_pose_offsets_deg", {}
    ).items():
        joint = mujoco.mj_name2id(
            model,
            mujoco.mjtObj.mjOBJ_JOINT,
            joint_name,
        )
        if joint >= 0:
            data.qpos[int(model.jnt_qposadr[joint])] = np.deg2rad(
                float(offset)
            )
    mujoco.mj_forward(model, data)
    body = mujoco.mj_name2id(
        model,
        mujoco.mjtObj.mjOBJ_BODY,
        "scapula",
    )
    return (
        data.xmat[body].reshape(3, 3).copy(),
        data.xpos[body].copy(),
    )


def build_shoulder_model(roller_count: int = 1) -> ET.ElementTree:
    if roller_count not in (1, 2):
        raise ValueError("roller_count must be 1 or 2")
    tree = ET.parse(V6_MODEL)
    root = tree.getroot()
    root.set("model", f"v7_shoulder_standalone_{roller_count}roller")
    remove_body(root, "gh_flex")
    clear_children(root, "tendon")
    clear_children(root, "actuator")
    clear_children(root, "sensor")
    for parent in root.iter():
        for child in list(parent):
            name = child.get("name", "")
            if (
                child.tag in ("geom", "site")
                and (
                    name == "V6_roller"
                    or name.startswith("v6_shoulder_roller")
                    or name.startswith("V2_guide_")
                    or name.startswith("V2_side_")
                    or name in ("V6_E", "V6_W")
                )
            ):
                parent.remove(child)
    scapula = next(
        (
            body
            for body in root.iter("body")
            if body.get("name") == "scapula"
        ),
        None,
    )
    if scapula is None:
        raise RuntimeError("scapula body not found")
    config = json.loads(V6_CONFIG.read_text(encoding="utf-8"))
    scapula_rotation, scapula_position = neutral_scapula_transform(config)
    outward_mj = design_vector_to_mj(OUTWARD_DESIGN)
    local_outward_mj = normalize(scapula_rotation.T @ outward_mj)
    spring_ref = -5.0 / SPRING_STIFFNESS_N_PER_M
    for index, centre_world in enumerate(roller_centres(roller_count), start=1):
        centre_mj = design_to_mj(centre_world)
        local = scapula_rotation.T @ (
            centre_mj - scapula_position
        )
        suffix = f"_{index}" if roller_count == 2 else ""
        carriage = ET.SubElement(
            scapula,
            "body",
            {
                "name": f"v7_roller_carriage{suffix}",
                "pos": vec(local),
            },
        )
        ET.SubElement(
            carriage,
            "joint",
            {
                "name": f"roller_normal{suffix}",
                "type": "slide",
                "axis": vec(local_outward_mj),
                "range": f"{-TRAVEL_M:.6f} {TRAVEL_M:.6f}",
                "stiffness": str(SPRING_STIFFNESS_N_PER_M),
                "springref": f"{spring_ref:.9f}",
                "damping": "2.0",
            },
        )
        ET.SubElement(
            carriage,
            "geom",
            {
                "name": f"v7_shoulder_roller{suffix}",
                "type": "cylinder",
                "size": (
                    f"{ROLLER_RADIUS_MM * 0.001:.9f} "
                    f"{ROLLER_WIDTH_MM * 0.0005:.9f}"
                ),
                "material": "v6_roller_mat",
                "contype": "2",
                "conaffinity": "1",
                "condim": "3",
                "friction": "0.10 0.005 0.0001",
                "solref": "0.002 1",
                "solimp": "0.9 0.95 0.001",
                "density": "1400",
                "group": "1",
            },
        )
        ET.SubElement(
            carriage,
            "site",
            {
                "name": f"v7_roller_marker{suffix}",
                "pos": "0 0 0",
                "size": "0.003",
                "rgba": "0.1 0.85 1 0.9",
                "group": "3",
            },
        )
    ET.indent(root, space="  ")
    return tree


def scalar_joint(model: mujoco.MjModel, name: str) -> int:
    joint = mujoco.mj_name2id(
        model,
        mujoco.mjtObj.mjOBJ_JOINT,
        name,
    )
    if joint < 0:
        raise KeyError(name)
    return joint


def interpolate_correction(
    config: dict,
    action_id: str,
    fraction: float,
) -> dict[str, float]:
    corrections = config.get("path_corrections", {}).get(action_id, {})
    if not corrections:
        return {}
    knots = sorted(float(value) for value in corrections)
    values = np.asarray(
        [
            corrections[str(knot)]
            if str(knot) in corrections
            else corrections[knot]
            for knot in knots
        ],
        dtype=float,
    )
    result = {
        name: float(
            np.interp(float(fraction), knots, values[:, index])
        )
        for index, name in enumerate(
            ("sc_elevation", "sc_protraction", "ac_ext")
        )
    }
    return result


def add_correction(
    q_target: dict[str, float],
    correction: dict[str, float],
    scale: float,
) -> dict[str, float]:
    result = dict(q_target)
    for name, value in correction.items():
        result[name] = float(result.get(name, 0.0)) + value * scale
    return result


def pose_targets(config: dict) -> list[dict]:
    k6 = {
        "sc_elevation": 42.26192461124019,
        "sc_protraction": -8.289513955493888,
        "ac_ext": 15.58573377932018,
    }
    targets = [{"id": "ZERO", "q_target_deg": {}}]
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        q_target = {
            name: value * fraction for name, value in k6.items()
        }
        correction = interpolate_correction(config, "K6", fraction)
        q_target = add_correction(
            q_target,
            correction,
            fraction,
        )
        targets.append(
            {
                "id": f"K6_{fraction:.2f}",
                "q_target_deg": q_target,
            }
        )
    for value in (-5.0, -2.5, 0.0, 2.5, 5.0):
        q_target = {"sc_protraction": value}
        if value < 0.0:
            correction = interpolate_correction(
                config,
                "SLIDE_BACK",
                1.0,
            )
            q_target = add_correction(
                q_target,
                correction,
                abs(value) / 5.0,
            )
        targets.append(
            {
                "id": f"SLIDE_{value:+.1f}",
                "q_target_deg": q_target,
            }
        )
    for value in (0.0, 5.0, 10.0, 15.0, 20.0):
        targets.append(
            {
                "id": f"SCAP_EXT_{value:.1f}",
                "q_target_deg": {"ac_ext": value},
            }
        )
    for value in (0.0, -1.25, -2.5, -3.75, -5.0):
        correction = interpolate_correction(
            config,
            "SCAP_INT",
            1.0,
        )
        targets.append(
            {
                "id": f"SCAP_INT_{value:.2f}",
                "q_target_deg": add_correction(
                    {"ac_ext": value},
                    correction,
                    abs(value) / 5.0,
                ),
            }
        )
    return targets


def contact_force(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    roller_geoms: list[int],
    band_geom: int,
) -> tuple[float, int]:
    total = 0.0
    count = 0
    force = np.zeros(6)
    for index in range(int(data.ncon)):
        contact = data.contact[index]
        pair = {int(contact.geom1), int(contact.geom2)}
        if band_geom in pair and any(geom in pair for geom in roller_geoms):
            mujoco.mj_contactForce(model, data, index, force)
            total += abs(float(force[0]))
            count += 1
    return total, count


def audit_pose(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    evaluator: ShoulderGeometry,
    config: dict,
    target: dict,
    preload_n: float,
) -> dict:
    shoulder_names = ("sc_elevation", "sc_protraction", "ac_ext")
    joint_ids = {
        name: scalar_joint(model, name) for name in shoulder_names
    }
    qpos_indices = {
        name: int(model.jnt_qposadr[joint])
        for name, joint in joint_ids.items()
    }
    dof_indices = {
        name: int(model.jnt_dofadr[joint])
        for name, joint in joint_ids.items()
    }
    zero_offsets = config.get("zero_pose_offsets_deg", {})
    mujoco.mj_resetData(model, data)
    for name in shoulder_names:
        data.qpos[qpos_indices[name]] = np.deg2rad(
            float(zero_offsets.get(name, 0.0))
            + float(target["q_target_deg"].get(name, 0.0))
        )
    roller_joints = [
        joint
        for joint in range(model.njnt)
        if (
            mujoco.mj_id2name(
                model,
                mujoco.mjtObj.mjOBJ_JOINT,
                joint,
            )
            or ""
        ).startswith("roller_normal")
    ]
    for joint in roller_joints:
        model.jnt_stiffness[joint] = SPRING_STIFFNESS_N_PER_M
        qpos_index = int(model.jnt_qposadr[joint])
        model.qpos_spring[qpos_index] = (
            -preload_n / SPRING_STIFFNESS_N_PER_M
        )
    mujoco.mj_forward(model, data)
    for _ in range(1000):
        mujoco.mj_step(model, data)
        for name in shoulder_names:
            data.qpos[qpos_indices[name]] = np.deg2rad(
                float(zero_offsets.get(name, 0.0))
                + float(target["q_target_deg"].get(name, 0.0))
            )
            data.qvel[dof_indices[name]] = 0.0
    mujoco.mj_forward(model, data)
    marker_ids = [
        mujoco.mj_name2id(
            model,
            mujoco.mjtObj.mjOBJ_SITE,
            f"v7_roller_marker{suffix}",
        )
        for suffix in (("",) if len(roller_joints) == 1 else ("_1", "_2"))
    ]
    band = r400_center(config)
    band_radius = float(config["contact"]["thorax_band"]["radius_mm"])
    gaps = []
    for site_id in marker_ids:
        design = mj_to_design(data.site_xpos[site_id])
        radial = float(
            np.linalg.norm(design[[0, 2]] - band)
            - band_radius
            - ROLLER_RADIUS_MM
        )
        gaps.append(radial)
    roller_geoms = [
        mujoco.mj_name2id(
            model,
            mujoco.mjtObj.mjOBJ_GEOM,
            f"v7_shoulder_roller{suffix}",
        )
        for suffix in (("",) if len(roller_joints) == 1 else ("_1", "_2"))
    ]
    band_geom = evaluator.contact_band_geom
    force, contacts = contact_force(
        model,
        data,
        roller_geoms,
        band_geom,
    )
    roller_q = [
        float(data.qpos[int(model.jnt_qposadr[joint])])
        for joint in roller_joints
    ]
    return {
        "id": target["id"],
        "q_target_deg": {
            name: float(zero_offsets.get(name, 0.0))
            + float(target["q_target_deg"].get(name, 0.0))
            for name in shoulder_names
        },
        "visual_clearance_mm": evaluator.visual_clearance(data),
        "pad_gap_mm": evaluator.pad_gap(data),
        "roller_gap_mm": gaps,
        "roller_contact_force_N": force,
        "roller_contact_count": contacts,
        "roller_joint_position_mm": [value * 1000.0 for value in roller_q],
    }


def evaluate_shoulder_model(
    model: mujoco.MjModel,
    config: dict,
    roller_count: int,
    preload_n: float = 5.0,
) -> dict:
    local_config = copy.deepcopy(config)
    local_config["routes"] = []
    evaluator = ShoulderGeometry(model, local_config)
    data = mujoco.MjData(model)
    rows = [
        audit_pose(
            model,
            data,
            evaluator,
            local_config,
            target,
            preload_n,
        )
        for target in pose_targets(config)
    ]
    max_detach = max(
        max(max(row["roller_gap_mm"]), 0.0) for row in rows
    )
    min_penetration = min(
        min(min(row["roller_gap_mm"]), 0.0) for row in rows
    )
    min_visual = min(row["visual_clearance_mm"] for row in rows)
    max_contact = max(row["roller_contact_force_N"] for row in rows)
    min_contact = min(
        (
            row["roller_contact_force_N"]
            for row in rows
            if row["roller_contact_count"] > 0
        ),
        default=0.0,
    )
    unused = max_contact + min_contact
    del unused
    return {
        "roller_count": roller_count,
        "preload_N": float(preload_n),
        "max_roller_detachment_mm": float(max_detach),
        "min_roller_penetration_mm": float(min_penetration),
        "min_visual_clearance_mm": float(min_visual),
        "max_contact_force_N": float(max_contact),
        "min_nonzero_contact_force_N": float(min_contact),
        "pass": bool(
            min_visual >= -0.5
            and min_penetration >= -0.5
            and max_detach <= 1.0
            and 3.0 <= max_contact <= 10.0
        ),
        "poses": rows,
    }


def load_config() -> dict:
    return json.loads(V6_CONFIG.read_text(encoding="utf-8"))
