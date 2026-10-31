"""Build and evaluate the isolated V7 elbow sector-spool model."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import json
import math
from project_runtime import WorkspacePath as Path
import sys
import xml.etree.ElementTree as ET

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v7"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from evaluate_v2 import dense_tendon_jacobian  # noqa: E402
from v2_common import quat_from_z_axis, vec  # noqa: E402


DEFAULT_PARAMETERS = {
    "radius_mm": 18.0,
    "entry_radius_mm": 10.0,
    "entry_width_mm": 12.0,
    "humerus_length_mm": 290.0,
    "ulna_length_mm": 237.0,
    "bone_radius_mm": 12.0,
    "elbow_range_deg": [0.0, 116.0],
    "external_anchor_mm": [0.0, 260.0, 65.0],
    "entry_guide_mm": [0.0, 95.0, 105.0],
    "anchor_angle_deg": 0.0,
    "entry_side_angle_deg": 0.0,
    "sector_side_angle_deg": 45.0,
    "timestep_s": 0.0005,
    "density_kg_m3": 1100.0,
}


def design_to_mj(point_mm: list[float] | np.ndarray) -> np.ndarray:
    point = np.asarray(point_mm, dtype=float)
    return np.array([point[0], -point[2], point[1]]) * 0.001


def radial_point(radius_m, angle_deg: float) -> np.ndarray:
    angle = math.radians(angle_deg)
    return np.array([0.0, radius_m * math.cos(angle), radius_m * math.sin(angle)])


def cylinder_x_quat() -> tuple[float, float, float, float]:
    return quat_from_z_axis(np.array([1.0, 0.0, 0.0]))


def mirror_design_y(point_mm: list[float]) -> list[float]:
    return [float(point_mm[0]), -float(point_mm[1]), float(point_mm[2])]


def build_elbow_model(parameters: dict | None = None) -> ET.ElementTree:
    params = {**DEFAULT_PARAMETERS, **(parameters or {})}
    radius = float(params["radius_mm"]) * 0.001
    entry_radius = float(params["entry_radius_mm"]) * 0.001
    humerus_length = float(params["humerus_length_mm"]) * 0.001
    ulna_length = float(params["ulna_length_mm"]) * 0.001
    bone_radius = float(params["bone_radius_mm"]) * 0.001
    quat = cylinder_x_quat()

    root = ET.Element("mujoco", {"model": "v7_elbow_standalone"})
    ET.SubElement(
        root,
        "compiler",
        {"angle": "degree", "autolimits": "true"},
    )
    ET.SubElement(
        root,
        "option",
        {
            "timestep": str(params["timestep_s"]),
            "gravity": "0 0 -9.81",
            "integrator": "implicitfast",
        },
    )
    visual = ET.SubElement(root, "visual")
    ET.SubElement(visual, "global", {"azimuth": "45", "elevation": "-18"})
    asset = ET.SubElement(root, "asset")
    ET.SubElement(
        asset,
        "material",
        {"name": "bone_mat", "rgba": "0.82 0.80 0.72 1"},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "sector_mat", "rgba": "0.20 0.62 0.84 1"},
    )
    ET.SubElement(
        asset,
        "material",
        {"name": "guide_mat", "rgba": "0.92 0.35 0.18 1"},
    )
    world = ET.SubElement(root, "worldbody")
    ET.SubElement(
        world,
        "geom",
        {
            "name": "floor",
            "type": "plane",
            "size": "1 1 0.02",
            "pos": "0 0 -0.30",
            "rgba": "0.3 0.32 0.34 0.2",
        },
    )
    humerus = ET.SubElement(world, "body", {"name": "humerus", "pos": "0 0 0"})
    ET.SubElement(
        humerus,
        "geom",
        {
            "name": "humerus_link",
            "type": "cylinder",
            "size": f"{bone_radius:.6f} {humerus_length * 0.5:.6f}",
            "pos": f"0 0 {humerus_length * 0.5:.6f}",
            "material": "bone_mat",
            "density": str(params["density_kg_m3"]),
        },
    )
    entry_centres = {
        "E01": design_to_mj(
            params.get(
                "entry_guide_E01_mm",
                mirror_design_y(params["entry_guide_mm"]),
            )
        ),
        "E02": design_to_mj(
            params.get(
                "entry_guide_E02_mm",
                params["entry_guide_mm"],
            )
        ),
    }
    external_anchors = {
        "E01": design_to_mj(
            params.get(
                "external_anchor_E01_mm",
                mirror_design_y(params["external_anchor_mm"]),
            )
        ),
        "E02": design_to_mj(
            params.get(
                "external_anchor_E02_mm",
                params["external_anchor_mm"],
            )
        ),
    }
    for route_id, centre in entry_centres.items():
        ET.SubElement(
            humerus,
            "geom",
            {
                "name": f"entry_{route_id}",
                "type": "cylinder",
                "size": f"{entry_radius:.6f} {float(params['entry_width_mm']) * 0.0005:.6f}",
                "pos": vec(centre),
                "quat": vec(np.asarray(quat)),
                "material": "guide_mat",
                "contype": "0",
                "conaffinity": "0",
                "density": "0",
            },
        )
        side_angle = float(
            params.get(
                f"entry_side_angle_{route_id}_deg",
                (
                    float(params["entry_side_angle_deg"])
                    if route_id == "E01"
                    else 180.0
                ),
            )
        )
        side = radial_point(entry_radius * 1.6, side_angle)
        ET.SubElement(
            humerus,
            "site",
            {
                "name": f"entry_side_{route_id}",
                "pos": vec(centre + side),
                "size": "0.0015",
                "rgba": "0.95 0.2 0.2 0.7",
            },
        )
        outgoing = external_anchors[route_id] - centre
        outgoing /= max(float(np.linalg.norm(outgoing)), 1e-9)
        ET.SubElement(
            humerus,
            "site",
            {
                "name": f"bridge_{route_id}",
                "pos": vec(centre + outgoing * 0.025),
                "size": "0.0015",
                "rgba": "0.25 0.8 0.9 0.7",
            },
        )
    ulna = ET.SubElement(humerus, "body", {"name": "ulna", "pos": "0 0 0"})
    ET.SubElement(
        ulna,
        "joint",
        {
            "name": "elbow",
            "type": "hinge",
            "axis": "1 0 0",
            "range": vec(np.asarray(params["elbow_range_deg"])),
            "damping": "0.2",
            "armature": "0.005",
        },
    )
    ET.SubElement(
        ulna,
        "geom",
        {
            "name": "ulna_link",
            "type": "cylinder",
            "size": f"{bone_radius * 0.82:.6f} {ulna_length * 0.5:.6f}",
            "pos": f"0 0 {-ulna_length * 0.5:.6f}",
            "material": "bone_mat",
            "density": str(params["density_kg_m3"]),
        },
    )
    anchor_angle = float(params["anchor_angle_deg"])
    sector_side_angle = float(params["sector_side_angle_deg"])
    for route_id in ("E01", "E02"):
        ET.SubElement(
            ulna,
            "geom",
            {
                "name": f"sector_{route_id}",
                "type": "cylinder",
                "size": f"{radius:.6f} 0.007",
                "pos": "0 0 0",
                "quat": vec(np.asarray(quat)),
                "material": "sector_mat",
                "contype": "0",
                "conaffinity": "0",
                "density": "0",
            },
        )
        angle = float(
            params.get(
                f"anchor_angle_{route_id}_deg",
                (
                    anchor_angle
                    if route_id == "E01"
                    else 180.0
                ),
            )
        )
        anchor = radial_point(radius * 1.02, angle)
        ET.SubElement(
            ulna,
            "site",
            {
                "name": f"anchor_{route_id}",
                "pos": vec(anchor),
                "size": "0.0025",
                "rgba": "0.98 0.62 0.15 1",
            },
        )
        side_angle = float(
            params.get(
                f"sector_side_angle_{route_id}_deg",
                (
                    sector_side_angle
                    if route_id == "E01"
                    else 135.0
                ),
            )
        )
        side = radial_point(radius * 1.55, side_angle)
        ET.SubElement(
            ulna,
            "site",
            {
                "name": f"sector_side_{route_id}",
                "pos": vec(side),
                "size": "0.0015",
                "rgba": "0.95 0.2 0.2 0.7",
            },
        )
    for route_id, anchor in external_anchors.items():
        ET.SubElement(
            world,
            "site",
            {
                "name": f"external_{route_id}",
                "pos": vec(anchor),
                "size": "0.004",
                "rgba": "0.1 0.8 0.5 1",
            },
        )
    tendon_root = ET.SubElement(root, "tendon")
    for route_id in ("E01", "E02"):
        tendon = ET.SubElement(
            tendon_root,
            "spatial",
            {
                "name": f"cable_{route_id}",
                "width": "0.0018",
                "rgba": "0.95 0.45 0.12 0.95",
            },
        )
        ET.SubElement(tendon, "site", {"site": f"external_{route_id}"})
        ET.SubElement(
            tendon,
            "geom",
            {
                "geom": f"entry_{route_id}",
                "sidesite": f"entry_side_{route_id}",
            },
        )
        ET.SubElement(tendon, "site", {"site": f"bridge_{route_id}"})
        ET.SubElement(
            tendon,
            "geom",
            {
                "geom": f"sector_{route_id}",
                "sidesite": f"sector_side_{route_id}",
            },
        )
        ET.SubElement(tendon, "site", {"site": f"anchor_{route_id}"})
    actuator_root = ET.SubElement(root, "actuator")
    for route_id in ("E01", "E02"):
        ET.SubElement(
            actuator_root,
            "position",
            {
                "name": f"pull_{route_id}",
                "tendon": f"cable_{route_id}",
                "kp": "1200",
                "kv": "80",
                "ctrllimited": "true",
                "ctrlrange": "0 1",
                "forcelimited": "true",
                "forcerange": "-300 0",
            },
        )
    sensor_root = ET.SubElement(root, "sensor")
    for route_id in ("E01", "E02"):
        ET.SubElement(
            sensor_root,
            "tendonpos",
            {"name": f"length_{route_id}", "tendon": f"cable_{route_id}"},
        )
    ET.SubElement(sensor_root, "jointpos", {"name": "elbow_deg", "joint": "elbow"})
    ET.indent(root, space="  ")
    return ET.ElementTree(root)


def evaluate_elbow(
    model: mujoco.MjModel,
    radius_mm: float,
    samples: int = 21,
) -> dict:
    data = mujoco.MjData(model)
    joint = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "elbow")
    dof = int(model.jnt_dofadr[joint])
    qpos_index = int(model.jnt_qposadr[joint])
    tendon_ids = {
        route_id: mujoco.mj_name2id(
            model,
            mujoco.mjtObj.mjOBJ_TENDON,
            f"cable_{route_id}",
        )
        for route_id in ("E01", "E02")
    }
    neutral_lengths = {}
    data.qpos[qpos_index] = 0.0
    mujoco.mj_forward(model, data)
    for route_id, tendon_id in tendon_ids.items():
        neutral_lengths[route_id] = float(data.ten_length[tendon_id])
    rows = []
    for fraction in np.linspace(0.0, 1.0, samples):
        angle = math.radians(116.0 * float(fraction))
        data.qpos[qpos_index] = angle
        data.qvel[:] = 0.0
        mujoco.mj_forward(model, data)
        jacobian = dense_tendon_jacobian(model, data)
        columns = [int(model.tendon_adr[tendon_ids[route_id]]) for route_id in ()]
        del columns
        moments = []
        wrap = []
        lengths = []
        for route_id in ("E01", "E02"):
            tendon_id = tendon_ids[route_id]
            moments.append(-float(jacobian[tendon_id, dof]))
            wrap.append(int(data.ten_wrapnum[tendon_id]))
            lengths.append(float(data.ten_length[tendon_id]))
        rows.append(
            {
                "fraction": float(fraction),
                "angle_deg": float(np.degrees(angle)),
                "moment_m_per_n": moments,
                "effective_moment_arm_mm": [
                    abs(moments[0]) * 1000.0,
                    abs(moments[1]) * 1000.0,
                ],
                "wrap_points": wrap,
                "length_mm": [length * 1000.0 for length in lengths],
            }
        )
    opposite = all(
        row["moment_m_per_n"][0] * row["moment_m_per_n"][1] < 0.0
        for row in rows
    )
    min_arm = min(
        min(row["effective_moment_arm_mm"]) for row in rows
    )
    final_lengths = rows[-1]["length_mm"]
    mean_moment_e01 = float(
        np.mean([row["moment_m_per_n"][0] for row in rows])
    )
    mean_moment_e02 = float(
        np.mean([row["moment_m_per_n"][1] for row in rows])
    )
    if mean_moment_e01 > 0.0:
        contraction_e01 = (
            neutral_lengths["E01"] * 1000.0 - final_lengths[0]
        ) / (neutral_lengths["E01"] * 1000.0) * 100.0
    else:
        contraction_e01 = (
            final_lengths[0] - neutral_lengths["E01"] * 1000.0
        ) / final_lengths[0] * 100.0
    if mean_moment_e02 > 0.0:
        contraction_e02 = (
            neutral_lengths["E02"] * 1000.0 - final_lengths[1]
        ) / (neutral_lengths["E02"] * 1000.0) * 100.0
    else:
        contraction_e02 = (
            final_lengths[1] - neutral_lengths["E02"] * 1000.0
        ) / final_lengths[1] * 100.0
    wrap_ok = all(
        all(value >= 4 for value in row["wrap_points"])
        for row in rows
    )
    return {
        "radius_mm": float(radius_mm),
        "opposite": bool(opposite),
        "wrap_ok": bool(wrap_ok),
        "minimum_arm_mm": float(min_arm),
        "contraction_percent": {
            "E01": float(contraction_e01),
            "E02": float(contraction_e02),
            "max": float(max(contraction_e01, contraction_e02)),
        },
        "pass": bool(
            opposite
            and wrap_ok
            and min_arm >= 8.0
            and max(contraction_e01, contraction_e02) <= 18.0
        ),
        "samples": rows,
    }


def load_parameters(path: Path | None = None) -> dict:
    if path is None:
        return dict(DEFAULT_PARAMETERS)
    return json.loads(path.read_text(encoding="utf-8"))
