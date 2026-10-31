"""Shared interfaces for the cable-design V2 simulation."""

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

import numpy as np


ROOT = Path(__file__).resolve().parent
MUJOCO_ROOT = ROOT.parent
PROJECT_ROOT = MUJOCO_ROOT.parent
CONFIG_PATH = ROOT / "cable_design_v2.json"
ACTION_MATRIX_PATH = ROOT / "action_cable_matrix.json"
MODEL_PATH = ROOT / "arm_cable_design_v2.xml"
EXPORT_ROOT = MUJOCO_ROOT / "exports" / "cable_design_v2"

POINTS_MM = {
    "SC": np.array([40.0, 235.0, -106.0]),
    "AC": np.array([156.86, 233.24, -11.97]),
    "GH": np.array([159.06, 187.54, 7.3]),
    "E": np.array([184.34, -101.35, 7.3]),
    "P": np.array([62.18, 210.10, 91.39]),
    "SUP": np.array([72.76, 268.45, 82.21]),
    "INF": np.array([44.55, 112.87, 106.69]),
}

BODY_ORIGIN_MM = {
    "T": np.zeros(3),
    "C": POINTS_MM["SC"],
    "S": POINTS_MM["AC"],
    "H": POINTS_MM["GH"],
    "F": POINTS_MM["E"],
    "X": np.zeros(3),
}

BODY_NAME = {
    "T": "thorax",
    "C": "clavicle_body",
    "S": "scapula",
    "H": "humerus",
    "F": "ulna",
    "X": "world",
}

GROUP_BY_CABLE = {
    "E01": "elbow",
    "E02": "elbow",
    "G01": "gh",
    "G02": "gh",
    "G03": "gh",
    "G04": "gh",
    "G05": "gh",
    "G06": "gh",
    "S01": "scapula",
    "S02": "scapula",
    "S03": "scapula",
    "S04": "scapula",
    "S05": "scapula",
    "S06": "scapula",
    "P01": "pectoral",
    "P02": "pectoral",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    length = float(np.linalg.norm(values))
    if length < 1e-12:
        raise ValueError("Cannot normalize a zero vector")
    return values / length


def design_to_mj(point_mm: np.ndarray, relative_to_mm: np.ndarray | None = None) -> np.ndarray:
    point = np.asarray(point_mm, dtype=float)
    if relative_to_mm is not None:
        point = point - np.asarray(relative_to_mm, dtype=float)
    return np.array([point[0], -point[2], point[1]], dtype=float) * 0.001


def mj_to_design(point_m: np.ndarray) -> np.ndarray:
    point = np.asarray(point_m, dtype=float)
    return np.array([point[0], point[2], -point[1]], dtype=float) * 1000.0


def vec(values: np.ndarray | list[float]) -> str:
    return " ".join(f"{float(value):.9f}" for value in values)


def quat_from_z_axis(axis: np.ndarray) -> tuple[float, float, float, float]:
    axis = normalize(axis)
    z_axis = np.array([0.0, 0.0, 1.0])
    cross = np.cross(z_axis, axis)
    dot = float(np.clip(np.dot(z_axis, axis), -1.0, 1.0))
    if np.linalg.norm(cross) < 1e-10:
        return (1.0, 0.0, 0.0, 0.0) if dot > 0.0 else (0.0, 1.0, 0.0, 0.0)
    angle = math.acos(dot)
    cross = normalize(cross)
    return (
        math.cos(angle / 2.0),
        float(cross[0] * math.sin(angle / 2.0)),
        float(cross[1] * math.sin(angle / 2.0)),
        float(cross[2] * math.sin(angle / 2.0)),
    )


def action_definitions() -> list[dict]:
    k6_pose = {
        "sc_elevation": 42.26192461124019,
        "sc_protraction": -8.289513955493888,
        "ac_ext": 15.58573377932018,
    }
    return [
        {
            "id": "K1",
            "label": "肘屈曲",
            "group": "elbow",
            "q_target_deg": {"elbow": 116.0},
            "target_deg": 116.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": ["E01", "E02"],
            "metric": "joint_abs",
            "metric_joint": "elbow",
        },
        {
            "id": "K2",
            "label": "GH屈曲",
            "group": "gh",
            "q_target_deg": {"gh_flex": 96.0},
            "target_deg": 96.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": [
                "G01",
                "G02",
                "G03",
                "G04",
                "G05",
                "G06",
                "P01",
                "P02",
            ],
            "metric": "joint_abs",
            "metric_joint": "gh_flex",
        },
        {
            "id": "K3",
            "label": "GH外展",
            "group": "gh",
            "q_target_deg": {"gh_abd": 96.0},
            "target_deg": 96.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": ["G01", "G02", "G03", "G04", "G05", "G06"],
            "metric": "joint_abs",
            "metric_joint": "gh_abd",
        },
        {
            "id": "K4",
            "label": "GH内旋",
            "group": "gh",
            "q_target_deg": {"gh_rot": 56.0},
            "target_deg": 56.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": [
                "G01",
                "G02",
                "G03",
                "G04",
                "G05",
                "G06",
                "P01",
                "P02",
            ],
            "metric": "joint_abs",
            "metric_joint": "gh_rot",
        },
        {
            "id": "K5",
            "label": "GH外旋",
            "group": "gh",
            "q_target_deg": {"gh_rot": -72.0},
            "target_deg": 72.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": [
                "G01",
                "G02",
                "G03",
                "G04",
                "G05",
                "G06",
                "P01",
                "P02",
            ],
            "metric": "joint_abs",
            "metric_joint": "gh_rot",
        },
        {
            "id": "K6",
            "label": "肩胛上回旋",
            "group": "scapula",
            "q_target_deg": k6_pose,
            "target_deg": 30.0,
            "acceptance_deg": 24.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": ["S01", "S02", "S03", "S04", "S05", "S06"],
            "mandatory_cables": ["S01", "S02", "S03", "S04"],
            "minimum_mandatory_tension_N": 2.0,
            "metric": "scapula_up",
        },
        {
            "id": "K7",
            "label": "组合抬升",
            "group": "combined",
            "q_target_deg": {"gh_flex": 96.0, **k6_pose},
            "target_deg": 150.0,
            "acceptance_deg": 120.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": [
                "G01",
                "G02",
                "G03",
                "S01",
                "S02",
                "S03",
                "S04",
                "S05",
                "S06",
                "P01",
                "P02",
            ],
            "mandatory_cables": ["S01", "S02", "S03", "S04", "G01"],
            "minimum_mandatory_tension_N": 2.0,
            "metric": "combined_lift",
        },
        {
            "id": "SLIDE_FWD",
            "label": "肩胛前滑",
            "group": "scapula",
            "q_target_deg": {"sc_protraction": 5.0},
            "target_mm": 15.0,
            "minimum_visible_mm": 5.0,
            "allowed_cables": ["S03", "S04", "S01", "S02"],
            "metric": "scapula_slide",
        },
        {
            "id": "SLIDE_BACK",
            "label": "肩胛后滑",
            "group": "scapula",
            "q_target_deg": {"sc_protraction": -5.0},
            "target_mm": -15.0,
            "minimum_visible_mm": 5.0,
            "allowed_cables": ["S03", "S04", "S01", "S02"],
            "metric": "scapula_slide",
        },
        {
            "id": "SCAP_EXT",
            "label": "肩胛外旋",
            "group": "scapula",
            "q_target_deg": {"ac_ext": 20.0},
            "target_deg": 20.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": ["S01", "S02", "S03", "S04", "S05", "S06"],
            "metric": "scapula_external",
        },
        {
            "id": "SCAP_INT",
            "label": "肩胛内旋",
            "group": "scapula",
            "q_target_deg": {"ac_ext": -20.0},
            "target_deg": 20.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": ["S01", "S02", "S03", "S04", "S05", "S06"],
            "metric": "scapula_external",
        },
        {
            "id": "PECTORAL",
            "label": "胸大肌抱合",
            "group": "pectoral",
            "q_target_deg": {
                "gh_flex": 25.0,
                "gh_abd": 8.0,
                "gh_rot": -20.0,
            },
            "target_deg": 20.0,
            "minimum_visible_deg": 5.0,
            "allowed_cables": ["P01", "P02", "G05", "G06"],
            "metric": "pectoral",
        },
    ]


def guide_geometry_from_bend(
    a_mm: np.ndarray, g_mm: np.ndarray, b_mm: np.ndarray, radius_mm: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    incoming = g_mm - a_mm
    outgoing = b_mm - g_mm
    u1 = normalize(-incoming)
    u2 = normalize(outgoing)
    bend_angle = math.acos(float(np.clip(np.dot(u1, u2), -1.0, 1.0)))
    if bend_angle < 1e-4:
        raise ValueError("Guide bend angle is zero")
    axis = np.cross(incoming, outgoing)
    if np.linalg.norm(axis) < 1e-8:
        axis = np.array([0.6556, 0.0, 0.7552])
    axis = normalize(axis)
    bisector = normalize(u1 + u2)
    centre_offset = radius_mm / math.sin(bend_angle / 2.0)
    centre = g_mm + centre_offset * bisector
    side = -bisector
    wrap_angle = math.pi - bend_angle
    return centre, axis, normalize(side), math.degrees(wrap_angle)
