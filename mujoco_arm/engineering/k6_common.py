"""Shared geometry and configuration helpers for the K6/S01 baseline."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
MUJOCO_ROOT = ROOT.parent
PROJECT_ROOT = MUJOCO_ROOT.parent
DEFAULT_CONFIG = ROOT / "k6_baseline.json"
DEFAULT_ASSETS = ROOT / "assets"
DEFAULT_XML = ROOT / "k6_baseline.xml"
DEFAULT_REPORT = MUJOCO_ROOT / "exports" / "k6_baseline" / "report.json"


def load_config(path: Path = DEFAULT_CONFIG) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def to_mujoco_design(point_mm: np.ndarray) -> np.ndarray:
    """Map design x-right/y-up/z-back to MuJoCo x/y/z-up, in metres."""

    point_mm = np.asarray(point_mm, dtype=float)
    return np.array(
        [point_mm[0], -point_mm[2], point_mm[1]], dtype=float
    ) * 0.001


def vec(values: np.ndarray | list[float]) -> str:
    return " ".join(f"{float(value):.9f}" for value in values)


def normalize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    length = float(np.linalg.norm(values))
    if length < 1e-12:
        raise ValueError("Cannot normalize a zero vector")
    return values / length


def quat_from_z_axis(axis: np.ndarray) -> tuple[float, float, float, float]:
    axis = normalize(axis)
    z_axis = np.array([0.0, 0.0, 1.0])
    cross = np.cross(z_axis, axis)
    dot = float(np.clip(np.dot(z_axis, axis), -1.0, 1.0))
    if np.linalg.norm(cross) < 1e-10:
        if dot > 0.0:
            return (1.0, 0.0, 0.0, 0.0)
        return (0.0, 1.0, 0.0, 0.0)
    angle = math.acos(dot)
    cross = normalize(cross)
    return (
        math.cos(angle / 2.0),
        float(cross[0] * math.sin(angle / 2.0)),
        float(cross[1] * math.sin(angle / 2.0)),
        float(cross[2] * math.sin(angle / 2.0)),
    )


def r400_center(config: dict) -> np.ndarray:
    """Return the inner R400 arc centre in the design x-z plane."""

    band = config["contact"]["thorax_band"]
    x0 = float(band["x_min_mm"])
    x1 = float(band["x_max_mm"])
    semi_x = float(config["thorax"]["semi_x_mm"])
    semi_z = float(config["thorax"]["semi_z_mm"])
    radius = float(band["radius_mm"])

    def ellipse_z(x: float) -> float:
        return semi_z * math.sqrt(max(0.0, 1.0 - (x / semi_x) ** 2))

    p0 = np.array([x0, ellipse_z(x0)])
    p1 = np.array([x1, ellipse_z(x1)])
    midpoint = 0.5 * (p0 + p1)
    chord = p1 - p0
    chord_length = float(np.linalg.norm(chord))
    perpendicular = np.array([-chord[1], chord[0]]) / chord_length
    height = math.sqrt(radius**2 - (chord_length / 2.0) ** 2)
    candidate_a = midpoint + height * perpendicular
    candidate_b = midpoint - height * perpendicular
    return candidate_a if np.linalg.norm(candidate_a) < np.linalg.norm(candidate_b) else candidate_b


def band_surface(x_mm: np.ndarray, config: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return R400 surface points and outward x-z normals for x coordinates."""

    x_mm = np.asarray(x_mm, dtype=float)
    centre = r400_center(config)
    radius = float(config["contact"]["thorax_band"]["radius_mm"])
    z_mm = centre[1] + np.sqrt(
        np.maximum(0.0, radius**2 - (x_mm - centre[0]) ** 2)
    )
    normals = np.column_stack((x_mm - centre[0], z_mm - centre[1]))
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    return np.column_stack((x_mm, z_mm)), normals


def surface_point(
    x_mm: float, y_mm: float, outward_offset_mm: float, config: dict
) -> np.ndarray:
    surface_xz, normal_xz = band_surface(np.asarray([x_mm]), config)
    point = np.asarray([surface_xz[0, 0], y_mm, surface_xz[0, 1]], dtype=float)
    normal = np.asarray([normal_xz[0, 0], 0.0, normal_xz[0, 1]], dtype=float)
    return point + outward_offset_mm * normal


def scapula_basis(config: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    points = config["points_mm"]
    p = np.asarray(points["P"], dtype=float)
    gh = np.asarray(points["GH"], dtype=float)
    sup = np.asarray(points["SUP"], dtype=float)
    inf = np.asarray(points["INF"], dtype=float)
    u_axis = normalize(sup - inf)
    lateral = gh - p
    v_axis = normalize(lateral - np.dot(lateral, u_axis) * u_axis)
    w_axis = normalize(np.asarray(config["scapula"]["normal"], dtype=float))
    if abs(np.dot(w_axis, u_axis)) > 1e-3 or abs(np.dot(w_axis, v_axis)) > 1e-3:
        raise ValueError("Scapula W axis is not perpendicular to U and V")
    return u_axis, v_axis, w_axis


def scapula_local_to_design(
    local_mm: np.ndarray, config: dict
) -> np.ndarray:
    local_mm = np.asarray(local_mm, dtype=float)
    p = np.asarray(config["points_mm"]["P"], dtype=float)
    u_axis, v_axis, w_axis = scapula_basis(config)
    return (
        p
        + local_mm[..., 0, None] * u_axis
        + local_mm[..., 1, None] * v_axis
        + local_mm[..., 2, None] * w_axis
    )


def resolve_s01_route(config: dict, route: dict | None = None) -> dict:
    source = route or config["s01"].get("selected") or config["s01"]["initial"]
    a = surface_point(
        float(source["A"]["x_mm"]),
        float(source["A"]["y_mm"]),
        float(source["A"]["outward_offset_mm"]),
        config,
    )
    g = surface_point(
        float(source["G"]["x_mm"]),
        float(source["G"]["y_mm"]),
        float(source["G"]["outward_offset_mm"]),
        config,
    )
    local_b = np.asarray(
        [
            source["B"]["u_mm"],
            source["B"]["v_mm"],
            source["B"]["w_mm"],
        ],
        dtype=float,
    )
    b = scapula_local_to_design(local_b, config)
    return {
        "A_mm": a.tolist(),
        "G_mm": g.tolist(),
        "B_mm": b.tolist(),
        "B_local_uvw_mm": local_b.tolist(),
        "sidesite_sign": float(
            source.get("sidesite_sign", config["s01"]["sidesite_sign"])
        ),
        "input": source,
    }


def rounded_rectangle_polygon(
    half_short_mm: float,
    half_long_mm: float,
    radius_mm: float,
    arc_segments: int = 4,
) -> np.ndarray:
    """Return a counter-clockwise rounded rectangle in (short, long) axes."""

    if radius_mm <= 0.0:
        raise ValueError("Corner radius must be positive")
    if radius_mm > min(half_short_mm, half_long_mm):
        raise ValueError("Corner radius is too large")
    sx = half_short_mm - radius_mm
    sy = half_long_mm - radius_mm
    centres = (
        (sx, sy, 0.0),
        (-sx, sy, 90.0),
        (-sx, -sy, 180.0),
        (sx, -sy, 270.0),
    )
    points = []
    for cx, cy, start_deg in centres:
        for index in range(arc_segments + 1):
            angle = math.radians(start_deg + 90.0 * index / arc_segments)
            points.append(
                (cx + radius_mm * math.cos(angle), cy + radius_mm * math.sin(angle))
            )
    return np.asarray(points, dtype=float)
