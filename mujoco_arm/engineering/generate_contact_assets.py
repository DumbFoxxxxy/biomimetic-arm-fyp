"""Generate the analytic R400 band and scapula contact pad as binary STLs."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

from k6_common import (
    DEFAULT_ASSETS,
    DEFAULT_CONFIG,
    band_surface,
    load_config,
    rounded_rectangle_polygon,
    scapula_basis,
    scapula_local_to_design,
    to_mujoco_design,
)


MUJOCO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MUJOCO_ROOT))
from stl_io import write_binary_stl  # noqa: E402


def _grid_triangles(rows: int, columns: int) -> list[tuple[int, int, int]]:
    triangles = []
    for row in range(rows - 1):
        for column in range(columns - 1):
            a = row * columns + column
            b = a + 1
            c = a + columns
            d = c + 1
            triangles.append((a, b, d))
            triangles.append((a, d, c))
    return triangles


def generate_thorax_band(config: dict) -> tuple[np.ndarray, dict]:
    band = config["contact"]["thorax_band"]
    x_values = np.linspace(
        float(band["x_min_mm"]),
        float(band["x_max_mm"]),
        int(band["x_segments"]) + 1,
    )
    y_values = np.linspace(
        float(band["y_min_mm"]),
        float(band["y_max_mm"]),
        int(band["y_segments"]) + 1,
    )
    surface_xz, normals = band_surface(x_values, config)
    radius = float(band["radius_mm"])
    thickness = float(band["thickness_mm"])
    centre = surface_xz - radius * normals
    inner_xz = centre + (radius - thickness) * normals

    outer_design = []
    inner_design = []
    for y_value in y_values:
        for point, inner in zip(surface_xz, inner_xz):
            outer_design.append([point[0], y_value, point[1]])
            inner_design.append([inner[0], y_value, inner[1]])
    outer_design = np.asarray(outer_design)
    inner_design = np.asarray(inner_design)

    outer_mj = np.asarray([to_mujoco_design(point) for point in outer_design])
    inner_mj = np.asarray([to_mujoco_design(point) for point in inner_design])
    rows = len(y_values)
    columns = len(x_values)
    offset = rows * columns

    triangles = []
    for tri in _grid_triangles(rows, columns):
        triangles.append(tuple(outer_mj[index] for index in tri))
    for tri in _grid_triangles(rows, columns):
        triangles.append(tuple(inner_mj[index] for index in reversed(tri)))

    def quad(a: int, b: int, c: int, d: int) -> None:
        triangles.append((outer_mj[a], outer_mj[b], inner_mj[b - offset]))
        triangles.append((outer_mj[a], inner_mj[b - offset], inner_mj[a - offset]))

    for row in range(rows - 1):
        left_outer = row * columns
        right_outer = left_outer + columns - 1
        triangles.append((outer_mj[left_outer], inner_mj[left_outer], inner_mj[(row + 1) * columns]))
        triangles.append((outer_mj[left_outer], inner_mj[(row + 1) * columns], outer_mj[(row + 1) * columns]))
        triangles.append((outer_mj[right_outer], outer_mj[(row + 1) * columns + columns - 1], inner_mj[(row + 1) * columns + columns - 1]))
        triangles.append((outer_mj[right_outer], inner_mj[(row + 1) * columns + columns - 1], inner_mj[right_outer]))

    for column in range(columns - 1):
        bottom_outer = column
        bottom_next_outer = column + 1
        top_outer = (rows - 1) * columns + column
        top_next_outer = top_outer + 1
        triangles.append((outer_mj[bottom_outer], inner_mj[bottom_outer], inner_mj[bottom_next_outer]))
        triangles.append((outer_mj[bottom_outer], inner_mj[bottom_next_outer], outer_mj[bottom_next_outer]))
        triangles.append((outer_mj[top_outer], outer_mj[top_next_outer], inner_mj[top_next_outer]))
        triangles.append((outer_mj[top_outer], inner_mj[top_next_outer], inner_mj[top_outer]))

    metadata = {
        "surface_radius_mm": radius,
        "thickness_mm": thickness,
        "x_range_mm": [float(x_values[0]), float(x_values[-1])],
        "y_range_mm": [float(y_values[0]), float(y_values[-1])],
        "centre_xz_mm": centre[0].tolist(),
    }
    return np.asarray(triangles), metadata


def generate_scapula_pad(config: dict) -> tuple[np.ndarray, dict]:
    pad = config["contact"]["scapula_pad"]
    half_short = float(pad["short_mm"]) / 2.0
    half_long = float(pad["long_mm"]) / 2.0
    polygon = rounded_rectangle_polygon(
        half_short,
        half_long,
        float(pad["corner_radius_mm"]),
        int(pad["corner_segments"]),
    )
    gap = float(pad["zero_gap_mm"])
    thickness = float(pad["thickness_mm"])
    bottom_local = np.column_stack(
        (
            polygon[:, 1],
            polygon[:, 0],
            np.full(len(polygon), gap),
        )
    )
    top_local = np.column_stack(
        (
            polygon[:, 1],
            polygon[:, 0],
            np.full(len(polygon), gap + thickness),
        )
    )
    bottom_design = scapula_local_to_design(bottom_local, config)
    top_design = scapula_local_to_design(top_local, config)
    ac = np.asarray(config["points_mm"]["AC"], dtype=float)
    bottom_mj = np.asarray([to_mujoco_design(point - ac) for point in bottom_design])
    top_mj = np.asarray([to_mujoco_design(point - ac) for point in top_design])
    centre_bottom = bottom_mj.mean(axis=0)
    centre_top = top_mj.mean(axis=0)

    triangles = []
    count = len(polygon)
    for index in range(count):
        nxt = (index + 1) % count
        triangles.append((centre_bottom, bottom_mj[nxt], bottom_mj[index]))
        triangles.append((centre_top, top_mj[index], top_mj[nxt]))
        triangles.append((bottom_mj[index], bottom_mj[nxt], top_mj[nxt]))
        triangles.append((bottom_mj[index], top_mj[nxt], top_mj[index]))

    metadata = {
        "short_mm": float(pad["short_mm"]),
        "long_mm": float(pad["long_mm"]),
        "thickness_mm": thickness,
        "corner_radius_mm": float(pad["corner_radius_mm"]),
        "zero_gap_mm": gap,
        "boundary_points_local_uv_mm": polygon.tolist(),
    }
    return np.asarray(triangles), metadata


def main() -> None:
    config = load_config(DEFAULT_CONFIG)
    DEFAULT_ASSETS.mkdir(parents=True, exist_ok=True)
    band_triangles, band_metadata = generate_thorax_band(config)
    pad_triangles, pad_metadata = generate_scapula_pad(config)
    band_path = DEFAULT_ASSETS / "thorax_contact_band.stl"
    pad_path = DEFAULT_ASSETS / "scapula_contact_pad.stl"
    write_binary_stl(band_path, band_triangles, "K6 analytic R400 thorax band")
    write_binary_stl(pad_path, pad_triangles, "K6 analytic scapula contact pad")

    u_axis, v_axis, w_axis = scapula_basis(config)
    metadata = {
        "thorax_band": {
            **band_metadata,
            "file": band_path.name,
            "triangle_count": len(band_triangles),
        },
        "scapula_pad": {
            **pad_metadata,
            "file": pad_path.name,
            "triangle_count": len(pad_triangles),
            "basis_uvw": {
                "U": u_axis.tolist(),
                "V": v_axis.tolist(),
                "W": w_axis.tolist(),
            },
        },
    }
    metadata_path = DEFAULT_ASSETS / "contact_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {band_path} ({len(band_triangles)} triangles)")
    print(f"wrote {pad_path} ({len(pad_triangles)} triangles)")
    print(f"wrote {metadata_path}")


if __name__ == "__main__":
    main()
