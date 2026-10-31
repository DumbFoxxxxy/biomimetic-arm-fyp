"""Small binary-STL helpers used by the tendon demo asset pipeline."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def read_binary_stl(path: Path) -> np.ndarray:
    data = path.read_bytes()
    if len(data) < 84:
        raise ValueError(f"{path} is too small to be a binary STL.")
    triangle_count = int(np.frombuffer(data[80:84], dtype="<u4")[0])
    expected = 84 + triangle_count * 50
    if expected > len(data):
        raise ValueError(f"{path} is not a complete binary STL.")
    records = np.frombuffer(
        data[84:expected],
        dtype=np.dtype(
            [
                ("normal", "<f4", (3,)),
                ("vertices", "<f4", (3, 3)),
                ("attribute", "<u2"),
            ]
        ),
    )
    return records["vertices"].astype(float)


def write_binary_stl(path: Path, triangles: np.ndarray, header: str) -> None:
    triangles = np.asarray(triangles, dtype="<f4")
    if triangles.shape[1:] != (3, 3):
        raise ValueError("triangles must have shape (n, 3, 3)")
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0],
        triangles[:, 2] - triangles[:, 0],
    )
    lengths = np.linalg.norm(normals, axis=1)
    lengths[lengths < 1e-12] = 1.0
    normals = normals / lengths[:, None]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        handle.write(header.encode("ascii", errors="replace")[:80].ljust(80, b"\0"))
        handle.write(np.asarray([len(triangles)], dtype="<u4").tobytes())
        for normal, triangle in zip(normals, triangles):
            handle.write(np.asarray(normal, dtype="<f4").tobytes())
            handle.write(np.asarray(triangle, dtype="<f4").tobytes())
            handle.write(np.asarray([0], dtype="<u2").tobytes())
