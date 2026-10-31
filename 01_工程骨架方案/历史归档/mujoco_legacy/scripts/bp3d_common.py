"""Shared BodyParts3D access, geometry helpers and the muscle fusion table.

Every BodyParts3D element is one organ mesh in a single body-wide coordinate
frame measured in millimetres.  The right-hand limb elements are the ones
whose file id does not end in ``M``; the ``M`` variants are the mirrored left
side.

Nothing in this module touches MuJoCo.  It is the geometric foundation that
asset collection, the muscle manifest, the ligament manifest and the ROM
derivation all share, so that a joint centre or a muscle centre-line is
computed exactly once.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

BP3D_DIR = Path(r"D:\XJTLU\isa_BP3D_4.0_obj_99")
BP3D_INDEX = BP3D_DIR / "isa_element_parts.txt"

ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "bp3d"
CONFIG_DIR = ROOT / "config"
REPORT_DIR = ROOT / "exports" / "bp3d"

FRAMES_PATH = ASSET_DIR / "bp3d_frames.json"
MUSCLE_MANIFEST_PATH = CONFIG_DIR / "bp3d_muscle_manifest.json"
LIGAMENT_MANIFEST_PATH = CONFIG_DIR / "bp3d_ligament_manifest.json"
ROM_PATH = CONFIG_DIR / "bp3d_joint_rom.json"
MODEL_PATH = ROOT / "arm_bp3d.xml"
CONTRACT_PATH = CONFIG_DIR / "bp3d_arm.json"

MM = 0.001

# Anatomical axis convention of the BodyParts3D standing pose.  Verified from
# the element bounding boxes: +Z is superior, -Y is anterior, -X is the right
# side of the body.
SUPERIOR = np.array([0.0, 0.0, 1.0])
ANTERIOR = np.array([0.0, -1.0, 0.0])
RIGHT = np.array([-1.0, 0.0, 0.0])

BONE_SPECS = (
    {
        "name": "thorax",
        "patterns": (
            r"^right \w+ rib$",
            r"^manubrium$",
            r"^body of sternum$",
            r"^\w+ cervical vertebra$",
            r"^\w+ thoracic vertebra$",
        ),
        "density": 1500.0,
        "role": "base",
    },
    {"name": "clavicle", "patterns": (r"^right clavicle$",), "density": 1900.0},
    {"name": "scapula", "patterns": (r"^right scapula$",), "density": 1900.0},
    {"name": "humerus", "patterns": (r"^right humerus$",), "density": 1900.0},
    {"name": "radius", "patterns": (r"^right radius$",), "density": 1900.0},
    {"name": "ulna", "patterns": (r"^right ulna$",), "density": 1900.0},
)

BONE_NAMES = tuple(spec["name"] for spec in BONE_SPECS)

# Fused drive channels.  Each channel owns one to four BodyParts3D muscle
# elements; the elements are merged into a single wire path because the wire
# runs down the centre of the printed X-lattice muscle group.  The lattice
# parameters come from the Yoshimura 2026 parametric framework: low-density
# cells for the deformable belly and high-density cells for the stiff ends
# that behave like tendon.
LATTICE_MUSCLE = {"x_mm": 14.0, "y_mm": 12.0, "z_mm": 30.0, "p_mm": 2.0}
LATTICE_TENDON = {"x_mm": 6.0, "y_mm": 8.0, "z_mm": 25.0, "p_mm": 2.0}

CHANNELS = (
    {
        "id": "supraspinatus",
        "group": "rotator_cuff",
        "elements": (r"^supraspinatus$",),
        "origin_bone": "scapula",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_abduction",
        "channels": 1,
        "max_force_n": 220.0,
        "travel_m": 0.06,
    },
    {
        "id": "infraspinatus_teres_minor",
        "group": "rotator_cuff",
        "elements": (r"^infraspinatus$", r"^teres minor$"),
        "origin_bone": "scapula",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_external_rotation",
        "channels": 1,
        "max_force_n": 260.0,
        "travel_m": 0.06,
    },
    {
        "id": "subscapularis",
        "group": "rotator_cuff",
        "elements": (r"^subscapularis$",),
        "origin_bone": "scapula",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_internal_rotation",
        "channels": 1,
        "max_force_n": 300.0,
        "travel_m": 0.06,
    },
    {
        "id": "teres_major",
        "group": "adductor",
        "elements": (r"^teres major$",),
        "origin_bone": "scapula",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_adduction",
        "channels": 1,
        "max_force_n": 200.0,
        "travel_m": 0.06,
    },
    {
        "id": "deltoid_anterior",
        "group": "deltoid",
        "elements": (r"^clavicular part of deltoid$", r"^clavicular part of pectoralis major$"),
        "origin_bone": "clavicle",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_flexion",
        "channels": 2,
        "max_force_n": 320.0,
        "travel_m": 0.08,
    },
    {
        "id": "deltoid_middle",
        "group": "deltoid",
        "elements": (r"^acromial part of deltoid$",),
        "origin_bone": "scapula",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_abduction",
        "channels": 2,
        "max_force_n": 380.0,
        "travel_m": 0.08,
    },
    {
        "id": "deltoid_posterior",
        "group": "deltoid",
        "elements": (r"^spinal part of deltoid$",),
        "origin_bone": "scapula",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_extension",
        "channels": 1,
        "max_force_n": 300.0,
        "travel_m": 0.08,
    },
    {
        "id": "pectoralis_major",
        "group": "adductor",
        "elements": (r"^sternocostal part of pectoralis major$", r"^abdominal part of pectoralis major$"),
        "origin_bone": "thorax",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_adduction",
        "channels": 2,
        "max_force_n": 380.0,
        "travel_m": 0.08,
    },
    {
        "id": "coracobrachialis",
        "group": "adductor",
        "elements": (r"^coracobrachialis$",),
        "origin_bone": "scapula",
        "insertion_bone": "humerus",
        "target_joint": "glenohumeral_flexion",
        "channels": 1,
        "max_force_n": 160.0,
        "travel_m": 0.06,
    },
    {
        "id": "trapezius_levator",
        "group": "scapulothoracic",
        "elements": (
            r"^descending part of trapezius$",
            r"^transverse part of trapezius$",
            r"^levator scapulae$",
        ),
        "origin_bone": "thorax",
        "insertion_bone": "scapula",
        "target_joint": "scapulothoracic_upward_rotation",
        "channels": 2,
        "max_force_n": 320.0,
        "travel_m": 0.10,
    },
    {
        "id": "rhomboid_lower_trap",
        "group": "scapulothoracic",
        "elements": (
            r"^ascending part of trapezius$",
            r"^rhomboid major$",
            r"^rhomboid minor$",
        ),
        "origin_bone": "thorax",
        "insertion_bone": "scapula",
        "target_joint": "scapulothoracic_downward_rotation",
        "channels": 2,
        "max_force_n": 300.0,
        "travel_m": 0.10,
    },
    {
        "id": "serratus_anterior",
        "group": "scapulothoracic",
        "elements": (r"^serratus anterior$",),
        "origin_bone": "thorax",
        "insertion_bone": "scapula",
        "target_joint": "scapulothoracic_upward_rotation",
        "channels": 3,
        "max_force_n": 400.0,
        "travel_m": 0.10,
    },
    {
        "id": "elbow_flexors",
        "group": "elbow",
        "elements": (
            r"^long head of biceps brachii$",
            r"^short head of biceps brachii$",
            r"^brachialis$",
        ),
        "origin_bone": "humerus",
        "insertion_bone": "radius",
        "target_joint": "elbow_flexion",
        "channels": 2,
        "max_force_n": 650.0,
        "travel_m": 0.12,
    },
    {
        "id": "triceps",
        "group": "elbow",
        "elements": (
            r"^long head of triceps brachii$",
            r"^lateral head of triceps brachii$",
            r"^medial head of triceps brachii$",
        ),
        "origin_bone": "humerus",
        "insertion_bone": "ulna",
        "target_joint": "elbow_extension",
        "channels": 2,
        "max_force_n": 700.0,
        "travel_m": 0.12,
    },
)

OPTIONAL_CHANNELS = (
    {
        "id": "subclavius",
        "group": "clavicle",
        "elements": (r"^subclavius$",),
        "origin_bone": "thorax",
        "insertion_bone": "clavicle",
        "target_joint": "sternoclavicular_depression",
        "channels": 1,
        "max_force_n": 100.0,
        "travel_m": 0.04,
    },
    {
        "id": "pectoralis_minor",
        "group": "scapulothoracic",
        "elements": (r"^pectoralis minor$",),
        "origin_bone": "thorax",
        "insertion_bone": "scapula",
        "target_joint": "scapulothoracic_downward_rotation",
        "channels": 1,
        "max_force_n": 160.0,
        "travel_m": 0.08,
    },
)

# Channel -> (joint, anatomical axis, action sign).  The action sign says which
# way a positive pull has to turn the joint, so it fixes the sign that dL/dq
# must have: tension drives the axis the other way, therefore a wire is only
# useful when sign(dL/dq) == -action sign.
CHANNEL_TARGETS: dict[str, tuple[str, str, float]] = {
    "supraspinatus": ("glenohumeral", "abduction", +1.0),
    "infraspinatus_teres_minor": ("glenohumeral", "internal_rotation", -1.0),
    "subscapularis": ("glenohumeral", "internal_rotation", +1.0),
    "teres_major": ("glenohumeral", "abduction", -1.0),
    "deltoid_anterior": ("glenohumeral", "flexion", +1.0),
    "deltoid_middle": ("glenohumeral", "abduction", +1.0),
    "deltoid_posterior": ("glenohumeral", "flexion", -1.0),
    "pectoralis_major": ("glenohumeral", "abduction", -1.0),
    "coracobrachialis": ("glenohumeral", "flexion", +1.0),
    "trapezius_levator": ("acromioclavicular", "upward_rotation", +1.0),
    "rhomboid_lower_trap": ("acromioclavicular", "upward_rotation", -1.0),
    "serratus_anterior": ("acromioclavicular", "upward_rotation", +1.0),
    "elbow_flexors": ("elbow", "flexion", +1.0),
    "triceps": ("elbow", "flexion", -1.0),
    "subclavius": ("sternoclavicular", "elevation", -1.0),
    "pectoralis_minor": ("acromioclavicular", "upward_rotation", -1.0),
}


# --------------------------------------------------------------------------
# BodyParts3D element index
# --------------------------------------------------------------------------


def load_index(path: Path = BP3D_INDEX) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8-sig") as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 3:
                rows.append(
                    {
                        "fma": parts[0].strip(),
                        "name": parts[1].strip(),
                        "file": parts[2].strip(),
                    }
                )
    if not rows:
        raise ValueError(f"{path} contained no usable rows.")
    return rows


def element_path(file_id: str, directory: Path = BP3D_DIR) -> Path:
    return directory / f"{file_id}.obj"


def resolve_elements(
    rows: list[dict[str, str]],
    patterns,
    require_files: bool = True,
) -> list[str]:
    """Right-side element ids whose name matches any pattern.

    File ids ending in ``M`` are the mirrored left side and are skipped, as
    are duplicates that appear under several concept names.
    """

    if isinstance(patterns, str):
        patterns = (patterns,)
    compiled = [re.compile(pattern) for pattern in patterns]
    found: list[str] = []
    for row in rows:
        name = row["name"]
        file_id = row["file"]
        if file_id.endswith("M"):
            continue
        if not any(pattern.search(name) for pattern in compiled):
            continue
        if file_id in found:
            continue
        if require_files and not element_path(file_id).exists():
            continue
        found.append(file_id)
    return found


# --------------------------------------------------------------------------
# OBJ input and output
# --------------------------------------------------------------------------


def read_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
    vertices: list[list[float]] = []
    faces: list[list[int]] = []
    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if line.startswith("v "):
                _, x, y, z = line.split()[:4]
                vertices.append([float(x), float(y), float(z)])
            elif line.startswith("f "):
                index = [int(token.split("/")[0]) - 1 for token in line.split()[1:]]
                for k in range(1, len(index) - 1):
                    faces.append([index[0], index[k], index[k + 1]])
    if not vertices or not faces:
        raise ValueError(f"{path} has no usable triangles.")
    return np.asarray(vertices, dtype=float), np.asarray(faces, dtype=np.int64)


def write_obj(path: Path, vertices: np.ndarray, faces: np.ndarray, header: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {header}" if header else "# bp3d arm asset"]
    lines += [f"v {v[0]:.9g} {v[1]:.9g} {v[2]:.9g}" for v in vertices]
    lines += [
        f"f {int(f[0]) + 1} {int(f[1]) + 1} {int(f[2]) + 1}" for f in faces
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def load_element_triangles(
    file_ids,
    directory: Path = BP3D_DIR,
    scale: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate the meshes of several elements into one vertex/face pair."""

    vertex_blocks: list[np.ndarray] = []
    face_blocks: list[np.ndarray] = []
    offset = 0
    for file_id in file_ids:
        vertices, faces = read_obj(element_path(file_id, directory))
        vertex_blocks.append(vertices * scale)
        face_blocks.append(faces + offset)
        offset += len(vertices)
    if not vertex_blocks:
        raise ValueError("No element ids were supplied.")
    return np.vstack(vertex_blocks), np.vstack(face_blocks)


def triangles_of(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    return vertices[faces]


# --------------------------------------------------------------------------
# Rigid-body mass properties
# --------------------------------------------------------------------------


def mesh_mass_properties(
    triangles: np.ndarray,
    density: float,
) -> dict[str, object]:
    """Volume, centre of mass and inertia of a closed triangle mesh.

    The integral of ``x x^T`` over one tetrahedron spanned by the origin and a
    triangle ``(a, b, c)`` is ``(V / 10) (a a^T + b b^T + c c^T + s s^T)``
    with ``s = a + b + c``, which is the compact form of the standard
    polyhedron mass-property sum.
    """

    a = triangles[:, 0]
    b = triangles[:, 1]
    c = triangles[:, 2]
    volume_six = np.einsum("ij,ij->i", a, np.cross(b, c))
    total_six = float(volume_six.sum())
    if abs(total_six) < 1e-12:
        raise ValueError("Mesh has zero volume.")
    # A consistently wound mesh can still integrate to a negative signed
    # volume depending on its winding convention, so the sum is flipped once
    # rather than taking absolute values per tetrahedron.
    if total_six < 0.0:
        volume_six = -volume_six
        total_six = -total_six
    volume = total_six / 6.0

    sum_six = total_six
    centre = (volume_six[:, None] * (a + b + c)).sum(axis=0) / (4.0 * sum_six)
    s = a + b + c
    second = (
        a[:, :, None] * a[:, None, :]
        + b[:, :, None] * b[:, None, :]
        + c[:, :, None] * c[:, None, :]
        + s[:, :, None] * s[:, None, :]
    )
    # vol_i = volume_six / 6 and the tetrahedron integral is vol_i/20, so the
    # per-triangle weight is volume_six / 120.  The density has to multiply
    # the whole integral, otherwise the parallel-axis step below would mix a
    # unit-density inertia with a density-scaled mass.
    moment = density * (volume_six[:, None, None] / 120.0 * second).sum(axis=0)

    mass = volume * density
    inertia_origin = np.trace(moment) * np.eye(3) - moment
    offset = centre
    inertia_centre = inertia_origin - mass * (
        float(offset @ offset) * np.eye(3) - np.outer(offset, offset)
    )
    inertia_centre = 0.5 * (inertia_centre + inertia_centre.T)
    values, vectors = np.linalg.eigh(inertia_centre)

    hull_volume = 0.0
    try:
        hull = ConvexHull(_decimate(triangles.reshape(-1, 3), 6000))
        hull_volume = float(hull.volume)
    except Exception:  # pragma: no cover - degenerate meshes only
        hull_volume = 0.0

    return {
        "volume_m3": volume,
        "hull_volume_m3": hull_volume,
        "watertight_ratio": (
            volume / hull_volume if hull_volume > 1e-12 else None
        ),
        "mass_kg": mass,
        "centre_m": centre,
        "inertia_diagonal": values,
        "inertia_quat": vectors,
    }


def _decimate(points: np.ndarray, limit: int) -> np.ndarray:
    if len(points) <= limit:
        return points
    step = int(np.ceil(len(points) / limit))
    return points[::step]


def convex_hull_mesh(
    vertices: np.ndarray,
    limit: int = 4000,
) -> tuple[np.ndarray, np.ndarray]:
    points = _decimate(np.unique(np.round(vertices, 6), axis=0), limit)
    hull = ConvexHull(points)
    # ConvexHull.simplices index into ``points``, so the vertex list has to be
    # remapped explicitly; trimming to hull.vertices without remapping would
    # leave every face pointing at a vertex that no longer exists.
    used = np.unique(hull.simplices)
    remap = np.zeros(len(points), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return points[used], remap[hull.simplices]


# --------------------------------------------------------------------------
# Surface fitting
# --------------------------------------------------------------------------


def sphere_fit(points: np.ndarray) -> dict[str, object]:
    """Least-squares sphere through a point set."""

    guess = points.mean(axis=0)
    centred = points - guess
    matrix = np.hstack([2.0 * centred, np.ones((len(centred), 1))])
    rhs = np.einsum("ij,ij->i", centred, centred)
    solution, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
    offset = solution[:3]
    centre = guess + offset
    radius = float(np.sqrt(max(solution[3] + float(offset @ offset), 0.0)))
    distances = np.linalg.norm(points - centre, axis=1)
    radius = 0.5 * (radius + float(distances.mean()))
    residual = float(np.sqrt(np.mean((distances - radius) ** 2)))
    return {"centre": centre, "radius": radius, "rms": residual}


def cylinder_fit(triangles: np.ndarray) -> dict[str, object]:
    """Axis, point and radius of the cylinder that best fits a surface band.

    The face normals of a cylinder lie in the plane perpendicular to its axis,
    so the axis is the direction of least normal variance.  The axis position
    then follows from an algebraic circle fit in the perpendicular plane.
    """

    a = triangles[:, 0]
    b = triangles[:, 1]
    c = triangles[:, 2]
    normals = np.cross(b - a, c - a)
    lengths = np.linalg.norm(normals, axis=1)
    keep = lengths > 1e-12
    normals = normals[keep] / lengths[keep, None]
    values, vectors = np.linalg.eigh(normals.T @ normals / len(normals))
    axis = vectors[:, 0]

    points = (a + b + c)[keep] / 3.0
    origin = points.mean(axis=0)
    u, v = _plane_basis(axis)
    local = points - origin
    x = local @ u
    y = local @ v
    matrix = np.column_stack([2.0 * x, 2.0 * y, np.ones(len(x))])
    rhs = x * x + y * y
    solution, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
    centre_2d = solution[:2]
    radius = float(np.sqrt(max(solution[2] + centre_2d @ centre_2d, 0.0)))
    axis_point = origin + centre_2d[0] * u + centre_2d[1] * v
    distances = np.linalg.norm(
        local - np.outer(local @ axis, axis), axis=1
    )
    residual = float(np.sqrt(np.mean((distances - radius) ** 2)))
    return {
        "point": axis_point,
        "axis": axis,
        "radius": radius,
        "rms": residual,
    }


def _plane_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    normal = normal / np.linalg.norm(normal)
    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(normal @ helper)) > 0.9:
        helper = np.array([1.0, 0.0, 0.0])
    u = np.cross(normal, helper)
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    return u, v / np.linalg.norm(v)


def principal_frame(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centre = points.mean(axis=0)
    centred = points - centre
    values, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    order = np.argsort(values)[::-1]
    return centre, vectors[:, order]


def plane_normal_fit(points: np.ndarray) -> np.ndarray:
    """Unit normal of the plane that best fits a flat point cloud."""

    centred = points - points.mean(axis=0)
    values, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    return vectors[:, 0]


def closest_pair(
    first: np.ndarray,
    second: np.ndarray,
    sample: int = 20000,
) -> tuple[np.ndarray, np.ndarray, float]:
    base = _decimate(first, sample)
    tree = cKDTree(second)
    distance, index = tree.query(base, workers=-1)
    best = int(np.argmin(distance))
    return base[best], second[index[best]], float(distance[best])


def ray_hits(
    origin: np.ndarray,
    direction: np.ndarray,
    triangles: np.ndarray,
) -> np.ndarray:
    """Positive ray parameters where the ray crosses the mesh (Moller-Trumbore)."""

    direction = direction / np.linalg.norm(direction)
    a = triangles[:, 0]
    b = triangles[:, 1]
    c = triangles[:, 2]
    edge1 = b - a
    edge2 = c - a
    pvec = np.cross(direction, edge2)
    det = np.einsum("ij,ij->i", edge1, pvec)
    ok = np.abs(det) > 1e-12
    inv = np.zeros_like(det)
    inv[ok] = 1.0 / det[ok]
    tvec = origin - a
    u = np.einsum("ij,ij->i", tvec, pvec) * inv
    qvec = np.cross(tvec, edge1)
    v = qvec @ direction * inv
    distance = np.einsum("ij,ij->i", edge2, qvec) * inv
    good = ok & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1.0 + 1e-9) & (distance > 1e-6)
    return distance[good]


# --------------------------------------------------------------------------
# Centre lines
# --------------------------------------------------------------------------


def plane_section_points(
    triangles: np.ndarray,
    point: np.ndarray,
    normal: np.ndarray,
) -> np.ndarray:
    """Intersection points of a plane with a triangle mesh."""

    normal = normal / np.linalg.norm(normal)
    distances = (triangles - point) @ normal
    segments: list[np.ndarray] = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        d0 = distances[:, i]
        d1 = distances[:, j]
        crossing = (d0 * d1) < 0.0
        if not crossing.any():
            continue
        weight = (d0[crossing] / (d0[crossing] - d1[crossing]))[:, None]
        edge = triangles[crossing, j] - triangles[crossing, i]
        segments.append(triangles[crossing, i] + weight * edge)
    if not segments:
        return np.empty((0, 3))
    return np.vstack(segments)


def largest_cluster_centroid(points: np.ndarray, tolerance: float) -> np.ndarray | None:
    """Centroid of the largest connected group of nearby points."""

    if len(points) == 0:
        return None
    tree = cKDTree(points)
    pairs = tree.query_pairs(tolerance, output_type="ndarray")
    parent = np.arange(len(points))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for left, right in pairs:
        root_left, root_right = find(int(left)), find(int(right))
        if root_left != root_right:
            parent[root_right] = root_left
    roots = np.array([find(i) for i in range(len(points))])
    counts = np.bincount(roots)
    best = int(np.argmax(counts))
    return points[roots == best].mean(axis=0)


def centre_line(
    triangles: np.ndarray,
    samples: int = 20,
    shrink: float = 0.06,
) -> np.ndarray:
    """Medial line of a roughly tubular muscle volume.

    The mesh is sliced perpendicular to its principal axis and the centroid of
    each cross-section becomes a point on the centre-line.  The two extreme
    slices are dropped because a slicing plane that only grazes the ends gives
    an unstable centroid.
    """

    points = triangles.reshape(-1, 3)
    centre, axes = principal_frame(points)
    axis = axes[:, 0]
    projection = (points - centre) @ axis
    low, high = np.percentile(projection, [shrink * 100.0, (1.0 - shrink) * 100.0])
    span = float(high - low)
    if span <= 0.0:
        raise ValueError("Muscle mesh has no extent along its principal axis.")
    spacing = span / samples
    tolerance = max(spacing * 0.75, 1e-4)

    line: list[np.ndarray] = []
    for index in range(samples + 1):
        offset = low + span * index / samples
        section = plane_section_points(triangles, centre + axis * offset, axis)
        centroid = largest_cluster_centroid(section, tolerance)
        if centroid is None:
            continue
        line.append(centroid)
    if len(line) < 2:
        raise ValueError("Could not extract a centre line from the muscle mesh.")
    return np.asarray(line)


# --------------------------------------------------------------------------
# Frames, joints and units
# --------------------------------------------------------------------------


def unit(vector: np.ndarray, label: str) -> np.ndarray:
    length = float(np.linalg.norm(vector))
    if length < 1e-12:
        raise ValueError(f"{label} has zero length.")
    return np.asarray(vector, dtype=float) / length


def perpendicular_to(axis: np.ndarray, reference: np.ndarray) -> np.ndarray:
    axis = unit(axis, "axis")
    reference = np.asarray(reference, dtype=float)
    projected = reference - float(reference @ axis) * axis
    if np.linalg.norm(projected) < 1e-9:
        u, _ = _plane_basis(axis)
        return u
    return unit(projected, "perpendicular")


def transform_points(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    return points @ matrix[:3, :3].T + matrix[:3, 3]


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, default=_json_default), encoding="utf-8"
    )


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    raise TypeError(f"Object of type {type(value)} is not JSON serializable")


def vector3(value) -> np.ndarray:
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.shape != (3,):
        raise ValueError(f"Expected three numbers, got {value!r}.")
    return array


def apply_pull_only_tension(model, data) -> int:
    """Cancel the compressive half of every ``assist_`` tendon spring.

    Sign convention, verified on a minimal MuJoCo model: when the tendon is
    shorter than its ``springlength`` the spring produces a POSITIVE
    generalized force that pushes the attachment points apart.  A cable or a
    ligament must never push, so the assist actuator has to supply the exact
    opposite, which is a NEGATIVE force.  The actuator therefore needs a
    ctrlrange and forcerange that include negative values, and the control
    written here is ``min(0, k (L - rest))``.

    Returns the number of assisted tendons.
    """

    import mujoco

    applied = 0
    for actuator in range(model.nu):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, actuator)
        if not name or not name.startswith("assist_"):
            continue
        if model.actuator_trntype[actuator] != mujoco.mjtTrn.mjTRN_TENDON:
            continue
        tendon = int(model.actuator_trnid[actuator, 0])
        stiffness = float(model.tendon_stiffness[tendon])
        rest = float(model.tendon_lengthspring[tendon, 0])
        length = float(data.ten_length[tendon])
        data.ctrl[actuator] = min(0.0, stiffness * (length - rest))
        applied += 1
    return applied


def bone_subtree(contract: dict, root: str) -> list[str]:
    """The bone and every bone downstream of it in the joint chain."""

    names = [root]
    changed = True
    while changed:
        changed = False
        for entry in contract["joints"].values():
            if entry["parent_bone"] in names and entry["child_bone"] not in names:
                names.append(entry["child_bone"])
                changed = True
    return names


def rotate_points(
    points: np.ndarray, centre: np.ndarray, axis: np.ndarray, degrees: float
) -> np.ndarray:
    """Rodrigues rotation of a point array about an axis through ``centre``."""

    axis = unit(axis, "rotation axis")
    angle = np.radians(degrees)
    local = np.asarray(points, dtype=float) - centre
    cos, sin = np.cos(angle), np.sin(angle)
    along = np.outer(local @ axis, axis)
    return centre + local * cos + np.cross(axis, local) * sin + along * (1.0 - cos)


def path_moment_arm(
    sites,
    contract: dict,
    joint: str,
    axis_name: str,
    delta_deg: float = 0.5,
    base_deg: float = 0.0,
) -> float:
    """dL/dq of a wire about one joint axis, in metres per radian.

    Rotating the joint moves the wire attachment points on the child bone and
    on everything downstream of it, exactly like the real articulated chain,
    so this reproduces what MuJoCo reports without needing the compiled model.
    """

    entry = contract["joints"][joint]
    centre = np.asarray(entry["origin_m"], dtype=float)
    axis = np.asarray(entry["axes"][axis_name], dtype=float)
    moving = set(bone_subtree(contract, entry["child_bone"]))
    bodies = [str(site["body"]) for site in sites]
    points = np.asarray([site["pos_m"] for site in sites], dtype=float)
    mask = np.array([body in moving for body in bodies])

    def length(degrees: float) -> float:
        moved = points.copy()
        moved[mask] = rotate_points(points[mask], centre, axis, base_deg + degrees)
        return float(np.linalg.norm(np.diff(moved, axis=0), axis=1).sum())

    return (length(delta_deg) - length(-delta_deg)) / np.radians(2.0 * delta_deg)


def wire_moment_arm_profile(
    sites,
    contract: dict,
    joint: str,
    axis_name: str,
    samples: int = 9,
) -> list[float]:
    """dL/dq sampled across the joint's range, to catch a sign flip."""

    low, high = (float(v) for v in contract["joints"][joint]["range_deg"])
    return [
        path_moment_arm(sites, contract, joint, axis_name, 0.5, angle)
        for angle in np.linspace(low, high, samples)
    ]


def insert_joint_pulley(
    sites,
    contract: dict,
    joint: str,
    axis_name: str,
    radius: float,
    wanted_sign: float,
) -> tuple[list[dict], dict | None]:
    """Add a wrapping point on the parent bone at the joint cylinder.

    A straight wire from a muscle belly to its insertion runs through the
    joint centre once the joint bends far enough, and its moment arm changes
    sign: past that angle the muscle pushes the joint the wrong way instead of
    driving it.  Real tendons do not have that problem because they wrap over
    the bone, so a wrapping point is added on the parent bone on whichever
    side of the axis keeps the moment arm the right way round.
    """

    entry = contract["joints"][joint]
    centre = np.asarray(entry["origin_m"], dtype=float)
    axis = unit(np.asarray(entry["axes"][axis_name], dtype=float), "pulley axis")
    moving = set(bone_subtree(contract, entry["child_bone"]))
    parent_bone = entry["parent_bone"]
    fixed = [index for index, site in enumerate(sites) if site["body"] not in moving]
    if not fixed or fixed[-1] + 1 >= len(sites):
        return sites, None
    # The wrapping point is bolted to the joint's parent bone even when the
    # muscle's own attachment is further up the chain, which is the case for
    # the deltoid (clavicle to humerus across a scapula-held joint).
    fixed_index = fixed[-1]
    moving_index = fixed_index + 1

    def score(profile) -> tuple[int, float]:
        aligned = sum(1 for value in profile if value * wanted_sign > 0.0)
        worst = min(abs(value) for value in profile)
        return aligned, worst

    baseline = wire_moment_arm_profile(sites, contract, joint, axis_name, samples=9)
    best = (score(baseline), None, None, baseline)
    u, v = _plane_basis(axis)
    for scale in (1.0, 0.65):
        for step in range(24):
            theta = 2.0 * np.pi * step / 24.0
            point = centre + radius * scale * (np.cos(theta) * u + np.sin(theta) * v)
            trial = list(sites)
            trial.insert(
                moving_index,
                {"body": parent_bone, "pos_m": point, "role": "pulley"},
            )
            profile = wire_moment_arm_profile(
                trial, contract, joint, axis_name, samples=9
            )
            candidate = score(profile)
            if candidate > best[0]:
                best = (candidate, point, trial, profile)
    if best[1] is None:
        return sites, None
    _, point, trial, profile = best
    return trial, {
        "joint": joint,
        "axis": axis_name,
        "point_m": point,
        "moment_arm_profile_mm_per_rad": [
            round(value * 1000.0, 2) for value in profile
        ],
    }
