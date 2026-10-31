"""Shared data contracts and constants for the cable-driven arm model."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parent
MODEL_ROOT = PROJECT_ROOT / "model"
DEFAULT_EXPORT_DIR = MODEL_ROOT / "exports"
DEFAULT_MANIFEST = ROOT / "config" / "cable_drive_manifest.json"
DEFAULT_OUTPUT_XML = ROOT / "arm_cable_drive.xml"
DEFAULT_REPORT_DIR = ROOT / "exports" / "cable_drive"
DEFAULT_BONE_SCALE_JSON = DEFAULT_EXPORT_DIR / "bone_scale.json"

BONE_NAMES = ("thorax", "clavicle", "scapula", "humerus", "ulna")
EXPECTED_JOINT_FRAMES = ("SC_J", "ACST_J", "GH_J", "ELBOW_J")

FRAME_AXIS_NAMES = {
    "SC_J": ("abduction_adduction", "elevation_depression", "axial_rotation"),
    "ACST_J": ("upward_rotation", "posterior_tilt", "external_rotation"),
    "GH_J": ("rotation", "flexion", "abduction"),
    "ELBOW_J": ("flexion_extension",),
}

# Root-to-leaf chain. Each entry names the joint that couples a child bone to
# its parent bone and the coordinate system that describes that joint.
JOINT_CHAIN = (
    {"joint": "SC_J", "parent": "thorax", "child": "clavicle"},
    {"joint": "ACST_J", "parent": "clavicle", "child": "scapula"},
    {"joint": "GH_J", "parent": "scapula", "child": "humerus"},
    {"joint": "ELBOW_J", "parent": "humerus", "child": "ulna"},
)

JOINT_SPECS = (
    {
        "body": "clavicle_protraction",
        "parent": "torso",
        "joint": "clavicle_protraction",
        "frame": "SC_J",
        "axis_index": 0,
        "mesh": None,
        "mass_kg": 0.01,
        "range_deg": (-30.0, 40.0),
    },
    {
        "body": "clavicle_elevation",
        "parent": "clavicle_protraction",
        "joint": "clavicle_elevation",
        "frame": "SC_J",
        "axis_index": 1,
        "mesh": None,
        "mass_kg": 0.01,
        "range_deg": (-20.0, 40.0),
    },
    {
        "body": "clavicle_axial_rotation",
        "parent": "clavicle_elevation",
        "joint": "clavicle_axial_rotation",
        "frame": "SC_J",
        "axis_index": 2,
        "mesh": "clavicle",
        "mass_kg": 0.15,
        "range_deg": (-30.0, 30.0),
    },
    {
        "body": "scapula_upward_rotation",
        "parent": "clavicle_axial_rotation",
        "joint": "scapula_upward_rotation",
        "frame": "ACST_J",
        "axis_index": 0,
        "mesh": None,
        "mass_kg": 0.01,
        "range_deg": (-20.0, 50.0),
    },
    {
        "body": "scapula_posterior_tilt",
        "parent": "scapula_upward_rotation",
        "joint": "scapula_posterior_tilt",
        "frame": "ACST_J",
        "axis_index": 1,
        "mesh": None,
        "mass_kg": 0.01,
        "range_deg": (-15.0, 35.0),
    },
    {
        "body": "scapula_external_rotation",
        "parent": "scapula_posterior_tilt",
        "joint": "scapula_external_rotation",
        "frame": "ACST_J",
        "axis_index": 2,
        "mesh": "scapula",
        "mass_kg": 0.15,
        "range_deg": (-20.0, 30.0),
    },
    {
        "body": "shoulder_rotation",
        "parent": "scapula_external_rotation",
        "joint": "shoulder_rotation",
        "frame": "GH_J",
        "axis_index": 0,
        "mesh": None,
        "mass_kg": 0.01,
        "range_deg": (-90.0, 90.0),
    },
    {
        "body": "shoulder_flexion",
        "parent": "shoulder_rotation",
        "joint": "shoulder_flexion",
        "frame": "GH_J",
        "axis_index": 1,
        "mesh": None,
        "mass_kg": 0.01,
        "range_deg": (-45.0, 180.0),
    },
    {
        "body": "shoulder_abduction",
        "parent": "shoulder_flexion",
        "joint": "shoulder_abduction",
        "frame": "GH_J",
        "axis_index": 2,
        "mesh": "humerus",
        "mass_kg": 1.8,
        "range_deg": (-45.0, 180.0),
    },
    {
        "body": "forearm_ulna",
        "parent": "shoulder_abduction",
        "joint": "elbow_flexion",
        "frame": "ELBOW_J",
        "axis_index": 0,
        "mesh": "ulna",
        "mass_kg": 0.45,
        "range_deg": (0.0, 145.0),
    },
)

SKELETON_JOINT_NAMES = tuple(spec["joint"] for spec in JOINT_SPECS)
SLIDER_JOINT_PREFIX = "slider_"
SLIDER_ACTUATOR_PREFIX = "pull_"

RHYTHM_ACTUATOR_TARGETS = {
    "clavicle_protraction_act": 0.10,
    "clavicle_elevation_act": 0.20,
    "clavicle_axial_rotation_act": 0.10,
    "scapula_upward_rotation_act": 0.30,
    "scapula_posterior_tilt_act": 0.15,
    "scapula_external_rotation_act": 0.10,
}

TASK_ENVELOPE_JOINTS = (
    "shoulder_flexion",
    "shoulder_abduction",
    "shoulder_rotation",
    "elbow_flexion",
)

RHYTHM_JOINT_GAINS = {
    "clavicle_protraction": 0.10,
    "clavicle_elevation": 0.20,
    "clavicle_axial_rotation": 0.10,
    "scapula_upward_rotation": 0.30,
    "scapula_posterior_tilt": 0.15,
    "scapula_external_rotation": 0.10,
}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def vec3(value: Any, label: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (3,):
        raise ValueError(f"{label} must contain exactly three values.")
    if not np.isfinite(array).all():
        raise ValueError(f"{label} contains a non-finite value.")
    return array


def normalize(vector: np.ndarray, label: str) -> np.ndarray:
    length = float(np.linalg.norm(vector))
    if length < 1e-12:
        raise ValueError(f"{label} has zero length.")
    return vector / length


def frame_axes(frame: dict[str, Any], frame_name: str) -> list[np.ndarray]:
    raw_axes = frame.get("axes")
    if not isinstance(raw_axes, list) or not raw_axes:
        raise ValueError(f"{frame_name} has no axes.")
    axes = []
    for index, axis in enumerate(raw_axes):
        if isinstance(axis, dict):
            direction = axis.get("direction")
            label = f"{frame_name}.axes[{index}].direction"
        else:
            direction = axis
            label = f"{frame_name}.axes[{index}]"
        axes.append(normalize(vec3(direction, label), label))
    return axes


def load_joint_frames(path: Path) -> dict[str, dict[str, Any]]:
    payload = load_json(path)
    frames = payload.get("frames", payload)
    missing = [name for name in EXPECTED_JOINT_FRAMES if name not in frames]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"joint_frames.json is missing: {joined}")

    result: dict[str, dict[str, Any]] = {}
    for name in EXPECTED_JOINT_FRAMES:
        frame = frames[name]
        origin = vec3(frame.get("origin_m"), f"{name}.origin_m")
        axes = frame_axes(frame, name)
        expected_axes = FRAME_AXIS_NAMES[name]
        if len(axes) < len(expected_axes):
            raise ValueError(
                f"{name} needs at least {len(expected_axes)} axes."
            )
        result[name] = {
            "origin_m": origin,
            "axes": axes,
            "range_deg": frame.get("range_deg"),
        }
    return result


def load_cable_manifest(path: Path) -> dict[str, Any]:
    payload = load_json(path)
    cables = payload.get("cables")
    if not isinstance(cables, list) or not cables:
        raise ValueError(f"{path} contains no cables.")
    enabled = [cable for cable in cables if cable.get("enabled", True)]
    if not enabled:
        raise ValueError(f"{path} contains no enabled cables.")
    return payload


def cable_assist_pairs(model) -> list[tuple[int, int]]:
    """(assist actuator id, tendon id) for every one-sided cable channel."""

    import mujoco

    pairs: list[tuple[int, int]] = []
    for actuator in range(model.nu):
        name = mujoco.mj_id2name(
            model, mujoco.mjtObj.mjOBJ_ACTUATOR, actuator
        )
        if not name or not name.startswith("assist_"):
            continue
        if model.actuator_trntype[actuator] != mujoco.mjtTrn.mjTRN_TENDON:
            continue
        pairs.append((actuator, int(model.actuator_trnid[actuator, 0])))
    return pairs


def apply_one_sided_tension(model, data, pairs) -> None:
    """Make every cable pull-only.

    The tendon springs are two-sided, so a cable that is shorter than its rest
    length would push.  Writing ``max(0, -k (L - rest))`` into the assist
    actuator cancels exactly that half and leaves ``max(0, k (L - rest))``.
    """

    for actuator, tendon in pairs:
        stiffness = float(model.tendon_stiffness[tendon])
        rest = float(model.tendon_lengthspring[tendon, 0])
        push = -stiffness * (float(data.ten_length[tendon]) - rest)
        data.ctrl[actuator] = push if push > 0.0 else 0.0


def load_bone_transforms(path: Path) -> dict[str, np.ndarray]:
    """part-to-world transform for every bone, as written by the exporter."""

    if not path.exists():
        return {}
    payload = load_json(path)
    raw = payload.get("transforms", payload)
    result: dict[str, np.ndarray] = {}
    for bone in BONE_NAMES:
        value = raw.get(bone)
        if value is None:
            continue
        matrix = np.asarray(value, dtype=float)
        if matrix.shape != (4, 4):
            raise ValueError(f"bone_transforms.json: {bone} is not a 4x4 matrix.")
        result[bone] = matrix
    return result


def load_bone_scales(path: Path = DEFAULT_BONE_SCALE_JSON) -> dict[str, float]:
    """Uniform part scale per bone, applied about the part origin.

    A part whose origin sits at its proximal joint centre can be lengthened by
    scaling the mesh and its part coordinate systems about that origin: the
    joint centre stays put and the bone grows distally.
    """

    scales = {bone: 1.0 for bone in BONE_NAMES}
    if not path.exists():
        return scales
    payload = load_json(path)
    for bone, value in (payload.get("scales") or {}).items():
        if bone not in scales:
            raise ValueError(f"{path.name}: unknown bone {bone}")
        factor = float(value)
        if factor <= 0.0:
            raise ValueError(f"{path.name}: {bone} has a non-positive scale")
        scales[bone] = factor
    return scales


def cable_guides(muscle: dict[str, Any]) -> list[dict[str, Any]]:
    """Return rigid bone-mounted wire guides in path order."""


    raw_guides = muscle.get("guides")
    if raw_guides is None:
        raw_guides = [
            {
                "body": muscle["via_body"],
                "local_m": muscle["via_local_m"],
                "radius_m": muscle.get("guide_radius_m", 0.006),
            }
        ]
    if not isinstance(raw_guides, list) or not raw_guides:
        raise ValueError(f"{muscle.get('id', 'cable')} has no wire guides.")

    result = []
    for index, guide in enumerate(raw_guides):
        part_point = guide.get("part_m")
        result.append(
            {
                "body": guide["body"],
                "local_m": vec3(
                    guide.get("local_m", guide.get("local")),
                    f"{muscle.get('id', 'cable')}.guides[{index}].local_m",
                ),
                "radius_m": float(guide.get("radius_m", 0.006)),
                "part_m": (
                    None
                    if part_point is None
                    else vec3(
                        part_point,
                        f"{muscle.get('id', 'cable')}.guides[{index}].part_m",
                    )
                ),
            }
        )
    return result


def component_base_name(component_name: str) -> str:
    stem = Path(component_name).stem
    if "-" in stem:
        stem = stem.rsplit("-", 1)[0]
    return stem.lower()


def triangle_centroid_and_box(stl_path: Path) -> tuple[np.ndarray, np.ndarray]:
    triangles = read_binary_stl(stl_path)
    vertices = triangles.reshape(-1, 3)
    return vertices.mean(axis=0), np.ptp(vertices, axis=0)


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
        raise ValueError("triangles must have shape (n, 3, 3).")
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0],
        triangles[:, 2] - triangles[:, 0],
    )
    lengths = np.linalg.norm(normals, axis=1)
    lengths[lengths < 1e-12] = 1.0
    normals = normals / lengths[:, None]

    with path.open("wb") as handle:
        header_bytes = header.encode("ascii", errors="replace")[:80]
        handle.write(header_bytes.ljust(80, b"\0"))
        handle.write(np.asarray([len(triangles)], dtype="<u4").tobytes())
        for normal, triangle in zip(normals, triangles):
            handle.write(np.asarray(normal, dtype="<f4").tobytes())
            handle.write(np.asarray(triangle, dtype="<f4").tobytes())
            handle.write(np.asarray([0], dtype="<u2").tobytes())


def transform_from_array(values: list[float]) -> np.ndarray:
    if len(values) != 16:
        raise ValueError("A SolidWorks transform must contain 16 values.")
    data = np.asarray(values, dtype=float)
    transform = np.eye(4)
    transform[:3, :3] = data[:9].reshape(3, 3)
    transform[:3, 3] = data[9:12]
    return transform


def task_envelope_states(
    model,
    random_samples: int = 600,
    seed: int = 11,
) -> list[np.ndarray]:
    """Return qpos samples covering the intended shoulder/elbow motions."""

    import mujoco

    rng = np.random.default_rng(seed)
    addresses: dict[str, int] = {}
    for joint_name in (
        "shoulder_flexion",
        "shoulder_abduction",
        "shoulder_rotation",
        "elbow_flexion",
        "clavicle_protraction",
        "clavicle_elevation",
        "clavicle_axial_rotation",
        "scapula_upward_rotation",
        "scapula_posterior_tilt",
        "scapula_external_rotation",
    ):
        joint_id = mujoco.mj_name2id(
            model, mujoco.mjtObj.mjOBJ_JOINT, joint_name
        )
        if joint_id < 0:
            raise ValueError(f"Joint required for calibration is missing: {joint_name}")
        addresses[joint_name] = int(model.jnt_qposadr[joint_id])

    def state(
        flexion_deg: float,
        abduction_deg: float,
        rotation_deg: float,
        elbow_deg: float,
    ) -> np.ndarray:
        qpos = np.zeros(model.nq)
        flexion = np.radians(flexion_deg)
        abduction = np.radians(abduction_deg)
        rotation = np.radians(rotation_deg)
        elbow = np.radians(elbow_deg)
        elevation = float(np.hypot(max(flexion, 0.0), max(abduction, 0.0)))
        depression = min(flexion, abduction, 0.0)
        qpos[addresses["shoulder_flexion"]] = flexion
        qpos[addresses["shoulder_abduction"]] = abduction
        qpos[addresses["shoulder_rotation"]] = rotation
        qpos[addresses["elbow_flexion"]] = elbow
        for joint_name, gain in RHYTHM_JOINT_GAINS.items():
            qpos[addresses[joint_name]] = gain * elevation + 0.03 * depression
        return qpos

    states = [np.zeros(model.nq)]
    flexion_values = (0.0, 30.0, 60.0, 90.0, 120.0)
    abduction_values = (0.0, 30.0, 60.0, 90.0, 120.0)
    rotation_values = (-45.0, -22.5, 0.0, 22.5, 45.0)
    elbow_values = (0.0, 45.0, 90.0, 145.0)
    for flexion in flexion_values:
        for abduction in abduction_values:
            for elbow in elbow_values:
                states.append(state(flexion, abduction, 0.0, elbow))
    for _ in range(random_samples):
        states.append(
            state(
                float(rng.uniform(0.0, 130.0)),
                float(rng.uniform(0.0, 130.0)),
                float(rng.uniform(-60.0, 60.0)),
                float(rng.uniform(0.0, 145.0)),
            )
        )
    return states
