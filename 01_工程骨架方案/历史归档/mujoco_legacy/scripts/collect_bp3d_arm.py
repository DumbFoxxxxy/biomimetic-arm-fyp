"""Collect the BodyParts3D right-arm assets and derive the joint frames.

This replaces the whole SolidWorks export chain: every BodyParts3D element is
already in one body-wide frame, so the bones land exactly where they belong
and no per-bone coordinate system has to be built by hand.  Joint centres come
from geometry:

* glenohumeral and proximal radioulnar: least-squares sphere through the
  articular head,
* elbow: least-squares cylinder through the distal humerus,
* sternoclavicular and acromioclavicular: midpoint of the closest vertex pair
  of the two bones that meet there.

Everything inside this module is in metres; millimetres appear only for
display.  Outputs land in ``mujoco_arm/bp3d`` and
``mujoco_arm/config/bp3d_arm.json``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import bp3d_common as B

# Anatomical mass of each dry bone.  The BodyParts3D meshes describe the outer
# surface of the whole bone, so applying cortical-bone density to their volume
# over-estimates the mass of a real, porous bone by two to three times.
BONE_MASS_KG = {
    "thorax": 5.0,
    "clavicle": 0.03,
    "scapula": 0.10,
    "humerus": 0.17,
    "radius": 0.05,
    "ulna": 0.07,
}

JOINT_FRAME_SPECS = (
    {
        "name": "sternoclavicular",
        "type": "ball",
        "parent_bone": "thorax",
        "child_bone": "clavicle",
        "axes": (
            ("elevation", "anterior_positive"),
            ("protraction", "superior"),
            ("axial_rotation", "clavicle_axis"),
        ),
        "primary_axis": "elevation",
        "range_deg": (-25.0, 30.0),
        "measure_axes": {
            "protraction": (-25.0, 25.0),
            "axial_rotation": (-30.0, 30.0),
        },
    },
    {
        "name": "acromioclavicular",
        "type": "ball",
        "parent_bone": "clavicle",
        "child_bone": "scapula",
        "axes": (
            ("upward_rotation", "anterior_positive"),
            ("posterior_tilt", "right"),
            ("external_rotation", "scapula_normal"),
        ),
        "primary_axis": "upward_rotation",
        "range_deg": (-10.0, 55.0),
        "measure_axes": {
            "posterior_tilt": (-25.0, 25.0),
            "external_rotation": (-25.0, 25.0),
        },
    },
    {
        "name": "glenohumeral",
        "type": "ball",
        "parent_bone": "scapula",
        "child_bone": "humerus",
        "axes": (
            ("abduction", "anterior_positive"),
            ("flexion", "right"),
            ("internal_rotation", "humerus_axis_proximal"),
        ),
        "primary_axis": "abduction",
        "range_deg": (-20.0, 115.0),
        "measure_axes": {
            "flexion": (-40.0, 115.0),
            "internal_rotation": (-80.0, 80.0),
        },
    },
    {
        "name": "elbow",
        "type": "hinge",
        "parent_bone": "humerus",
        "child_bone": "ulna",
        "axes": (("flexion", "elbow_axis"),),
        "primary_axis": "flexion",
        "range_deg": (0.0, 150.0),
        "measure_axes": {},
    },
    {
        "name": "proximal_radioulnar",
        "type": "ball",
        "parent_bone": "ulna",
        "child_bone": "radius",
        "axes": (
            ("pronation", "radius_axis"),
            ("valgus", "anterior_positive"),
            ("varus_twist", "right"),
        ),
        "primary_axis": "valgus",
        "range_deg": (-8.0, 8.0),
        "measure_axes": {
            "pronation": (-15.0, 15.0),
            "varus_twist": (-8.0, 8.0),
        },
    },
)


def load_bone(
    rows: list[dict[str, str]],
    name: str,
    directory: Path = B.BP3D_DIR,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Bone mesh in metres, still expressed in the BodyParts3D global frame."""

    spec = next(item for item in B.BONE_SPECS if item["name"] == name)
    elements = B.resolve_elements(rows, spec["patterns"])
    if not elements:
        raise FileNotFoundError(f"No BodyParts3D element matched bone {name}")
    vertices, faces = B.load_element_triangles(elements, directory, scale=B.MM)
    return elements, vertices, faces


def articular_band(
    vertices: np.ndarray,
    other: np.ndarray,
    radius: float,
) -> np.ndarray:
    """Vertices within ``radius`` of the closest approach to the other bone."""

    contact, _, _ = B.closest_pair(vertices, other)
    distance = np.linalg.norm(vertices - contact, axis=1)
    band = vertices[distance <= radius]
    if len(band) < 20:
        band = vertices[np.argsort(distance)[:20]]
    return band


def sphere_joint(
    child_vertices: np.ndarray,
    parent_vertices: np.ndarray,
    band_radius: float,
) -> dict[str, object]:
    band = articular_band(child_vertices, parent_vertices, band_radius)
    fit = B.sphere_fit(band)
    return {
        "origin": np.asarray(fit["centre"], dtype=float),
        "head_radius": float(fit["radius"]),
        "fit_rms_mm": float(fit["rms"]) * 1000.0,
        "band_vertices": int(len(band)),
    }


def cylinder_joint(
    child_vertices: np.ndarray,
    child_faces: np.ndarray,
    parent_vertices: np.ndarray,
    band_radius: float,
) -> dict[str, object]:
    contact, _, _ = B.closest_pair(child_vertices, parent_vertices)
    distance = np.linalg.norm(child_vertices - contact, axis=1)
    selected = distance <= band_radius
    if int(selected.sum()) < 30:
        selected = distance <= np.percentile(distance, 25.0)
    triangle_mask = selected[child_faces].all(axis=1)
    triangles = child_vertices[child_faces[triangle_mask]]
    if len(triangles) < 20:
        raise ValueError("Not enough triangles to fit an axis.")
    fit = B.cylinder_fit(triangles)
    return {
        "origin": np.asarray(fit["point"], dtype=float),
        "axis": np.asarray(fit["axis"], dtype=float),
        "cylinder_radius": float(fit["radius"]),
        "fit_rms_mm": float(fit["rms"]) * 1000.0,
        "band_triangles": int(len(triangles)),
    }


def design_radius(
    centre: np.ndarray,
    axis: np.ndarray,
    vertices: np.ndarray,
    band: float = 0.06,
    percentile: float = 60.0,
) -> float:
    """Radius of the bone surface a ligament collar would wrap around.

    Only vertices close to the joint centre are used, because the shoulder
    joints sit against a thorax that is far bigger than the joint itself.
    """

    axis = B.unit(axis, "joint axis")
    local = vertices - centre
    distance = np.linalg.norm(local, axis=1)
    near = distance <= band
    if int(near.sum()) < 30:
        near = distance <= np.percentile(distance, 5.0)
    local = local[near]
    radial = np.linalg.norm(local - np.outer(local @ axis, axis), axis=1)
    return float(np.percentile(radial, percentile))


def rib_gliding_spheres(
    rows: list[dict[str, str]],
    count: int = 20,
    radius_pad: float = 0.002,
) -> list[dict[str, object]]:
    """Sphere chain along the right ribs: the scapulothoracic gliding surface."""

    spheres: list[dict[str, object]] = []
    for file_id in B.resolve_elements(rows, r"^right \w+ rib$"):
        vertices, faces = B.read_obj(B.element_path(file_id))
        triangles = vertices[faces]
        line = B.centre_line(triangles, samples=30)
        length = float(np.linalg.norm(np.diff(line, axis=0), axis=1).sum())
        volume = float(B.mesh_mass_properties(triangles, 1.0)["volume_m3"])
        radius = float(np.sqrt(max(volume / (np.pi * max(length, 1e-6)), 1e-8)))
        radius = float(np.clip(radius, 3.0, 12.0)) + radius_pad
        segments = np.linalg.norm(np.diff(line, axis=0), axis=1)
        cumulative = np.concatenate([[0.0], np.cumsum(segments)])
        for target in np.linspace(0.0, cumulative[-1], count):
            index = int(
                np.clip(np.searchsorted(cumulative, target) - 1, 0, len(line) - 2)
            )
            span = max(cumulative[index + 1] - cumulative[index], 1e-9)
            weight = (target - cumulative[index]) / span
            point = line[index] + weight * (line[index + 1] - line[index])
            spheres.append(
                {
                    "centre_m": point * B.MM,
                    "radius_m": radius * B.MM,
                    "rib": file_id,
                }
            )
    return spheres


def stable_elbow_axis(
    humerus_vertices: np.ndarray,
    humerus_faces: np.ndarray,
    ulna_vertices: np.ndarray,
    bands_mm: tuple[float, ...] = (18.0, 22.0, 26.0),
) -> tuple[np.ndarray, float, float]:
    """Elbow axis that is stable against the width of the fitted band.

    A wide band also contains the humeral shaft, whose normals then win the
    least-variance fit and return a vertical axis.  The axis is therefore
    fitted at several band widths, aligned and averaged, and the worst angular
    disagreement between them is reported as a quality number.
    """

    axes: list[np.ndarray] = []
    radius = 0.0
    rms = 0.0
    for band in bands_mm:
        fit = cylinder_joint(
            humerus_vertices, humerus_faces, ulna_vertices, band * B.MM
        )
        axis = np.asarray(fit["axis"], dtype=float)
        if axes and float(axis @ axes[0]) < 0.0:
            axis = -axis
        axes.append(axis)
        radius = float(fit["cylinder_radius"])
        rms = float(fit["fit_rms_mm"])
    mean = B.unit(np.mean(axes, axis=0), "elbow axis")
    spread = float(
        np.degrees(
            max(
                np.arccos(np.clip(abs(float(axis @ mean)), -1.0, 1.0))
                for axis in axes
            )
        )
    )
    return mean, radius, spread, rms



def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assets", type=Path, default=B.ASSET_DIR)
    parser.add_argument("--contract", type=Path, default=B.CONTRACT_PATH)
    parser.add_argument("--head-band-mm", type=float, default=42.0)
    parser.add_argument("--elbow-band-mm", type=float, default=22.0)
    parser.add_argument("--rib-spheres", type=int, default=20)
    args = parser.parse_args()

    head_band = args.head_band_mm * B.MM
    elbow_band = args.elbow_band_mm * B.MM
    elbow_bands = (
        args.elbow_band_mm - 4.0,
        args.elbow_band_mm,
        args.elbow_band_mm + 4.0,
    )

    rows = B.load_index()
    print(f"BodyParts3D index rows: {len(rows)}")
    args.assets.mkdir(parents=True, exist_ok=True)

    bones: dict[str, dict[str, object]] = {}
    meshes: dict[str, tuple[np.ndarray, np.ndarray, list[str]]] = {}
    for spec in B.BONE_SPECS:
        name = spec["name"]
        elements, vertices, faces = load_bone(rows, name)
        meshes[name] = (vertices, faces, elements)
        properties = B.mesh_mass_properties(vertices[faces], 1000.0)
        mass = BONE_MASS_KG[name]
        inertia_scale = mass / max(float(properties["mass_kg"]), 1e-12)
        B.write_obj(
            args.assets / f"{name}.obj",
            vertices,
            faces,
            header=f"BodyParts3D {name} in metres, global frame",
        )
        extent_mm = np.ptp(vertices, axis=0) * 1000.0
        bones[name] = {
            "elements": elements,
            "vertices": int(len(vertices)),
            "triangles": int(len(faces)),
            "volume_cm3": round(float(properties["volume_m3"]) * 1e6, 1),
            "watertight_ratio": (
                None
                if properties["watertight_ratio"] is None
                else round(float(properties["watertight_ratio"]), 3)
            ),
            "mass_kg": mass,
            "mass_source": "anatomical dry bone mass",
            "centre_m": np.asarray(properties["centre_m"], dtype=float),
            "inertia_diagonal": np.asarray(
                properties["inertia_diagonal"], dtype=float
            )
            * inertia_scale,
            "inertia_quat_wxyz": _quat_wxyz(
                np.asarray(properties["inertia_quat"], dtype=float)
            ),
            "extent_mm": [round(float(v), 1) for v in extent_mm],
            "mesh_file": f"{name}.obj",
        }
        print(
            f"  {name:9s} elements={len(elements):3d} verts={len(vertices):6d} "
            f"vol={bones[name]['volume_cm3']:9.1f} cm3  mass={mass:5.2f} kg  "
            f"extent={np.round(extent_mm, 1)}"
        )

    clavicle_v = meshes["clavicle"][0]
    scapula_v = meshes["scapula"][0]
    humerus_v = meshes["humerus"][0]
    humerus_f = meshes["humerus"][1]
    radius_v = meshes["radius"][0]
    ulna_v = meshes["ulna"][0]

    sternum_ids = B.resolve_elements(rows, (r"^manubrium$", r"^body of sternum$"))
    sternum_v, _ = B.load_element_triangles(sternum_ids, scale=B.MM)

    # --- joint centres -----------------------------------------------------
    sc_point, sc_other, sc_gap = B.closest_pair(clavicle_v, sternum_v)
    ac_point, ac_other, ac_gap = B.closest_pair(clavicle_v, scapula_v)
    gh_fit = sphere_joint(humerus_v, scapula_v, head_band)
    elbow_fit = cylinder_joint(humerus_v, humerus_f, ulna_v, elbow_band)
    elbow_axis, elbow_radius, elbow_spread, elbow_rms = stable_elbow_axis(
        humerus_v, humerus_f, ulna_v, elbow_bands
    )
    elbow_fit["axis"] = elbow_axis
    elbow_fit["cylinder_radius"] = elbow_radius
    elbow_fit["fit_rms_mm"] = elbow_rms
    pruj_fit = sphere_joint(radius_v, ulna_v, head_band)
    druj_point, _, druj_gap = B.closest_pair(radius_v, ulna_v)
    elbow_gap = B.closest_pair(humerus_v, ulna_v)[2]

    sc_centre = 0.5 * (sc_point + sc_other)
    ac_centre = 0.5 * (ac_point + ac_other)

    scapula_normal = B.plane_normal_fit(scapula_v)
    if float(scapula_normal @ B.ANTERIOR) < 0.0:
        scapula_normal = -scapula_normal

    clavicle_axis = B.unit(ac_centre - sc_centre, "clavicle axis")
    humerus_axis = B.unit(elbow_fit["origin"] - gh_fit["origin"], "humerus axis")
    radius_axis = B.unit(druj_point - pruj_fit["origin"], "radius axis")

    axis_map = {
        "anterior": B.ANTERIOR,
        # A muscle's action axis: turning the right arm about the anterior
        # direction raises it medially, so abduction and clavicle / scapula
        # elevation use the opposite sign.  Keeping the sign equal to the
        # named action is what makes "a positive cable pull produces this
        # motion" true throughout the model.
        "anterior_positive": -B.ANTERIOR,
        "superior": B.SUPERIOR,
        "right": B.RIGHT,
        "clavicle_axis": clavicle_axis,
        "humerus_axis": humerus_axis,
        # Internal rotation of a hanging right arm is a negative turn about the
        # distal-pointing humeral axis, so the action axis points proximally.
        "humerus_axis_proximal": -humerus_axis,
        "radius_axis": radius_axis,
        "scapula_normal": scapula_normal,
        "elbow_axis": np.asarray(elbow_fit["axis"], dtype=float),
    }

    origins = {
        "sternoclavicular": sc_centre,
        "acromioclavicular": ac_centre,
        "glenohumeral": np.asarray(gh_fit["origin"], dtype=float),
        "elbow": np.asarray(elbow_fit["origin"], dtype=float),
        "proximal_radioulnar": np.asarray(pruj_fit["origin"], dtype=float),
    }

    # Flexion has to swing the ulna anteriorly, so the fitted cylinder axis is
    # oriented such that a positive angle is a real elbow flexion.
    ulna_distal = ulna_v[
        int(np.argmax(np.linalg.norm(ulna_v - origins["elbow"], axis=1)))
    ]
    lever = np.cross(axis_map["elbow_axis"], ulna_distal - origins["elbow"])
    if float(lever @ B.ANTERIOR) < 0.0:
        axis_map["elbow_axis"] = -axis_map["elbow_axis"]

    gaps_mm = {
        "sternoclavicular": sc_gap * 1000.0,
        "acromioclavicular": ac_gap * 1000.0,
        "glenohumeral": 0.0,
        "elbow": elbow_gap * 1000.0,
        "proximal_radioulnar": druj_gap * 1000.0,
    }

    joints: dict[str, dict[str, object]] = {}
    for spec in JOINT_FRAME_SPECS:
        name = spec["name"]
        origin = origins[name]
        axes = {
            axis_name: B.unit(axis_map[key], f"{name}.{axis_name}")
            for axis_name, key in spec["axes"]
        }
        radius = design_radius(
            origin,
            axes[spec["primary_axis"]],
            np.vstack(
                [meshes[spec["child_bone"]][0], meshes[spec["parent_bone"]][0]]
            ),
        )
        joints[name] = {
            "type": spec["type"],
            "parent_bone": spec["parent_bone"],
            "child_bone": spec["child_bone"],
            "origin_m": origin,
            "axes": axes,
            "primary_axis": spec["primary_axis"],
            "range_deg": list(spec["range_deg"]),
            "measure_axes": {
                axis: list(bounds) for axis, bounds in spec["measure_axes"].items()
            },
            "design_radius_m": radius,
            "gap_mm": round(float(gaps_mm[name]), 2),
        }
        if name == "glenohumeral":
            joints[name]["head_radius_m"] = float(gh_fit["head_radius"])
            joints[name]["fit_rms_mm"] = round(float(gh_fit["fit_rms_mm"]), 2)
        elif name == "elbow":
            joints[name]["fit_rms_mm"] = round(float(elbow_fit["fit_rms_mm"]), 2)
            joints[name]["cylinder_radius_m"] = float(elbow_fit["cylinder_radius"])
            joints[name]["axis_band_spread_deg"] = round(float(elbow_spread), 2)
        elif name == "proximal_radioulnar":
            joints[name]["fit_rms_mm"] = round(float(pruj_fit["fit_rms_mm"]), 2)

    print("\n== derived joint frames (millimetres, BodyParts3D frame) ==")
    for name, entry in joints.items():
        extra = (
            f"  fit rms={entry['fit_rms_mm']:.2f} mm" if "fit_rms_mm" in entry else ""
        )
        print(
            f"  {name:20s} {entry['type']:6s} "
            f"centre={np.round(np.asarray(entry['origin_m']) * 1000.0, 1)}"
            f"  design r={entry['design_radius_m'] * 1000.0:5.1f} mm"
            f"  gap={entry['gap_mm']:5.2f} mm{extra}"
        )
    print(f"  distal radioulnar contact gap = {druj_gap * 1000.0:.2f} mm")

    # --- collision proxies --------------------------------------------------
    collision: dict[str, object] = {}
    for name in ("clavicle", "scapula", "humerus", "radius", "ulna"):
        hull_vertices, hull_faces = B.convex_hull_mesh(meshes[name][0])
        B.write_obj(
            args.assets / f"{name}_hull.obj",
            hull_vertices,
            hull_faces,
            header=f"convex collision hull of {name}",
        )
        collision[f"{name}_hull"] = {
            "file": f"{name}_hull.obj",
            "vertices": int(len(hull_vertices)),
            "body": name,
        }
    for name, ids in (
        ("sternum", sternum_ids),
        (
            "spine",
            B.resolve_elements(
                rows, (r"^\w+ cervical vertebra$", r"^\w+ thoracic vertebra$")
            ),
        ),
    ):
        vertices, _ = B.load_element_triangles(ids, scale=B.MM)
        hull_vertices, hull_faces = B.convex_hull_mesh(vertices, limit=2500)
        B.write_obj(
            args.assets / f"{name}_hull.obj",
            hull_vertices,
            hull_faces,
            header=f"convex collision hull of {name}",
        )
        collision[f"{name}_hull"] = {
            "file": f"{name}_hull.obj",
            "vertices": int(len(hull_vertices)),
            "body": "thorax",
        }

    collision["rib_gliding"] = rib_gliding_spheres(rows, count=args.rib_spheres)
    collision["humeral_head"] = {
        "centre_m": joints["glenohumeral"]["origin_m"],
        "radius_m": joints["glenohumeral"]["head_radius_m"],
        "body": "humerus",
    }
    collision["elbow_cylinder"] = {
        "point_m": joints["elbow"]["origin_m"],
        "axis": joints["elbow"]["axes"]["flexion"],
        "radius_m": joints["elbow"]["cylinder_radius_m"],
        "body": "humerus",
    }

    contract = {
        "schema_version": 1,
        "source": str(B.BP3D_DIR),
        "units": "meters",
        "frame_convention": {
            "superior": B.SUPERIOR,
            "anterior": B.ANTERIOR,
            "right": B.RIGHT,
        },
        "bones": bones,
        "joints": joints,
        "collision": collision,
        "distal_radioulnar": {
            "contact_m": druj_point,
            "gap_mm": round(float(druj_gap * 1000.0), 2),
        },
    }
    B.save_json(args.contract, contract)
    B.save_json(
        B.FRAMES_PATH,
        {
            "units": "meters",
            "bones": {name: entry["mesh_file"] for name, entry in bones.items()},
            "joints": joints,
            "collision": collision,
        },
    )
    print(f"\ncontract: {args.contract}")
    print(f"frames:   {B.FRAMES_PATH}")


def _quat_wxyz(matrix: np.ndarray) -> list[float]:
    """Quaternion of the rotation whose columns are the inertia axes."""

    rotation = np.array(matrix, dtype=float)
    if np.linalg.det(rotation) < 0.0:
        rotation[:, 0] = -rotation[:, 0]
    trace = float(np.trace(rotation))
    if trace > 0.0:
        scale = np.sqrt(trace + 1.0) * 2.0
        w = 0.25 * scale
        x = (rotation[2, 1] - rotation[1, 2]) / scale
        y = (rotation[0, 2] - rotation[2, 0]) / scale
        z = (rotation[1, 0] - rotation[0, 1]) / scale
    elif rotation[0, 0] >= rotation[1, 1] and rotation[0, 0] >= rotation[2, 2]:
        scale = np.sqrt(1.0 + rotation[0, 0] - rotation[1, 1] - rotation[2, 2]) * 2.0
        w = (rotation[2, 1] - rotation[1, 2]) / scale
        x = 0.25 * scale
        y = (rotation[0, 1] + rotation[1, 0]) / scale
        z = (rotation[0, 2] + rotation[2, 0]) / scale
    elif rotation[1, 1] >= rotation[2, 2]:
        scale = np.sqrt(1.0 + rotation[1, 1] - rotation[0, 0] - rotation[2, 2]) * 2.0
        w = (rotation[0, 2] - rotation[2, 0]) / scale
        x = (rotation[0, 1] + rotation[1, 0]) / scale
        y = 0.25 * scale
        z = (rotation[1, 2] + rotation[2, 1]) / scale
    else:
        scale = np.sqrt(1.0 + rotation[2, 2] - rotation[0, 0] - rotation[1, 1]) * 2.0
        w = (rotation[1, 0] - rotation[0, 1]) / scale
        x = (rotation[0, 2] + rotation[2, 0]) / scale
        y = (rotation[1, 2] + rotation[2, 1]) / scale
        z = 0.25 * scale
    quaternion = np.array([w, x, y, z], dtype=float)
    return quaternion / np.linalg.norm(quaternion)


if __name__ == "__main__":
    main()
