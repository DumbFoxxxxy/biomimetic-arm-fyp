"""Adopt the exported right-arm STLs and build a CAD-based reference model."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
import struct
import xml.etree.ElementTree as ET

import numpy as np


ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parent
BONE_DIR = PROJECT_ROOT / "骨骼"
CAD_ASSET_DIR = ROOT / "assets_cad"
BASE_XML = ROOT / "arm_reference.xml"
CAD_XML = ROOT / "arm_cad_reference.xml"
ARM_CONFIG = ROOT / "config" / "arm_model.json"
CAD_ARM_CONFIG = ROOT / "config" / "arm_model_cad.json"
MUSCLE_CONFIG = ROOT / "config" / "muscle_manifest.json"
CAD_MUSCLE_CONFIG = ROOT / "config" / "muscle_manifest_cad.json"
REGISTRATION_PATH = ROOT / "config" / "cad_registration.json"


OBJ_FILES = {
    "clavicle": BONE_DIR / "锁骨.obj",
    "scapula": BONE_DIR / "肩胛骨.obj",
    "humerus": BONE_DIR / "肱骨.obj",
    "radius": BONE_DIR / "桡骨.obj",
    "ulna": BONE_DIR / "尺骨.obj",
}

STL_FILES = {
    "clavicle": BONE_DIR / "part - 锁骨-2.STL",
    "scapula": BONE_DIR / "part - 肩胛骨-1.STL",
    "humerus": BONE_DIR / "part - 肱骨-1.STL",
    "radius": BONE_DIR / "part - 桡骨-1.STL",
    "ulna": BONE_DIR / "part - 尺骨-1.STL",
}

# Three corresponding landmarks per bone, in millimetres. These values are
# close to the anatomical landmarks used by the reference model.
LANDMARKS = {
    "clavicle": {
        "old": [
            [-11.42, -129.58, 1319.03],
            [-149.46, -70.63, 1349.59],
            [-80.00, -100.00, 1345.00],
        ],
        "new": [
            [270.86, 18.37, 1319.49],
            [147.17, 81.49, 1344.38],
            [205.00, 60.00, 1354.00],
        ],
    },
    "scapula": {
        "old": [
            [-149.46, -70.63, 1349.59],
            [-166.87, -76.59, 1332.51],
            [-94.77, -4.30, 1190.73],
        ],
        "new": [
            [147.17, 81.49, 1344.38],
            [137.13, 83.33, 1306.56],
            [187.56, 144.12, 1190.75],
        ],
    },
    "humerus": {
        "old": [
            [-166.87, -76.59, 1332.51],
            [-202.31, -70.81, 1033.51],
            [-237.32, -66.95, 1050.09],
        ],
        "new": [
            [137.13, 83.33, 1306.56],
            [80.54, 77.17, 1033.41],
            [44.55, 81.29, 1049.85],
        ],
    },
    "radius": {
        "old": [
            [-233.68, -71.29, 1036.22],
            [-274.33, -116.87, 806.70],
            [-232.25, -76.89, 1033.64],
        ],
        "new": [
            [53.92, 78.86, 1035.94],
            [7.72, 31.13, 806.95],
            [49.94, 70.75, 1033.82],
        ],
    },
    "ulna": {
        "old": [
            [-207.81, -58.64, 1051.39],
            [-236.63, -112.33, 806.50],
            [-202.31, -70.81, 1033.51],
        ],
        "new": [
            [74.77, 90.06, 1051.41],
            [45.62, 36.23, 807.83],
            [80.54, 77.17, 1033.41],
        ],
    },
}

OLD_PIVOTS = {
    "torso": np.array([0.0, 0.0, 0.0]),
    "clavicle_protraction": np.array([-11.42, -129.58, 1319.03]),
    "clavicle_elevation": np.array([-11.42, -129.58, 1319.03]),
    "clavicle_axial_rotation": np.array([-11.42, -129.58, 1319.03]),
    "scapula_upward_rotation": np.array([-149.46, -70.63, 1349.59]),
    "scapula_posterior_tilt": np.array([-149.46, -70.63, 1349.59]),
    "scapula_external_rotation": np.array([-149.46, -70.63, 1349.59]),
    "shoulder_rotation": np.array([-166.87, -76.59, 1332.51]),
    "shoulder_flexion": np.array([-166.87, -76.59, 1332.51]),
    "shoulder_abduction": np.array([-166.87, -76.59, 1332.51]),
    "forearm_ulna": np.array([-214.12, -68.78, 1039.50]),
    "radius": np.array([-233.68, -71.29, 1036.22]),
}

NEW_PIVOTS = {
    "torso": np.array([274.432, 18.062, 1322.786]),
    "clavicle_protraction": np.array([274.432, 18.062, 1322.786]),
    "clavicle_elevation": np.array([274.432, 18.062, 1322.786]),
    "clavicle_axial_rotation": np.array([274.432, 18.062, 1322.786]),
    "scapula_upward_rotation": np.array([147.170, 81.490, 1344.380]),
    "scapula_posterior_tilt": np.array([147.170, 81.490, 1344.380]),
    "scapula_external_rotation": np.array([147.170, 81.490, 1344.380]),
    "shoulder_rotation": np.array([137.128, 83.325, 1306.558]),
    "shoulder_flexion": np.array([137.128, 83.325, 1306.558]),
    "shoulder_abduction": np.array([137.128, 83.325, 1306.558]),
    "forearm_ulna": np.array([74.831, 88.674, 1043.236]),
    "radius": np.array([53.919, 78.856, 1035.942]),
}

BODY_TO_BONE = {
    "torso": "clavicle",
    "clavicle_protraction": "clavicle",
    "clavicle_elevation": "clavicle",
    "clavicle_axial_rotation": "clavicle",
    "scapula_upward_rotation": "scapula",
    "scapula_posterior_tilt": "scapula",
    "scapula_external_rotation": "scapula",
    "shoulder_rotation": "humerus",
    "shoulder_flexion": "humerus",
    "shoulder_abduction": "humerus",
    "forearm_ulna": "ulna",
    "radius": "radius",
}


def read_obj(path: Path) -> np.ndarray:
    points = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("v "):
            points.append([float(x) for x in line.split()[1:4]])
    return np.asarray(points, dtype=float)


def read_stl(path: Path) -> np.ndarray:
    data = path.read_bytes()
    count = struct.unpack("<I", data[80:84])[0]
    raw = np.frombuffer(
        data[84 : 84 + count * 50], dtype=np.uint8
    ).reshape(count, 50)
    triangles = raw[:, 12:48].copy().view("<f4").reshape(-1, 3, 3)
    return np.unique(triangles.reshape(-1, 3).astype(float), axis=0)


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    translation = target_center - rotation @ source_center
    return rotation, translation


def bone_transforms() -> dict:
    result = {}
    for bone, landmarks in LANDMARKS.items():
        source = np.asarray(landmarks["old"], dtype=float)
        target = np.asarray(landmarks["new"], dtype=float)
        rotation, translation = kabsch(source, target)
        source_points = read_obj(OBJ_FILES[bone])
        target_points = read_stl(STL_FILES[bone])
        transformed = source_points @ rotation.T + translation
        # A coarse nearest-point residual is enough to detect a bad landmark set.
        sample = transformed[:: max(1, len(transformed) // 1200)]
        target_sample = target_points[:: max(1, len(target_points) // 1200)]
        distance = np.linalg.norm(
            sample[:, None, :] - target_sample[None, :, :], axis=2
        ).min(axis=1)
        result[bone] = {
            "rotation": rotation.tolist(),
            "translation_mm": translation.tolist(),
            "sample_rmse_mm": float(np.sqrt(np.mean(distance**2))),
        }
    return result


def transform_point(
    transforms: dict,
    body_name: str,
    local_point_mm,
) -> np.ndarray:
    bone = BODY_TO_BONE[body_name]
    old_global = OLD_PIVOTS[body_name] + np.asarray(
        local_point_mm, dtype=float
    )
    rotation = np.asarray(transforms[bone]["rotation"], dtype=float)
    translation = np.asarray(transforms[bone]["translation_mm"], dtype=float)
    new_global = rotation @ old_global + translation
    return (new_global - NEW_PIVOTS[body_name]) * 0.001


def set_pos(element: ET.Element, point_mm) -> None:
    point_m = np.asarray(point_mm, dtype=float) * 0.001
    element.set("pos", " ".join(f"{x:.9g}" for x in point_m))


def geom_by_name(root: ET.Element, name: str) -> ET.Element:
    for geom in root.iter("geom"):
        if geom.get("name") == name:
            return geom
    raise ValueError(f"Geometry not found: {name}")


def site_by_name(root: ET.Element, name: str) -> ET.Element:
    for site in root.iter("site"):
        if site.get("name") == name:
            return site
    raise ValueError(f"Site not found: {name}")


def body_by_name(root: ET.Element, name: str) -> ET.Element:
    for body in root.iter("body"):
        if body.get("name") == name:
            return body
    raise ValueError(f"Body not found: {name}")


def main() -> None:
    CAD_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for bone, source in STL_FILES.items():
        shutil.copyfile(source, CAD_ASSET_DIR / f"{bone}.stl")

    transforms = bone_transforms()
    REGISTRATION_PATH.write_text(
        json.dumps(
            {
                "source": "SolidWorks assembly STL exports",
                "units": "millimetres",
                "transforms": transforms,
                "old_pivots_mm": {
                    key: value.tolist() for key, value in OLD_PIVOTS.items()
                },
                "new_pivots_mm": {
                    key: value.tolist() for key, value in NEW_PIVOTS.items()
                },
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    tree = ET.parse(BASE_XML)
    root = tree.getroot()
    compiler = root.find("compiler")
    if compiler is not None:
        compiler.set("meshdir", CAD_ASSET_DIR.resolve().as_posix())

    for mesh in root.find("asset").iter("mesh"):
        name = mesh.get("name")
        if name in STL_FILES:
            mesh.set("file", f"{name}.stl")

    for body_name, position in NEW_PIVOTS.items():
        body = body_by_name(root, body_name)
        if body_name == "torso":
            set_pos(body, position)
        elif body_name in {
            "clavicle_protraction",
            "clavicle_elevation",
            "clavicle_axial_rotation",
            "scapula_posterior_tilt",
            "scapula_external_rotation",
            "shoulder_flexion",
            "shoulder_abduction",
        }:
            set_pos(body, [0.0, 0.0, 0.0])
        else:
            parent_name = body.get("name")
            parent_map = {
                "scapula_upward_rotation": "torso",
                "shoulder_rotation": "scapula_external_rotation",
                "forearm_ulna": "shoulder_abduction",
                "radius": "forearm_ulna",
            }
            set_pos(body, position - NEW_PIVOTS[parent_map[parent_name]])

    geom_positions = {
        "clavicle_visual": NEW_PIVOTS["clavicle_axial_rotation"],
        "scapula_visual": NEW_PIVOTS["scapula_external_rotation"],
        "humerus_visual": NEW_PIVOTS["shoulder_abduction"],
        "ulna_visual": NEW_PIVOTS["forearm_ulna"],
        "radius_visual": NEW_PIVOTS["radius"],
    }
    for name, pivot in geom_positions.items():
        set_pos(geom_by_name(root, name), -pivot)

    # Transform the reference thorax capsules into the new assembly frame via
    # the clavicle registration. The torso body origin is the new SC centre.
    torso_transform = transforms["clavicle"]
    rotation = np.asarray(torso_transform["rotation"], dtype=float)
    translation = np.asarray(torso_transform["translation_mm"], dtype=float)
    for geom in body_by_name(root, "torso").findall("geom"):
        fromto = geom.get("fromto")
        if not fromto:
            continue
        start, end = np.array(
            [float(x) for x in fromto.split()], dtype=float
        ).reshape(2, 3)
        start_new = rotation @ start + translation - NEW_PIVOTS["torso"]
        end_new = rotation @ end + translation - NEW_PIVOTS["torso"]
        geom.set(
            "fromto",
            " ".join(f"{x:.9g}" for x in np.concatenate((start_new, end_new)) * 0.001),
        )
    set_pos(site_by_name(root, "sternoclavicular_ref"), [0.0, 0.0, 0.0])
    set_pos(site_by_name(root, "acromioclavicular_ref"), [0.0, 0.0, 0.0])

    elbow_local = NEW_PIVOTS["forearm_ulna"] - NEW_PIVOTS["shoulder_abduction"]
    set_pos(site_by_name(root, "elbow_ref"), elbow_local)
    contact_local = NEW_PIVOTS["radius"] - NEW_PIVOTS["forearm_ulna"]
    distal_ulna_local = np.array([45.62, 36.23, 807.83]) - NEW_PIVOTS[
        "forearm_ulna"
    ]
    wrist_local = np.array([26.64, 33.61, 807.09]) - NEW_PIVOTS[
        "forearm_ulna"
    ]
    set_pos(site_by_name(root, "radiocapitellar_contact_ref"), contact_local)
    set_pos(site_by_name(root, "distal_ulnar_head_ref"), distal_ulna_local)
    set_pos(site_by_name(root, "wrist_ref"), wrist_local)

    axis_geom = geom_by_name(root, "pronation_axis_visual")
    set_pos(axis_geom, [0.0, 0.0, 0.0])
    axis_geom.attrib.pop("pos", None)
    axis_geom.set(
        "fromto",
        " ".join(
            f"{x:.9g}"
            for x in np.concatenate((contact_local, distal_ulna_local)) * 0.001
        ),
    )
    pronation_axis = distal_ulna_local - contact_local
    pronation_axis = pronation_axis / np.linalg.norm(pronation_axis)
    for joint in root.iter("joint"):
        if joint.get("name") == "forearm_pronation":
            joint.set("axis", " ".join(f"{x:.9g}" for x in pronation_axis))

    ET.indent(tree, space="  ")
    tree.write(CAD_XML, encoding="utf-8", xml_declaration=True)

    arm_config = json.loads(ARM_CONFIG.read_text(encoding="utf-8"))
    for ligament in arm_config["ligaments"]:
        for prefix, point_key in (
            ("origin", "origin_local"),
            ("insertion", "insertion_local"),
        ):
            body_name = ligament[f"{prefix}_body"]
            transformed = transform_point(
                transforms, body_name, ligament[point_key]
            )
            ligament[point_key] = transformed.tolist()
    CAD_ARM_CONFIG.write_text(
        json.dumps(
            {**arm_config, "base_xml": CAD_XML.name, "mesh_asset_dir": "assets_cad",
             "output_xml": "arm_psm_generated_cad.xml"},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    muscle_config = json.loads(MUSCLE_CONFIG.read_text(encoding="utf-8"))
    for muscle in muscle_config["muscles"]:
        for prefix, point_key in (
            ("origin", "origin_local"),
            ("insertion", "insertion_local"),
        ):
            body_name = muscle[f"{prefix}_body"]
            transformed = transform_point(
                transforms, body_name, muscle[point_key]
            )
            muscle[point_key] = transformed.tolist()
        for via in muscle.get("via", []):
            transformed = transform_point(transforms, via["body"], via["local"])
            via["local"] = transformed.tolist()
    CAD_MUSCLE_CONFIG.write_text(
        json.dumps(muscle_config, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(f"CAD assets: {CAD_ASSET_DIR}")
    print(f"CAD reference: {CAD_XML}")
    print(f"CAD arm config: {CAD_ARM_CONFIG}")
    print(f"CAD muscle manifest: {CAD_MUSCLE_CONFIG}")
    for bone, data in transforms.items():
        print(f"{bone:10s} sample_rmse={data['sample_rmse_mm']:.3f} mm")


if __name__ == "__main__":
    main()
