"""Build the multi-cable force-aware V2 MuJoCo model."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import json
from project_runtime import WorkspacePath as Path
import xml.etree.ElementTree as ET

import numpy as np

from v2_common import (
    BODY_NAME,
    CONFIG_PATH,
    MODEL_PATH,
    MUJOCO_ROOT,
    load_json,
    normalize,
    quat_from_z_axis,
    vec,
)


BASE_MODEL = MUJOCO_ROOT / "arm_tendon_demo.xml"


def find_body(root: ET.Element, name: str) -> ET.Element:
    for body in root.iter("body"):
        if body.get("name") == name:
            return body
    raise KeyError(name)


def find_spatial(root: ET.Element, name: str) -> ET.Element:
    for tendon in root.iter("spatial"):
        if tendon.get("name") == name:
            return tendon
    raise KeyError(name)


def ensure_spatial(root: ET.Element, name: str) -> ET.Element:
    for tendon in root.iter("spatial"):
        if tendon.get("name") == name:
            return tendon
    tendon_root = root.find("tendon")
    if tendon_root is None:
        tendon_root = ET.SubElement(root, "tendon")
    return ET.SubElement(
        tendon_root,
        "spatial",
        {"name": name, "width": "0.0018", "rgba": "0.95 0.45 0.12 0.95"},
    )


def ensure_asset_mesh(asset: ET.Element, name: str, path: Path) -> None:
    for mesh in asset.findall("mesh"):
        if mesh.get("name") == name:
            return
    ET.SubElement(
        asset,
        "mesh",
        {"name": name, "file": str(path.resolve()).replace("\\", "/")},
    )


def ensure_asset_material(asset: ET.Element, name: str, rgba: str) -> None:
    for material in asset.findall("material"):
        if material.get("name") == name:
            return
    ET.SubElement(asset, "material", {"name": name, "rgba": rgba})


def add_site(
    parent: ET.Element,
    name: str,
    local_position_mj: np.ndarray,
    rgba: str,
    size: float = 0.003,
) -> None:
    ET.SubElement(
        parent,
        "site",
        {
            "name": name,
            "pos": vec(local_position_mj),
            "size": f"{size:.6f}",
            "rgba": rgba,
            "group": "3",
        },
    )


def build(config: dict) -> tuple[ET.ElementTree, dict]:
    tree = ET.parse(BASE_MODEL)
    root = tree.getroot()
    root.set("model", config["model_name"])
    option = root.find("option")
    if option is None:
        option = ET.SubElement(root, "option")
    option.set("gravity", f"0 0 {-float(config['load_case']['gravity_m_s2']):.6f}")
    option.set("timestep", str(config["simulation"]["timestep_s"]))

    asset = root.find("asset")
    if asset is None:
        raise RuntimeError("Base model has no asset section")
    band_path = Path(config["contact"]["band_mesh"])
    pad_path = Path(config["contact"]["pad_mesh"])
    ensure_asset_mesh(asset, "v2_thorax_contact_band", band_path)
    ensure_asset_mesh(asset, "v2_scapula_contact_pad", pad_path)
    ensure_asset_material(asset, "v2_band_mat", "0.10 0.75 0.78 0.42")
    ensure_asset_material(asset, "v2_pad_mat", "0.95 0.48 0.12 0.75")
    ensure_asset_material(asset, "v2_guide_mat", "0.94 0.26 0.16 1")

    # Use geometry mass for the actual bone meshes.
    density = float(config["materials"]["polymer_density_kg_m3"])
    for name in ("clavicle_visual", "scapula_mesh", "humerus_shaft", "ulna_shaft"):
        for geom in root.iter("geom"):
            if geom.get("name") == name:
                geom.set("density", f"{density:.1f}")

    # Full bone meshes are visual-only in V2.
    for name in ("thorax_collision", "clavicle_collision", "scapula_collision"):
        for geom in root.iter("geom"):
            if geom.get("name") == name:
                geom.set("contype", "0")
                geom.set("conaffinity", "0")

    contact = config["contact"]
    thorax = find_body(root, "thorax")
    ET.SubElement(
        thorax,
        "geom",
        {
            "name": "v2_thorax_contact_band",
            "type": "mesh",
            "mesh": "v2_thorax_contact_band",
            "material": "v2_band_mat",
            "contype": "1",
            "conaffinity": "2",
            "condim": str(contact["condim"]),
            "friction": vec(contact["friction"]),
            "solref": vec(contact["solref"]),
            "solimp": vec(contact["solimp"]),
            "density": "0",
            "group": "1",
        },
    )
    scapula = find_body(root, "scapula")
    ET.SubElement(
        scapula,
        "geom",
        {
            "name": "v2_scapula_contact_pad",
            "type": "mesh",
            "mesh": "v2_scapula_contact_pad",
            "material": "v2_pad_mat",
            "contype": "2",
            "conaffinity": "1",
            "condim": str(contact["condim"]),
            "friction": vec(contact["friction"]),
            "solref": vec(contact["solref"]),
            "solimp": vec(contact["solimp"]),
            "density": "0",
            "group": "1",
        },
    )

    for joint in root.iter("joint"):
        name = joint.get("name")
        if name in config["joint_limits_deg"]:
            joint.set("range", vec(config["joint_limits_deg"][name]))

    route_metadata = {}
    for route in config["routes"]:
        route_id = route["id"]
        tendon = ensure_spatial(root, f"cable_{route_id}")
        for child in list(tendon):
            tendon.remove(child)
        metadata_path = []
        anchor_number = 0
        for element in route["path"]:
            if element["kind"] == "anchor":
                owner = element["owner"]
                body_name = BODY_NAME[owner]
                site_name = f"V2_{'AB'[anchor_number]}_{route_id}"
                anchor_number += 1
                if owner == "X":
                    parent = root.find("worldbody")
                    local = np.asarray(element["local_mm"], dtype=float)
                else:
                    parent = find_body(root, body_name)
                    local = np.asarray(element["local_mm"], dtype=float)
                add_site(
                    parent,
                    site_name,
                    local,
                    "0.98 0.62 0.15 0.95",
                )
                ET.SubElement(tendon, "site", {"site": site_name})
                metadata_path.append(
                    {
                        "kind": "anchor",
                        "site": site_name,
                        "owner": owner,
                        "body": body_name,
                        "local_mm": local.tolist(),
                    }
                )
                continue

            if element["kind"] != "guide":
                if element["kind"] != "pass_through":
                    raise ValueError(element["kind"])
                owner = element["owner"]
                body_name = BODY_NAME[owner]
                if owner == "X":
                    parent = root.find("worldbody")
                else:
                    parent = find_body(root, body_name)
                site_name = f"V2_pass_{route_id}_{len(metadata_path)}"
                local = np.asarray(element["local_mm"], dtype=float)
                add_site(
                    parent,
                    site_name,
                    local,
                    "0.35 0.85 0.95 0.9",
                    size=0.002,
                )
                ET.SubElement(tendon, "site", {"site": site_name})
                metadata_path.append(
                    {
                        "kind": "pass_through",
                        "site": site_name,
                        "owner": owner,
                        "body": body_name,
                        "local_mm": local.tolist(),
                    }
                )
                continue
            owner = element["owner"]
            body_name = BODY_NAME[owner]
            if owner == "X":
                raise ValueError("Guides on external world anchors are unsupported")
            body = find_body(root, body_name)
            guide_name = f"V2_guide_{route_id}_{element['name']}"
            local = np.asarray(element["local_mm"], dtype=float)
            axis_design = normalize(np.asarray(element["axis_design"], dtype=float))
            side_design = normalize(np.asarray(element["side_design"], dtype=float))
            axis_mj = normalize(
                np.array([axis_design[0], -axis_design[2], axis_design[1]])
            )
            side_mj = normalize(
                np.array([side_design[0], -side_design[2], side_design[1]])
            )
            quat = quat_from_z_axis(axis_mj)
            radius = float(element["radius_mm"]) * 0.001
            ET.SubElement(
                body,
                "geom",
                {
                    "name": guide_name,
                    "type": "cylinder",
                    "size": f"{radius:.9f} 0.008",
                    "pos": vec(local),
                    "quat": vec(quat),
                    "material": "v2_guide_mat",
                    "contype": "0",
                    "conaffinity": "0",
                    "density": "0",
                    "group": "1",
                },
            )
            side_name = f"V2_side_{route_id}"
            side_local = local + side_mj * radius * 1.5
            add_site(
                body,
                side_name,
                side_local,
                "0.95 0.2 0.2 0.7",
                size=0.0015,
            )
            ET.SubElement(
                tendon,
                "geom",
                {"geom": guide_name, "sidesite": side_name},
            )
            metadata_path.append(
                {
                    "kind": "guide",
                    "geom": guide_name,
                    "sidesite": side_name,
                    "owner": owner,
                    "body": body_name,
                    "local_mm": local.tolist(),
                    "radius_mm": float(element["radius_mm"]),
                    "axis_design": axis_design.tolist(),
                    "wrap_angle_deg": float(element.get("wrap_angle_deg", 0.0)),
                }
            )
        route_metadata[route_id] = {
            "path": metadata_path,
            "guide_count": sum(
                element["kind"] == "guide" for element in route["path"]
            ),
        }

    actuator_root = root.find("actuator")
    sensor_root = root.find("sensor")
    if actuator_root is None:
        actuator_root = ET.SubElement(root, "actuator")
    if sensor_root is None:
        sensor_root = ET.SubElement(root, "sensor")
    for route in config["routes"]:
        route_id = route["id"]
        if any(
            actuator.get("name") == f"pull_{route_id}"
            for actuator in actuator_root.findall("position")
        ):
            continue
        ET.SubElement(
            actuator_root,
            "position",
            {
                "name": f"pull_{route_id}",
                "tendon": f"cable_{route_id}",
                "kp": str(config["simulation"].get("kp", 1200.0)),
                "kv": str(config["simulation"].get("kv", 80.0)),
                "ctrlrange": "0 1",
                "forcerange": "-160 0",
                "forcelimited": "true",
                "group": "2",
            },
        )
        ET.SubElement(
            sensor_root,
            "tendonpos",
            {"name": f"length_{route_id}", "tendon": f"cable_{route_id}"},
        )

    return tree, route_metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument("--out", type=Path, default=MODEL_PATH)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=CONFIG_PATH.parent / "route_metadata_v2.json",
    )
    args = parser.parse_args()
    config = load_json(args.config)
    tree, metadata = build(config)
    ET.indent(tree, space="  ")
    args.out.write_text(
        ET.tostring(tree.getroot(), encoding="unicode"),
        encoding="utf-8",
    )
    args.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out}")
    print(f"wrote {args.metadata}")


if __name__ == "__main__":
    main()
