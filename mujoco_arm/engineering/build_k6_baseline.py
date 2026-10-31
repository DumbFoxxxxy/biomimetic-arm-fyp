"""Build the isolated K6/S01 quasi-static MuJoCo baseline."""

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

from k6_common import (
    DEFAULT_ASSETS,
    DEFAULT_CONFIG,
    DEFAULT_XML,
    MUJOCO_ROOT,
    load_config,
    normalize,
    quat_from_z_axis,
    resolve_s01_route,
    scapula_basis,
    to_mujoco_design,
    vec,
)


VISUAL_DIR = MUJOCO_ROOT / "assets_tendon_demo"


def add_mesh(asset: ET.Element, name: str, path: Path) -> None:
    ET.SubElement(
        asset,
        "mesh",
        {"name": name, "file": str(path.resolve()).replace("\\", "/")},
    )


def add_site(
    parent: ET.Element,
    name: str,
    position_mj: np.ndarray,
    rgba: str,
    size: float = 0.0035,
) -> None:
    ET.SubElement(
        parent,
        "site",
        {
            "name": name,
            "pos": vec(position_mj),
            "size": f"{size:.6f}",
            "rgba": rgba,
            "group": "3",
        },
    )


def add_visual_geom(
    parent: ET.Element,
    name: str,
    mesh: str,
    material: str,
    position_mj: np.ndarray,
) -> None:
    ET.SubElement(
        parent,
        "geom",
        {
            "name": name,
            "type": "mesh",
            "mesh": mesh,
            "pos": vec(position_mj),
            "material": material,
            "contype": "0",
            "conaffinity": "0",
            "density": "0",
        },
    )


def build_model(config: dict, route_override: dict | None = None) -> tuple[ET.ElementTree, dict]:
    points = {key: np.asarray(value, dtype=float) for key, value in config["points_mm"].items()}
    route = resolve_s01_route(config, route_override)
    a_design = np.asarray(route["A_mm"])
    g_design = np.asarray(route["G_mm"])
    b_design = np.asarray(route["B_mm"])
    sc = points["SC"]
    ac = points["AC"]
    gh = points["GH"]
    p = points["P"]
    sup = points["SUP"]
    inf = points["INF"]
    u_axis, _, _ = scapula_basis(config)

    incoming = g_design - a_design
    outgoing = b_design - g_design
    axis_design = normalize(np.cross(incoming, outgoing))
    side_design = normalize(
        np.cross(axis_design, outgoing + incoming)
    )
    sign = float(route["sidesite_sign"])
    radius_mm = float(config["s01"]["guide_radius_mm"])
    side_design = sign * side_design

    root = ET.Element("mujoco", {"model": config["model_name"]})
    ET.SubElement(
        root,
        "compiler",
        {"angle": "degree", "autolimits": "true"},
    )
    ET.SubElement(
        root,
        "option",
        {
            "timestep": str(config["simulation"]["timestep_s"]),
            "gravity": "0 0 0",
            "integrator": "implicitfast",
        },
    )
    visual = ET.SubElement(root, "visual")
    ET.SubElement(
        visual,
        "global",
        {"azimuth": "140", "elevation": "-12"},
    )
    ET.SubElement(
        visual,
        "headlight",
        {"diffuse": "0.85 0.85 0.85", "ambient": "0.35 0.35 0.35"},
    )

    asset = ET.SubElement(root, "asset")
    for name in ("thorax", "clavicle", "scapula", "humerus"):
        add_mesh(asset, f"{name}_visual", VISUAL_DIR / f"{name}.stl")
    add_mesh(asset, "thorax_contact_band", DEFAULT_ASSETS / "thorax_contact_band.stl")
    add_mesh(asset, "scapula_contact_pad", DEFAULT_ASSETS / "scapula_contact_pad.stl")
    for name, rgba in (
        ("thorax_mat", "0.46 0.56 0.66 0.28"),
        ("band_mat", "0.20 0.70 0.75 0.45"),
        ("pad_mat", "0.95 0.50 0.16 0.72"),
        ("scapula_mat", "0.28 0.68 0.46 1"),
        ("clavicle_mat", "0.80 0.72 0.58 1"),
        ("humerus_mat", "0.86 0.73 0.30 1"),
        ("guide_mat", "0.92 0.30 0.18 1"),
    ):
        ET.SubElement(asset, "material", {"name": name, "rgba": rgba})

    worldbody = ET.SubElement(root, "worldbody")
    ET.SubElement(
        worldbody,
        "light",
        {
            "name": "key",
            "pos": "-0.4 -0.5 1.5",
            "dir": "0.4 0.3 -1",
            "directional": "true",
        },
    )
    ET.SubElement(
        worldbody,
        "geom",
        {
            "name": "floor",
            "type": "plane",
            "size": "1.5 1.5 0.1",
            "pos": "0 0 -0.45",
            "rgba": "0.12 0.14 0.17 1",
            "contype": "0",
            "conaffinity": "0",
        },
    )

    thorax = ET.SubElement(worldbody, "body", {"name": "thorax", "pos": "0 0 0"})
    add_visual_geom(thorax, "thorax_visual", "thorax_visual", "thorax_mat", np.zeros(3))
    ET.SubElement(
        thorax,
        "geom",
        {
            "name": "thorax_contact_band",
            "type": "mesh",
            "mesh": "thorax_contact_band",
            "material": "band_mat",
            "contype": "1",
            "conaffinity": "2",
            "condim": str(config["contact"]["condim"]),
            "friction": vec(config["contact"]["friction"]),
            "solref": vec(config["contact"]["solref"]),
            "solimp": vec(config["contact"]["solimp"]),
            "density": "0",
            "group": "1",
        },
    )

    guide_position_mj = to_mujoco_design(g_design)
    guide_axis_mj = normalize(
        np.array([axis_design[0], -axis_design[2], axis_design[1]])
    )
    guide_quat = quat_from_z_axis(guide_axis_mj)
    ET.SubElement(
        thorax,
        "geom",
        {
            "name": "s01_guide",
            "type": "cylinder",
            "size": (
                f"{radius_mm * 0.001:.9f} "
                f"{float(config['s01']['guide_half_length_mm']) * 0.001:.9f}"
            ),
            "pos": vec(guide_position_mj),
            "quat": vec(guide_quat),
            "material": "guide_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "0",
            "group": "1",
        },
    )
    side_point = g_design + side_design * radius_mm * 1.5
    add_site(
        thorax,
        "guide_side_S01",
        to_mujoco_design(side_point),
        "1 0.2 0.2 0.8",
        size=0.002,
    )
    add_site(
        thorax,
        "A_S01",
        to_mujoco_design(a_design),
        "0.95 0.55 0.15 0.9",
    )

    sc_local = to_mujoco_design(sc)
    ac_from_sc = to_mujoco_design(ac) - sc_local
    clavicle_mesh_from_sc = to_mujoco_design(p) - sc_local

    clavicle_elevation = ET.SubElement(
        thorax, "body", {"name": "clavicle_elevation", "pos": vec(sc_local)}
    )
    ET.SubElement(
        clavicle_elevation,
        "joint",
        {
            "name": "sc_elevation",
            "type": "hinge",
            "axis": vec(to_mujoco_design(np.array([0.0, 0.0, 1.0]))),
            "range": vec(config["joint_limits_deg"]["sc_elevation"]),
            "damping": "0.6",
            "armature": "0.005",
        },
    )
    ET.SubElement(
        clavicle_elevation,
        "inertial",
        {"pos": "0 0 0", "mass": "0.03", "diaginertia": "2e-5 2e-5 2e-5"},
    )
    clavicle_protraction = ET.SubElement(
        clavicle_elevation, "body", {"name": "clavicle_protraction", "pos": "0 0 0"}
    )
    ET.SubElement(
        clavicle_protraction,
        "joint",
        {
            "name": "sc_protraction",
            "type": "hinge",
            "axis": vec(to_mujoco_design(np.array([0.0, 1.0, 0.0]))),
            "range": vec(config["joint_limits_deg"]["sc_protraction"]),
            "damping": "0.6",
            "armature": "0.005",
        },
    )
    ET.SubElement(
        clavicle_protraction,
        "inertial",
        {"pos": "0 0 0", "mass": "0.03", "diaginertia": "2e-5 2e-5 2e-5"},
    )
    clavicle_body = ET.SubElement(
        clavicle_protraction, "body", {"name": "clavicle_body", "pos": "0 0 0"}
    )
    add_visual_geom(
        clavicle_body,
        "clavicle_visual",
        "clavicle_visual",
        "clavicle_mat",
        clavicle_mesh_from_sc,
    )

    scapula = ET.SubElement(
        clavicle_body, "body", {"name": "scapula", "pos": vec(ac_from_sc)}
    )
    ET.SubElement(
        scapula,
        "joint",
        {
            "name": "ac_ext",
            "type": "hinge",
            "axis": vec(to_mujoco_design(u_axis)),
            "range": vec(config["joint_limits_deg"]["ac_ext"]),
            "damping": "0.6",
            "armature": "0.01",
        },
    )
    ET.SubElement(
        scapula,
        "inertial",
        {"pos": "0 0 0", "mass": "0.08", "diaginertia": "8e-5 8e-5 8e-5"},
    )
    scapula_mesh_from_ac = to_mujoco_design(p) - to_mujoco_design(ac)
    add_visual_geom(
        scapula,
        "scapula_visual",
        "scapula_visual",
        "scapula_mat",
        scapula_mesh_from_ac,
    )
    ET.SubElement(
        scapula,
        "geom",
        {
            "name": "scapula_contact_pad",
            "type": "mesh",
            "mesh": "scapula_contact_pad",
            "material": "pad_mat",
            "contype": "2",
            "conaffinity": "1",
            "condim": str(config["contact"]["condim"]),
            "friction": vec(config["contact"]["friction"]),
            "solref": vec(config["contact"]["solref"]),
            "solimp": vec(config["contact"]["solimp"]),
            "density": "0",
            "group": "1",
        },
    )
    add_site(
        scapula,
        "B_S01",
        to_mujoco_design(b_design) - to_mujoco_design(ac),
        "0.95 0.55 0.15 0.9",
    )

    gh_local = to_mujoco_design(gh) - to_mujoco_design(ac)
    gh_flex = ET.SubElement(scapula, "body", {"name": "gh_flex", "pos": vec(gh_local)})
    ET.SubElement(
        gh_flex,
        "joint",
        {
            "name": "gh_flex",
            "type": "hinge",
            "axis": vec(to_mujoco_design(np.array([0.996194698, 0.087155743, 0.0]))),
            "range": "-120 120",
            "damping": "0.6",
        },
    )
    ET.SubElement(
        gh_flex,
        "inertial",
        {"pos": "0 0 0", "mass": "0.22", "diaginertia": "2e-4 2e-4 2e-4"},
    )
    add_visual_geom(gh_flex, "humerus_visual", "humerus_visual", "humerus_mat", np.zeros(3))

    tendon = ET.SubElement(root, "tendon")
    spatial = ET.SubElement(
        tendon,
        "spatial",
        {
            "name": "cable_S01",
            "width": "0.0012",
            "rgba": "0.95 0.45 0.12 1",
        },
    )
    ET.SubElement(spatial, "site", {"site": "A_S01"})
    ET.SubElement(
        spatial,
        "geom",
        {"geom": "s01_guide", "sidesite": "guide_side_S01"},
    )
    ET.SubElement(spatial, "site", {"site": "B_S01"})

    actuator = ET.SubElement(root, "actuator")
    ET.SubElement(
        actuator,
        "position",
        {
            "name": "pull_S01",
            "tendon": "cable_S01",
            "kp": str(config["simulation"]["kp"]),
            "kv": str(config["simulation"]["kv"]),
            "ctrlrange": "0 1",
            "forcerange": "-250 0",
            "forcelimited": "true",
        },
    )

    sensor = ET.SubElement(root, "sensor")
    ET.SubElement(sensor, "tendonpos", {"name": "length_S01", "tendon": "cable_S01"})
    for joint_name in ("sc_elevation", "sc_protraction", "ac_ext"):
        ET.SubElement(
            sensor,
            "jointpos",
            {"name": f"angle_{joint_name}", "joint": joint_name},
        )

    route_metadata = {
        "A_mm": a_design.tolist(),
        "G_mm": g_design.tolist(),
        "B_mm": b_design.tolist(),
        "B_local_uvw_mm": route["B_local_uvw_mm"],
        "guide_axis_design": axis_design.tolist(),
        "guide_side_design": side_design.tolist(),
        "guide_radius_mm": radius_mm,
        "sidesite_sign": sign,
        "input": route["input"],
    }
    return ET.ElementTree(root), route_metadata


def write_model(
    config_path: Path,
    output_path: Path,
    route_override_path: Path | None,
    route_output: Path,
) -> None:
    config = load_config(config_path)
    route_override = None
    if route_override_path is not None:
        route_override = json.loads(route_override_path.read_text(encoding="utf-8"))
    tree, route_metadata = build_model(config, route_override)
    ET.indent(tree, space="  ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    route_output.write_text(
        json.dumps(route_metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {output_path}")
    print(f"wrote {route_output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--out", type=Path, default=DEFAULT_XML)
    parser.add_argument("--route", type=Path, default=None)
    parser.add_argument(
        "--route-out",
        type=Path,
        default=DEFAULT_ASSETS / "route_metadata.json",
    )
    args = parser.parse_args()
    write_model(args.config, args.out, args.route, args.route_out)


if __name__ == "__main__":
    main()
