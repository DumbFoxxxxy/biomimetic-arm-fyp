"""Build the interactive 7-DOF, 16-tendon demonstration model."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parent
DEFAULT_TSV = PROJECT_ROOT / "tmp" / "cable_anchors_final.tsv"
DEFAULT_JSON = ROOT / "config" / "tendon_demo_16.json"
DEFAULT_XML = ROOT / "arm_tendon_demo.xml"

GH_MM = np.array([159.06, 187.54, 7.30])
P_MM = np.array([62.18, 210.10, 91.39])
E_MM = np.array([184.34, -101.35, 7.30])
SC_MM = np.array([40.0, 235.0, -106.0])
AC_MM = np.array([156.86, 233.24, -11.97])
SUP_MM = np.array([72.76, 268.45, 82.21])
INF_MM = np.array([44.55, 112.87, 106.69])

A_MM = np.array([0.087155743, -0.996194698, 0.0])
B_MM = np.array([0.996194698, 0.087155743, 0.0])
C_MM = np.array([0.0, 0.0, 1.0])
U_MM = SUP_MM - INF_MM
U_MM /= np.linalg.norm(U_MM)
V0_MM = GH_MM - P_MM
V_MM = V0_MM - np.dot(V0_MM, U_MM) * U_MM
V_MM /= np.linalg.norm(V_MM)
W_MM = np.array([0.6556, 0.0, 0.7552])
E_SLIDE_MM = np.array([0.7552, 0.0, -0.6556])
E_SLIDE_MM /= np.linalg.norm(E_SLIDE_MM)


def mj(point_mm: np.ndarray) -> np.ndarray:
    """Map design x-right/y-up/z-back to MuJoCo x/y/z-up."""

    return np.array([point_mm[0], -point_mm[2], point_mm[1]]) * 0.001


def vec(values: np.ndarray) -> str:
    return " ".join(f"{value:.9f}" for value in values)


def quat_from_z_axis(axis: np.ndarray) -> tuple[float, float, float, float]:
    axis = axis / np.linalg.norm(axis)
    z_axis = np.array([0.0, 0.0, 1.0])
    cross = np.cross(z_axis, axis)
    dot = float(np.clip(np.dot(z_axis, axis), -1.0, 1.0))
    if np.linalg.norm(cross) < 1e-10:
        if dot > 0.0:
            return (1.0, 0.0, 0.0, 0.0)
        return (0.0, 1.0, 0.0, 0.0)
    angle = math.acos(dot)
    cross /= np.linalg.norm(cross)
    return (
        math.cos(angle / 2.0),
        float(cross[0] * math.sin(angle / 2.0)),
        float(cross[1] * math.sin(angle / 2.0)),
        float(cross[2] * math.sin(angle / 2.0)),
    )


def fmt3(values: np.ndarray) -> tuple[float, float, float]:
    return tuple(float(value) for value in values)


OWNERS = {
    "E01": ("H", "H", "F"),
    "E02": ("H", "H", "F"),
    "G01": ("S", "S", "H"),
    "G02": ("S", "S", "H"),
    "G03": ("S", "S", "H"),
    "G04": ("S", "S", "H"),
    "G05": ("S", "S", "H"),
    "G06": ("S", "S", "H"),
    "S01": ("T", None, "S"),
    "S02": ("T", None, "S"),
    "S03": ("T", None, "S"),
    "S04": ("T", None, "S"),
    "S05": ("X", None, "S"),
    "S06": ("X", None, "S"),
    "P01": ("C", "T", "H"),
    "P02": ("T", "T", "H"),
}


def load_routes(path: Path) -> list[dict]:
    routes = []
    lines = path.read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:
        if not line.strip():
            continue
        values = line.split("\t")
        route_id = values[0]
        a = np.array([float(values[1]), float(values[2]), float(values[3])])
        g = np.array([float(values[4]), float(values[5]), float(values[6])])
        b = np.array([float(values[7]), float(values[8]), float(values[9])])
        has_guide = bool(int(values[10]))
        owners = OWNERS[route_id]
        points = [a]
        point_owners = [owners[0]]
        if has_guide:
            points.append(g)
            point_owners.append(owners[1])
        points.append(b)
        point_owners.append(owners[2])
        routes.append(
            {
                "id": route_id,
                "owners": point_owners,
                "points_mm": [point.tolist() for point in points],
                "has_guide": has_guide,
                "neutral_length_m": float(
                    sum(
                        np.linalg.norm(points[index + 1] - points[index])
                        for index in range(len(points) - 1)
                    )
                    * 0.001
                ),
            }
        )
    if len(routes) != 16:
        raise ValueError(f"Expected 16 routes, found {len(routes)}")
    return routes


def local_point(owner: str, point_mm: np.ndarray) -> np.ndarray:
    world = mj(point_mm)
    if owner == "T":
        return world
    if owner == "X":
        return world
    if owner == "S":
        return world - mj(AC_MM)
    if owner == "C":
        return world - mj(SC_MM)
    if owner == "H":
        return world - mj(GH_MM)
    if owner == "F":
        return world - mj(E_MM)
    raise ValueError(owner)


def add_site(parent: ET.Element, name: str, pos: np.ndarray, rgba: str) -> None:
    ET.SubElement(
        parent,
        "site",
        {
            "name": name,
            "pos": vec(pos),
            "size": "0.0035",
            "rgba": rgba,
            "group": "3",
        },
    )


def build_model(routes: list[dict], max_pull_m: float) -> ET.ElementTree:
    root = ET.Element("mujoco", {"model": "tendon_demo_16"})
    ET.SubElement(
        root,
        "compiler",
        {
            "angle": "degree",
            "autolimits": "true",
            "meshdir": "D:/XJTLU/PSP/mujoco_arm/assets_tendon_demo",
        },
    )
    ET.SubElement(
        root,
        "option",
        {
            "timestep": "0.001",
            "gravity": "0 0 0",
            "integrator": "implicitfast",
        },
    )
    visual = ET.SubElement(root, "visual")
    ET.SubElement(
        visual,
        "global",
        {
            "azimuth": "45",
            "elevation": "-16",
        },
    )
    ET.SubElement(
        visual,
        "headlight",
        {"diffuse": "0.8 0.8 0.8", "ambient": "0.3 0.3 0.3"},
    )

    asset = ET.SubElement(root, "asset")
    for mesh_name in ("thorax", "clavicle", "scapula", "humerus", "ulna"):
        ET.SubElement(
            asset,
            "mesh",
            {"name": mesh_name, "file": f"{mesh_name}.stl", "scale": "1 1 1"},
        )
    for name, rgba in (
        ("thorax_mat", "0.45 0.55 0.65 0.55"),
        ("scapula_mat", "0.30 0.72 0.52 1"),
        ("humerus_mat", "0.85 0.72 0.30 1"),
        ("ulna_mat", "0.35 0.58 0.88 1"),
        ("joint_mat", "0.95 0.35 0.18 1"),
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
    ET.SubElement(
        thorax,
        "geom",
        {
            "name": "thorax_geom",
            "type": "mesh",
            "mesh": "thorax",
            "material": "thorax_mat",
            "contype": "0",
            "conaffinity": "0",
        },
    )
    ET.SubElement(
        thorax,
        "geom",
        {
            "name": "thorax_collision",
            "type": "ellipsoid",
            "pos": "0 0 0.14",
            "size": "0.135 0.100 0.135",
            "rgba": "0.5 0.5 0.5 0",
            "contype": "1",
            "conaffinity": "6",
            "group": "0",
            "density": "0",
        },
    )

    # One fixed site per external machine anchor.
    for route in routes:
        for index, owner in enumerate(route["owners"]):
            if owner == "X":
                add_site(
                    worldbody,
                    f"{'AGB'[index]}_{route['id']}",
                    mj(np.asarray(route["points_mm"][index])),
                    "0.95 0.55 0.15 0.9",
                )

    sc_local = mj(SC_MM)
    ac_from_sc = mj(AC_MM) - sc_local
    p_from_sc = mj(P_MM) - sc_local
    p_from_ac = mj(P_MM) - mj(AC_MM)

    clavicle_elev = ET.SubElement(
        thorax, "body", {"name": "clavicle_elevation", "pos": vec(sc_local)}
    )
    ET.SubElement(
        clavicle_elev,
        "joint",
        {
            "name": "sc_elevation",
            "type": "hinge",
            "axis": vec(mj(C_MM)),
            "range": "-5 30",
            "damping": "0.4",
            "armature": "0.005",
        },
    )
    ET.SubElement(
        clavicle_elev,
        "inertial",
        {"pos": "0 0 0", "mass": "0.005", "diaginertia": "1e-6 1e-6 1e-6"},
    )
    clavicle_protract = ET.SubElement(
        clavicle_elev, "body", {"name": "clavicle_protraction", "pos": "0 0 0"}
    )
    ET.SubElement(
        clavicle_protract,
        "joint",
        {
            "name": "sc_protraction",
            "type": "hinge",
            "axis": "0 0 1",
            "range": "-15 15",
            "damping": "0.4",
            "armature": "0.005",
        },
    )
    ET.SubElement(
        clavicle_protract,
        "inertial",
        {"pos": "0 0 0", "mass": "0.005", "diaginertia": "1e-6 1e-6 1e-6"},
    )
    clavicle_body = ET.SubElement(
        clavicle_protract, "body", {"name": "clavicle_body", "pos": "0 0 0"}
    )
    ET.SubElement(
        clavicle_body,
        "geom",
        {
            "name": "clavicle_visual",
            "type": "mesh",
            "mesh": "clavicle",
            "pos": vec(p_from_sc),
            "material": "scapula_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "700",
        },
    )
    ET.SubElement(
        clavicle_body,
        "geom",
        {
            "name": "clavicle_collision",
            "type": "mesh",
            "mesh": "clavicle",
            "pos": vec(p_from_sc),
            "rgba": "0.8 0.2 0.2 0",
            "contype": "2",
            "conaffinity": "1",
            "group": "0",
            "density": "0",
        },
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
            "axis": vec(mj(U_MM)),
            "range": "-15 15",
            "damping": "0.5",
            "armature": "0.01",
        },
    )
    ET.SubElement(
        scapula,
        "geom",
        {
            "name": "scapula_mesh",
            "type": "mesh",
            "mesh": "scapula",
            "pos": vec(p_from_ac),
            "material": "scapula_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "700",
        },
    )
    ET.SubElement(
        scapula,
        "geom",
        {
            "name": "scapula_collision",
            "type": "mesh",
            "mesh": "scapula",
            "pos": vec(p_from_ac),
            "rgba": "0.8 0.2 0.2 0",
            "contype": "4",
            "conaffinity": "1",
            "group": "0",
            "density": "0",
        },
    )
    gh_local = mj(GH_MM) - mj(AC_MM)

    gh_flex = ET.SubElement(
        scapula, "body", {"name": "gh_flex", "pos": vec(gh_local)}
    )
    ET.SubElement(
        gh_flex,
        "joint",
        {
            "name": "gh_flex",
            "type": "hinge",
            "axis": vec(mj(B_MM)),
            "range": "-150 150",
            "damping": "0.5",
            "armature": "0.008",
        },
    )
    ET.SubElement(
        gh_flex,
        "inertial",
        {"pos": "0 0 0", "mass": "0.004", "diaginertia": "1e-6 1e-6 1e-6"},
    )
    gh_abd = ET.SubElement(gh_flex, "body", {"name": "gh_abd", "pos": "0 0 0"})
    ET.SubElement(
        gh_abd,
        "joint",
        {
            "name": "gh_abd",
            "type": "hinge",
            "axis": vec(mj(C_MM)),
            "range": "-150 150",
            "damping": "0.5",
            "armature": "0.008",
        },
    )
    ET.SubElement(
        gh_abd,
        "inertial",
        {"pos": "0 0 0", "mass": "0.004", "diaginertia": "1e-6 1e-6 1e-6"},
    )
    gh_rot = ET.SubElement(gh_abd, "body", {"name": "gh_rot", "pos": "0 0 0"})
    ET.SubElement(
        gh_rot,
        "joint",
        {
            "name": "gh_rot",
            "type": "hinge",
            "axis": vec(mj(A_MM)),
            "range": "-120 120",
            "damping": "0.5",
            "armature": "0.008",
        },
    )
    ET.SubElement(
        gh_rot,
        "inertial",
        {"pos": "0 0 0", "mass": "0.004", "diaginertia": "1e-6 1e-6 1e-6"},
    )
    humerus = ET.SubElement(gh_rot, "body", {"name": "humerus", "pos": "0 0 0"})
    ET.SubElement(
        humerus,
        "geom",
        {
            "name": "humerus_shaft",
            "type": "mesh",
            "mesh": "humerus",
            "material": "humerus_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "900",
        },
    )
    ET.SubElement(
        humerus,
        "geom",
        {
            "name": "gh_ball",
            "type": "sphere",
            "size": "0.014",
            "material": "joint_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "0",
        },
    )

    ulna = ET.SubElement(
        humerus, "body", {"name": "ulna", "pos": vec(mj(E_MM) - mj(GH_MM))}
    )
    ET.SubElement(
        ulna,
        "joint",
        {
            "name": "elbow",
            "type": "hinge",
            "axis": vec(mj(B_MM)),
            "range": "0 145",
            "damping": "0.2",
            "armature": "0.005",
        },
    )
    ET.SubElement(
        ulna,
        "geom",
        {
            "name": "ulna_shaft",
            "type": "mesh",
            "mesh": "ulna",
            "material": "ulna_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "800",
        },
    )
    ET.SubElement(
        ulna,
        "geom",
        {
            "name": "elbow_joint",
            "type": "sphere",
            "size": "0.012",
            "material": "joint_mat",
            "contype": "0",
            "conaffinity": "0",
            "density": "0",
        },
    )

    body_by_owner = {
        "T": thorax,
        "S": scapula,
        "C": clavicle_body,
        "H": humerus,
        "F": ulna,
    }
    for route in routes:
        labels = ("A", "G", "B") if route["has_guide"] else ("A", "B")
        for index, (owner, point_mm) in enumerate(
            zip(route["owners"], route["points_mm"])
        ):
            if owner == "X":
                continue
            site_name = f"{labels[index]}_{route['id']}"
            add_site(
                body_by_owner[owner],
                site_name,
                local_point(owner, np.asarray(point_mm)),
                "0.95 0.55 0.15 0.9",
            )

    tendon_root = ET.SubElement(root, "tendon")
    actuator_root = ET.SubElement(root, "actuator")
    sensor_root = ET.SubElement(root, "sensor")
    neutral = {}
    for route in routes:
        route_id = route["id"]
        labels = ("A", "G", "B") if route["has_guide"] else ("A", "B")
        point_names = [f"{labels[index]}_{route_id}" for index in range(len(route["owners"]))]
        if route["has_guide"]:
            path_names = [point_names[0], point_names[1], point_names[2]]
        else:
            path_names = [point_names[0], point_names[1]]
        spatial = ET.SubElement(
            tendon_root,
            "spatial",
            {
                "name": f"cable_{route_id}",
                "width": "0.0018",
                "rgba": "0.95 0.45 0.12 0.95",
            },
        )
        for site_name in path_names:
            ET.SubElement(spatial, "site", {"site": site_name})
        neutral_length = float(route["neutral_length_m"])
        neutral[route_id] = neutral_length
        ET.SubElement(
            actuator_root,
            "position",
            {
                "name": f"pull_{route_id}",
                "tendon": f"cable_{route_id}",
                "kp": "1200",
                "kv": "80",
                "ctrlrange": f"{max(0.0, neutral_length - max_pull_m):.6f} {neutral_length + 0.01:.6f}",
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
    for joint_name in (
        "sc_elevation",
        "sc_protraction",
        "gh_flex",
        "gh_abd",
        "gh_rot",
        "elbow",
    ):
        ET.SubElement(
            sensor_root,
            "jointpos",
            {"name": f"angle_{joint_name}", "joint": joint_name},
        )

    return ET.ElementTree(root)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", type=Path, default=DEFAULT_TSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--xml", type=Path, default=DEFAULT_XML)
    parser.add_argument("--max-pull-mm", type=float, default=80.0)
    args = parser.parse_args()

    routes = load_routes(args.tsv)
    tree = build_model(routes, args.max_pull_mm * 0.001)
    ET.indent(tree, space="  ")
    tree.write(args.xml, encoding="utf-8", xml_declaration=True)
    args.json.write_text(
        json.dumps(
            {
                "routes": routes,
                "joints": {
                    "sc_elevation": {"range_deg": [-5.0, 30.0]},
                    "sc_protraction": {"range_deg": [-15.0, 15.0]},
                    "ac_ext": {"range_deg": [-15.0, 15.0]},
                    "gh_flex": {"range_deg": [-150.0, 150.0]},
                    "gh_abd": {"range_deg": [-150.0, 150.0]},
                    "gh_rot": {"range_deg": [-120.0, 120.0]},
                    "elbow": {"range_deg": [0.0, 145.0]},
                },
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.xml}")
    print(f"wrote {args.json}")
    for route in routes:
        print(
            f"{route['id']:3s} L0={route['neutral_length_m']*1000:7.2f} mm "
            f"path={'A-G-B' if route['has_guide'] else 'A-B'}"
        )


if __name__ == "__main__":
    main()
