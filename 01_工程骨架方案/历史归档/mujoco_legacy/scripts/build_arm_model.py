"""Build the right-arm MuJoCo model from the JSON design contracts."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np


ROOT = Path(__file__).resolve().parent
CONFIG_DIR = ROOT / "config"
DEFAULT_ARM_CONFIG = CONFIG_DIR / "arm_model.json"
DEFAULT_MUSCLE_CONFIG = CONFIG_DIR / "muscle_manifest.json"


def vec3(value) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (3,):
        raise ValueError(f"Expected a three-element vector, got {value!r}")
    return array


def fmt(value: float) -> str:
    return f"{value:.9g}"


def fmt_vec(value) -> str:
    return " ".join(fmt(float(x)) for x in value)


def quat_from_x_axis(axis: np.ndarray) -> np.ndarray:
    """Return a quaternion whose local +X axis equals the supplied vector."""

    x_axis = axis / np.linalg.norm(axis)
    reference = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(x_axis, reference))) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    y_axis = np.cross(reference, x_axis)
    y_axis /= np.linalg.norm(y_axis)
    z_axis = np.cross(x_axis, y_axis)
    rotation = np.column_stack((x_axis, y_axis, z_axis))
    quat = np.zeros(4, dtype=float)
    mujoco.mju_mat2Quat(quat, rotation.reshape(-1))
    return quat


class ArmBuilder:
    def __init__(self, arm_config_path: Path, muscle_config_path: Path) -> None:
        self.arm_config_path = arm_config_path
        self.muscle_config_path = muscle_config_path
        self.arm_config = json.loads(arm_config_path.read_text(encoding="utf-8"))
        self.muscle_config = json.loads(
            muscle_config_path.read_text(encoding="utf-8")
        )
        self.base_path = ROOT / self.arm_config["base_xml"]
        self.output_path = ROOT / self.arm_config["output_xml"]
        self.tree = ET.parse(self.base_path)
        self.root = self.tree.getroot()
        self.worldbody = self.root.find("worldbody")
        if self.worldbody is None:
            raise ValueError("The base XML has no <worldbody>.")

        self.base_model = mujoco.MjModel.from_xml_path(str(self.base_path))
        self.base_data = mujoco.MjData(self.base_model)
        mujoco.mj_forward(self.base_model, self.base_data)

        self.body_elements = self._index_body_elements()
        self.added_sites: set[str] = set()
        self.added_tendons: list[dict] = []
        self.added_actuators: list[str] = []

    def _index_body_elements(self) -> dict[str, ET.Element]:
        result = {}
        for body in self.worldbody.iter("body"):
            name = body.get("name")
            if name:
                result[name] = body
        return result

    def body_id(self, name: str) -> int:
        body_id = mujoco.mj_name2id(
            self.base_model, mujoco.mjtObj.mjOBJ_BODY, name
        )
        if body_id < 0:
            raise ValueError(f"Body not found in base XML: {name}")
        return int(body_id)

    def body_to_world(self, body_name: str, local: np.ndarray) -> np.ndarray:
        body_id = self.body_id(body_name)
        position = self.base_data.xpos[body_id]
        rotation = self.base_data.xmat[body_id].reshape(3, 3)
        return position + rotation @ local

    def world_to_body(self, body_name: str, world: np.ndarray) -> np.ndarray:
        body_id = self.body_id(body_name)
        position = self.base_data.xpos[body_id]
        rotation = self.base_data.xmat[body_id].reshape(3, 3)
        return rotation.T @ (world - position)

    def add_site(
        self,
        body_name: str,
        site_name: str,
        local: np.ndarray,
        rgba: tuple[float, float, float, float],
    ) -> None:
        if site_name in self.added_sites:
            raise ValueError(f"Duplicate site name: {site_name}")
        body = self.body_elements.get(body_name)
        if body is None:
            raise ValueError(f"No XML body element named {body_name}")
        ET.SubElement(
            body,
            "site",
            {
                "name": site_name,
                "pos": fmt_vec(local),
                "size": "0.004",
                "rgba": " ".join(fmt(x) for x in rgba),
                "group": "3",
            },
        )
        self.added_sites.add(site_name)

    def configure_root(self) -> None:
        compiler = self.root.find("compiler")
        if compiler is None:
            compiler = ET.SubElement(self.root, "compiler")
        compiler.set("angle", "degree")
        compiler.set("coordinate", "local")
        compiler.set("autolimits", "true")
        mesh_asset_dir = self.arm_config.get("mesh_asset_dir", "assets")
        compiler.set("meshdir", (ROOT / mesh_asset_dir).resolve().as_posix())

        option = self.root.find("option")
        if option is None:
            option = ET.SubElement(self.root, "option")
        option.set("timestep", fmt(float(self.arm_config["timestep"])))
        option.set("gravity", fmt_vec(self.arm_config["gravity"]))
        option.set("integrator", self.arm_config["integrator"])

        for tag in ("actuator", "tendon", "sensor"):
            existing = self.root.find(tag)
            if existing is not None:
                self.root.remove(existing)

    def add_ligaments(self) -> None:
        tendon_root = ET.SubElement(self.root, "tendon")
        for ligament in self.arm_config["ligaments"]:
            ligament_id = ligament["id"]
            origin_body = ligament["origin_body"]
            insertion_body = ligament["insertion_body"]
            origin_local = vec3(ligament["origin_local"])
            insertion_local = vec3(ligament["insertion_local"])
            self.add_site(
                origin_body,
                f"lig_{ligament_id}_o",
                origin_local,
                (0.20, 0.80, 0.95, 0.75),
            )
            self.add_site(
                insertion_body,
                f"lig_{ligament_id}_i",
                insertion_local,
                (0.20, 0.80, 0.95, 0.75),
            )
            origin_world = self.body_to_world(origin_body, origin_local)
            insertion_world = self.body_to_world(
                insertion_body, insertion_local
            )
            rest_length = float(np.linalg.norm(insertion_world - origin_world))
            spatial = ET.SubElement(
                tendon_root,
                "spatial",
                {
                    "name": f"ligament_{ligament_id}",
                    "width": fmt(ligament["width"]),
                    "stiffness": fmt(ligament["stiffness"]),
                    "damping": fmt(ligament["damping"]),
                    "springlength": fmt(rest_length),
                    "rgba": "0.2 0.8 0.95 0.65",
                },
            )
            ET.SubElement(spatial, "site", {"site": f"lig_{ligament_id}_o"})
            ET.SubElement(spatial, "site", {"site": f"lig_{ligament_id}_i"})
            self.added_tendons.append(
                {
                    "name": f"ligament_{ligament_id}",
                    "kind": "ligament",
                    "rest_length": rest_length,
                }
            )

    def add_muscles(self) -> None:
        tendon_root = self.root.find("tendon")
        if tendon_root is None:
            tendon_root = ET.SubElement(self.root, "tendon")
        actuator_root = ET.SubElement(self.root, "actuator")

        defaults = self.muscle_config["defaults"]
        for muscle in self.muscle_config["muscles"]:
            if not muscle.get("enabled", False):
                continue

            muscle_id = muscle["id"]
            origin_body = muscle["origin_body"]
            insertion_body = muscle["insertion_body"]
            origin_local = vec3(muscle["origin_local"])
            insertion_local = vec3(muscle["insertion_local"])
            origin_world = self.body_to_world(origin_body, origin_local)
            insertion_world = self.body_to_world(
                insertion_body, insertion_local
            )
            axis_vector = insertion_world - origin_world
            length = float(np.linalg.norm(axis_vector))
            if length < 0.005:
                raise ValueError(
                    f"Muscle {muscle_id} is too short: {length:.6f} m"
                )
            axis = axis_vector / length
            quat = quat_from_x_axis(axis)

            lattice = muscle.get("lattice", {})
            width = (
                float(lattice.get("y_mm", defaults["lattice_y_mm"]))
                * 1e-3
                * (int(defaults["transverse_nodes"]) - 1)
            )
            longitudinal_nodes = int(
                lattice.get("longitudinal_nodes", defaults["longitudinal_nodes"])
            )
            transverse_nodes = int(
                lattice.get("transverse_nodes", defaults["transverse_nodes"])
            )
            flex_name = f"psm_{muscle_id}"
            flex = ET.SubElement(
                self.worldbody,
                "flexcomp",
                {
                    "name": flex_name,
                    "type": "grid",
                    "dim": "2",
                    "count": (
                        f"{longitudinal_nodes} "
                        f"{transverse_nodes} 1"
                    ),
                    "spacing": (
                        f"{fmt(length / (longitudinal_nodes - 1))} "
                        f"{fmt(width / (transverse_nodes - 1))} 0.01"
                    ),
                    "radius": fmt(
                        float(
                            lattice.get(
                                "branch_diameter_mm",
                                defaults["branch_diameter_mm"],
                            )
                        )
                        * 0.0005
                    ),
                    "mass": fmt(
                        float(lattice.get("mass_kg", defaults["mass_kg"]))
                    ),
                    "pos": fmt_vec(origin_world),
                    "quat": fmt_vec(quat),
                    "rgba": "0.82 0.86 0.78 0.85",
                },
            )
            ET.SubElement(
                flex,
                "elasticity",
                {
                    "young": fmt(float(defaults["young_pa"])),
                    "poisson": fmt(
                        float(lattice.get("poisson", defaults["poisson"]))
                    ),
                    "damping": fmt(
                        float(lattice.get("damping", defaults["damping"]))
                    ),
                },
            )
            ET.SubElement(
                flex,
                "contact",
                {"contype": "0", "conaffinity": "0"},
            )

            # Pin both end faces to their bone bodies. The body ids are stable
            # because flex bodies are appended after all skeleton bodies.
            origin_body_id = self.body_id(origin_body)
            insertion_body_id = self.body_id(insertion_body)
            for transverse in range(transverse_nodes):
                ET.SubElement(
                    flex,
                    "pin",
                    {
                        "id": str(origin_body_id),
                        "grid": f"0 {transverse}",
                    },
                )
                ET.SubElement(
                    flex,
                    "pin",
                    {
                        "id": str(insertion_body_id),
                        "grid": f"{longitudinal_nodes - 1} {transverse}",
                    },
                )

            # Cable paths are generated parallel to the flexible body. The
            # equivalent elastic coupling comes from the pinned flex endpoints.
            perpendicular = np.cross(axis, np.array([0.0, 0.0, 1.0]))
            if np.linalg.norm(perpendicular) < 1e-8:
                perpendicular = np.cross(axis, np.array([0.0, 1.0, 0.0]))
            perpendicular /= np.linalg.norm(perpendicular)

            cable_count = int(muscle.get("cables", 1))
            for cable_index in range(cable_count):
                if cable_count == 1:
                    offset = 0.0
                else:
                    offset = width * 0.22 * (
                        cable_index - (cable_count - 1) / 2.0
                    )
                world_offset = perpendicular * offset
                site_names = []
                for suffix, body_name, local, world_point in (
                    (
                        "origin",
                        origin_body,
                        origin_local,
                        origin_world + world_offset,
                    ),
                    (
                        "insertion",
                        insertion_body,
                        insertion_local,
                        insertion_world + world_offset,
                    ),
                ):
                    local_point = self.world_to_body(body_name, world_point)
                    site_name = (
                        f"cable_{muscle_id}_{cable_index}_{suffix}"
                    )
                    self.add_site(
                        body_name,
                        site_name,
                        local_point,
                        (0.95, 0.38, 0.12, 0.85),
                    )
                    site_names.append(site_name)

                # Insert via sites between the endpoints in manifest order.
                path_items: list[tuple[str, str]] = []
                for via_index, via in enumerate(muscle.get("via", [])):
                    body_name = via["body"]
                    local = vec3(via["local"])
                    world_point = self.body_to_world(body_name, local)
                    local_point = self.world_to_body(
                        body_name, world_point + world_offset
                    )
                    site_name = (
                        f"cable_{muscle_id}_{cable_index}_via_{via_index}"
                    )
                    self.add_site(
                        body_name,
                        site_name,
                        local_point,
                        (0.95, 0.38, 0.12, 0.85),
                    )
                    path_items.append(("site", site_name))

                cable_name = f"cable_{muscle_id}_{cable_index}"
                spatial = ET.SubElement(
                    tendon_root,
                    "spatial",
                    {
                        "name": cable_name,
                        "width": "0.0015",
                        "rgba": "0.95 0.38 0.12 0.9",
                    },
                )
                ET.SubElement(spatial, "site", {"site": site_names[0]})
                for _, site_name in path_items:
                    ET.SubElement(spatial, "site", {"site": site_name})
                ET.SubElement(spatial, "site", {"site": site_names[1]})

                actuator_name = f"motor_{cable_name}"
                ET.SubElement(
                    actuator_root,
                    "motor",
                    {
                        "name": actuator_name,
                        "tendon": cable_name,
                        "gear": "1",
                        "ctrllimited": "true",
                        "ctrlrange": f"0 {fmt(float(muscle['max_force_n']))}",
                        "forcelimited": "true",
                        "forcerange": f"0 {fmt(float(muscle['max_force_n']))}",
                        "group": "2",
                    },
                )
                self.added_actuators.append(actuator_name)
                self.added_tendons.append(
                    {
                        "name": cable_name,
                        "kind": "muscle_cable",
                        "muscle": muscle_id,
                        "channel": cable_index,
                    }
                )

    def add_sensors(self) -> None:
        sensor_root = ET.SubElement(self.root, "sensor")
        for joint_id in range(self.base_model.njnt):
            name = mujoco.mj_id2name(
                self.base_model, mujoco.mjtObj.mjOBJ_JOINT, joint_id
            )
            if name:
                ET.SubElement(
                    sensor_root,
                    "jointpos",
                    {"name": f"sensor_{name}", "joint": name},
                )
        for tendon in self.added_tendons:
            ET.SubElement(
                sensor_root,
                "tendonpos",
                {"name": f"sensor_{tendon['name']}", "tendon": tendon["name"]},
            )
        for actuator_name in self.added_actuators:
            ET.SubElement(
                sensor_root,
                "actuatorfrc",
                {"name": f"sensor_{actuator_name}", "actuator": actuator_name},
            )

    def build(self) -> Path:
        self.configure_root()
        self.add_ligaments()
        self.add_muscles()
        self.add_sensors()
        ET.indent(self.tree, space="  ")
        self.tree.write(
            self.output_path,
            encoding="utf-8",
            xml_declaration=True,
        )
        return self.output_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--arm-config", type=Path, default=DEFAULT_ARM_CONFIG
    )
    parser.add_argument(
        "--muscle-config", type=Path, default=DEFAULT_MUSCLE_CONFIG
    )
    args = parser.parse_args()

    output = ArmBuilder(args.arm_config, args.muscle_config).build()
    model = mujoco.MjModel.from_xml_path(str(output))
    print(f"generated: {output}")
    print(f"bodies: {model.nbody}")
    print(f"joints: {model.njnt}")
    print(f"flex bodies: {model.nflex}")
    print(f"tendons: {model.ntendon}")
    print(f"actuators: {model.nu}")
    print(f"DOF: {model.nv}")


if __name__ == "__main__":
    main()
