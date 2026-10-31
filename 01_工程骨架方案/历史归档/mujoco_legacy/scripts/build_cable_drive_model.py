"""Build the five-bone shoulder/elbow cable-drive MuJoCo model."""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from cable_drive_common import (
    BONE_NAMES,
    DEFAULT_EXPORT_DIR,
    DEFAULT_MANIFEST,
    DEFAULT_OUTPUT_XML,
    JOINT_SPECS,
    RHYTHM_ACTUATOR_TARGETS,
    SKELETON_JOINT_NAMES,
    SLIDER_ACTUATOR_PREFIX,
    SLIDER_JOINT_PREFIX,
    load_cable_manifest,
    load_bone_transforms,
    load_joint_frames,
    cable_guides,
    normalize,
    read_binary_stl,
    task_envelope_states,
    triangle_centroid_and_box,
    vec3,
)


def fmt(value: float) -> str:
    return f"{value:.9g}"


def fmt_vec(values) -> str:
    return " ".join(fmt(float(value)) for value in values)


def diagonal_inertia(mass: float, box: np.ndarray) -> np.ndarray:
    x, y, z = (float(value) for value in box)
    return mass / 12.0 * np.array(
        [y * y + z * z, x * x + z * z, x * x + y * y]
    )


def detect_mesh_scale(path: Path) -> float:
    triangles = read_binary_stl(path)
    extent = np.ptp(triangles.reshape(-1, 3), axis=0)
    return 0.001 if float(np.max(extent)) > 5.0 else 1.0


def make_default(compiler_angle: str = "degree") -> ET.Element:
    default = ET.Element("default")
    visual = ET.SubElement(default, "default", {"class": "visual_mesh"})
    ET.SubElement(
        visual,
        "geom",
        {
            "type": "mesh",
            "contype": "0",
            "conaffinity": "0",
            "group": "1",
            "density": "0",
            "rgba": "0.80 0.78 0.70 1",
        },
    )
    arm_joint = ET.SubElement(default, "default", {"class": "arm_joint"})
    ET.SubElement(
        arm_joint,
        "joint",
        {
            "type": "hinge",
            "limited": "true",
            "damping": "0.25",
            "armature": "0.02",
        },
    )
    slider_joint = ET.SubElement(default, "default", {"class": "slider_joint"})
    ET.SubElement(
        slider_joint,
        "joint",
        {
            "type": "slide",
            "limited": "true",
            "damping": "0.05",
            "armature": "0.001",
        },
    )
    slider_actuator = ET.SubElement(
        default, "default", {"class": "slider_actuator"}
    )
    ET.SubElement(slider_actuator, "position", {})
    rhythm_actuator = ET.SubElement(
        default, "default", {"class": "rhythm_actuator"}
    )
    ET.SubElement(rhythm_actuator, "position", {"kp": "5000", "kv": "200"})
    return default


class CableDriveBuilder:
    def __init__(
        self,
        export_dir: Path,
        joint_frames_path: Path,
        manifest_path: Path,
        output_path: Path,
    ) -> None:
        self.export_dir = export_dir.resolve()
        self.bone_dir = self.export_dir / "bones"
        self.manifest = load_cable_manifest(manifest_path)
        self.frames = load_joint_frames(joint_frames_path)
        self.bone_transforms = load_bone_transforms(
            export_dir / "bone_transforms.json"
        )
        self.body_bone = {"torso": "thorax"}
        for spec in JOINT_SPECS:
            if spec["mesh"]:
                self.body_bone[spec["body"]] = spec["mesh"]
        self.output_path = output_path.resolve()
        self.defaults = self.manifest["defaults"]
        self.root = ET.Element("mujoco", {"model": self.manifest["model_name"]})
        self.body_elements: dict[str, ET.Element] = {}
        self.body_origin_global: dict[str, np.ndarray] = {}
        self.body_mass: dict[str, float] = {}
        self.cable_channels: list[dict] = []
        self.tendon_elements: dict[str, ET.Element] = {}
        self.tendon_stiffness: dict[str, float] = {}
        self.joint_stiffness: dict[str, float] = {
            str(name): float(value)
            for name, value in (
                self.defaults.get("joint_stiffness_nm_per_rad") or {}
            ).items()
        }

    def build(self) -> Path:
        self._add_top_level()
        self._add_assets()
        self._add_worldbody()
        self._add_cables()
        self._add_sensors()
        ET.indent(self.root, space="  ")
        ET.ElementTree(self.root).write(
            self.output_path, encoding="utf-8", xml_declaration=True
        )
        self._calibrate_springlengths()
        return self.output_path

    def _add_top_level(self) -> None:
        ET.SubElement(
            self.root,
            "compiler",
            {
                "angle": "degree",
                "coordinate": "local",
                "autolimits": "true",
                "meshdir": self.bone_dir.as_posix(),
            },
        )
        ET.SubElement(
            self.root,
            "option",
            {
                "timestep": fmt(float(self.manifest.get("timestep", 0.001))),
                "gravity": fmt_vec(self.manifest.get("gravity", [0, 0, -9.81])),
                "integrator": self.manifest.get("integrator", "implicitfast"),
            },
        )
        visual = ET.SubElement(self.root, "visual")
        ET.SubElement(
            visual,
            "headlight",
            {"diffuse": "0.8 0.8 0.8", "ambient": "0.25 0.25 0.25"},
        )
        ET.SubElement(visual, "rgba", {"haze": "0.15 0.2 0.25 1"})
        ET.SubElement(
            visual,
            "global",
            {"azimuth": "135", "elevation": "-18"},
        )
        self.root.append(make_default())

    def _add_assets(self) -> None:
        asset = ET.SubElement(self.root, "asset")
        for bone_name in BONE_NAMES:
            path = self.bone_dir / f"{bone_name}.stl"
            if not path.exists():
                raise FileNotFoundError(
                    f"Missing global bone STL: {path}. Run "
                    "prepare_cable_drive_exports.py first."
                )
            ET.SubElement(
                asset,
                "mesh",
                {
                    "name": bone_name,
                    "file": f"{bone_name}.stl",
                    "scale": fmt_vec([detect_mesh_scale(path)] * 3),
                },
            )
        materials = {
            "clavicle_mat": "0.78 0.76 0.68 1",
            "scapula_mat": "0.72 0.75 0.69 1",
            "humerus_mat": "0.82 0.80 0.72 1",
            "ulna_mat": "0.70 0.72 0.67 1",
            "torso_mat": "0.48 0.52 0.55 1",
            "slider_mat": "0.85 0.35 0.18 1",
            "joint_mat": "0.95 0.55 0.12 1",
        }
        for name, rgba in materials.items():
            ET.SubElement(asset, "material", {"name": name, "rgba": rgba})

    def _add_worldbody(self) -> None:
        worldbody = ET.SubElement(self.root, "worldbody")
        ET.SubElement(
            worldbody,
            "light",
            {
                "name": "key_light",
                "pos": "-0.35 -0.45 1.80",
                "dir": "0.35 0.25 -1",
                "directional": "true",
                "diffuse": "0.9 0.9 0.9",
                "specular": "0.2 0.2 0.2",
            },
        )
        ET.SubElement(
            worldbody,
            "geom",
            {
                "name": "floor",
                "type": "plane",
                "size": "2 2 0.1",
                "pos": "0 0 0",
                "material": "torso_mat",
                "contype": "0",
                "conaffinity": "0",
                "group": "0",
            },
        )

        sc_origin = self.frames["SC_J"]["origin_m"]
        torso = ET.SubElement(
            worldbody,
            "body",
            {"name": "torso", "pos": fmt_vec(sc_origin)},
        )
        self.body_elements["torso"] = torso
        self.body_origin_global["torso"] = sc_origin
        self._add_mesh_geom(torso, "thorax", sc_origin, "torso_mat", 5.0)
        self.body_mass["torso"] = 5.0

        for spec in JOINT_SPECS:
            parent = self.body_elements[spec["parent"]]
            frame = self.frames[spec["frame"]]
            origin = frame["origin_m"]
            parent_origin = self.body_origin_global[spec["parent"]]
            body = ET.SubElement(
                parent,
                "body",
                {"name": spec["body"], "pos": fmt_vec(origin - parent_origin)},
            )
            self.body_elements[spec["body"]] = body
            self.body_origin_global[spec["body"]] = origin

            axis = frame["axes"][int(spec["axis_index"])]
            lower, upper = spec["range_deg"]
            attributes = {
                "name": spec["joint"],
                "class": "arm_joint",
                "axis": fmt_vec(axis),
                "range": f"{fmt(lower)} {fmt(upper)}",
            }
            stiffness = self.joint_stiffness.get(spec["joint"], 0.0)
            if stiffness:
                # Passive capsule / ligament stand-in: a wire that also torques
                # this joint no longer drags it to its limit for free.
                attributes["stiffness"] = fmt(stiffness)
                attributes["springref"] = "0"
            ET.SubElement(
                body,
                "joint",
                attributes,
            )
            if spec["mesh"]:
                self._add_mesh_geom(
                    body,
                    spec["mesh"],
                    origin,
                    f"{spec['mesh']}_mat",
                    float(spec["mass_kg"]),
                )
            else:
                self._add_massive_body_inertial(body, float(spec["mass_kg"]))
            self.body_mass[spec["body"]] = float(spec["mass_kg"])

        actuator = ET.SubElement(self.root, "actuator")
        for spec in JOINT_SPECS:
            if spec["body"] in (
                "clavicle_protraction",
                "clavicle_elevation",
                "clavicle_axial_rotation",
                "scapula_upward_rotation",
                "scapula_posterior_tilt",
                "scapula_external_rotation",
            ):
                lower, upper = spec["range_deg"]
                ET.SubElement(
                    actuator,
                    "position",
                    {
                        "name": f"{spec['joint']}_act",
                        "class": "rhythm_actuator",
                        "joint": spec["joint"],
                        "ctrlrange": f"{fmt(lower)} {fmt(upper)}",
                    },
                )
        self.rhythm_actuator_root = actuator

    def _add_mesh_geom(
        self,
        body: ET.Element,
        mesh_name: str,
        body_origin: np.ndarray,
        material: str,
        mass: float,
    ) -> None:
        path = self.bone_dir / f"{mesh_name}.stl"
        triangles = read_binary_stl(path)
        scale = detect_mesh_scale(path)
        vertices_m = triangles.reshape(-1, 3) * scale
        centroid_global = vertices_m.mean(axis=0)
        extent = np.ptp(vertices_m, axis=0)
        ET.SubElement(
            body,
            "inertial",
            {
                "pos": fmt_vec(centroid_global - body_origin),
                "mass": fmt(mass),
                "diaginertia": fmt_vec(diagonal_inertia(mass, extent)),
            },
        )
        ET.SubElement(
            body,
            "geom",
            {
                "name": f"{mesh_name}_visual",
                "class": "visual_mesh",
                "mesh": mesh_name,
                "material": material,
                "pos": fmt_vec(-body_origin),
            },
        )

    def _add_massive_body_inertial(self, body: ET.Element, mass: float) -> None:
        ET.SubElement(
            body,
            "inertial",
            {
                "pos": "0 0 0",
                "mass": fmt(mass),
                "diaginertia": "1e-5 1e-5 1e-5",
            },
        )

    def _add_cables(self) -> None:
        tendon_root = ET.SubElement(self.root, "tendon")
        worldbody = self.root.find("worldbody")
        if worldbody is None:
            raise ValueError("worldbody is missing.")
        slider_offset = float(self.defaults["slider_offset_m"])
        slider_radius = float(self.defaults["slider_radius_m"])
        # A slider may be retracted below its neutral position to let the wire
        # go slack; without that the other cables would hold the arm.
        release = float(self.defaults.get("release_m", 0.08))

        for muscle in self.manifest["cables"]:
            if not muscle.get("enabled", True):
                continue
            insertion_body = muscle["insertion_body"]
            if insertion_body not in self.body_elements:
                raise ValueError(f"Unknown insertion body: {insertion_body}")

            insertion_local = self._attach_local(
                insertion_body,
                muscle.get("insertion_part_m"),
                muscle["insertion_local_m"],
                "insertion_local_m",
            )
            insertion_global = self.body_origin_global[insertion_body] + insertion_local
            guides = self._resolve_guides(muscle, insertion_global)
            guides = [
                {
                    "body": guide["body"],
                    "local_m": self._attach_local(
                        guide["body"],
                        guide.get("part_m"),
                        guide["local_m"],
                        f"{muscle['id']}.guides.local_m",
                    ),
                    "radius_m": guide["radius_m"],
                }
                for guide in guides
            ]
            for guide in guides:
                if guide["body"] not in self.body_elements:
                    raise ValueError(f"Unknown guide body: {guide['body']}")
            guide_globals = [
                self.body_origin_global[guide["body"]] + guide["local_m"]
                for guide in guides
            ]
            last_guide_global = guide_globals[-1]
            if len(guide_globals) >= 2:
                direction = normalize(
                    guide_globals[-1] - guide_globals[-2],
                    f"{muscle['id']}.direction",
                )
            else:
                direction = normalize(
                    last_guide_global - insertion_global,
                    f"{muscle['id']}.direction",
                )
            perpendicular = np.cross(direction, np.array([0.0, 0.0, 1.0]))
            if np.linalg.norm(perpendicular) < 1e-8:
                perpendicular = np.cross(direction, np.array([0.0, 1.0, 0.0]))
            perpendicular = normalize(perpendicular, "cable_perpendicular")

            channels = int(muscle.get("channels", 1))
            channel_offset = float(muscle.get("channel_offset_m", 0.0))
            for channel in range(channels):
                offset = 0.0
                if channels > 1:
                    offset = channel_offset * (channel - (channels - 1) / 2.0)
                offset_vector = perpendicular * offset
                cable_id = f"{muscle['id']}_{channel}"
                insertion_site = f"cable_{cable_id}_insertion"
                guide_sites = [
                    f"cable_{cable_id}_via_{index}"
                    for index in range(len(guides))
                ]
                slider_body_name = f"{SLIDER_JOINT_PREFIX}{cable_id}"
                slider_site = f"slider_{cable_id}_site"
                cable_name = f"cable_{cable_id}"

                self._add_site(
                    self.body_elements[insertion_body],
                    insertion_site,
                    insertion_local + offset_vector,
                    (0.95, 0.38, 0.12, 0.9),
                )
                for guide, guide_site in zip(guides, guide_sites):
                    self._add_site(
                        self.body_elements[guide["body"]],
                        guide_site,
                        guide["local_m"] + offset_vector,
                        (0.95, 0.38, 0.12, 0.9),
                    )
                    ET.SubElement(
                        self.body_elements[guide["body"]],
                        "geom",
                        {
                            "name": f"{guide_site}_fixture",
                            "type": "sphere",
                            "pos": fmt_vec(guide["local_m"] + offset_vector),
                            "size": fmt(guide["radius_m"]),
                            "material": "slider_mat",
                            "contype": "0",
                            "conaffinity": "0",
                            "group": "2",
                        },
                    )

                slider_origin = (
                    last_guide_global + direction * slider_offset + offset_vector
                )
                slider = ET.SubElement(
                    worldbody,
                    "body",
                    {"name": slider_body_name, "pos": fmt_vec(slider_origin)},
                )
                ET.SubElement(
                    slider,
                    "joint",
                    {
                        "name": slider_body_name,
                        "class": "slider_joint",
                        "axis": fmt_vec(direction),
                        "range": f"{fmt(-release)} {fmt(self._travel_for(muscle))}",
                    },
                )
                ET.SubElement(
                    slider,
                    "inertial",
                    {
                        "pos": "0 0 0",
                        "mass": "0.02",
                        "diaginertia": "1e-6 1e-6 1e-6",
                    },
                )
                ET.SubElement(
                    slider,
                    "geom",
                    {
                        "name": f"{slider_body_name}_geom",
                        "type": "sphere",
                        "size": fmt(slider_radius),
                        "material": "slider_mat",
                        "contype": "0",
                        "conaffinity": "0",
                    },
                )
                ET.SubElement(
                    slider,
                    "site",
                    {
                        "name": slider_site,
                        "pos": "0 0 0",
                        "size": "0.004",
                        "rgba": "0.95 0.38 0.12 0.9",
                        "group": "3",
                    },
                )

                stiffness = self._stiffness_for(muscle)
                effective_mass = self._effective_mass_for(muscle)
                damping = self._damping_for(effective_mass, stiffness)
                tendon = ET.SubElement(
                    tendon_root,
                    "spatial",
                    {
                        "name": cable_name,
                        "width": fmt(float(self.defaults["cable_width_m"])),
                        "stiffness": fmt(stiffness),
                        "damping": fmt(damping),
                        "springlength": "0",
                        "rgba": "0.95 0.38 0.12 0.9",
                    },
                )
                ET.SubElement(tendon, "site", {"site": insertion_site})
                for guide_site in guide_sites:
                    ET.SubElement(tendon, "site", {"site": guide_site})
                ET.SubElement(tendon, "site", {"site": slider_site})
                self.tendon_elements[cable_name] = tendon
                self.tendon_stiffness[cable_name] = stiffness
                actuator_element = self._add_slider_actuator(
                    muscle, cable_id, slider_body_name
                )
                self._add_assist_actuator(muscle, cable_id, cable_name)
                self.cable_channels.append(
                    {
                        "id": cable_id,
                        "muscle_id": muscle["id"],
                        "group": muscle.get("group", muscle["id"]),
                        "cable": cable_name,
                        "slider_joint": slider_body_name,
                        "actuator": f"{SLIDER_ACTUATOR_PREFIX}{cable_id}",
                        "max_force_n": float(muscle["max_force_n"]),
                        "travel_m": self._travel_for(muscle),
                        "pretension_n": float(self.defaults["pretension_n"]),
                        "damping_ratio": float(self.defaults["damping_ratio"]),
                        "effective_mass_kg": effective_mass,
                        "guides": guide_sites,
                        "guide_bodies": [guide["body"] for guide in guides],
                        "tendon": tendon,
                        "actuator_element": actuator_element,
                    }
                )

    def _resolve_guides(
        self, muscle: dict, insertion_global: np.ndarray
    ) -> list[dict]:
        guide_mode = muscle.get("guide_mode", "constant_moment_arm")
        if guide_mode not in ("constant_moment_arm", "pulley"):
            return cable_guides(muscle)

        target_joint = muscle["target_joint"]
        spec = next(
            (item for item in JOINT_SPECS if item["joint"] == target_joint),
            None,
        )
        if spec is None:
            raise ValueError(
                f"Unknown target joint for {muscle['id']}: {target_joint}"
            )
        guide_body = muscle["guide_body"]
        if guide_body not in self.body_elements:
            raise ValueError(f"Unknown guide body: {guide_body}")

        axis = self.frames[spec["frame"]]["axes"][int(spec["axis_index"])]
        joint_origin = self.body_origin_global[spec["body"]]
        insertion_vector = insertion_global - joint_origin
        along = float(np.dot(insertion_vector, axis))
        perpendicular_raw = insertion_vector - axis * along
        if np.linalg.norm(perpendicular_raw) < 1e-8:
            reference = np.array([0.0, 0.0, 1.0])
            if abs(float(np.dot(axis, reference))) > 0.9:
                reference = np.array([0.0, 1.0, 0.0])
            perpendicular_raw = np.cross(axis, reference)
        perpendicular = normalize(perpendicular_raw, "guide perpendicular")
        moment_sign = float(muscle.get("guide_moment_sign", 1.0))

        if guide_mode == "constant_moment_arm":
            tangent = normalize(
                np.cross(axis, perpendicular),
                "guide tangent",
            )
            if moment_sign < 0.0:
                tangent = -tangent
            guide_offset = float(muscle.get("guide_offset_m", 0.06))
            guide_global = (
                joint_origin
                + axis * along
                + perpendicular_raw
                + tangent * guide_offset
            )
            guide_local = guide_global - self.body_origin_global[guide_body]
            return [
                {
                    "body": guide_body,
                    "local_m": guide_local,
                    "radius_m": float(muscle.get("guide_radius_m", 0.006)),
                }
            ]

        e1 = perpendicular
        e2 = normalize(np.cross(axis, e1), "pulley tangent")
        # A positive moment sign means the cable must drive the joint in the
        # positive direction, which needs dL/dq < 0, so the wrap has to go the
        # other way round the pulley.
        if moment_sign > 0.0:
            # Move the arc to the other side of the axis; mirroring alone
            # would only reverse the same polyline.
            e1 = -e1
            e2 = -e2
        radius = float(muscle.get("pulley_radius_m", 0.025))
        wrap_angle = math.radians(float(muscle.get("wrap_angle_deg", 140.0)))
        point_count = max(2, int(muscle.get("pulley_points", 6)))
        # Enter the arc at the tangent point, otherwise the first cable
        # segment is radial and the joint has no moment arm at the neutral
        # pose, so it can never start moving.
        reach = float(np.linalg.norm(perpendicular_raw))
        if reach > radius:
            start = -math.acos(max(-1.0, min(1.0, radius / reach)))
        else:
            start = 0.0
        guides = []
        for index in range(point_count):
            ratio = index / (point_count - 1)
            theta = start + wrap_angle * ratio
            guide_global = (
                joint_origin
                + axis * along
                + radius * (math.cos(theta) * e1 + math.sin(theta) * e2)
            )
            guides.append(
                {
                    "body": guide_body,
                    "local_m": guide_global
                    - self.body_origin_global[guide_body],
                    "radius_m": float(muscle.get("guide_radius_m", 0.006)),
                }
            )
        return guides

    def _attach_local(
        self,
        body: str,
        part_point,
        local_point,
        label: str,
    ) -> np.ndarray:
        """Offset of an attach point from the body origin.

        A point stored in bone part coordinates is mapped through the bone
        transform, so it stays on the bone no matter how the parts are placed
        in the SolidWorks assembly.  Points stored as a plain offset keep the
        old behaviour.
        """

        if part_point is not None:
            bone = self.body_bone.get(body)
            matrix = self.bone_transforms.get(bone)
            if matrix is None:
                raise ValueError(
                    f"{label}: {body} has no bone transform in "
                    "bone_transforms.json."
                )
            point = np.asarray(part_point, dtype=float).reshape(3)
            world = matrix[:3, :3] @ point + matrix[:3, 3]
            return world - self.body_origin_global[body]
        return vec3(local_point, label)

    def _add_site(
        self,
        body: ET.Element,
        name: str,
        local: np.ndarray,
        rgba,
    ) -> None:
        ET.SubElement(
            body,
            "site",
            {
                "name": name,
                "pos": fmt_vec(local),
                "size": "0.004",
                "rgba": " ".join(fmt(value) for value in rgba),
                "group": "3",
            },
        )

    def _add_slider_actuator(
        self, muscle: dict, cable_id: str, slider_joint: str
    ) -> ET.Element:
        return ET.SubElement(
            self.rhythm_actuator_root,
            "position",
            {
                "name": f"{SLIDER_ACTUATOR_PREFIX}{cable_id}",
                "class": "slider_actuator",
                "joint": slider_joint,
                "kp": fmt(float(self.defaults["slider_kp"])),
                "kv": fmt(float(self.defaults["slider_kv"])),
                "ctrlrange": (
                    f"{fmt(-float(self.defaults.get('release_m', 0.08)))} "
                    f"{fmt(self._travel_for(muscle))}"
                ),
                "group": "2",
            },
        )

    def _add_assist_actuator(
        self, muscle: dict, cable_id: str, cable_name: str
    ) -> ET.Element:
        """One-sided tension source that cancels the spring's push half.

        MuJoCo tendon springs are two-sided, so a "slack" cable would push
        instead of releasing.  The control loop writes
        ``max(0, -k (L - rest))`` into this actuator, which leaves the net
        tendon force at ``max(0, k (L - rest))`` -- a wire that only pulls.
        """

        limit = float(muscle["max_force_n"])
        return ET.SubElement(
            self.rhythm_actuator_root,
            "general",
            {
                "name": f"assist_{cable_id}",
                "tendon": cable_name,
                "gainprm": "1 0 0",
                "biasprm": "0 0 0",
                "ctrlrange": f"0 {fmt(limit)}",
                "forcerange": f"0 {fmt(limit)}",
                "group": "3",
            },
        )

    def _stiffness_for(self, muscle: dict) -> float:
        max_force = float(muscle["max_force_n"])
        rated_stretch = float(
            muscle.get("rated_stretch_m", self.defaults["rated_stretch_m"])
        )
        return max_force / max(rated_stretch, 1e-6)

    def _travel_for(self, muscle: dict) -> float:
        requested = float(
            muscle.get("travel_m", self.defaults["travel_m"])
        )
        return max(0.001, requested)

    def _effective_mass_for(self, muscle: dict) -> float:
        insertion_body = muscle["insertion_body"]
        via_body = muscle["via_body"]
        return max(
            0.02,
            min(
                self.body_mass.get(insertion_body, 0.1),
                self.body_mass.get(via_body, 0.1),
                0.25,
            ),
        )

    def _damping_for(self, effective_mass: float, stiffness: float) -> float:
        zeta = float(self.defaults["damping_ratio"])
        return 2.0 * zeta * math.sqrt(stiffness * effective_mass)

    def _add_sensors(self) -> None:
        sensor = ET.SubElement(self.root, "sensor")
        for joint_name in SKELETON_JOINT_NAMES:
            ET.SubElement(
                sensor,
                "jointpos",
                {"name": f"sensor_{joint_name}", "joint": joint_name},
            )
        for channel in self.cable_channels:
            ET.SubElement(
                sensor,
                "jointpos",
                {
                    "name": f"sensor_{channel['slider_joint']}",
                    "joint": channel["slider_joint"],
                },
            )
            ET.SubElement(
                sensor,
                "tendonpos",
                {"name": f"sensor_{channel['cable']}", "tendon": channel["cable"]},
            )
            ET.SubElement(
                sensor,
                "tendonvel",
                {
                    "name": f"sensor_{channel['cable']}_velocity",
                    "tendon": channel["cable"],
                },
            )
            ET.SubElement(
                sensor,
                "actuatorfrc",
                {
                    "name": f"sensor_{channel['actuator']}",
                    "actuator": channel["actuator"],
                },
            )

    def _sample_states(self, model: mujoco.MjModel) -> list[np.ndarray]:
        return task_envelope_states(model, random_samples=600, seed=11)

    def _calibrate_springlengths(self) -> None:
        model = mujoco.MjModel.from_xml_path(str(self.output_path))
        data = mujoco.MjData(model)
        tendon_ids = {}
        for cable_name in self.tendon_elements:
            tendon_id = mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_TENDON, cable_name
            )
            if tendon_id < 0:
                raise RuntimeError(f"Tendon missing after build: {cable_name}")
            tendon_ids[cable_name] = int(tendon_id)

        minimum_lengths = {
            cable_name: float("inf") for cable_name in tendon_ids
        }
        maximum_lengths = {
            cable_name: float("-inf") for cable_name in tendon_ids
        }
        for state in self._sample_states(model):
            data.qpos[:] = state
            mujoco.mj_forward(model, data)
            for cable_name, tendon_id in tendon_ids.items():
                minimum_lengths[cable_name] = min(
                    minimum_lengths[cable_name],
                    float(data.ten_length[tendon_id]),
                )
                maximum_lengths[cable_name] = max(
                    maximum_lengths[cable_name],
                    float(data.ten_length[tendon_id]),
                )

        pretension = float(self.defaults["pretension_n"])
        rest_mode = str(self.defaults.get("spring_rest_mode", "reference"))
        reference_lengths: dict[str, float] = {}
        if rest_mode == "reference":
            data.qpos[:] = 0.0
            mujoco.mj_forward(model, data)
            reference_lengths = {
                cable_name: float(data.ten_length[tendon_id])
                for cable_name, tendon_id in tendon_ids.items()
            }
        for channel in self.cable_channels:
            cable_name = channel["cable"]
            tendon = channel["tendon"]
            minimum = minimum_lengths[cable_name]
            maximum = maximum_lengths[cable_name]
            length_span = max(0.0, maximum - minimum)
            if rest_mode == "reference":
                # Taut at the neutral pose and pulled by the slider.  A cable
                # that would shorten below this rest length just goes slack,
                # which the linear spring can only approximate.
                travel = max(channel["travel_m"], 1e-4)
                stiffness = max(
                    1.0, channel["max_force_n"] - pretension
                ) / travel
            else:
                stiffness = (
                    max(1.0, channel["max_force_n"] - pretension)
                    / max(length_span + channel["travel_m"], 1e-4)
                )
            damping = self._damping_for(channel["effective_mass_kg"], stiffness)
            safety_margin = float(
                self.defaults.get("spring_safety_margin_m", 0.001)
            )
            if rest_mode == "reference":
                reference = reference_lengths[cable_name]
                # Leave a little slack at the neutral pose.  With the
                # one-sided assist actuator a slack cable pulls nothing at
                # all, so the arm sits still instead of drifting on the
                # balance of eight taut springs.
                slack = float(
                    self.defaults.get("spring_reference_slack_m", 0.005)
                )
                rest = reference + slack
            else:
                rest = max(
                    0.0, minimum - max(pretension / stiffness, safety_margin)
                )
            tendon.set("stiffness", fmt(stiffness))
            tendon.set("damping", fmt(damping))
            tendon.set("springlength", fmt(rest))
            self.tendon_stiffness[cable_name] = stiffness
            actuator = channel.get("actuator_element")
            if actuator is not None:
                limit = channel["max_force_n"]
                actuator.set("forcelimited", "true")
                actuator.set("forcerange", f"{-limit:g} {limit:g}")

        ET.indent(self.root, space="  ")
        ET.ElementTree(self.root).write(
            self.output_path, encoding="utf-8", xml_declaration=True
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT_DIR)
    parser.add_argument(
        "--joint-frames",
        type=Path,
        default=DEFAULT_EXPORT_DIR / "joint_frames.json",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_XML)
    args = parser.parse_args()

    builder = CableDriveBuilder(
        args.export_dir,
        args.joint_frames,
        args.manifest,
        args.output,
    )
    output = builder.build()
    model = mujoco.MjModel.from_xml_path(str(output))
    print(f"generated: {output}")
    print(f"bodies: {model.nbody}")
    print(f"joints: {model.njnt}")
    print(f"tendons: {model.ntendon}")
    print(f"actuators: {model.nu}")
    print(f"DOF: {model.nv}")
    print(f"cable channels: {len(builder.cable_channels)}")


if __name__ == "__main__":
    main()
