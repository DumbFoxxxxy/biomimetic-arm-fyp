"""Build the BodyParts3D five-bone cable-drive MuJoCo model.

Body tree (every body origin is the proximal joint centre of that bone, so a
global point maps to a local one by subtracting the body origin)::

    thorax -> clavicle -> scapula -> humerus -> ulna -> radius

The thorax is the fixed base.  Joints are ball joints except the elbow, which
is a hinge built on the cylinder fitted through the distal humerus.  Each wire
is a spatial tendon anchored on the insertion bone, wrapped over the origin
bone and finished on a world-fixed slider, exactly like a muscle that is
attached at both ends and pulled from outside the body.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

import bp3d_common as B
SLIDER_ACTUATOR_PREFIX = "pull_"
SLIDER_JOINT_PREFIX = "slider_"
WIRE_PREFIX = "wire_"

# Passive stiffness of the printed TPU collar that connects two bones.  A ball
# joint applies it along the rotation vector, which is the isotropic lattice
# case the source framework describes, and the ligament strips add the
# directional limits on top.  The girdle values are set so that the weight of
# the arm (about 0.55 N.m about the sternoclavicular joint) does not collapse
# the shoulder between cable commands: a real arm holds its pose with muscle
# tone, and with every cable slack the collar is the only thing left.
DEFAULT_JOINT_STIFFNESS = {
    "sternoclavicular": 8.0,
    "acromioclavicular": 8.0,
    # The girdle is stiff because it has to carry the arm's weight with every
    # cable slack, but the humerus and the elbow must stay free or the collar
    # would eat the cable torque instead of letting the joint travel.
    "glenohumeral": 2.5,
    "elbow": 1.5,
    "proximal_radioulnar": 8.0,
}

# Which bone pairs are allowed to touch.  A convex hull fills in the
# concavities that let the real clavicle and scapula lie on the rib cage, so
# the girdle hulls intrude up to 14 mm into the rib proxies at the neutral
# pose.  Disabling that pair is better than letting the penetration shove the
# arm around; the scapulothoracic interface needs a concave decomposition
# before it can be switched on.  The pairs that stay enabled are the ones
# whose proxies are geometrically faithful: the humeral head sphere against
# the scapula hull, and the two forearm bones against each other.
CONTACT_BITS = {
    "thorax": (1, 0),
    "clavicle": (8, 0),
    "scapula": (16, 2),
    "humerus": (2, 16 | 64),
    "ulna": (64, 2 | 128),
    "radius": (128, 64),
}


def fmt(value: float) -> str:
    return f"{value:.9g}"


def fmt_vec(values) -> str:
    return " ".join(fmt(float(value)) for value in values)


class Bp3dArmBuilder:
    def __init__(
        self,
        contract_path: Path,
        muscle_path: Path,
        output_path: Path,
        collision: bool = True,
        ligament_path: Path | None = None,
        rom_path: Path | None = None,
    ) -> None:
        self.contract = B.load_json(contract_path)
        self.muscles = B.load_json(muscle_path)
        self.ligaments = (
            B.load_json(ligament_path)
            if ligament_path is not None and ligament_path.exists()
            else None
        )
        self.rom = (
            B.load_json(rom_path)
            if rom_path is not None and rom_path.exists()
            else None
        )
        self.output_path = output_path
        self.collision = collision
        self.asset_dir = B.ASSET_DIR

        self.joints = self.contract["joints"]
        self.bones = self.contract["bones"]
        self.body_origin: dict[str, np.ndarray] = {
            "thorax": np.zeros(3),
            **{
                entry["child_bone"]: np.asarray(entry["origin_m"], dtype=float)
                for entry in self.joints.values()
            },
        }
        self.parent_joint: dict[str, str] = {
            entry["child_bone"]: name for name, entry in self.joints.items()
        }
        self.child_of: dict[str, str | None] = {"thorax": None}
        for name, entry in self.joints.items():
            self.child_of[entry["child_bone"]] = entry["parent_bone"]

        self.root = ET.Element("mujoco", {"model": self.muscles["model_name"]})
        self.body_elements: dict[str, ET.Element] = {}
        self.tendon_elements: dict[str, ET.Element] = {}
        self.wires: list[dict[str, object]] = []
        self.ligament_tendons: list[dict[str, object]] = []
        self.tendon_root: ET.Element | None = None
        self.actuator_root: ET.Element | None = None

    # -- top level ---------------------------------------------------------

    def build(self) -> Path:
        self._add_top_level()
        self._add_assets()
        self._add_worldbody()
        # MuJoCo requires tendon and actuator after worldbody and in that
        # order, so both containers are created here rather than on demand.
        self.tendon_root = ET.SubElement(self.root, "tendon")
        self.actuator_root = ET.SubElement(self.root, "actuator")
        self._add_ligaments()
        self._add_wires()
        self._add_sensors()
        ET.indent(self.root, space="  ")
        ET.ElementTree(self.root).write(
            self.output_path, encoding="utf-8", xml_declaration=True
        )
        self._calibrate()
        return self.output_path

    def _add_top_level(self) -> None:
        ET.SubElement(
            self.root,
            "compiler",
            {
                "angle": "degree",
                "autolimits": "true",
                "meshdir": self.asset_dir.as_posix(),
            },
        )
        ET.SubElement(
            self.root,
            "option",
            {
                "timestep": fmt(float(self.muscles.get("timestep", 0.001))),
                "gravity": fmt_vec(self.muscles.get("gravity", [0.0, 0.0, -9.81])),
                "integrator": self.muscles.get("integrator", "implicitfast"),
                "cone": "elliptic",
            },
        )
        visual = ET.SubElement(self.root, "visual")
        ET.SubElement(
            visual,
            "headlight",
            {"diffuse": "0.8 0.8 0.8", "ambient": "0.3 0.3 0.3"},
        )
        ET.SubElement(visual, "rgba", {"haze": "0.15 0.2 0.25 1"})
        ET.SubElement(
            visual,
            "global",
            {
                "azimuth": "150",
                "elevation": "-12",
                # A larger offscreen buffer so the model can be rendered for
                # visual QA without a window.
                "offwidth": "1920",
                "offheight": "1080",
            },
        )

        default = ET.SubElement(self.root, "default")
        visual_mesh = ET.SubElement(default, "default", {"class": "visual_mesh"})
        ET.SubElement(
            visual_mesh,
            "geom",
            {
                "type": "mesh",
                "contype": "0",
                "conaffinity": "0",
                "density": "0",
                "group": "2",
                "rgba": "0.82 0.80 0.72 1",
            },
        )
        bone = ET.SubElement(default, "default", {"class": "bone"})
        ET.SubElement(
            bone,
            "geom",
            {
                "contype": "1",
                "conaffinity": "1",
                "group": "0",
                "rgba": "0.55 0.62 0.72 1",
                "solref": "0.02 1",
                "margin": "0.0005",
            },
        )
        slider_joint = ET.SubElement(default, "default", {"class": "slider_joint"})
        ET.SubElement(
            slider_joint,
            "joint",
            {"type": "slide", "limited": "true", "damping": "0.05", "armature": "0.001"},
        )
        slider_actuator = ET.SubElement(default, "default", {"class": "slider_actuator"})
        ET.SubElement(slider_actuator, "position", {})

    def _add_assets(self) -> None:
        asset = ET.SubElement(self.root, "asset")
        for name in self.bones:
            ET.SubElement(asset, "mesh", {"name": name, "file": f"{name}.obj"})
        if self.collision:
            for entry in self.contract["collision"].values():
                if isinstance(entry, dict) and "file" in entry:
                    mesh_name = Path(entry["file"]).stem
                    ET.SubElement(
                        asset, "mesh", {"name": mesh_name, "file": entry["file"]}
                    )
        for name, rgba in (
            ("bone_mat", "0.85 0.83 0.76 1"),
            ("slider_mat", "0.95 0.42 0.15 1"),
            ("ligament_mat", "0.30 0.65 0.85 1"),
        ):
            ET.SubElement(
                asset, "material", {"name": name, "rgba": rgba, "specular": "0.25"}
            )

    # -- skeleton ----------------------------------------------------------

    def _add_worldbody(self) -> None:
        worldbody = ET.SubElement(self.root, "worldbody")
        ET.SubElement(
            worldbody,
            "light",
            {
                "name": "key_light",
                "pos": "-0.4 -0.6 1.6",
                "dir": "0.4 0.4 -1",
                "directional": "true",
                "diffuse": "0.9 0.9 0.9",
                "specular": "0.2 0.2 0.2",
            },
        )

        order = self._body_order()
        for name in order:
            origin = self.body_origin[name]
            parent = self.child_of[name]
            if parent is None:
                body = ET.SubElement(
                    worldbody, "body", {"name": name, "pos": fmt_vec(origin)}
                )
                body.set("pos", fmt_vec(origin))
            else:
                body = ET.SubElement(
                    self.body_elements[parent],
                    "body",
                    {"name": name, "pos": fmt_vec(origin - self.body_origin[parent])},
                )
            self.body_elements[name] = body

            joint_name = self.parent_joint.get(name)
            if joint_name is not None:
                self._add_joint(body, joint_name)

            self._add_inertial(body, name, origin)
            ET.SubElement(
                body,
                "geom",
                {
                    "name": f"{name}_visual",
                    "class": "visual_mesh",
                    "mesh": name,
                    "material": "bone_mat",
                    "pos": fmt_vec(-origin),
                },
            )
            if self.collision:
                self._add_collision(body, name, origin)


    def _body_order(self) -> list[str]:
        order: list[str] = []
        pending = ["thorax"]
        while pending:
            name = pending.pop(0)
            order.append(name)
            pending.extend(
                entry["child_bone"]
                for entry in self.joints.values()
                if entry["parent_bone"] == name
            )
        return order

    def _add_joint(self, body: ET.Element, joint_name: str) -> None:
        entry = self.joints[joint_name]
        stiffness = DEFAULT_JOINT_STIFFNESS.get(joint_name, 0.0)
        if entry["type"] == "ball":
            limit = max(abs(float(v)) for v in entry["range_deg"])
            ET.SubElement(
                body,
                "joint",
                {
                    "name": joint_name,
                    "type": "ball",
                    "pos": "0 0 0",
                    "damping": fmt(0.35),
                    "armature": fmt(0.02),
                    "stiffness": fmt(stiffness),
                    "springref": "0",
                    "limited": "true",
                    "range": f"0 {fmt(limit)}",
                },
            )
            return
        axis = np.asarray(entry["axes"][entry["primary_axis"]], dtype=float)
        low, high = (float(v) for v in entry["range_deg"])
        ET.SubElement(
            body,
            "joint",
            {
                "name": joint_name,
                "type": "hinge",
                "pos": "0 0 0",
                "axis": fmt_vec(axis),
                "damping": fmt(0.2),
                "armature": fmt(0.02),
                "stiffness": fmt(stiffness),
                "springref": "0",
                "limited": "true",
                "range": f"{fmt(low)} {fmt(high)}",
            },
        )

    def _add_inertial(self, body: ET.Element, name: str, origin: np.ndarray) -> None:
        entry = self.bones[name]
        ET.SubElement(
            body,
            "inertial",
            {
                "pos": fmt_vec(np.asarray(entry["centre_m"], dtype=float) - origin),
                "quat": fmt_vec(entry["inertia_quat_wxyz"]),
                "mass": fmt(float(entry["mass_kg"])),
                "diaginertia": fmt_vec(entry["inertia_diagonal"]),
            },
        )

    def _add_collision(self, body: ET.Element, name: str, origin: np.ndarray) -> None:
        collision = self.contract["collision"]
        contype, conaffinity = CONTACT_BITS[name]
        mask = {"contype": str(contype), "conaffinity": str(conaffinity)}
        if name == "thorax":
            for key in ("sternum_hull", "spine_hull"):
                entry = collision[key]
                ET.SubElement(
                    body,
                    "geom",
                    {
                        "name": key,
                        "class": "bone",
                        "type": "mesh",
                        "mesh": Path(entry["file"]).stem,
                        "pos": "0 0 0",
                        "material": "slider_mat",
                        **mask,
                    },
                )
            for index, sphere in enumerate(collision["rib_gliding"]):
                ET.SubElement(
                    body,
                    "geom",
                    {
                        "name": f"rib_{index}",
                        "class": "bone",
                        "type": "sphere",
                        "pos": fmt_vec(sphere["centre_m"]),
                        "size": fmt(float(sphere["radius_m"])),
                        "rgba": "0.55 0.62 0.72 1",
                        **mask,
                    },
                )
            return
        hull = collision[f"{name}_hull"]
        ET.SubElement(
            body,
            "geom",
            {
                "name": f"{name}_hull",
                "class": "bone",
                "type": "mesh",
                "mesh": Path(hull["file"]).stem,
                "pos": fmt_vec(-origin),
                **mask,
            },
        )
        if name == "humerus":
            head = collision["humeral_head"]
            ET.SubElement(
                body,
                "geom",
                {
                    "name": "humeral_head",
                    "class": "bone",
                    "type": "sphere",
                    "pos": fmt_vec(np.asarray(head["centre_m"], dtype=float) - origin),
                    "size": fmt(float(head["radius_m"])),
                    "rgba": "0.70 0.55 0.35 1",
                    **mask,
                },
            )

    # -- wires -------------------------------------------------------------

    def _add_wires(self) -> None:
        tendon_root = self.tendon_root
        assert tendon_root is not None
        worldbody = self.root.find("worldbody")
        defaults = self.muscles.setdefault("defaults", {})
        slider_offset = float(defaults.get("slider_offset_m", 0.10))
        slider_radius = float(defaults.get("slider_radius_m", 0.008))
        release = float(defaults.get("release_m", 0.08))
        width = float(defaults.get("wire_width_m", 0.0015))

        for channel in self.muscles["channels"]:
            if not channel.get("enabled", True):
                continue
            for path in channel["paths"]:
                sites = list(path["sites"])
                # The wire is anchored on the insertion bone, runs back along
                # the muscle centre line and leaves the body at the origin
                # end, where the slider pulls it.
                ordered = list(reversed(sites))
                anchor = ordered[0]
                tail = ordered[-1]
                pivot = self._slider_pivot(str(tail["body"]), str(anchor["body"]))
                outward = np.asarray(tail["pos_m"], dtype=float) - pivot
                if float(np.linalg.norm(outward)) < 1e-6:
                    outward = np.asarray(ordered[-1]["pos_m"], dtype=float) - np.asarray(
                        ordered[-2]["pos_m"], dtype=float
                    )
                outward = B.unit(outward, f"{path['id']} slider direction")

                wire_id = str(path["id"])
                site_names: list[str] = []
                for index, site in enumerate(ordered):
                    body = str(site["body"])
                    local = np.asarray(site["pos_m"], dtype=float) - self.body_origin[body]
                    name = f"w_{wire_id}_{index}"
                    self._add_site(self.body_elements[body], name, local)
                    site_names.append(name)

                slider_body = f"{SLIDER_JOINT_PREFIX}{wire_id}"
                slider_site = f"{SLIDER_JOINT_PREFIX}{wire_id}_site"
                # The slider rides on the tail bone, not on the world.  That is
                # the Bowden-tube case: the tube end is clamped to the limb, so
                # the wire's outlet moves with the body and the last segment
                # contributes no moment about any joint.  A world-fixed slider
                # would give every wire a huge lever arm about the most
                # proximal joint and the arm would simply pivot from the
                # sternoclavicular joint whenever any channel was pulled.
                tail_body = str(tail["body"])
                tail_local = (
                    np.asarray(tail["pos_m"], dtype=float) - self.body_origin[tail_body]
                )
                slider_pos = tail_local + outward * slider_offset
                slider = ET.SubElement(
                    self.body_elements[tail_body],
                    "body",
                    {"name": slider_body, "pos": fmt_vec(slider_pos)},
                )
                ET.SubElement(
                    slider,
                    "joint",
                    {
                        "name": slider_body,
                        "class": "slider_joint",
                        "axis": fmt_vec(outward),
                        "range": f"{fmt(-release)} {fmt(float(channel['travel_m']))}",
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
                        "name": f"{slider_body}_geom",
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
                        "rgba": "0.95 0.42 0.15 0.9",
                        "group": "3",
                    },
                )

                tendon = ET.SubElement(
                    tendon_root,
                    "spatial",
                    {
                        "name": f"{WIRE_PREFIX}{wire_id}",
                        "width": fmt(width),
                        "stiffness": "1000",
                        "damping": "1",
                        "springlength": "0",
                        "rgba": "0.95 0.42 0.15 0.9",
                    },
                )
                for name in site_names:
                    ET.SubElement(tendon, "site", {"site": name})
                ET.SubElement(tendon, "site", {"site": slider_site})
                self.tendon_elements[f"{WIRE_PREFIX}{wire_id}"] = tendon

                self._add_slider_actuator(channel, wire_id, slider_body, release)
                self._add_assist_actuator(channel, wire_id)
                self.wires.append(
                    {
                        "id": wire_id,
                        "channel": channel["id"],
                        "group": channel["group"],
                        "tendon": f"{WIRE_PREFIX}{wire_id}",
                        "tendon_element": tendon,
                        "slider_joint": slider_body,
                        "actuator": f"{SLIDER_ACTUATOR_PREFIX}{wire_id}",
                        "max_force_n": float(channel["max_force_n"]),
                        "travel_m": float(channel["travel_m"]),
                        "anchor_body": str(anchor["body"]),
                        "tail_body": str(tail["body"]),
                        "series_stiffness_n_per_m": float(
                            channel.get("series_stiffness_n_per_m", 0.0)
                        ),
                    }
                )

    def _slider_pivot(self, tail_bone: str, anchor_bone: str) -> np.ndarray:
        """Joint centre the wire tail is furthest from, so the slider is radial.

        A final wire segment that runs along the radius of a joint contributes
        no moment about that joint, which keeps the drive moment coming from
        the muscle path rather than from the exit hole.
        """

        for bone in (tail_bone, anchor_bone):
            joint_name = self.parent_joint.get(bone)
            if joint_name is not None:
                return np.asarray(self.joints[joint_name]["origin_m"], dtype=float)
        return np.zeros(3)

    def _add_site(self, body: ET.Element, name: str, local: np.ndarray) -> None:
        ET.SubElement(
            body,
            "site",
            {
                "name": name,
                "pos": fmt_vec(local),
                "size": "0.003",
                "rgba": "0.95 0.42 0.15 0.9",
                "group": "3",
            },
        )

    def _add_slider_actuator(
        self, channel: dict, wire_id: str, slider_joint: str, release: float
    ) -> None:
        assert self.actuator_root is not None
        defaults = self.muscles["defaults"]
        ET.SubElement(
            self.actuator_root,
            "position",
            {
                "name": f"{SLIDER_ACTUATOR_PREFIX}{wire_id}",
                "class": "slider_actuator",
                "joint": slider_joint,
                "kp": fmt(float(defaults.get("slider_kp", 50000.0))),
                "kv": fmt(float(defaults.get("slider_kv", 500.0))),
                "ctrlrange": f"{fmt(-release)} {fmt(float(channel['travel_m']))}",
                "group": "2",
            },
        )

    def _add_assist_actuator(self, channel: dict, wire_id: str) -> None:
        """Pull-only correction for a two-sided tendon spring.

        The range must include negative values: neutralizing the spring's
        compressive half means driving the actuator backwards.
        """

        assert self.actuator_root is not None
        limit = float(channel["max_force_n"])
        ET.SubElement(
            self.actuator_root,
            "general",
            {
                "name": f"assist_{wire_id}",
                "tendon": f"{WIRE_PREFIX}{wire_id}",
                "gainprm": "1 0 0",
                "biasprm": "0 0 0",
                "ctrlrange": f"{fmt(-limit)} {fmt(limit)}",
                "forcerange": f"{fmt(-limit)} {fmt(limit)}",
                "group": "3",
            },
        )

    # -- ligaments ---------------------------------------------------------

    def _add_ligaments(self) -> None:
        if not self.ligaments:
            return
        tendon_root = self.tendon_root
        assert tendon_root is not None
        assert self.actuator_root is not None
        for strip in self.ligaments["strips"]:
            site_names = []
            for index, site in enumerate(strip["sites"]):
                body = str(site["body"])
                local = np.asarray(site["pos_m"], dtype=float) - self.body_origin[body]
                name = f"lig_{strip['id']}_{index}"
                ET.SubElement(
                    self.body_elements[body],
                    "site",
                    {
                        "name": name,
                        "pos": fmt_vec(local),
                        "size": "0.002",
                        "rgba": "0.30 0.65 0.85 0.9",
                        "group": "3",
                    },
                )
                site_names.append(name)
            tendon = ET.SubElement(
                tendon_root,
                "spatial",
                {
                    "name": f"lig_{strip['id']}",
                    "width": fmt(float(strip["width_m"])),
                    "stiffness": fmt(float(strip["stiffness_n_per_m"])),
                    "damping": fmt(float(strip.get("damping", 0.05))),
                    "springlength": fmt(float(strip["springlength_m"])),
                    "rgba": "0.30 0.65 0.85 0.9",
                },
            )
            for name in site_names:
                ET.SubElement(tendon, "site", {"site": name})
            limit = float(strip["max_force_n"])
            ET.SubElement(
                self.actuator_root,
                "general",
                {
                    "name": f"assist_lig_{strip['id']}",
                    "tendon": f"lig_{strip['id']}",
                    "gainprm": "1 0 0",
                    "biasprm": "0 0 0",
                    "ctrlrange": f"{fmt(-limit)} {fmt(limit)}",
                    "forcerange": f"{fmt(-limit)} {fmt(limit)}",
                    "group": "3",
                },
            )
            self.ligament_tendons.append(
                {
                    "id": strip["id"],
                    "tendon": f"lig_{strip['id']}",
                    "joint": strip.get("joint"),
                }
            )

    # -- sensors -----------------------------------------------------------

    def _add_sensors(self) -> None:
        sensor = ET.SubElement(self.root, "sensor")
        for name in self.body_elements:
            if name == "thorax":
                continue
            ET.SubElement(
                sensor,
                "framequat",
                {"name": f"sensor_{name}_quat", "objtype": "body", "objname": name},
            )
        ET.SubElement(
            sensor,
            "jointpos",
            {"name": "sensor_elbow_flexion", "joint": "elbow"},
        )
        for wire in self.wires:
            ET.SubElement(
                sensor,
                "jointpos",
                {"name": f"sensor_{wire['slider_joint']}", "joint": wire["slider_joint"]},
            )
            ET.SubElement(
                sensor,
                "tendonpos",
                {"name": f"sensor_{wire['tendon']}", "tendon": wire["tendon"]},
            )
            ET.SubElement(
                sensor,
                "tendonvel",
                {
                    "name": f"sensor_{wire['tendon']}_velocity",
                    "tendon": wire["tendon"],
                },
            )
        for ligament in self.ligament_tendons:
            ET.SubElement(
                sensor,
                "tendonpos",
                {"name": f"sensor_{ligament['tendon']}", "tendon": ligament["tendon"]},
            )

    # -- spring rest lengths ----------------------------------------------

    def _sample_states(self, model: mujoco.MjModel) -> list[np.ndarray]:
        rng = np.random.default_rng(7)
        states = [np.zeros(model.nq)]
        joint_info: list[tuple[int, int, str]] = []
        for name in self.parent_joint.values():
            joint_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
            joint_info.append(
                (joint_id, int(model.jnt_qposadr[joint_id]), self.joints[name]["type"])
            )

        def write(qpos: np.ndarray, values: dict[str, float | np.ndarray]) -> None:
            for joint_id, address, kind in joint_info:
                value = values.get(mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, joint_id))
                if value is None:
                    continue
                if kind == "ball":
                    vector = np.asarray(value, dtype=float)
                    angle = float(np.linalg.norm(vector))
                    if angle < 1e-9:
                        qpos[address : address + 4] = [1.0, 0.0, 0.0, 0.0]
                    else:
                        axis = vector / angle
                        half = angle * 0.5
                        qpos[address : address + 4] = np.concatenate(
                            [[np.cos(half)], axis * np.sin(half)]
                        )
                else:
                    qpos[address] = float(value)

        # Structured sweeps: one joint at a time, everywhere in its range.
        for name, entry in self.joints.items():
            low, high = (float(v) for v in entry["range_deg"])
            axis = (
                np.zeros(3)
                if entry["type"] == "ball"
                else np.asarray(entry["axes"][entry["primary_axis"]], dtype=float)
            )
            for step in np.linspace(low, high, 7):
                qpos = np.zeros(model.nq)
                if entry["type"] == "ball":
                    direction = np.asarray(
                        entry["axes"][entry["primary_axis"]], dtype=float
                    )
                    write(qpos, {name: direction * np.radians(step)})
                else:
                    write(qpos, {name: np.radians(step)})
                states.append(qpos)

        # Random combinations so the rest-length calibration sees the whole
        # envelope rather than one joint at a time.
        for _ in range(240):
            qpos = np.zeros(model.nq)
            values: dict[str, np.ndarray] = {}
            for name, entry in self.joints.items():
                low, high = (float(v) for v in entry["range_deg"])
                angle = np.radians(rng.uniform(low, high))
                if entry["type"] == "ball":
                    direction = rng.normal(size=3)
                    values[name] = B.unit(direction, "random axis") * angle
                else:
                    values[name] = np.array([angle])
            for name, value in values.items():
                if self.joints[name]["type"] == "ball":
                    write(qpos, {name: value})
                else:
                    write(qpos, {name: float(value[0])})
            states.append(qpos)
        return states

    def _calibrate(self) -> None:
        model = mujoco.MjModel.from_xml_path(str(self.output_path))
        data = mujoco.MjData(model)
        defaults = self.muscles.setdefault("defaults", {})
        pretension = float(defaults.get("pretension_n", 0.2))
        slack = float(defaults.get("reference_slack_m", 0.004))
        damping_ratio = float(defaults.get("damping_ratio", 0.7))

        tendon_ids = {
            name: mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_TENDON, name)
            for name in self.tendon_elements
        }
        minimum = {name: float("inf") for name in tendon_ids}
        for state in self._sample_states(model):
            data.qpos[:] = state
            mujoco.mj_forward(model, data)
            for name, tendon_id in tendon_ids.items():
                minimum[name] = min(minimum[name], float(data.ten_length[tendon_id]))

        data.qpos[:] = 0.0
        mujoco.mj_forward(model, data)
        neutral = {
            name: float(data.ten_length[tendon_id])
            for name, tendon_id in tendon_ids.items()
        }

        for wire in self.wires:
            tendon_name = str(wire["tendon"])
            travel = max(float(wire["travel_m"]), 1e-4)
            stiffness = max(1.0, float(wire["max_force_n"]) - pretension) / travel
            mass = self._effective_mass(str(wire["anchor_body"]))
            damping = 2.0 * damping_ratio * math.sqrt(stiffness * mass)
            springlength = neutral[tendon_name] + slack
            tendon = wire["tendon_element"]
            assert isinstance(tendon, ET.Element)
            tendon.set("stiffness", fmt(stiffness))
            tendon.set("damping", fmt(damping))
            tendon.set("springlength", fmt(springlength))
            wire["stiffness_n_per_m"] = stiffness
            wire["springlength_m"] = springlength
            wire["neutral_length_m"] = neutral[tendon_name]
            wire["min_length_m"] = minimum[tendon_name]
            wire["slack_m"] = springlength - minimum[tendon_name]

        ET.indent(self.root, space="  ")
        ET.ElementTree(self.root).write(
            self.output_path, encoding="utf-8", xml_declaration=True
        )
        B.save_json(
            B.REPORT_DIR / "model_wires.json",
            {
                "wires": [
                    {key: value for key, value in wire.items() if key != "tendon_element"}
                    for wire in self.wires
                ],
                "ligaments": self.ligament_tendons,
            },
        )

    def _effective_mass(self, body: str) -> float:
        entry = self.bones.get(body)
        if entry is None:
            return 0.1
        return max(0.02, min(float(entry["mass_kg"]), 0.25))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=B.CONTRACT_PATH)
    parser.add_argument("--muscles", type=Path, default=B.MUSCLE_MANIFEST_PATH)
    parser.add_argument("--ligaments", type=Path, default=B.LIGAMENT_MANIFEST_PATH)
    parser.add_argument("--rom", type=Path, default=B.ROM_PATH)
    parser.add_argument("--output", type=Path, default=B.MODEL_PATH)
    parser.add_argument("--no-collision", action="store_true")
    args = parser.parse_args()

    builder = Bp3dArmBuilder(
        args.contract,
        args.muscles,
        args.output,
        collision=not args.no_collision,
        ligament_path=args.ligaments,
        rom_path=args.rom,
    )
    output = builder.build()
    model = mujoco.MjModel.from_xml_path(str(output))
    print(f"generated: {output}")
    print(f"bodies:        {model.nbody}")
    print(f"joints:        {model.njnt}")
    print(f"DOF:           {model.nv}")
    print(f"tendons:       {model.ntendon}")
    print(f"actuators:     {model.nu}")
    print(f"geoms:         {model.ngeom}")
    print(f"cable wires:   {len(builder.wires)}")
    print(f"ligament ties: {len(builder.ligament_tendons)}")
    slider_actuators = [
        mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, index)
        for index in range(model.nu)
    ]
    print(
        "sliders:       "
        f"{sum(1 for name in slider_actuators if name and name.startswith(SLIDER_ACTUATOR_PREFIX))}"
    )


if __name__ == "__main__":
    main()
