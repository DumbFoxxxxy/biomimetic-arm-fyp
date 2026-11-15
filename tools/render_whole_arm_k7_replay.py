"""Render the published K7 geometry audit as a MuJoCo animation."""

from copy import deepcopy
import struct
import xml.etree.ElementTree as ET

import imageio_ffmpeg
import mujoco
import numpy as np
from PIL import Image, ImageDraw

from build_channel_first_inventory import ROOT, provenance
from build_fyp_showcase_media import font
from build_s1_teres_origin_mount import save
from project_runtime import load_mujoco_xml
from study_s1_local_force_geometry import load_reference


def xml_vector(values):
    return " ".join(map(str, values))


class K7Replay:
    def __init__(self, folder, combo, tpu_sources, cuff_ref, cuffs, spring_trial=None):
        self.folder = folder
        self.with_spring = spring_trial is not None
        audit_name = "whole_arm_k7_spring_relief.json" if self.with_spring else "whole_arm_k7_tpu_moving_cuff_audit.json"
        self.audit_ref = provenance(folder / audit_name)
        self.audit = load_reference(self.audit_ref)
        if (self.audit["combined"] != combo or self.audit["tpu_sources"] != tpu_sources
                or self.audit["moving_cuff_reference"] != cuff_ref
                or len(self.audit["samples"]) != len(cuffs["frames"])):
            raise ValueError("Replay input does not match the published combined audit")
        if self.with_spring:
            manifest = load_reference(provenance(folder / "whole_arm_first_pass_assembly.json"))
            if manifest["k1"] != spring_trial or manifest["k7"] != self.audit_ref:
                raise ValueError("Selected assembly and reported spring check differ")
            self.assembly = manifest["installed"]
            if provenance(ROOT / self.assembly["file"]) != self.assembly:
                raise ValueError("Spring relief assembly has changed")
        else:
            self.assembly = provenance(folder / "whole_arm_latest_neutral.xml")
        self.xml = ET.parse(ROOT / self.assembly["file"]).getroot()
        assets = self.xml.find("asset")
        world = self.xml.find("worldbody")
        body = world.find(".//body[@name='neutral_TPU_reference']")
        self.muscles = {}
        for key, source in tpu_sources.items():
            name = key.replace(":", "_")
            asset = assets.find(f"mesh[@name='{name}']")
            if asset.get("file") != source["file"]:
                raise ValueError(f"Assembly contains a different TPU: {key}")
            file = name + ".msh"
            asset.set("file", file)
            geom = body.find(f"geom[@name='{name}']")
            body.remove(geom)
            geom.set("pos", "0 0 0")
            world.append(geom)
            self.muscles[key] = name, file
        self.cuffs = cuffs
        self.hardware = {}
        for index, (key, source) in enumerate(cuffs["moving_sources"].items()):
            if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
                raise ValueError(f"Cuff hardware changed: {key}")
            name = f"moving_cuff_{index}"
            matches = [mesh for mesh in assets.findall("mesh") if mesh.get("file") == source["file"]]
            for mesh in matches:
                # Replace any neutral copy of the same component with its posed version.
                for parent in world.iter():
                    for geom in parent.findall(f"geom[@mesh='{mesh.get('name')}']"):
                        parent.remove(geom)
                assets.remove(mesh)
            ET.SubElement(assets, "mesh", name=name, file=source["file"], scale=".001 .001 .001")
            ET.SubElement(world, "geom", name=name, type="mesh", mesh=name,
                          rgba=".62 .70 .77 1", conaffinity="0", contype="0", group="1")
            self.hardware[key] = name, np.asarray(source["translation_source_mm"])
        self.width, self.height = 600, 640
        ET.SubElement(ET.SubElement(self.xml, "visual"), "global",
                      offwidth=str(self.width), offheight=str(self.height))
        self.frames = []

    def render(self, index, pose, muscles, paths):
        audit_frame = self.audit["samples"][index]
        hardware_frame = self.cuffs["frames"][index]
        if (pose["input_deg"] != audit_frame["input_deg"] or pose["input_deg"] != hardware_frame["input_deg"]
                or audit_frame["invalid_meshes"] or muscles.keys() != self.muscles.keys()
                or hardware_frame["matrices"].keys() != self.hardware.keys() or len(paths) != 10):
            raise ValueError(f"Incomplete or mismatched replay frame {index}")
        xml = deepcopy(self.xml)
        world = xml.find("worldbody")
        assets = {}
        for key, (vertices, faces) in muscles.items():
            _, filename = self.muscles[key]
            if not np.isfinite(vertices).all():
                raise ValueError(f"Nonfinite TPU vertices: {key}")
            assets[filename] = (struct.pack("<4i", len(vertices), 0, 0, len(faces))
                                + vertices.astype("<f4").tobytes() + faces.astype("<i4").tobytes())
        for key, (name, offset) in self.hardware.items():
            basis = np.asarray(hardware_frame["matrices"][key])
            rotation = np.ascontiguousarray(basis[:, :3])
            np.testing.assert_allclose(rotation @ rotation.T, np.eye(3), atol=1e-8, rtol=0.)
            quaternion = np.empty(4)
            mujoco.mju_mat2Quat(quaternion, rotation.ravel())
            geom = world.find(f"geom[@name='{name}']")
            geom.set("pos", xml_vector(.001 * (rotation @ offset + basis[:, 3])))
            geom.set("quat", xml_vector(quaternion))
        model = load_mujoco_xml(ET.tostring(xml, encoding="unicode"), assets, self.folder)
        model.vis.headlight.ambient[:] = [.50, .50, .50]
        model.vis.headlight.diffuse[:] = [.72, .72, .72]
        data = mujoco.MjData(model)
        data.qpos[:] = pose["qpos"]
        mujoco.mj_forward(model, data)
        for key, (vertices, _) in muscles.items():
            name, _ = self.muscles[key]
            geom, mesh = model.geom(name), model.mesh(name)
            start = mesh.vertadr.item()
            local = model.mesh_vert[start:start + mesh.vertnum.item()]
            xyz = local @ data.geom_xmat[geom.id].reshape(3, 3).T + data.geom_xpos[geom.id]
            for extent in (np.min, np.max):
                np.testing.assert_allclose(extent(xyz, axis=0), extent(vertices, axis=0) * .001, atol=5e-7, rtol=0.)
        rendered = []
        visible_meshes = {model.geom(name).id for name, _ in (*self.muscles.values(), *self.hardware.values())}
        if self.with_spring:
            visible_meshes.add(model.geom("FJ1477_return_spring_neutral").id)
        options = mujoco.MjvOption()
        options.geomgroup[3] = 1
        with mujoco.Renderer(model, self.height, self.width) as renderer:
            for azimuth in (125., 10.):
                view = mujoco.MjvCamera()
                view.lookat[:] = [-.10, -.09, 1.10]
                view.distance, view.azimuth, view.elevation = 1.04, azimuth, 6.
                renderer.update_scene(data, camera=view, scene_option=options)
                scene = renderer.scene
                shown = {geom.objid for geom in scene.geoms[:scene.ngeom]
                         if geom.objtype == mujoco.mjtObj.mjOBJ_GEOM}
                if not visible_meshes <= shown:
                    raise ValueError("Rendered scene omitted a TPU or cuff mesh")
                scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = False
                for key, points in paths.items():
                    color = [.59, .96, .35, 1.] if key.startswith(("FJ1500", "FJ1504")) else [.32, .82, 1., 1.]
                    for a, b in zip(points[:-1], points[1:]):
                        if scene.ngeom >= scene.maxgeom:
                            raise ValueError("No room for the wire overlay")
                        geom = scene.geoms[scene.ngeom]
                        mujoco.mjv_initGeom(geom, mujoco.mjtGeom.mjGEOM_CAPSULE, np.zeros(3),
                                            np.zeros(3), np.eye(3).ravel(), np.array(color, dtype=np.float32))
                        mujoco.mjv_connector(geom, mujoco.mjtGeom.mjGEOM_CAPSULE, .0003, a * .001, b * .001)
                        scene.ngeom += 1
                rendered.append(Image.fromarray(renderer.render()))
        image = Image.new("RGB", (2 * self.width, self.height + 84), "#131e25")
        for panel, pixels in enumerate(rendered):
            image.paste(pixels, (panel * self.width, 52))
        draw = ImageDraw.Draw(image)
        draw.text((20, 13), "K7  /  COMBINED ARM ELEVATION", font=font(23, True), fill="#e6f2f2")
        draw.text((image.width - 20, 17), f"Humerus {pose['humerus_direction_change_deg']:04.1f}° from neutral",
                  font=font(19), fill="#8edfff", anchor="ra")
        draw.text((20, self.height + 56), "MuJoCo geometric replay  /  TPU shape interpolation  /  10 N cuff pose",
                  font=font(18), fill="#c4d2d6")
        self.frames.append(image)
        print(f"Rendered K7 pose {index + 1}/{len(self.audit['samples'])}", flush=True)

    def finish(self):
        if len(self.frames) != len(self.audit["samples"]) or len(self.frames) < 2:
            raise ValueError("K7 animation must include the entire audited range")
        stem = "whole_arm_k7_first_pass" if self.with_spring else "whole_arm_k7_combined"
        movie = self.folder / f"{stem}.mp4"
        poster = self.folder / f"{stem}.png"
        self.frames[-1].save(poster, optimize=True)
        timeline = self.frames + self.frames[-2::-1]
        writer = imageio_ffmpeg.write_frames(str(movie), self.frames[0].size, fps=6, quality=8,
                                             codec="libx264", macro_block_size=2,
                                             output_params=["-movflags", "+faststart"])
        writer.send(None)
        try:
            for frame in timeline:
                writer.send(frame.tobytes())
        finally:
            writer.close()
        save(self.folder / f"{stem}_media.json", dict(
            assembly=self.assembly, audit=self.audit_ref,
            files=[provenance(file) for file in (movie, poster)],
            forward_samples=len(self.frames), frames=len(timeline), playback_fps=6,
            scope="Prescribed K7 geometry playback with eight interpolated TPU meshes, ten wire overlays, and both "
                  "cuff assemblies in their local 10 N probe poses. "
                  + ("Includes the checked neutral spring shape fixed to the moving humerus. " if self.with_spring else "")
                  + "No force-driven integration."))
