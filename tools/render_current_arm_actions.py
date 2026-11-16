"""Display the audited whole-arm actions in MuJoCo or export a two-view replay."""

import argparse
import time
import xml.etree.ElementTree as ET

import imageio_ffmpeg
import mujoco
import numpy as np
from PIL import Image, ImageDraw

from render_scapular_routing import (
    BG, FOLDER, GUIDE, HEIGHT, INK, NEW, RESULT, ROOT, SUBTLE, WIDTH,
    arm_points, camera, font, forearm_reference, marker, panel, project, provenance, read, setup_colors, tube,
)
from project_runtime import load_mujoco_xml


DEST = FOLDER / "action_showcase"
FPS = 15
ACCENT = "#68edce"
AXIS_MARKER = (1., .73, .29, 1.)
PULL_LINES = {
    "trapezius_superior": ((.46, .76, 1., 1.), "#82c5ff", "上斜方肌"),
    "FJ1536_g2": ((.87, .63, 1., 1.), "#db9eff", "大菱形肌"),
}


class Replay:
    def __init__(self):
        record = read(DEST / "current_actions.json")
        self.actions = {key: action for key, action in record["actions"].items() if action["playback"]}
        xml_path = ROOT / record["assembly"]["file"]
        cap = record["cap"]
        cap_path = ROOT / cap["file"]
        for path, source in [(xml_path, record["assembly"]), (cap_path, cap), (RESULT, record["routing"])]:
            if provenance(path)["sha256"] != source["sha256"]:
                raise ValueError(f"Replay asset changed after validation: {path}")
        relief = read(FOLDER / "whole_arm_k1_cap_relief.json")["cap_relief"]
        xml = ET.parse(xml_path).getroot()
        xml.find(f"./asset/mesh[@file='{relief['original']['file']}']").set("file", cap["file"])
        self.model = load_mujoco_xml(ET.tostring(xml, encoding="unicode"), base_dir=xml_path.parent)
        self.model.vis.global_.offwidth = max(self.model.vis.global_.offwidth, WIDTH)
        self.model.vis.global_.offheight = max(self.model.vis.global_.offheight, HEIGHT)
        setup_colors(self.model, xml_path, whole_arm=True)
        self.data = mujoco.MjData(self.model)
        mujoco.mj_forward(self.model, self.data)
        self.distal = forearm_reference(self.model, self.data, xml_path)
        self.neutral = arm_points(self.data, self.distal)
        axis = self.neutral[1] - self.neutral[0]
        axis /= np.linalg.norm(axis)
        elbow = self.model.joint("elbow")
        tangent = self.data.body("forearm").xmat.reshape(3, 3) @ elbow.axis
        radial = np.cross(tangent, axis)
        self.axial_marker = self.data.body("humerus").xmat.reshape(3, 3).T @ (
            .053 * radial / np.linalg.norm(radial)
        )
        self.options = mujoco.MjvOption()
        self.options.geomgroup[2:4] = 0
        self.options.sitegroup[2] = 0
        self.options.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = 0
        origins = read(FOLDER / "whole_arm_shared_scene_check.json")["neutral_body_origins_mm"]
        routing = read(RESULT)
        sites = routing["candidate"]["sites"]
        source = ROOT / routing["source"]["file"]
        if provenance(source) != routing["source"]:
            raise ValueError("Shoulder pull lines have changed since validation")
        saved_routes = read(source)["routes"]

        def local_sites(items):
            return [(site["owner"], (np.asarray(site["port_mm"]) - origins[site["owner"]]) * .001)
                    for site in items]

        self.sites = local_sites(sites)
        self.pull_sites = {key: local_sites(saved_routes[key]["sites"]) for key in PULL_LINES}
        self.support = (sites[2]["owner"], (np.asarray(sites[2]["base_mm"]) - origins[sites[2]["owner"]]) * .001)
        self.cameras = [camera(True, azimuth) for azimuth in (125., 218.)]

    def set_frame(self, frame):
        self.data.qpos[:] = frame["qpos"]
        mujoco.mj_forward(self.model, self.data)

    def world(self, site):
        owner, point = site
        body = self.data.body(owner)
        return body.xmat.reshape(3, 3) @ point + body.xpos

    def points(self):
        return [self.world(site) for site in self.sites], self.world(self.support)

    def pulling_lines(self, action):
        keys = {"K6": ("trapezius_superior",), "retraction": tuple(PULL_LINES)}.get(action, ())
        return [([self.world(site) for site in self.pull_sites[key]], *PULL_LINES[key]) for key in keys]

    def rotation_marker(self):
        elbow = self.data.body("forearm").xpos
        radial = self.data.body("humerus").xmat.reshape(3, 3) @ self.axial_marker
        return elbow, elbow + radial

    def frames(self, key):
        frames = self.actions[key]["playback"]
        if self.actions[key]["status"] == "active_return":
            return frames
        # Visual timing only. Preserve the peak with a short pause, and play the
        # same geometric path in reverse.
        peak = len(frames) // 2
        outbound = [frames[round(position)] for position in np.linspace(0, peak, 22)]
        return [*outbound, *[outbound[-1]] * 7, *outbound[-2::-1]]

    def display(self, key):
        from mujoco import viewer

        keys = list(self.actions)
        state = {"index": keys.index(key), "epoch": time.monotonic(), "paused_at": None}

        def press(code):
            now = time.monotonic()
            if ord("1") <= code <= ord("9") and code - ord("1") < len(keys):
                state.update(index=code - ord("1"), epoch=now, paused_at=None)
                print(self.actions[keys[state["index"]]]["label"], flush=True)
            if code == ord("P"):
                paused_at = state["paused_at"]
                if paused_at is None:
                    state["paused_at"] = now
                else:
                    state.update(epoch=state["epoch"] + now - paused_at, paused_at=None)

        print("数字键切换动作，P 暂停，鼠标旋转和缩放。", flush=True)
        for index, entry in enumerate(keys, start=1):
            print(index, entry, self.actions[entry]["label"], flush=True)
        traces = {entry: self.frames(entry) for entry in keys}
        with viewer.launch_passive(self.model, self.data, key_callback=press) as view:
            with view.lock():
                view.cam.lookat[:] = self.cameras[0].lookat
                view.cam.distance = self.cameras[0].distance
                view.cam.azimuth = self.cameras[0].azimuth
                view.cam.elevation = self.cameras[0].elevation
                view.opt.geomgroup[:] = self.options.geomgroup
                view.opt.sitegroup[:] = self.options.sitegroup
                view.opt.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = 0
            while view.is_running():
                entry = keys[state["index"]]
                frames = traces[entry]
                elapsed = (state["paused_at"] or time.monotonic()) - state["epoch"]
                duration = frames[-1].get("time_s", (len(frames) - 1) / FPS)
                index = round((len(frames) - 1) * (max(0., elapsed) % (duration + .6)) / duration)
                with view.lock():
                    self.set_frame(frames[min(index, len(frames) - 1)])
                    points, support = self.points()
                    scene = view.user_scn
                    scene.ngeom = 0
                    for a, b in zip(points, points[1:]):
                        tube(scene, a, b, NEW)
                    for path, color, _, _ in self.pulling_lines(entry):
                        for a, b in zip(path, path[1:]):
                            tube(scene, a, b, color)
                    tube(scene, support, points[2], GUIDE, .0009)
                    marker(scene, points[2], GUIDE)
                    if entry in {"K4", "K5"}:
                        a, b = self.rotation_marker()
                        tube(scene, a, b, AXIS_MARKER)
                        marker(scene, b, AXIS_MARKER)
                view.sync()
                time.sleep(.02)

    def picture(self, renderer, key, frame):
        action = self.actions[key]
        self.set_frame(frame)
        points, support = self.points()
        pulling = self.pulling_lines(key)
        current = arm_points(self.data, self.distal)
        image = Image.new("RGB", (2 * WIDTH + 24, HEIGHT + 144), BG)
        draw = ImageDraw.Draw(image)
        tag = f"{key}  ·  " if key.startswith("K") else ""
        draw.text((20, 13), tag + action["label"], font=font(28, True), fill=INK)
        reading = f"{frame['input_value']:.2f} {action['units']}"
        draw.text((image.width - 20, 15), reading, font=font(27, True), fill=ACCENT, anchor="ra")
        for index, viewpoint in enumerate(self.cameras):
            x = 8 + index * (WIDTH + 8)
            picture = panel(renderer, self.data, viewpoint, self.options, points, NEW,
                            support=support, arm_reference=(self.neutral, current))
            overlay = ImageDraw.Draw(picture)
            for path, _, color, _ in pulling:
                overlay.line([project(renderer.scene, point) for point in path], fill=color, width=2)
            if key in {"K4", "K5"}:
                a, b = [project(renderer.scene, point) for point in self.rotation_marker()]
                overlay.line([a, b], fill="#ffc062", width=3)
                bx, by = b
                overlay.ellipse((bx - 4, by - 4, bx + 4, by + 4), fill="#ffc062")
            image.paste(picture, (x, 73))
            draw.text((x + 12, 53), ("斜前视角", "后侧视角")[index], font=font(14), fill=SUBTLE)
        active = action["status"] == "active_return"
        label = "MuJoCo 钢丝牵引并回位 · 理想重力预平衡" if active else "规定关节姿态 · 几何回放"
        if active:
            label += f" · 最大线张力 {max(frame['tension_N']):.1f} N"
        draw.text((20, HEIGHT + 89), label, font=font(19, True), fill=INK)
        legend = "青绿 前锯肌走线     黄绿 新导向     白色虚线 中立位手臂"
        if key in {"K4", "K5"}:
            legend += "     橙色 肱骨转动标记"
        draw.text((20, HEIGHT + 119), legend, font=font(14), fill=SUBTLE)
        for index, (_, _, color, name) in enumerate(pulling):
            draw.text((490 + 95 * index, HEIGHT + 119), name, font=font(14), fill=color)
        draw.text((image.width - 20, HEIGHT + 119),
                  f"本段峰值 {action['screened_extent']:.2f} {action['units']}",
                  font=font(14), fill=SUBTLE, anchor="ra")
        return image

    def export(self, keys, poster_only):
        with mujoco.Renderer(self.model, HEIGHT, WIDTH) as renderer:
            cards = []
            writer = None
            try:
                for key in keys:
                    action = self.actions[key]
                    peak = max(action["playback"], key=lambda frame: frame["input_value"])
                    card = self.picture(renderer, key, peak)
                    cards.append(card)
                    if poster_only:
                        continue
                    if writer is None:
                        name = "current_actions" if len(keys) > 1 else key
                        writer = imageio_ffmpeg.write_frames(str(DEST / f"{name}.mp4"), card.size, fps=FPS,
                                                            quality=8, macro_block_size=2,
                                                            output_params=["-movflags", "+faststart"])
                        writer.send(None)
                    frames = self.frames(key)
                    if action["status"] == "active_return":
                        times = [frame["time_s"] for frame in frames]
                        offsets = np.linspace(0., times[-1], round(FPS * times[-1]) + 1)
                        frames = [frames[int(np.argmin(np.abs(np.asarray(times) - offset)))] for offset in offsets]
                    previous, pixels = None, None
                    for frame in frames:
                        if frame is not previous:
                            pixels = (card if frame is peak else self.picture(renderer, key, frame)).tobytes()
                            previous = frame
                        writer.send(pixels)
                    print(f"Rendered {key}", flush=True)
            finally:
                if writer is not None:
                    writer.close()

        # Arrange the peak frames for a single glance at all replayable actions.
        width, height = 448, 222
        columns = min(len(cards), 3)
        rows = (len(cards) + columns - 1) // columns
        sheet = Image.new("RGB", (columns * (width + 8) + 8, rows * (height + 8) + 8), "#0b1518")
        for index, card in enumerate(cards):
            column, row = index % columns, index // columns
            sheet.paste(card.resize((width, height), Image.Resampling.LANCZOS),
                        (8 + column * (width + 8), 8 + row * (height + 8)))
        name = "current_actions" if len(keys) > 1 else keys[0]
        sheet.save(DEST / f"{name}_overview.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", help="Initial action in the viewer or single export, e.g. K6")
    parser.add_argument("--view", action="store_true", help="Open the rotatable MuJoCo viewer")
    parser.add_argument("--poster", action="store_true", help="Export only the overview image")
    arguments = parser.parse_args()
    replay = Replay()
    if arguments.action is not None and arguments.action not in replay.actions:
        parser.error(f"Available actions: {', '.join(replay.actions)}")
    if arguments.view:
        replay.display(arguments.action or next(iter(replay.actions)))
    else:
        replay.export([arguments.action] if arguments.action else list(replay.actions), arguments.poster)
