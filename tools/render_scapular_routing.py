"""Render the audited K6 wire comparison as a MuJoCo video and peak frame."""

import argparse
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import imageio_ffmpeg
import mujoco
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from build_channel_first_inventory import DEST, ROOT, provenance, read

sys.path.insert(0, str(ROOT))
from project_runtime import load_mujoco_model


WIDTH, HEIGHT = 576, 440
FOLDER = DEST.parent / "whole_arm_first_pass"
OUTPUT = FOLDER / "scapular_control"
RESULT = OUTPUT / "serratus_guide_refinement.json"
FPS = 12
BG, INK, SUBTLE = "#142126", "#edf7f4", "#aabfc0"
OLD = (1., .34, .23, 1.)
NEW = (.22, .95, .76, 1.)
GUIDE = (.72, 1., .25, 1.)


def font(size, bold=False):
    path = Path("C:/Windows/Fonts") / ("msyhbd.ttc" if bold else "msyh.ttc")
    if path.exists():
        return ImageFont.truetype(str(path), size)
    return ImageFont.load_default(size)


def tube(scene, start, end, rgba, radius=.00125):
    geom = scene.geoms[scene.ngeom]
    mujoco.mjv_initGeom(geom, mujoco.mjtGeom.mjGEOM_CAPSULE, np.zeros(3),
                        np.zeros(3), np.eye(3).ravel(), np.array(rgba, dtype=np.float32))
    mujoco.mjv_connector(geom, mujoco.mjtGeom.mjGEOM_CAPSULE, radius, start, end)
    scene.ngeom += 1


def marker(scene, point, rgba, radius=.003):
    geom = scene.geoms[scene.ngeom]
    mujoco.mjv_initGeom(geom, mujoco.mjtGeom.mjGEOM_SPHERE, np.full(3, radius),
                        point, np.eye(3).ravel(), np.array(rgba, dtype=np.float32))
    scene.ngeom += 1


def setup_colors(model, xml, whole_arm=False):
    rib_names = {}
    for mesh in ET.parse(xml).findall("./asset/mesh"):
        filename = mesh.get("file", "").replace("\\", "/").rsplit("/", 1)[-1]
        if filename in {"FJ3342.obj", "FJ3344.obj"}:
            rib_names[model.mesh(mesh.get("name")).id] = filename
    if set(rib_names.values()) != {"FJ3342.obj", "FJ3344.obj"}:
        raise ValueError("Replay does not contain the two reference ribs")
    for geom in range(model.ngeom):
        rgb = model.geom_rgba[geom].copy()
        if model.geom_type[geom] == mujoco.mjtGeom.mjGEOM_MESH:
            rib = rib_names.get(model.geom_dataid[geom])
            if rib == "FJ3344.obj":
                rgb = [.92, .70, .39, .88]
            elif rib:
                rgb = [.69, .63, .51, .76]
            else:
                owner = model.body(model.geom_bodyid[geom]).name
                if owner == "thorax":
                    rgb = [.58, .65, .64, .57 if whole_arm else .36]
                elif owner == "scapula":
                    rgb = [.79, .85, .81, .93 if whole_arm else .72]
                else:
                    rgb = [.61, .68, .71, .82]
        model.geom_rgba[geom] = rgb
    model.vis.headlight.ambient[:] = [.47, .47, .47]
    model.vis.headlight.diffuse[:] = [.68, .68, .68]


def camera(whole_arm=False, azimuth=218.):
    cam = mujoco.MjvCamera()
    cam.lookat[:] = [-.096, -.106, 1.085] if whole_arm else [-.108, -.058, 1.249]
    cam.distance, cam.azimuth, cam.elevation = (.93 if whole_arm else .34), azimuth, (7. if whole_arm else 14.)
    return cam


def project(scene, point):
    left, right = scene.camera
    delta = point - (left.pos + right.pos) / 2.
    lateral = np.cross(left.forward, left.up)
    scale = left.frustum_near / (delta @ left.forward)
    height = left.frustum_top - left.frustum_bottom
    width = height * WIDTH / HEIGHT
    x = WIDTH * (.5 + (scale * (delta @ lateral) - left.frustum_center) / width)
    y = HEIGHT * (left.frustum_top - scale * (delta @ left.up)) / height
    return x, y


def panel(renderer, data, cam, options, points, color, support=None, arm_reference=None):
    renderer.update_scene(data, camera=cam, scene_option=options)
    scene = renderer.scene
    scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = False
    for a, b in zip(points, points[1:]):
        tube(scene, a, b, color)
    if support is not None:
        tube(scene, support, points[2], GUIDE, .0009)
    for index, point in enumerate(points):
        changed = support is not None and index == 2
        marker(scene, point, GUIDE if changed else color, .0042 if changed else .0019)
    pixels = renderer.render().astype(np.float32)
    background = np.array([17., 29., 35.], dtype=np.float32)
    image = Image.fromarray((background + pixels * (1. - background / 255.)).astype(np.uint8))
    projected = [project(scene, point) for point in points]
    labels = ImageDraw.Draw(image, "RGBA")
    labels.line(projected, fill=(*np.round(np.array(color[:3]) * 255).astype(int), 164), width=2)
    if support is not None:
        x, y = projected[2]
        accent = "#d4ff83"
        labels.ellipse((x - 6, y - 6, x + 6, y + 6), outline=accent, width=2)
        labels.line((x + 5, y - 4, x + 32, y - 34, x + 50, y - 34), fill=accent, width=1)
        labels.rounded_rectangle((x + 50, y - 46, x + 126, y - 22), radius=4, fill=(20, 33, 38, 230))
        labels.text((x + 56, y - 43), "新增导向", font=font(14, True), fill=accent)
    if arm_reference:
        neutral, current = ([np.array(project(scene, point)) for point in path] for path in arm_reference)
        for a, b in zip(neutral, neutral[1:]):
            vector = b - a
            spans = max(1, round(np.linalg.norm(vector) / 9.))
            for index in range(spans):
                start = a + vector * index / spans
                end = a + vector * (index + .52) / spans
                labels.line((*start, *end), fill=(218, 235, 231, 138), width=1)
        for point in current:
            x, y = point
            labels.ellipse((x - 3, y - 3, x + 3, y + 3), fill=(233, 252, 248, 231))
    return image


def compose(left, right, frame, trace, peak, whole_arm=False):
    image = Image.new("RGB", (2 * WIDTH + 24, HEIGHT + 170), BG)
    draw = ImageDraw.Draw(image)
    draw.text((16, 12), "整臂运动回放" if whole_arm else "前锯肌钢丝走线对照", font=font(25, True), fill=INK)
    draw.text((image.width - 16, 15), f"上回旋 {frame['sc_deg']:.2f}°", font=font(24, True),
              fill="#68edce", anchor="ra")
    draw.text((16, 50), "整臂 · 斜前视角" if whole_arm else "原走线 · 同一姿态的几何对照",
              font=font(15, True), fill="#edf7f4" if whole_arm else "#ff987f")
    draw.text((WIDTH + 20, 50), "整臂 · 后侧视角" if whole_arm else "新增一处导向 · MuJoCo 主动牵引回放",
              font=font(15, True), fill="#edf7f4" if whole_arm else "#68edce")
    image.paste(left, (8, 76))
    image.paste(right, (WIDTH + 16, 76))
    legend_y = HEIGHT + 84
    if whole_arm:
        legend = "青绿 前锯肌走线   黄绿 新增导向   白色虚线 中立位手臂参考"
        scope = f"盂肱与肘固定 · 肩关节中心上移 {frame['superior_mm']:.1f} mm"
    else:
        legend = "橙红 原钢丝   青绿 新钢丝   黄绿 新导向   琥珀色 右第六肋   路径叠加"
        scope = f"K6 单轴 · 理想重力预平衡 · 峰值 {peak:.2f}°"
    draw.text((16, legend_y), legend, font=font(14), fill=INK)
    draw.text((image.width - 16, legend_y), scope, font=font(14),
              fill=SUBTLE, anchor="ra")
    x0, x1 = 17, image.width - 17
    y0, y1 = HEIGHT + 119, HEIGHT + 151
    times = trace[-1]["time_s"]

    def chart(sample):
        return (x0 + (x1 - x0) * sample["time_s"] / times,
                y1 - (y1 - y0) * sample["sc_deg"] / peak)

    draw.line((x0, y1, x1, y1), fill="#465659", width=1)
    draw.line([chart(sample) for sample in trace], fill="#79c4b4", width=2)
    x, y = chart(frame)
    draw.line((x, y0, x, y1), fill="#9cb7b5")
    draw.ellipse((x - 3, y - 3, x + 3, y + 3), fill=INK)
    return image


def forearm_reference(model, data, xml_path):
    xml = ET.parse(xml_path)
    asset = next(mesh for mesh in xml.findall("./asset/mesh")
                 if mesh.get("file", "").endswith("direct_forearm_seated.obj"))
    node = xml.find(f".//geom[@mesh='{asset.get('name')}']")
    if node is None:
        raise ValueError("Forearm reference mesh is missing")
    geom = model.geom(node.get("name"))
    mesh = model.mesh(asset.get("name"))
    start = mesh.vertadr.item()
    vertices = model.mesh_vert[start:start + mesh.vertnum.item()]
    world = vertices @ data.geom_xmat[geom.id].reshape(3, 3).T + data.geom_xpos[geom.id]
    tip = world[world[:, 2] <= world[:, 2].min() + .012].mean(axis=0)
    forearm = data.body("forearm")
    return (tip - forearm.xpos) @ forearm.xmat.reshape(3, 3)


def arm_points(data, distal):
    forearm = data.body("forearm")
    return [data.body("humerus").xpos.copy(), forearm.xpos.copy(),
            forearm.xmat.reshape(3, 3) @ distal + forearm.xpos]


def main(preview_only=False, whole_arm=False):
    bundle = read(RESULT)
    model_path = ROOT / bundle["model"]["file"]
    if provenance(model_path) != bundle["model"]:
        raise ValueError("The replay model changed after trajectory validation")
    model = load_mujoco_model(model_path)
    model.vis.global_.offwidth = max(model.vis.global_.offwidth, WIDTH)
    model.vis.global_.offheight = max(model.vis.global_.offheight, HEIGHT)
    setup_colors(model, model_path, whole_arm)
    trace = bundle["best"]["frames"]
    peak_frame = max(trace, key=lambda frame: frame["sc_deg"])
    source_path = ROOT / bundle["source"]["file"]
    if provenance(source_path) != bundle["source"]:
        raise ValueError("The original routing changed after trajectory validation")
    source = read(source_path)
    old = source["routes"]["serratus_middle"]["sites"]
    origins = read(FOLDER / "whole_arm_shared_scene_check.json")["neutral_body_origins_mm"]
    old_local = [(item["owner"], (np.asarray(item["port_mm"]) - origins[item["owner"]]) * .001) for item in old]
    base = np.asarray(bundle["candidate"]["sites"][2]["base_mm"]) * .001
    site_names = [f"serratus_middle_{i}" for i in range(len(bundle["candidate"]["sites"]))]
    opts = mujoco.MjvOption()
    opts.geomgroup[2:4] = 0
    opts.sitegroup[2] = 0
    opts.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = 0
    cam = camera(whole_arm, 125. if whole_arm else 218.)
    side_cam = camera(whole_arm, 218.)
    data = mujoco.MjData(model)
    data.qpos[:] = trace[0]["qpos"]
    mujoco.mj_forward(model, data)
    if whole_arm:
        distal = forearm_reference(model, data, model_path)
        neutral = arm_points(data, distal)
    with mujoco.Renderer(model, HEIGHT, WIDTH) as renderer:
        def render(frame):
            data.qpos[:] = frame["qpos"]
            mujoco.mj_forward(model, data)
            if abs(data.ten_length[0] * 1000. - frame["length_mm"][0]) > 1e-6:
                raise ValueError("Cable geometry differs from the recorded simulation")
            refined = [data.site(name).xpos.copy() for name in site_names]
            if whole_arm:
                reference = (neutral, arm_points(data, distal))
                left, right = [panel(renderer, data, viewpoint, opts, refined, NEW,
                                     support=base, arm_reference=reference) for viewpoint in (cam, side_cam)]
                return compose(left, right, frame, trace, peak_frame["sc_deg"], whole_arm=True)
            baseline = []
            for owner, point in old_local:
                body = data.body(owner)
                baseline.append(body.xmat.reshape(3, 3) @ point + body.xpos)
            left = panel(renderer, data, cam, opts, baseline, OLD)
            right = panel(renderer, data, cam, opts, refined, NEW, support=base)
            return compose(left, right, frame, trace, peak_frame["sc_deg"])

        poster = render(peak_frame)
        stem = "serratus_whole_arm" if whole_arm else "serratus_routing"
        poster.save(OUTPUT / f"{stem}_peak.png")
        if preview_only:
            return
        destination = OUTPUT / f"{stem}_replay.mp4"
        writer = imageio_ffmpeg.write_frames(str(destination), poster.size, fps=FPS, quality=8, macro_block_size=2,
                                            output_params=["-movflags", "+faststart"])
        writer.send(None)
        try:
            times = np.asarray([frame["time_s"] for frame in trace])
            indices = [int(np.argmin(abs(times - t))) for t in np.linspace(0., times[-1], FPS * 4 + 1)]
            saved_peak = False
            for index in indices:
                frame = trace[index]
                if frame["time_s"] >= peak_frame["time_s"] and not saved_peak:
                    for _ in range(FPS):
                        writer.send(poster.tobytes())
                    saved_peak = True
                writer.send(render(frame).tobytes())
        finally:
            writer.close()
        print(f"Rendered {destination.name}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poster", action="store_true", help="Render only the peak-angle preview")
    parser.add_argument("--whole-arm", action="store_true", help="Show the new routing in two whole-arm views")
    args = parser.parse_args()
    main(args.poster, args.whole_arm)
