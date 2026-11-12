"""Render GitHub showcase media from the audited MuJoCo shoulder–elbow replay."""

import csv
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import imageio_ffmpeg
import mujoco
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from build_channel_first_inventory import DEST, ROOT, provenance, read
from study_s1_local_force_geometry import load_reference
from view_s1_series_elbow_start_pose import load_scene


BACKGROUND = np.array([11, 26, 34], dtype=float)
INK = "#e9f3f1"
MUTED = "#a4b8bd"
TEAL = "#5ce1bc"


def font(size, bold=False):
    windows = Path("C:/Windows/Fonts") / ("segoeuib.ttf" if bold else "segoeui.ttf")
    try:
        return ImageFont.truetype(str(windows if windows.exists() else
                                      "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"), size)
    except OSError:
        return ImageFont.load_default(size)


def camera(azimuth, height=1.112, distance=.74):
    view = mujoco.MjvCamera()
    view.lookat[:] = [-.16, -.053, height]
    view.distance = distance
    view.azimuth = azimuth
    view.elevation = -5.
    return view


def recolor(model, assembly):
    for i, (key, _) in enumerate(assembly["sources"].items()):
        geom = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, f"replay_{i}")
        if geom < 0:
            continue
        if key.startswith("thorax_reference:"):
            color = [.49, .60, .64, .57]
        elif key == "bone:forearm":
            color = [.29, .91, .72, 1.]
        elif key == "elbow:elbow_shoulder_pin":
            color = [.96, .61, .25, 1.]
        else:
            color = [.98, .88, .71, 1.]
        model.geom_rgba[geom] = color
    model.vis.headlight.ambient[:] = [.52, .52, .52]
    model.vis.headlight.diffuse[:] = [.72, .72, .72]


def render_panel(renderer, data, view):
    renderer.update_scene(data, camera=view)
    renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = False
    pixels = renderer.render().astype(float)
    tinted = BACKGROUND + pixels * (1. - BACKGROUND / 255.)
    return Image.fromarray(tinted.clip(0, 255).astype(np.uint8))


def compose(left, right, time, angle, chart):
    image = Image.new("RGB", (960, 576), "#10212b")
    draw = ImageDraw.Draw(image)
    draw.text((25, 17), "SHOULDER–ELBOW SIMULATION", font=font(19, True), fill=INK)
    draw.text((935, 18), f"{time:04.2f} s   ·   {max(0., angle):05.1f}° elbow", font=font(18, True),
              fill=TEAL, anchor="ra")
    image.paste(left, (24, 52))
    image.paste(right, (488, 52))
    for x, label in ((39, "Side view"), (503, "Oblique view")):
        draw.text((x, 65), label, font=font(16, True), fill=INK)
    draw.text((24, 425), "Elbow flexion", font=font(16, True), fill=INK)
    draw.text((936, 426), "MuJoCo replay  /  model prediction", font=font(14), fill=MUTED, anchor="ra")
    for level in (0., 50., 100.):
        y = 542 - level * .8
        draw.line((67, y, 936, y), fill="#314750")
        draw.text((59, y), f"{level:.0f}°", font=font(12), fill=MUTED, anchor="rm")
    draw.line(chart, fill=TEAL, width=3, joint="curve")
    cursor_x = 67 + 869 * time / 4.
    cursor_y = 542 - max(0., angle) * .8
    draw.line((cursor_x, 455, cursor_x, 542), fill="#527f80")
    draw.ellipse((cursor_x-4, cursor_y-4, cursor_x+4, cursor_y+4), fill="#f1fcf8")
    for t in range(5):
        draw.text((67 + 869 * t / 4., 549), f"{t} s", font=font(12), fill=MUTED, anchor="ma")
    return image


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    report_path = folder / "gravity_balanced_release_timing.json"
    report = read(report_path)
    case = report["selected"]
    trial = report["cases"][case]
    if not trial["fully_audited"] or trial["first_collision"] is not None:
        raise ValueError("Showcase requires the checked rigid trajectory")
    source = load_reference(report["baseline"])
    assembly = load_reference(source["geometry_assembly"])
    model, frames, _ = load_scene(report, case)
    recolor(model, assembly)
    trace = ROOT / trial["trace"]["file"]
    with trace.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(frames):
        raise ValueError("Replay and source trace lengths differ")
    times = np.array([float(row["time_s"]) for row in rows])
    angles = np.array([float(row["elbow_deg"]) for row in rows])
    chart = list(zip(67 + times[::4] * (869 / 4.), 542 - np.maximum(0., angles[::4]) * .8))
    peak = int(np.argmax(angles))
    fps = 24
    sample_indices = np.searchsorted(times, np.arange(round(times[-1] * fps) + 1) / fps).clip(0, len(frames)-1)
    out = ROOT / "docs/showcase/media"
    out.mkdir(parents=True, exist_ok=True)
    model.vis.global_.offwidth, model.vis.global_.offheight = 448, 360
    data = mujoco.MjData(model)
    views = (camera(8.), camera(139., height=1.143, distance=.80))
    scratch = ROOT / "tmp"
    scratch.mkdir(exist_ok=True)
    with TemporaryDirectory(prefix="fyp-showcase-", dir=scratch) as temporary, mujoco.Renderer(model, height=360, width=448) as renderer:
        staging = Path(temporary)

        def at(index):
            data.qpos[:] = frames[index]
            mujoco.mj_forward(model, data)
            images = [render_panel(renderer, data, view) for view in views]
            return compose(*images, times[index], angles[index], chart)

        at(peak).save(out / "motion-poster.png", optimize=True)
        for position, index in enumerate(sample_indices):
            at(index).save(staging / f"frame-{position:04d}.png")
        command = [imageio_ffmpeg.get_ffmpeg_exe(), "-hide_banner", "-loglevel", "error", "-y",
                   "-framerate", str(fps), "-i", str(staging / "frame-%04d.png")]
        subprocess.run([*command, "-c:v", "libx264", "-crf", "20", "-pix_fmt", "yuv420p",
                        "-movflags", "+faststart", str(out / "shoulder-elbow.mp4")], check=True)
        subprocess.run([*command, "-filter_complex",
                        "fps=16,scale=800:-2:flags=lanczos,split[a][b];[a]palettegen=max_colors=96[p];[b][p]paletteuse=dither=bayer",
                        "-loop", "0", str(out / "shoulder-elbow.gif")], check=True)
    (out / "source.json").write_text(json.dumps(dict(
        replay=provenance(report_path), trajectory=trial["trace"], frames=len(frames),
        peak_elbow_deg=trial["model_peak_elbow_deg"], peak_sample_time_s=times[peak],
        scope="Recorded MuJoCo rigid-body equivalent-drive replay with anatomical meshes. No physical footage or TPU deformation.",
        media=[provenance(out / name) for name in ("motion-poster.png", "shoulder-elbow.mp4", "shoulder-elbow.gif")]),
        indent=2) + "\n", encoding="utf-8")
    print(f"Rendered {len(sample_indices)} frames from {len(frames)} verified samples")


if __name__ == "__main__":
    main()
