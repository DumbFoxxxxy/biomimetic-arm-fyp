"""Chinese MuJoCo GUI for V6 trajectory playback and wrist-path comparison."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import argparse
import copy
import json
from project_runtime import WorkspacePath as Path
import sys
import time
import tkinter as tk
from tkinter import ttk

import mujoco
import mujoco.viewer
import numpy as np


ROOT = Path(__file__).resolve().parent
V5_ROOT = ROOT.parent / "cable_design_v5"
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V5_ROOT), str(V2_ROOT)]

from trajectory_v6 import TrajectoryRunner, load_inputs  # noqa: E402


class PathCanvas(tk.Canvas):
    def __init__(self, master, **kwargs):
        super().__init__(master, background="#f7f8fa", **kwargs)
        self.target = []
        self.actual = []

    def reset(self) -> None:
        self.target = []
        self.actual = []
        self.delete("all")

    def add(self, target: np.ndarray, actual: np.ndarray) -> None:
        self.target.append(np.asarray(target))
        self.actual.append(np.asarray(actual))
        self.redraw()

    def redraw(self) -> None:
        self.delete("all")
        width = max(int(self.winfo_width()), 200)
        height = max(int(self.winfo_height()), 200)
        self.create_text(
            10,
            10,
            text="腕端 X-Y 路径",
            anchor="nw",
            fill="#28313d",
            font=("Microsoft YaHei UI", 9, "bold"),
        )
        points = self.target + self.actual
        if len(points) < 2:
            return
        values = np.asarray(points)
        minimum = values[:, :2].min(axis=0)
        maximum = values[:, :2].max(axis=0)
        span = np.maximum(maximum - minimum, 20.0)
        margin = 24

        def convert(point):
            x = margin + (point[0] - minimum[0]) / span[0] * (
                width - 2 * margin
            )
            y = height - margin - (point[1] - minimum[1]) / span[1] * (
                height - 2 * margin
            )
            return x, y

        for values_set, color, label in (
            (self.target, "#147d64", "目标"),
            (self.actual, "#d66a1f", "实际"),
        ):
            if len(values_set) < 2:
                continue
            coordinates = []
            for point in values_set:
                coordinates.extend(convert(point))
            self.create_line(
                *coordinates,
                fill=color,
                width=2,
                smooth=True,
            )
            self.create_text(
                width - 60,
                10 + (0 if label == "目标" else 20),
                text=label,
                fill=color,
                anchor="nw",
                font=("Microsoft YaHei UI", 8),
            )


class ErrorCanvas(tk.Canvas):
    def __init__(self, master, **kwargs):
        super().__init__(master, background="#f7f8fa", **kwargs)
        self.values = []

    def reset(self) -> None:
        self.values = []
        self.delete("all")
        self.create_text(
            10,
            10,
            text="路径误差 (mm)",
            anchor="nw",
            fill="#28313d",
            font=("Microsoft YaHei UI", 9, "bold"),
        )

    def add(self, value: float) -> None:
        self.values.append(float(value))
        self.redraw()

    def redraw(self) -> None:
        width = max(int(self.winfo_width()), 200)
        height = max(int(self.winfo_height()), 120)
        self.delete("all")
        self.create_text(
            10,
            8,
            text="路径误差 (mm)",
            anchor="nw",
            fill="#28313d",
            font=("Microsoft YaHei UI", 8, "bold"),
        )
        if len(self.values) < 2:
            return
        maximum = max(25.0, max(self.values) * 1.1)
        coordinates = []
        for index, value in enumerate(self.values):
            x = 12 + index / max(1, len(self.values) - 1) * (width - 24)
            y = height - 12 - value / maximum * (height - 36)
            coordinates.extend((x, y))
        self.create_line(*coordinates, fill="#456990", width=2)
        self.create_text(
            width - 12,
            8,
            text=f"{self.values[-1]:.1f}",
            anchor="ne",
            fill="#456990",
            font=("Microsoft YaHei UI", 8),
        )


def static_action_trajectories(matrix: dict) -> list[dict]:
    result = []
    for action in matrix["actions"]:
        if action["id"] == "K1":
            duration = 1.6
        else:
            duration = 1.4
        result.append(
            {
                "id": f"STATIC_{action['id']}",
                "label": f"{action['id']} {action['label']}",
                "duration_s": duration,
                "return_to_neutral": False,
                "waypoints": [
                    {"time_s": 0.0, "q_target_deg": {}},
                    {
                        "time_s": duration * 0.75,
                        "q_target_deg": action.get("q_target_deg", {}),
                    },
                    {
                        "time_s": duration,
                        "q_target_deg": action.get("q_target_deg", {}),
                    },
                ],
                "active_cables": list(action["allowed_cables"]),
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config, matrix, library = load_inputs()
    library = copy.deepcopy(library)
    library["trajectories"].extend(static_action_trajectories(matrix))
    runner = TrajectoryRunner(
        ROOT / "arm_cable_design_v6.xml",
        config,
        matrix,
        library,
    )
    if args.check:
        print(
            f"V6 GUI ready: {runner.model.nq} joints, "
            f"{runner.model.ntendon} cables"
        )
        return

    current = {"trajectory": library["trajectories"][0], "time": 0.0}
    playing = {"value": False}
    root = tk.Tk()
    root.title("V6 连续轨迹仿生右臂")
    root.geometry("1080x640")
    style = ttk.Style(root)
    style.configure("TLabel", font=("Microsoft YaHei UI", 9))
    style.configure("TButton", font=("Microsoft YaHei UI", 9))

    left = ttk.Frame(root, padding=8)
    left.pack(side="left", fill="y")
    right = ttk.Frame(root, padding=8)
    right.pack(side="right", fill="both", expand=True)
    selection = tk.StringVar(
        value=f"{current['trajectory']['id']} | {current['trajectory']['label']}"
    )
    speed = tk.DoubleVar(value=1.0)
    status = tk.StringVar(value="就绪")

    ttk.Label(
        left,
        text="轨迹与动作",
        font=("Microsoft YaHei UI", 11, "bold"),
    ).pack(anchor="w", pady=(0, 6))
    choices = [
        f"{item['id']} | {item['label']}" for item in library["trajectories"]
    ]
    combo = ttk.Combobox(
        left,
        values=choices,
        textvariable=selection,
        state="readonly",
        width=32,
    )
    combo.pack(anchor="w", fill="x")
    ttk.Label(left, text="播放速度").pack(anchor="w", pady=(10, 2))
    ttk.Scale(
        left,
        from_=0.25,
        to=1.0,
        variable=speed,
        orient="horizontal",
        length=260,
    ).pack(anchor="w")
    speed_label = ttk.Label(left, text="1.00×")
    speed_label.pack(anchor="w")
    buttons = ttk.Frame(left)
    buttons.pack(anchor="w", pady=10)
    play_button = ttk.Button(buttons, text="播放")
    play_button.grid(row=0, column=0, padx=(0, 6))
    ttk.Button(buttons, text="暂停", command=lambda: playing.update(value=False)).grid(
        row=0, column=1, padx=(0, 6)
    )
    ttk.Button(buttons, text="复位", command=lambda: reset()).grid(
        row=0, column=2
    )
    ttk.Label(left, textvariable=status, wraplength=270).pack(
        anchor="w", pady=(8, 0)
    )
    path_canvas = PathCanvas(right, width=520, height=390)
    path_canvas.pack(fill="both", expand=True)
    error_canvas = ErrorCanvas(right, width=520, height=150)
    error_canvas.pack(fill="x", pady=(8, 0))

    def selected_trajectory() -> dict:
        selected = selection.get()
        for item in library["trajectories"]:
            if selected.startswith(item["id"] + " |"):
                return item
        return library["trajectories"][0]

    def reset() -> None:
        playing["value"] = False
        current["time"] = 0.0
        current["trajectory"] = selected_trajectory()
        runner.reset()
        path_canvas.reset()
        error_canvas.reset()
        status.set("已复位")

    def toggle_play() -> None:
        if selection.get() != (
            f"{current['trajectory']['id']} | {current['trajectory']['label']}"
        ):
            reset()
        if current["time"] >= float(current["trajectory"]["duration_s"]):
            reset()
        playing["value"] = not playing["value"]
        play_button.configure(text="暂停" if playing["value"] else "播放")

    play_button.configure(command=toggle_play)
    combo.bind("<<ComboboxSelected>>", lambda _event: reset())
    running = {"value": True}

    def on_close() -> None:
        running["value"] = False
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    with mujoco.viewer.launch_passive(
        runner.model,
        runner.data,
        show_right_ui=False,
    ) as viewer:
        last_wall = time.perf_counter()
        while running["value"] and viewer.is_running():
            root.update()
            speed_value = float(speed.get())
            speed_label.configure(text=f"{speed_value:.2f}×")
            if playing["value"]:
                now = time.perf_counter()
                wall_dt = min(0.05, max(0.0, now - last_wall))
                last_wall = now
                current["time"] = min(
                    float(current["trajectory"]["duration_s"]),
                    current["time"] + wall_dt * speed_value,
                )
                q_target, wrist_target, wrist_actual = runner.step_toward_time(
                    current["trajectory"],
                    current["time"],
                )
                error = float(np.linalg.norm(wrist_target - wrist_actual))
                path_canvas.add(wrist_target, wrist_actual)
                error_canvas.add(error)
                status.set(
                    f"时间 {current['time']:.2f}/"
                    f"{current['trajectory']['duration_s']:.2f} s，"
                    f"误差 {error:.1f} mm"
                )
                if current["time"] >= float(
                    current["trajectory"]["duration_s"]
                ):
                    playing["value"] = False
                    play_button.configure(text="播放")
            else:
                last_wall = time.perf_counter()
            viewer.sync()
            time.sleep(0.005)
        root.destroy()


if __name__ == "__main__":
    main()
