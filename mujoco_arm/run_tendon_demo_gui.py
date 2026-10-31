"""MuJoCo viewer plus a simple cable-pull control panel."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import json
from project_runtime import WorkspacePath as Path
import time
import tkinter as tk
from tkinter import ttk

import mujoco
import mujoco.viewer
import numpy as np


ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "arm_tendon_demo.xml"
CONFIG_PATH = ROOT / "config" / "tendon_demo_16.json"

GROUPS = (
    ("肘关节", ("E01", "E02")),
    ("盂肱关节屈伸", ("G01", "G02")),
    ("盂肱关节外展内收", ("G03", "G04")),
    ("盂肱关节内外旋", ("G05", "G06")),
    ("肩胛骨前后滑动", ("S03", "S04")),
    ("肩胛骨上下回旋", ("S01", "S02")),
    ("肩胛骨内外旋", ("S05", "S06")),
    ("胸大肌", ("P01", "P02")),
)

PRESETS = {
    "零位复位": {},
    "K1 屈肘": {"E01": 40.0},
    "K2 肩前屈": {"G01": 30.0},
    "K3 肩外展": {"G03": 40.0},
    "K4 肩内旋": {"G05": 30.0},
    "K5 肩外旋": {"G06": 30.0},
    "K6 肩胛上回旋": {"S01": 40.0},
    "K7 组合抬举": {"G01": 30.0, "S01": 40.0},
    "肩胛下压": {"S02": 40.0},
    "肩胛前移": {"S03": 30.0},
    "肩胛后移": {"S04": 30.0},
    "肩胛外旋": {"S06": 40.0},
    "肩胛内旋": {"S05": 30.0},
    "胸大肌抱合": {"P01": 30.0, "P02": 40.0},
}

CABLE_LABELS = {
    "E01": "屈肘",
    "E02": "伸肘",
    "G01": "肩前屈",
    "G02": "肩后伸",
    "G03": "肩外展",
    "G04": "肩内收",
    "G05": "肩内旋",
    "G06": "肩外旋",
    "S01": "肩胛上回旋",
    "S02": "肩胛下回旋",
    "S03": "肩胛前移",
    "S04": "肩胛后移",
    "S05": "肩胛内旋",
    "S06": "肩胛外旋",
    "P01": "胸大肌锁骨部",
    "P02": "胸大肌胸肋部",
}


def main() -> None:
    model = _load_project_model(str(MODEL_PATH))
    data = mujoco.MjData(model)
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    neutral = {
        route["id"]: float(route["neutral_length_m"])
        for route in config["routes"]
    }
    actuators = {}
    for index in range(model.nu):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, index) or ""
        if name.startswith("pull_"):
            actuators[name[5:]] = index

    data.eq_active[:] = 0
    data.ctrl[:] = 0.0
    for cable_id, actuator in actuators.items():
        data.ctrl[actuator] = neutral[cable_id]
    mujoco.mj_forward(model, data)

    root = tk.Tk()
    root.title("十四钢丝仿生手臂演示")
    root.resizable(False, False)
    style = ttk.Style(root)
    style.configure("TLabel", font=("Microsoft YaHei UI", 9))
    style.configure("TLabelframe.Label", font=("Microsoft YaHei UI", 9, "bold"))
    style.configure("TButton", font=("Microsoft YaHei UI", 9))
    values = {cable_id: tk.DoubleVar(value=0.0) for cable_id in actuators}

    left = ttk.Frame(root, padding=8)
    left.grid(row=0, column=0, sticky="nsew")
    right = ttk.Frame(root, padding=8)
    right.grid(row=0, column=1, sticky="nsew")

    ttk.Label(
        left,
        text="钢丝拉入量（mm）",
        font=("Microsoft YaHei UI", 10, "bold"),
    ).grid(
        row=0, column=0, columnspan=2, sticky="w", pady=(0, 6)
    )
    row = 1
    for group_name, cables in GROUPS:
        group = ttk.LabelFrame(left, text=group_name, padding=5)
        group.grid(row=row, column=0, columnspan=2, sticky="ew", pady=3)
        for index, cable_id in enumerate(cables):
            ttk.Label(
                group,
                text=f"{cable_id} {CABLE_LABELS[cable_id]}",
                width=15,
            ).grid(
                row=index, column=0, sticky="w"
            )
            ttk.Scale(
                group,
                from_=0.0,
                to=80.0,
                variable=values[cable_id],
                orient="horizontal",
                length=230,
            ).grid(row=index, column=1, sticky="ew", padx=4)
            ttk.Label(group, textvariable=values[cable_id], width=6).grid(
                row=index, column=2, sticky="e"
            )
        row += 1

    def reset() -> None:
        data.qpos[:] = 0.0
        data.qvel[:] = 0.0
        for variable in values.values():
            variable.set(0.0)
        status.set("当前状态：零位姿态")

    def apply_preset(pulls: dict[str, float]) -> None:
        reset()
        for cable_id, amount in pulls.items():
            values[cable_id].set(amount)
        if pulls:
            text = "、".join(
                f"{cable_id} {amount:.0f} mm"
                for cable_id, amount in pulls.items()
            )
            status.set(f"当前预设：{text}")
        else:
            status.set("当前状态：零位姿态")

    ttk.Label(
        right,
        text="姿态预设",
        font=("Microsoft YaHei UI", 10, "bold"),
    ).grid(row=0, column=0, sticky="w")
    for index, (name, pulls) in enumerate(PRESETS.items(), start=1):
        ttk.Button(
            right,
            text=name,
            width=18,
            command=lambda p=pulls: apply_preset(p),
        ).grid(row=index, column=0, sticky="ew", pady=2)

    info_index = len(PRESETS) + 1
    ttk.Separator(right).grid(
        row=info_index, column=0, sticky="ew", pady=8
    )
    ttk.Label(
        right,
        text=(
            "0 mm：钢丝松弛\n"
            "80 mm：最大拉入量\n\n"
            "拖动滑块会缩短对应钢丝的\n"
            "目标长度，钢丝拉力使关节运动。\n\n"
            "肩胛骨由 SC 两轴关节、AC\n"
            "铰链和胸廓接触共同约束。"
        ),
        justify="left",
    ).grid(row=info_index + 1, column=0, sticky="w")
    ttk.Button(right, text="恢复零位", command=reset).grid(
        row=info_index + 2, column=0, sticky="ew", pady=(10, 0)
    )
    status = tk.StringVar(value="当前状态：零位姿态")
    ttk.Label(
        right,
        textvariable=status,
        foreground="#444444",
        wraplength=220,
    ).grid(row=info_index + 3, column=0, sticky="w", pady=(8, 0))

    running = True

    def on_close() -> None:
        nonlocal running
        running = False
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)

    with mujoco.viewer.launch_passive(
        model,
        data,
        show_right_ui=False,
    ) as viewer:
        while running and viewer.is_running():
            root.update()
            for cable_id, actuator in actuators.items():
                data.ctrl[actuator] = (
                    neutral[cable_id] - float(values[cable_id].get()) * 0.001
                )
            mujoco.mj_step(model, data)
            viewer.sync()
            time.sleep(model.opt.timestep)
        root.destroy()


if __name__ == "__main__":
    main()
