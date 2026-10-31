"""Audit V7 single- and dual-roller scapulothoracic models."""

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
import sys
import xml.etree.ElementTree as ET

import mujoco


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v7"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from shoulder_roller_v7 import (  # noqa: E402
    build_shoulder_model,
    evaluate_shoulder_model,
    load_config,
)


MODEL = ROOT / "v7_shoulder_standalone.xml"


def build_model(count: int) -> mujoco.MjModel:
    tree = build_shoulder_model(count)
    xml = ET.tostring(tree.getroot(), encoding="unicode")
    return _load_project_xml(xml)


def main() -> None:
    config = load_config()
    nominal = []
    models = {}
    for count in (1, 2):
        model = build_model(count)
        models[count] = model
        nominal.append(
            evaluate_shoulder_model(model, config, count, preload_n=5.0)
        )
    passing = [item for item in nominal if item["pass"]]
    selected = passing[0] if passing else min(
        nominal,
        key=lambda item: (
            abs(item["max_roller_detachment_mm"] - 0.35)
            + abs(item["min_visual_clearance_mm"] + 0.1)
        ),
    )
    selected_count = int(selected["roller_count"])
    sensitivity = []
    for preload in (3.0, 5.0, 10.0):
        sensitivity.append(
            evaluate_shoulder_model(
                models[selected_count],
                config,
                selected_count,
                preload_n=preload,
            )
        )
    preload_stable = all(item["pass"] for item in sensitivity)
    report = {
        "stage": "shoulder",
        "single_roller_pass": nominal[0]["pass"],
        "dual_roller_pass": nominal[1]["pass"],
        "selected_roller_count": selected_count,
        "nominal_pass": bool(selected["pass"]),
        "preload_sensitivity_pass": bool(preload_stable),
        "pass": bool(selected["pass"] and preload_stable),
        "nominal": nominal,
        "preload_sensitivity": sensitivity,
        "failure_reasons": [],
    }
    if not selected["pass"]:
        if selected["min_visual_clearance_mm"] < -0.5:
            report["failure_reasons"].append("full_mesh_penetration")
        if selected["min_roller_penetration_mm"] < -0.5:
            report["failure_reasons"].append("roller_penetration")
        if selected["max_roller_detachment_mm"] > 1.0:
            report["failure_reasons"].append("roller_detachment")
        if not 3.0 <= selected["max_contact_force_N"] <= 10.0:
            report["failure_reasons"].append("roller_contact_force")
    if not preload_stable:
        report["failure_reasons"].append("preload_sensitivity")
    tree = build_shoulder_model(selected_count)
    tree.write(MODEL, encoding="utf-8", xml_declaration=True)
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "shoulder_report_v7.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "single_roller_pass": report["single_roller_pass"],
                "dual_roller_pass": report["dual_roller_pass"],
                "selected_roller_count": selected_count,
                "nominal_pass": report["nominal_pass"],
                "preload_sensitivity_pass": report["preload_sensitivity_pass"],
                "pass": report["pass"],
                "min_visual_clearance_mm": selected[
                    "min_visual_clearance_mm"
                ],
                "max_roller_detachment_mm": selected[
                    "max_roller_detachment_mm"
                ],
                "max_contact_force_N": selected["max_contact_force_N"],
                "failure_reasons": report["failure_reasons"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
