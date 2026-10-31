"""Export the V7 gate decision and diagnostic parameter summary."""

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


ROOT = Path(__file__).resolve().parent
EXPORT = ROOT.parent / "exports" / "cable_design_v7"
DESIGN = ROOT.parent.parent / "design"


def main() -> None:
    elbow = json.loads(
        (EXPORT / "elbow_report_v7.json").read_text(encoding="utf-8")
    )
    shoulder = json.loads(
        (EXPORT / "shoulder_report_v7.json").read_text(encoding="utf-8")
    )
    integration = json.loads(
        (EXPORT / "integration_report_v7.json").read_text(encoding="utf-8")
    )
    elbow_baseline = next(
        item
        for item in elbow["variants"]
        if float(item["radius_mm"]) == 18.0
    )
    shoulder_selected = next(
        item
        for item in shoulder["nominal"]
        if int(item["roller_count"]) == int(
            shoulder["selected_roller_count"]
        )
    )
    lines = [
        "# V7 subsystem gate decisions",
        "",
        "## Elbow sector spool",
        "",
        f"- Gate: `{'PASS' if elbow['pass'] else 'FAIL'}`",
        f"- Baseline radius: `18 mm`",
        f"- Minimum effective moment arm: "
        f"`{elbow_baseline['minimum_arm_mm']:.3f} mm`",
        f"- Maximum active-stroke contraction: "
        f"`{elbow_baseline['contraction_percent']['max']:.3f}%`",
        f"- Opposite moments: `{elbow_baseline['opposite']}`",
        f"- Real guide wrap: `{elbow_baseline['wrap_ok']}`",
        "",
        "## Shoulder R400 roller",
        "",
        f"- Gate: `{'PASS' if shoulder['pass'] else 'FAIL'}`",
        f"- Selected roller count: `{shoulder['selected_roller_count']}`",
        f"- Maximum roller detachment: "
        f"`{shoulder_selected['max_roller_detachment_mm']:.3f} mm`",
        f"- Minimum visual clearance: "
        f"`{shoulder_selected['min_visual_clearance_mm']:.3f} mm`",
        f"- Maximum roller contact force: "
        f"`{shoulder_selected['max_contact_force_N']:.3f} N`",
        f"- Preload sensitivity: `{shoulder['preload_sensitivity_pass']}`",
        "",
        "## Integration",
        "",
        f"- Status: `{integration['status']}`",
        f"- Gate pass: `{integration['gate_pass']}`",
        f"- Frozen hardpoints: `{integration['frozen_hardpoints']}`",
        "",
    ]
    if integration["failure_reasons"]:
        lines.append("## Blocking reasons")
        lines.append("")
        lines.extend(
            f"- `{reason}`" for reason in integration["failure_reasons"]
        )
    else:
        lines.append(
            "Both local mechanisms passed. K1/K2 integration may proceed."
        )
    lines.extend(
        [
            "",
            "No V7 SolidWorks hardpoint file was written because the "
            "integration gate is closed.",
        ]
    )
    DESIGN.mkdir(parents=True, exist_ok=True)
    path = DESIGN / "V7_SUBSYSTEM_GATE.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
