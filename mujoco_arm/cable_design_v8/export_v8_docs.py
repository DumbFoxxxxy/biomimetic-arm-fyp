"""Export the current V8 shoulder-gate decision document."""

from __future__ import annotations

import json

from v8_shoulder import EXPORT_ROOT, PROJECT_ROOT


DESIGN = PROJECT_ROOT / "design"


def main() -> None:
    kinematics = json.loads(
        (EXPORT_ROOT / "kinematics_report_v8.json").read_text(encoding="utf-8")
    )
    contact = json.loads(
        (EXPORT_ROOT / "contact_report_v8.json").read_text(encoding="utf-8")
    )
    integration = json.loads(
        (EXPORT_ROOT / "integration_report_v8.json").read_text(encoding="utf-8")
    )
    worst = kinematics.get("failure_diagnostics", {}).get(
        "worst_target",
        "none",
    )
    feasible = kinematics.get("failure_diagnostics", {}).get(
        "maximum_clearance_feasible_slide_back"
    )
    feasible_text = (
        f"`s={feasible['s_mm']:.1f} mm`, clearance "
        f"`{feasible['clearance_mm']:.3f} mm`"
        if feasible
        else "not found"
    )
    lines = [
        "# V8 five-DOF shoulder chain",
        "",
        "## Purpose",
        "",
        "V8 replaces the single AC hinge with an exact three-DOF AC ball joint",
        "and keeps the SC two-axis universal joint and the 150 mm clavicle.",
        "Shoulder motion is defined by the surface coordinates `s`, `t`, and",
        "`phi`; cable optimization is outside this gate.",
        "",
        "## Model",
        "",
        "- SC center: `(40, 235, -106) mm`",
        "- AC center: `(156.86, 233.24, -11.97) mm`",
        "- SC-AC design length: `150 mm`",
        "- MuJoCo joints: `sc_elevation`, `sc_protraction`, `ac_ball`",
        "- Scalar shoulder DOFs: `5`",
        "- Full scapula and clavicle meshes remain visual-only",
        "",
        "## Gate Results",
        "",
        "| Check | Result |",
        "| --- | --- |",
        f"| Kinematic gate | `{kinematics['pass']}` |",
        f"| Contact gate | `{contact['pass']}` |",
        f"| Integration gate | `{integration['gate_pass']}` |",
        f"| Frozen hardpoints | `{integration['frozen_hardpoints']}` |",
        f"| Maximum IK position residual | `{kinematics['max_position_residual_mm']:.6g} mm` |",
        f"| Maximum IK orientation residual | `{kinematics['max_orientation_residual_deg']:.6g} deg` |",
        f"| Minimum full-mesh clearance | `{kinematics['min_full_mesh_clearance_mm']:.3f} mm` |",
        f"| Minimum contact-pad clearance | `{kinematics['min_contact_pad_clearance_mm']:.3f} mm` |",
        f"| Worst target | `{worst}` |",
        "",
        "The contact gate reports:",
        "",
        "```text",
        ", ".join(contact["failure_reasons"]),
        "```",
        "",
        "The maximum currently clearance-feasible negative slide is",
        f"{feasible_text}.",
        "",
        "## Interpretation",
        "",
        "The five-DOF chain can follow the new surface targets to micrometre",
        "level. The remaining failure is geometric: at the requested",
        "`SLIDE_BACK=-15 mm`, the rigid scapula mesh intersects the current",
        "R400/ellipsoid thorax envelope even after the permitted surface path",
        "adjustment. Per the approved gate, the R8 roller contact stage was not",
        "used to hide this error and no hardpoint table was frozen.",
        "",
        "## Outputs",
        "",
        "- `mujoco_arm/cable_design_v8/v8_shoulder_open_chain.xml`",
        "- `mujoco_arm/cable_design_v8/v8_shoulder_contact.xml`",
        "- `mujoco_arm/exports/cable_design_v8/kinematics_report_v8.json`",
        "- `mujoco_arm/exports/cable_design_v8/contact_report_v8.json`",
        "- `mujoco_arm/exports/cable_design_v8/integration_report_v8.json`",
        "",
        "The V7 R18 elbow regression and the 16-cable demo regression both pass.",
        "V5/V6/V7 source and model files were not edited.",
        "",
    ]
    DESIGN.mkdir(parents=True, exist_ok=True)
    (DESIGN / "V8_SHOULDER_CHAIN.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )
    print(f"wrote {DESIGN / 'V8_SHOULDER_CHAIN.md'}")


if __name__ == "__main__":
    main()

