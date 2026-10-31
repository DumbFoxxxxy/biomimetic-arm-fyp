"""Generate V5 design, anchor and action documents from current results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from v2_common import PROJECT_ROOT, load_json  # noqa: E402


CONFIG = ROOT / "cable_design_v5.json"
MATRIX = ROOT / "action_cable_matrix_v5.json"
EXPORT = ROOT.parent / "exports" / "cable_design_v5"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    parser.add_argument("--evaluation", type=Path, default=EXPORT / "evaluation_v5.json")
    parser.add_argument("--check", type=Path, default=EXPORT / "check_report_v5.json")
    parser.add_argument("--sweep", type=Path, default=EXPORT / "channel_sweep_v5.json")
    args = parser.parse_args()
    config = load_json(args.config)
    matrix = load_json(args.matrix)
    evaluation = load_json(args.evaluation)
    check = load_json(args.check)
    sweep = load_json(args.sweep) if args.sweep.exists() else {"variants": []}

    anchors = [
        "# V5 cable anchors",
        "",
        f"Status: **{'FROZEN' if check['frozen_hardpoints'] else 'NOT FROZEN - AUDIT FAILED'}**",
        "",
        f"- Channels: `{check['channel_count']}`",
        f"- Joints: `{check['joint_count']}`",
        f"- Contact pad offset: `{config.get('contact_pad_offset_mm', 0.0)} mm`",
        "- Zero-pose offsets: "
        + ", ".join(
            f"`{name}={value:.3f} deg`"
            for name, value in config["zero_pose_offsets_deg"].items()
        ),
        "",
        "| Cable | Element | Body | World (x,y,z) mm | Radius mm | Axis (design) |",
        "| --- | --- | --- | --- | ---: | --- |",
    ]
    for route in config["routes"]:
        for index, element in enumerate(route["path"]):
            world = element.get("world_mm", [0.0, 0.0, 0.0])
            axis = element.get("axis_design", [0.0, 0.0, 0.0])
            anchors.append(
                f"| {route['id']} | {index}:{element['kind']} | "
                f"{element['body']} | ({world[0]:.3f}, {world[1]:.3f}, "
                f"{world[2]:.3f}) | {element.get('radius_mm', '')} | "
                f"({axis[0]:.4f}, {axis[1]:.4f}, {axis[2]:.4f}) |"
            )
    (PROJECT_ROOT / "design" / "CABLE_ANCHORS_V5.md").write_text(
        "\n".join(anchors) + "\n", encoding="utf-8"
    )

    lines = [
        "# V5 action-to-cable matrix",
        "",
        f"Status: **{'PASS' if check['full_pass'] else 'FAIL - hard constraints not satisfied'}**",
        "",
        "| Action | Target | Acceptance | Selected active cables | Tension | Span | Contraction | Clearance | Visual | Contact |",
        "| --- | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |",
    ]
    by_id = {action["id"]: action for action in evaluation["actions"]}
    for action in matrix["actions"]:
        result = by_id[action["id"]]
        lines.append(
            f"| {action['id']} | {action.get('target_deg', action.get('target_mm', ''))} | "
            f"{action.get('acceptance_deg', action.get('acceptance_mm', ''))} | "
            f"{', '.join(result['selected_cables'])} | "
            f"{'PASS' if result['pass_tension'] else 'FAIL'} | "
            f"{'PASS' if result['pass_positive_spanning'] else 'FAIL'} | "
            f"{'PASS' if result['pass_contraction'] else 'FAIL'} | "
            f"{'PASS' if result['pass_cable_clearance'] else 'FAIL'} | "
            f"{'PASS' if result['pass_visual_clearance'] else 'FAIL'} | "
            f"{'PASS' if result['pass_contact'] else 'FAIL'} |"
        )
    (PROJECT_ROOT / "design" / "ACTION_CABLE_MATRIX_V5.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    decisions = [
        "# V5 design decisions and audit",
        "",
        "## Scope",
        "",
        "- Seven controlled DOFs and 16-18 independent cable channels.",
        "- K7 target is 126 deg; the previous 150 deg target is retained only as legacy metadata.",
        "- SCAP_INT is formally -5 deg; the full -20 deg target is recorded as infeasible.",
        "- R400 is the manufacturing contact surface; action-specific ellipsoids are sensitivity only.",
        "",
        "## Implemented",
        "",
        "- Shoulder zero offsets are limited to +/-10 deg and contact pad offset to 0-10 mm.",
        "- K6/K7/SLIDE_BACK/SCAP_INT use five-knot shoulder-path corrections.",
        "- Full scapula vertices are checked at 21 poses against the thorax.",
        "- Tension solving uses exact nonnegative equilibrium, mandatory minimum tension, real guide wrap and hard positive-spanning checks.",
        "- A 16/17/18-channel sweep was generated without overwriting V2/V3/V4.",
        "",
        "## Current result",
        "",
        f"- Structural pass: `{check['structural_pass']}`",
        f"- Full pass: `{check['full_pass']}`",
        f"- Frozen hardpoints: `{check['frozen_hardpoints']}`",
        "",
        "### Channel sweep",
        "",
        "| Variant | Channels | Full pass |",
        "| --- | ---: | --- |",
    ]
    for item in sweep.get("variants", []):
        decisions.append(
            f"| {item['name']} | {item['channel_count']} | {item['full_pass']} |"
        )
    decisions += [
        "",
        "## Action status",
        "",
        "Passing the complete current V5 audit: "
        + ", ".join(
            action["id"]
            for action in evaluation["actions"]
            if action["pass_minimum_functional"]
        ),
        "",
        "Blocking actions: "
        + ", ".join(check["failure_actions"].keys()),
        "",
        "K7 achieved: "
        f"`{by_id['K7']['final']['achievement']['achieved']:.3f} deg`; "
        "SCAP_INT achieved: "
        f"`{by_id['SCAP_INT']['final']['achievement']['achieved']:.3f} deg`.",
        "",
        "## Blocking findings",
        "",
        "- The current E01/E02 elbow topology cannot produce opposing moments over the full 116 deg stroke.",
        "- The preferred and global cable sets still fail hard positive spanning in several actions.",
        "- Adding only G08 or S07 does not remove the K1 failure and does not make the global design pass.",
        "- Therefore this V5 package must not be used as a frozen manufacturing drawing.",
        "",
    ]
    (PROJECT_ROOT / "design" / "V5_DESIGN_DECISIONS.md").write_text(
        "\n".join(decisions), encoding="utf-8"
    )
    print("wrote CABLE_ANCHORS_V5.md, ACTION_CABLE_MATRIX_V5.md and V5_DESIGN_DECISIONS.md")


if __name__ == "__main__":
    main()
