"""Audit V7 elbow sector-spool geometry at R15, R18 and R20."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import mujoco


ROOT = Path(__file__).resolve().parent
V2_ROOT = ROOT.parent / "cable_design_v2"
EXPORT = ROOT.parent / "exports" / "cable_design_v7"
sys.path[:0] = [str(ROOT), str(V2_ROOT)]

from elbow_spool_v7 import (  # noqa: E402
    DEFAULT_PARAMETERS,
    build_elbow_model,
    evaluate_elbow,
)


def main() -> None:
    parameters = json.loads(
        (ROOT / "elbow_parameters_v7.json").read_text(encoding="utf-8")
    )
    variants = []
    for radius in (15.0, 18.0, 20.0):
        trial = {**parameters, **DEFAULT_PARAMETERS, "radius_mm": radius}
        tree = build_elbow_model(trial)
        xml = ET.tostring(tree.getroot(), encoding="unicode")
        model = mujoco.MjModel.from_xml_string(xml)
        result = evaluate_elbow(model, radius)
        result["structural_pass"] = bool(
            model.nq == 1
            and model.ntendon == 2
            and model.nu == 2
        )
        result["joint_count"] = int(model.nq)
        result["tendon_count"] = int(model.ntendon)
        result["actuator_count"] = int(model.nu)
        variants.append(result)
        tree.write(
            EXPORT / f"v7_elbow_r{radius:g}.xml",
            encoding="utf-8",
            xml_declaration=True,
        )
    baseline = next(item for item in variants if item["radius_mm"] == 18.0)
    report = {
        "stage": "elbow",
        "baseline_radius_mm": 18.0,
        "baseline_pass": bool(baseline["pass"]),
        "pass": bool(baseline["pass"]),
        "variants": variants,
        "failure_reasons": [],
    }
    if not baseline["opposite"]:
        report["failure_reasons"].append("E01_E02_moment_sign")
    if not baseline["wrap_ok"]:
        report["failure_reasons"].append("guide_wrap")
    if baseline["minimum_arm_mm"] < 8.0:
        report["failure_reasons"].append("effective_moment_arm")
    if baseline["contraction_percent"]["max"] > 18.0:
        report["failure_reasons"].append("cable_contraction")
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "elbow_report_v7.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "baseline_pass": report["baseline_pass"],
                "minimum_arm_mm": baseline["minimum_arm_mm"],
                "contraction_percent": baseline["contraction_percent"],
                "failure_reasons": report["failure_reasons"],
                "variants": {
                    item["radius_mm"]: item["pass"] for item in variants
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
