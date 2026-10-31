"""Apply the V7 integration gate. No assembly is built before both stages pass."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPORT = ROOT.parent / "exports" / "cable_design_v7"


def main() -> None:
    elbow = json.loads(
        (EXPORT / "elbow_report_v7.json").read_text(encoding="utf-8")
    )
    shoulder = json.loads(
        (EXPORT / "shoulder_report_v7.json").read_text(encoding="utf-8")
    )
    elbow_pass = bool(elbow["pass"])
    shoulder_pass = bool(shoulder["pass"])
    gate_pass = elbow_pass and shoulder_pass
    reasons = []
    if not elbow_pass:
        reasons.extend(
            f"elbow:{reason}" for reason in elbow["failure_reasons"]
        )
    if not shoulder_pass:
        reasons.extend(
            f"shoulder:{reason}" for reason in shoulder["failure_reasons"]
        )
    report = {
        "stage": "integration",
        "elbow_gate_pass": elbow_pass,
        "shoulder_gate_pass": shoulder_pass,
        "gate_pass": gate_pass,
        "status": "ready_for_integration" if gate_pass else "blocked",
        "full_pass": False,
        "frozen_hardpoints": False,
        "integrated_model_written": False,
        "k1_target_deg": 116.0,
        "k2_target_deg": 96.0,
        "failure_reasons": reasons,
    }
    if not gate_pass:
        report["blocking_statement"] = (
            "V7 integration was not started because at least one isolated "
            "mechanism gate failed. The acceptance limits were not relaxed."
        )
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "integration_report_v7.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
