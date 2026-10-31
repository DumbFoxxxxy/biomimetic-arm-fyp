"""Check whether the SolidWorks audit contains a testable arm motion chain."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


EXPECTED_COMPONENTS = {
    "thorax": ("thorax", "胸廓", "躯干", "torso"),
    "clavicle": ("clavicle", "锁骨"),
    "scapula": ("scapula", "肩胛"),
    "humerus": ("humerus", "肱骨"),
    "radius": ("radius", "桡骨"),
    "ulna": ("ulna", "尺骨"),
}

REQUIRED_JOINT_TOKENS = {
    "sternoclavicular": ("SC_", "STERNOCLAVICULAR"),
    "acromioclavicular_or_scapulothoracic": (
        "AC_",
        "ST_",
        "SCAPULOTHORACIC",
    ),
    "glenohumeral": ("GH_", "GLENOHUMERAL"),
    "elbow": ("ELBOW_", "HUMERORADIAL"),
    "proximal_radioulnar": ("PRUJ_", "PROXIMAL_RADIOULNAR"),
    "distal_radioulnar": ("DRUJ_", "DISTAL_RADIOULNAR"),
}


def read_pipe_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-16").splitlines()
    if not lines:
        return []
    headers = lines[0].split("|")
    rows = []
    for line in lines[1:]:
        values = line.split("|")
        rows.append(dict(zip(headers, values)))
    return rows


def read_audit(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    result = {}
    for line in path.read_text(encoding="utf-16").splitlines():
        parts = line.split("|", 1)
        if len(parts) == 2:
            result[parts[0]] = parts[1]
        elif parts:
            result[parts[0]] = ""
    return result


def normalized(value: str) -> str:
    return "".join(character for character in value.lower() if character.isalnum())


def matching_components(
    components: list[dict[str, str]], aliases: tuple[str, ...]
) -> list[dict[str, str]]:
    normalized_aliases = tuple(normalized(alias) for alias in aliases)
    matches = []
    for component in components:
        text = normalized(
            component.get("name", "") + component.get("path", "")
        )
        if any(alias in text for alias in normalized_aliases):
            matches.append(component)
    return matches


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--validation-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
    )
    args = parser.parse_args()

    validation_dir = args.validation_dir.resolve()
    audit = read_audit(validation_dir / "assembly_audit.txt")
    components = read_pipe_table(validation_dir / "components.txt")
    features = read_pipe_table(validation_dir / "features.txt")
    mates = read_pipe_table(validation_dir / "mates.txt")

    findings = []
    component_status = {}
    for logical_name, aliases in EXPECTED_COMPONENTS.items():
        matches = matching_components(components, aliases)
        active = [
            component
            for component in matches
            if component.get("suppressed", "").lower() != "true"
        ]
        component_status[logical_name] = {
            "count": len(matches),
            "active_count": len(active),
            "names": [component.get("name", "") for component in matches],
        }
        if not active:
            findings.append(
                {
                    "severity": "error",
                    "code": f"missing_{logical_name}",
                    "message": f"No active {logical_name} component was found.",
                }
            )

    non_reference_types = {
        "Reference",
        "RefPlane",
        "OriginProfileFeature",
        "DocsFolder",
        "FavoriteFolder",
        "HistoryFolder",
        "SelectionSetFolder",
        "SensorFolder",
    }
    feature_text = " ".join(
        row.get("name", "")
        for row in features
        if row.get("type", "") not in non_reference_types
    )
    normalized_feature_text = normalized(feature_text)
    joint_status = {}
    for joint_name, tokens in REQUIRED_JOINT_TOKENS.items():
        found_tokens = [
            token
            for token in tokens
            if normalized(token) in normalized_feature_text
        ]
        joint_status[joint_name] = found_tokens
        if not found_tokens:
            findings.append(
                {
                    "severity": "warning",
                    "code": f"joint_names_not_found_{joint_name}",
                    "message": (
                        f"No feature name identifies the {joint_name} joint. "
                        "Its axis cannot yet be mapped reliably."
                    ),
                }
            )

    mate_feature_count = sum(
        1
        for row in features
        if row.get("type", "").startswith("Mate")
    )
    reference_axis_count = sum(
        1 for row in features if row.get("type", "") == "RefAxis"
    )
    coordinate_system_count = sum(
        1 for row in features if row.get("type", "") == "CoordSys"
    )

    if not components:
        findings.append(
            {
                "severity": "error",
                "code": "audit_missing",
                "message": (
                    "Run InspectBiomimeticAssembly.swb in SolidWorks before "
                    "running this validator."
                ),
            }
        )
    elif mate_feature_count == 0:
        findings.append(
            {
                "severity": "warning",
                "code": "mates_not_reported",
                "message": (
                    "No assembly mate features were reported by the audit. "
                    "The model may be positioned but not kinematically "
                    "constrained."
                ),
            }
        )

    if reference_axis_count == 0:
        findings.append(
            {
                "severity": "warning",
                "code": "no_reference_axes",
                "message": (
                    "No reference axes were found. Define named axes for the "
                    "shoulder, elbow and radioulnar joints."
                ),
            }
        )

    report = {
        "source": str(validation_dir / "assembly_audit.txt"),
        "audit": audit,
        "component_status": component_status,
        "joint_name_status": joint_status,
        "counts": {
            "components": len(components),
            "features": len(features),
            "mates": len(mates),
            "mate_features": mate_feature_count,
            "reference_axes": reference_axis_count,
            "coordinate_systems": coordinate_system_count,
        },
        "findings": findings,
        "ready_for_mujoco_kinematic_export": not any(
            finding["severity"] == "error" for finding in findings
        ),
    }

    output_path = args.output or validation_dir / "motion_contract_report.json"
    output_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"Components: {len(components)}")
    print(f"Mate features: {mate_feature_count}")
    print(f"Reference axes: {reference_axis_count}")
    print(f"Coordinate systems: {coordinate_system_count}")
    for logical_name, status in component_status.items():
        state = "OK" if status["active_count"] else "MISSING"
        print(
            f"  {logical_name:12s} {state:7s} "
            f"{', '.join(status['names']) or '-'}"
        )
    if findings:
        print("Findings:")
        for finding in findings:
            print(
                f"  [{finding['severity'].upper()}] "
                f"{finding['code']}: {finding['message']}"
            )
    else:
        print("Findings: none")
    print(f"Report: {output_path}")

    return 0 if report["ready_for_mujoco_kinematic_export"] else 1


if __name__ == "__main__":
    sys.exit(main())
