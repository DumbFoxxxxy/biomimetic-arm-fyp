"""Explore two-anchor elbow tendon geometry before adding guide hardware."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parent
EXPORT = ROOT.parent / "exports" / "cable_design_v6"
E = np.array([184.34, -101.35, 7.30], dtype=float)
AXIS = np.array([0.9962, 0.0872, 0.0], dtype=float)
AXIS /= np.linalg.norm(AXIS)
FRACTIONS = np.linspace(0.0, 1.0, 41)


def rotate_x(point: np.ndarray, angle: float) -> np.ndarray:
    cosine = np.cos(angle)
    sine = np.sin(angle)
    y = cosine * point[1] - sine * point[2]
    z = sine * point[1] + cosine * point[2]
    return np.array([point[0], y, z])


def mirror_z(point: np.ndarray) -> np.ndarray:
    result = point.copy()
    result[2] = 2.0 * E[2] - result[2]
    return result


def evaluate(values: np.ndarray) -> dict:
    a = values[0:3]
    b = values[3:6]
    b_local = b - E
    a2 = mirror_z(a)
    b_local2 = mirror_z(b) - E
    samples = []
    for fraction in FRACTIONS:
        angle = np.deg2rad(116.0 * float(fraction))
        b_q = E + rotate_x(b_local, angle)
        b_q2 = E + rotate_x(b_local2, angle)
        length_1 = np.linalg.norm(b_q - a)
        length_2 = np.linalg.norm(b_q2 - a2)
        delta = np.deg2rad(0.05)
        b_plus = E + rotate_x(b_local, angle + delta)
        b_minus = E + rotate_x(b_local, angle - delta)
        b2_plus = E + rotate_x(b_local2, angle + delta)
        b2_minus = E + rotate_x(b_local2, angle - delta)
        d_l_1 = (
            np.linalg.norm(b_plus - a) - np.linalg.norm(b_minus - a)
        ) / (2.0 * delta)
        d_l_2 = (
            np.linalg.norm(b2_plus - a2) - np.linalg.norm(b2_minus - a2)
        ) / (2.0 * delta)
        samples.append(
            {
                "fraction": float(fraction),
                "moment_1_m": float(-d_l_1),
                "moment_2_m": float(-d_l_2),
                "length_1_mm": float(length_1),
                "length_2_mm": float(length_2),
            }
        )
    opposite = all(
        sample["moment_1_m"] * sample["moment_2_m"] < 0.0
        for sample in samples
    )
    arm_ok = all(
        min(abs(sample["moment_1_m"]), abs(sample["moment_2_m"])) >= 0.008
        for sample in samples
    )
    return {
        "values": values.tolist(),
        "opposite": opposite,
        "arm_ok": arm_ok,
        "pass": opposite and arm_ok,
        "samples": samples,
    }


def objective(values: np.ndarray) -> float:
    result = evaluate(values)
    penalty = 0.0
    for sample in result["samples"]:
        product = sample["moment_1_m"] * sample["moment_2_m"]
        if product >= 0.0:
            penalty += 20.0 + min(10.0, abs(product) * 500.0)
        minimum_arm = min(
            abs(sample["moment_1_m"]), abs(sample["moment_2_m"])
        )
        penalty += max(0.0, 0.008 - minimum_arm) * 5000.0
    return penalty


def main() -> None:
    bounds = [
        (160.0, 210.0),
        (0.0, 160.0),
        (-80.0, -15.0),
        (160.0, 210.0),
        (-230.0, -110.0),
        (-80.0, -15.0),
    ]
    result = differential_evolution(
        objective,
        bounds,
        seed=106,
        maxiter=180,
        popsize=18,
        polish=True,
        workers=1,
        updating="immediate",
    )
    evaluation = evaluate(np.asarray(result.x))
    output = {
        "optimizer_success": bool(result.success),
        "score": float(result.fun),
        "result": evaluation,
        "interpretation": (
            "Direct two-anchor geometry is sufficient."
            if evaluation["pass"]
            else "Direct two-anchor geometry cannot remain antagonistic; "
            "use a capstan/differential pulley topology."
        ),
    }
    EXPORT.mkdir(parents=True, exist_ok=True)
    (EXPORT / "elbow_direct_exploration_v6.json").write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "score": output["score"],
                "pass": evaluation["pass"],
                "values": evaluation["values"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
