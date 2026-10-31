"""Search muscle subsets and PSM parameters with Optuna NSGA-II.

The script falls back to a deterministic random search when Optuna is not
installed, so the contract and model build can still be tested.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import random
import tempfile

import mujoco
import numpy as np

from build_arm_model import ArmBuilder


ROOT = Path(__file__).resolve().parent
ARM_CONFIG = ROOT / "config" / "arm_model.json"
MUSCLE_CONFIG = ROOT / "config" / "muscle_manifest.json"
RESULTS_PATH = ROOT / "exports" / "pareto_designs.json"

SHAPES = ("straight_diamond", "fan", "y_branch", "wrap")
NODES_BY_SHAPE = {
    "straight_diamond": (7, 3),
    "fan": (7, 4),
    "y_branch": (7, 3),
    "wrap": (6, 4),
}


def sampled_parameters(suggest) -> dict:
    return {
        "x_mm": suggest("x_mm", 8.0, 20.0),
        "y_mm": suggest("y_mm", 5.0, 16.0),
        "z_mm": suggest("z_mm", 20.0, 40.0),
        "p_mm": suggest("p_mm", 1.5, 3.0),
    }


def build_candidate(base_manifest: dict, suggest) -> dict:
    manifest = copy.deepcopy(base_manifest)
    global_parameters = sampled_parameters(suggest)
    for muscle in manifest["muscles"]:
        muscle_id = muscle["id"]
        if muscle.get("required", False):
            muscle["enabled"] = True
        else:
            muscle["enabled"] = suggest(
                f"enable_{muscle_id}", [False, True]
            )
        muscle["cables"] = suggest(f"cables_{muscle_id}", [1, 2])
        muscle["shape_class"] = suggest(
            f"shape_{muscle_id}", list(SHAPES)
        )
        longitudinal, transverse = NODES_BY_SHAPE[muscle["shape_class"]]
        muscle["lattice"] = {
            **global_parameters,
            "longitudinal_nodes": longitudinal,
            "transverse_nodes": transverse,
        }
    return manifest


def evaluate_candidate(
    manifest: dict, trial_id: int, arm_config_path: Path
) -> tuple[float, ...]:
    with tempfile.TemporaryDirectory(prefix="psm_opt_") as temp_dir:
        temp = Path(temp_dir)
        muscle_path = temp / "muscle_manifest.json"
        arm_path = temp / "arm_model.json"
        muscle_path.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        arm = json.loads(arm_config_path.read_text(encoding="utf-8"))
        arm["output_xml"] = f"generated/trial_{trial_id:05d}.xml"
        arm_path.write_text(
            json.dumps(arm, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        try:
            builder = ArmBuilder(arm_path, muscle_path)
            output_path = builder.build()
            model = mujoco.MjModel.from_xml_path(str(output_path))
            data = mujoco.MjData(model)
            for _ in range(80):
                mujoco.mj_step(model, data)
            if not np.isfinite(data.qpos).all():
                raise RuntimeError("non-finite state")

            enabled = [
                muscle
                for muscle in manifest["muscles"]
                if muscle.get("enabled", False)
            ]
            actuator_count = sum(int(muscle["cables"]) for muscle in enabled)

            cable_length = 0.0
            lattice_volume = 0.0
            for muscle in enabled:
                origin = builder.body_to_world(
                    muscle["origin_body"],
                    np.asarray(muscle["origin_local"], dtype=float),
                )
                insertion = builder.body_to_world(
                    muscle["insertion_body"],
                    np.asarray(muscle["insertion_local"], dtype=float),
                )
                length = float(np.linalg.norm(insertion - origin))
                cable_length += length * int(muscle["cables"])
                lattice = muscle["lattice"]
                lattice_volume += (
                    length
                    * float(lattice["y_mm"])
                    * 1e-3
                    * float(lattice["z_mm"])
                    * 1e-3
                )
            optional_count = sum(
                1 for muscle in enabled if not muscle.get("required", False)
            )
            return (
                float(actuator_count),
                float(cable_length + 0.02 * lattice_volume),
                float(optional_count),
            )
        except Exception:
            return (1e6, 1e6, 1e6)
        finally:
            trial_path = ROOT / arm["output_xml"]
            trial_path.unlink(missing_ok=True)


def run_optuna(
    base_manifest: dict, n_trials: int, arm_config_path: Path
) -> tuple[list, str]:
    try:
        import optuna
    except ImportError:
        return run_random(base_manifest, n_trials, arm_config_path)

    sampler = optuna.samplers.NSGAIISampler(seed=42)
    study = optuna.create_study(
        study_name="right_arm_psm",
        directions=["minimize", "minimize", "minimize"],
        sampler=sampler,
    )

    def objective(trial):
        trial_index = len(study.trials)
        def suggest(name, *args):
            if len(args) == 1:
                return trial.suggest_categorical(name, list(args[0]))
            return trial.suggest_float(name, float(args[0]), float(args[1]))

        candidate = build_candidate(base_manifest, suggest)
        return evaluate_candidate(candidate, trial_index, arm_config_path)

    study.optimize(objective, n_trials=n_trials)
    pareto_trials = study.best_trials
    records = []
    for trial in pareto_trials:
        records.append(
            {
                "trial": trial.number,
                "values": list(trial.values),
                "params": trial.params,
            }
        )
    return records, "optuna_nsga2"


def run_random(
    base_manifest: dict, n_trials: int, arm_config_path: Path
) -> tuple[list, str]:
    generator = random.Random(42)
    records = []

    def suggest(name: str, *args):
        if len(args) == 1:
            return generator.choice(list(args[0]))
        return generator.uniform(float(args[0]), float(args[1]))

    for trial in range(n_trials):
        candidate = build_candidate(base_manifest, suggest)
        values = evaluate_candidate(candidate, trial, arm_config_path)
        params = {
            muscle["id"]: {
                "enabled": bool(muscle["enabled"]),
                "cables": int(muscle["cables"]),
                "shape_class": muscle["shape_class"],
                "lattice": muscle["lattice"],
            }
            for muscle in candidate["muscles"]
        }
        records.append(
            {
                "trial": trial,
                "values": list(values),
                "params": params,
            }
        )

    records.sort(
        key=lambda record: (
            record["values"][0]
            + 10.0 * record["values"][1]
            + 20.0 * record["values"][2]
        )
    )
    return records[:20], "random_fallback"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--arm-config", type=Path, default=ARM_CONFIG)
    parser.add_argument(
        "--muscle-config", type=Path, default=MUSCLE_CONFIG
    )
    parser.add_argument("--results", type=Path, default=RESULTS_PATH)
    args = parser.parse_args()

    base_manifest = json.loads(
        args.muscle_config.read_text(encoding="utf-8")
    )
    records, method = run_optuna(
        base_manifest, args.trials, args.arm_config
    )
    output = {
        "method": method,
        "objective_names": [
            "actuator_count",
            "cable_length_and_lattice_volume",
            "optional_muscle_count",
        ],
        "records": records,
    }
    args.results.parent.mkdir(parents=True, exist_ok=True)
    args.results.write_text(
        json.dumps(output, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"method: {method}")
    print(f"records: {len(records)}")
    print(f"results: {args.results}")


if __name__ == "__main__":
    main()
