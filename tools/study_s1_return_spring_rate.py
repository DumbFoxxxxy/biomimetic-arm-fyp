"""Evaluate the sized extension spring in the extended-rail force model."""

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_spring_integration import load_candidate
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from simulate_elbow_series_slider import SLIDE_DAMPING
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference
from study_s1_series_elbow_release_timing import run


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    _, extended, _, spring, sources, owners, revisions, retention = load_candidate(folder)
    size = load_reference(spring["size"])
    budget = load_reference(extended["budget"])
    timing = load_reference(budget["timing"])
    active = load_reference(timing["input_cycles"])
    start = load_reference(active["start"])
    axis = load_reference(load_reference(timing["baseline"])["shoulder_endpoints"])["axes_source"][0]
    schedule = timing["cases"][timing["selected"]]
    friction_file = folder / "extended_anterior_friction_integration.json"
    friction = read(friction_file)
    if friction["model"] != extended["model"]:
        raise ValueError("Spring and friction studies use different models")
    anterior = next(row for row in friction["cases"] if row["mu"] == .3)
    model_file = ROOT / extended["model"]["file"]
    if provenance(model_file) != extended["model"]:
        raise ValueError("Extended model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {key: data.body(key).xpos.copy() * 1000. for key in MOBILE}
    names = active["shoulder_channels"]
    actuator_ids = [model.actuator(name).id for name in names]
    front = names.index("FJ1468")
    elbow = model.actuator("FJ1486").id
    targets = bounded({"thorax_reference:FJ3330": load_component(sources["thorax_reference:FJ3330"])})
    obstacles = {key: load_component({**sources, **revisions}[key]) for key in
                 ("bone:humerus", "elbow:elbow_washer_left", "elbow:elbow_shoulder_pin")}
    header = ["time_s", "GH_delta_deg", "flexion_projection_deg", "elbow_deg", "slide_mm",
              *(f"{model.actuator(i).name}_tension_N" for i in range(model.nu)),
              *(f"qpos_{i}" for i in range(model.nq))]
    cases = []
    rates = spring["design"]["rate_N_per_mm"]
    for label, rate in zip(("low", "nominal", "high"), rates, strict=True):
        for command in (0., anterior["anterior_input_increment_N"]):
            delta = .6 * np.array(active["full_flexion_increment_N"])
            delta[front] += command
            transmission = np.ones(len(names))
            transmission[front] = anterior["transmission"]
            rows = run(model, np.array(start["candidate"]["quaternion_wxyz"]), actuator_ids,
                       np.array(active["preload_N"]), delta, schedule["shoulder_release_start_s"],
                       schedule["elbow_release_start_s"], np.array(axis),
                       shoulder_transmission=transmission, return_preload_N=size["installed_preload_N"],
                       return_stiffness_N_per_m=1000. * rate)
            rows[:, 5 + elbow] *= abs(model.actuator_gear[elbow, 0])
            first = None
            for row in rows:
                data.qpos[:] = row[-model.nq:]
                mujoco.mj_kinematics(model, data)
                frames = transforms(model, data, neutral)
                hits = {key: hit for key, solid in obstacles.items()
                        if (hit := overlaps(solid.transform(frames[owners[key]]), targets))}
                if hits:
                    first = dict(time_s=float(row[0]), elbow_deg=float(row[3]), hits_mm3=hits)
                    break
            trace = folder / f"spring_rate_{label}_front_{int(command)}N.csv"
            write_csv(trace, header, rows)
            case = dict(rate_N_per_mm=rate, modulus_MPa=size["assumed_shear_modulus_MPa"][len(cases) // 2],
                        anterior_command_peak_N=command,
                        anterior_delivered_peak_N=float(rows[:, 5 + actuator_ids[front]].max()),
                        elbow_peak_deg=float(rows[:, 3].max()), elbow_end_deg=float(rows[-1, 3]),
                        slide_peak_mm=float(rows[:, 4].max()), slide_end_mm=float(rows[-1, 4]),
                        remaining_rail_mm=float(extended["rail_stroke_mm"] - rows[:, 4].max()),
                        first_known_rib_interference=first, sample_count=len(rows), trace=provenance(trace))
            cases.append(case)
            print(label, command, case["elbow_peak_deg"], case["slide_peak_mm"], first, flush=True)
    save(folder / "return_spring_rate_dynamics.json", dict(
        generator=provenance(ROOT / "tools/study_s1_return_spring_rate.py"), retention=retention,
        model=extended["model"], spring_size=spring["size"], friction=provenance(friction_file),
        preload_N=size["installed_preload_N"], assumed_damping_Ns_per_m=SLIDE_DAMPING, cases=cases,
        scope="Force-driven GH/elbow model using sized spring rates and assumed damping/preload. "
              "The three known rib pairs are checked at 5 ms. Coil, TPU and other rigid pairs "
              "are not validated on these new trajectories.",
    ))


if __name__ == "__main__":
    main()
