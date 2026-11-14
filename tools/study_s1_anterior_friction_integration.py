"""Apply the FJ1468 forward-sliding friction sensitivity to the shoulder–elbow model."""

import argparse
from itertools import combinations

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from study_s1_local_force_geometry import load_reference
from study_s1_series_elbow_release_timing import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extended", action="store_true", help="Use the 30 N extended-rail candidate")
    extended = parser.parse_args().extended
    folder = DEST.parent / "integration_review/s1_series_elbow"
    timing_file = folder / "gravity_balanced_release_timing.json"
    timing = read(timing_file)
    active = load_reference(timing["input_cycles"])
    start = load_reference(active["start"])
    baseline = load_reference(active["baseline"])
    shoulder = load_reference(baseline["shoulder_endpoints"])
    geometry = load_reference(baseline["geometry_assembly"])
    model_ref = baseline["model"]
    sources = geometry["sources"]
    extension_ref = None
    if extended:
        extension_file = folder / "elbow_extended_rail_free_drive.json"
        extension_ref = provenance(extension_file)
        extension = read(extension_file)
        budget = load_reference(extension["budget"])
        geometry = load_reference(budget["assembly"])
        sources = {**geometry["sources"], **extension["shafts"]}
        model_ref = extension["model"]
    path_file = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief/FJ1468_revised_channel_motion.json"
    channel = read(path_file)
    wrap = max(row["near_port_shortcut_screen"]["near_port_turn_rad"] for row in channel["cases"]
               if row.get("near_port_shortcut_screen", {}).get("valid"))
    model_file = ROOT / model_ref["file"]
    if provenance(model_file) != model_ref:
        raise ValueError("Source MuJoCo model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    origins = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    keys = list(sources)
    owners = [geometry["owners"][key] for key in keys]
    solids = [load_component(sources[key]) for key in keys]
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                        if owners[i] != owners[j] and (owners[i] in MOBILE or owners[j] in MOBILE)])
    pin, rib = (keys.index(key) for key in ("elbow:elbow_shoulder_pin", "thorax_reference:FJ3330"))
    local_pairs = np.array([[0, 1]])
    channels = [model.actuator(name).id for name in active["shoulder_channels"]]
    anterior = active["shoulder_channels"].index("FJ1468")
    preload = np.asarray(active["preload_N"])
    delta = .6 * np.asarray(active["full_flexion_increment_N"])
    # Use the existing documented anterior pull case amplitude. Reuse the
    # selected release schedule for every trial so only transmission changes.
    peak = timing["cases"]["anterior_pull_2N"]["added_anterior_peak_N"]
    delta[anterior] += peak
    schedule = timing["cases"][timing["selected"]]
    names = ["time_s", "GH_delta_deg", "flexion_projection_deg", "elbow_deg", "slide_mm",
             *(f"{model.actuator(i).name}_delivered_N" for i in range(model.nu)),
             *(f"qpos_{i}" for i in range(model.nq))]
    results = []
    cases = [(mu, peak) for mu in (0., .1, .2, .3)]
    if extended:
        cases.insert(0, (0., 0.))
    for mu, command in cases:
        trial_delta = delta.copy()
        trial_delta[anterior] += command - peak
        transmission = np.ones(len(channels))
        transmission[anterior] = np.exp(-mu * wrap)
        rows = run(model, np.asarray(start["candidate"]["quaternion_wxyz"]), channels, preload, trial_delta,
                   schedule["shoulder_release_start_s"], schedule["elbow_release_start_s"],
                   np.asarray(shoulder["axes_source"][0]), sample_period_s=.005,
                   shoulder_transmission=transmission)
        if extended and command == 0.:
            np.testing.assert_allclose(rows[:, 3].max(), extension["diagnostic_peak_elbow_deg"], atol=1e-8)
            np.testing.assert_allclose(rows[:, 4].max(), extension["peak_slide_mm"], atol=1e-8)
        # Stored actuator columns normally hold motor control; the extended
        # elbow motor converts its 25 N control to 30 N of tendon tension.
        if extended:
            elbow = model.actuator("FJ1486").id
            rows[:, 5 + elbow] *= abs(model.actuator_gear[elbow, 0])
        collision = None
        pin_collision = None
        for index, row in enumerate(rows):
            data.qpos[:] = row[-model.nq:]
            mujoco.mj_kinematics(model, data)
            positions = transforms(model, data, origins)
            if pin_collision is None:
                hits = find_hits([solids[pin].transform(positions[owners[pin]]), solids[rib]],
                                 [keys[pin], keys[rib]], local_pairs)
                if hits:
                    pin_collision = dict(time_s=row[0], hits=hits)
            if index % 4 == 0 and collision is None:
                placed = [part.transform(positions[owner]) if owner in positions else part
                          for part, owner in zip(solids, owners, strict=True)]
                hits = find_hits(placed, keys, pairs)
                if hits:
                    collision = dict(time_s=row[0], hits=hits)
            if collision and pin_collision:
                break
        label = "reference" if command == 0. else f"mu_{round(mu * 100):02d}"
        prefix = "extended_anterior_friction" if extended else "anterior_friction"
        trace = folder / f"{prefix}_{label}.csv"
        write_csv(trace, names, rows)
        stops = [event["time_s"] for event in (collision, pin_collision) if event is not None]
        clear = rows[rows[:, 0] < min(stops)] if stops else rows
        result = dict(mu=mu, anterior_input_increment_N=command, transmission=float(transmission[anterior]),
                      anterior_delivered_peak_N=float(rows[:, 5 + channels[anterior]].max()),
                      shoulder_flexion_projection_peak_deg=float(rows[:, 2].max()),
                      elbow_peak_deg=float(rows[:, 3].max()),
                      elbow_peak_before_first_sampled_collision_deg=float(clear[:, 3].max()) if len(clear) else None,
                      peak_slide_mm=float(rows[:, 4].max()), end_GH_offset_deg=float(rows[-1, 1]),
                      end_elbow_deg=float(rows[-1, 3]), first_rigid_collision=collision,
                      first_pin_rib_collision=pin_collision, sample_count=len(rows), trace=provenance(trace))
        results.append(result)
        print(mu, result["anterior_delivered_peak_N"], result["shoulder_flexion_projection_peak_deg"],
              result["elbow_peak_deg"], result["peak_slide_mm"], collision, pin_collision, flush=True)
    filename = "extended_anterior_friction_integration.json" if extended else "anterior_friction_integration.json"
    save(folder / filename, dict(
        generator=provenance(ROOT / "tools/study_s1_anterior_friction_integration.py"),
        timing=provenance(timing_file), channel=provenance(path_file), model=model_ref, extension=extension_ref,
        anterior_input_increment_N=peak, near_port_wrap_rad=wrap, cases=results,
        scope="Forward full-slip capstan attenuation of FJ1468 only, applied as a fixed sensitivity "
              "factor across the command profile in the existing ideal GH/elbow force model. "
              "Wrap is a conservative sampled local proxy from a different prescribed SC/AC path. "
              "This is neither a stick-slip law nor a reverse-sliding prediction. Existing rigid "
              "assembly is checked at 20 ms intervals and the known pin-rib pair at 5 ms. "
              "Revised TPU and SC/AC motion are excluded.",
    ))


if __name__ == "__main__":
    main()
