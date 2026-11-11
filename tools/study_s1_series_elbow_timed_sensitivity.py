"""Screen timing, elbow stiffness and tensile input around the checked cycle."""

from itertools import combinations, product

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from simulate_elbow_series_slider import SLIDE
from study_s1_local_force_geometry import load_reference
from study_s1_series_elbow_release_timing import run


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    timing_file = folder / "gravity_balanced_release_timing.json"
    timing = read(timing_file)
    selected = timing["cases"][timing["selected"]]
    active = load_reference(timing["input_cycles"])
    start = load_reference(active["start"])
    baseline = load_reference(timing["baseline"])
    assembly = load_reference(baseline["geometry_assembly"])
    shoulder = load_reference(baseline["shoulder_endpoints"])
    model_file = ROOT / baseline["model"]["file"]
    if provenance(model_file) != baseline["model"]:
        raise ValueError("Source MuJoCo model changed")
    xml = model_file.read_text(encoding="utf-8")
    check_model = mujoco.MjModel.from_xml_string(xml)
    data = mujoco.MjData(check_model)
    mujoco.mj_kinematics(check_model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    keys = list(assembly["sources"])
    owners = assembly["owners"]
    solids = [load_component(assembly["sources"][key]) for key in keys]
    pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                      if owners[keys[i]] != owners[keys[j]]
                      and (owners[keys[i]] in MOBILE or owners[keys[j]] in MOBILE)])
    probe_keys = ["elbow:elbow_shoulder_pin", "thorax_reference:FJ3330", "bone:forearm", "stop:arm"]
    probe_solids = [solids[keys.index(key)] for key in probe_keys]
    channels = [check_model.actuator(name).id for name in active["shoulder_channels"]]
    preload = np.array(active["preload_N"])
    delta = .6 * np.array(active["full_flexion_increment_N"])
    initial = np.array(start["candidate"]["quaternion_wxyz"])
    axis = np.array(shoulder["axes_source"][0])
    slide_limit = check_model.joint(SLIDE).range[1] * 1000.
    flexor_id = check_model.actuator("FJ1486").id
    control_column = 5 + flexor_id
    columns = ["time_s", "GH_delta_deg", "flexion_projection_deg", "elbow_deg", "slide_mm",
               *(f"{check_model.actuator(i).name}_input_N" for i in range(check_model.nu)),
               *(f"qpos_{i}" for i in range(check_model.nq))]

    def first_hit(rows, names, parts, interactions):
        for row in rows:
            data.qpos[:] = row[-check_model.nq:]
            mujoco.mj_kinematics(check_model, data)
            positions = transforms(check_model, data, neutral)
            placed = [part.transform(positions[owners[key]]) if owners[key] in positions else part
                      for key, part in zip(names, parts, strict=True)]
            hits = find_hits(placed, names, interactions)
            if hits:
                return dict(time_s=float(row[0]), elbow_deg=float(row[3]), hits=hits)
        return None

    samples = [(force, stiffness, selected["shoulder_release_start_s"])
               for force, stiffness in product((25., 30., 35.), (160., 200., 240.))]
    samples.extend((25., 200., when) for when in (.75, .95, 1.05))
    checks = {}
    traces = {}
    for force, stiffness, release in samples:
        label = f"force_{force:g}_k_{stiffness:g}_shoulder_{release:.2f}"
        model = mujoco.MjModel.from_xml_string(xml)
        model.tendon_stiffness[model.tendon("FJ1486").id] = stiffness
        # The reused protocol commands up to 25 N. A tendon gear magnitude
        # other than one scales that to the chosen physical line force.
        model.actuator_gear[flexor_id, 0] = -force / 25.
        rows = run(model, initial, channels, preload, delta, release, selected["elbow_release_start_s"], axis)
        rows[:, control_column] *= force / 25.
        if abs(rows[:, control_column].max() - force) > 1e-9:
            raise RuntimeError("Force sensitivity command was not applied")
        if (force, stiffness, release) == (25., 200., selected["shoulder_release_start_s"]):
            np.testing.assert_allclose(rows[:, 3].max(), selected["model_peak_elbow_deg"], atol=1e-8, rtol=0.)
        hit = first_hit(rows, probe_keys, probe_solids, np.array([[0, 1], [2, 3]]))
        limit = np.flatnonzero(rows[:, 4] >= slide_limit)
        first_limit = float(rows[limit[0], 0]) if len(limit) else None
        checks[label] = dict(flexor_force_N=force, axial_flexor_stiffness_N_per_m=stiffness,
                             shoulder_release_s=release, model_peak_elbow_deg=float(rows[:, 3].max()),
                             model_peak_slide_mm=float(rows[:, 4].max()), model_final_elbow_deg=float(rows[-1, 3]),
                             first_slide_limit_sample_s=first_limit, first_probe_collision=hit)
        traces[label] = rows
        print(label, "peak", round(rows[:, 3].max(), 3), "slide", round(rows[:, 4].max(), 4),
              "two-pair check", None if hit is None else hit["time_s"], flush=True)
    candidates = sorted((name for name, row in checks.items() if row["first_probe_collision"] is None
                         and row["first_slide_limit_sample_s"] is None and row["model_peak_elbow_deg"] >= 116.),
                        key=lambda name: (checks[name]["flexor_force_N"],
                                          abs(checks[name]["axial_flexor_stiffness_N_per_m"] - 200.)))
    confirmed = None
    for name in candidates:
        rows = traces[name]
        hit = first_hit(rows, keys, solids, pairs)
        checks[name]["first_full_collision"] = hit
        checks[name]["full_audit"] = True
        if hit is None:
            confirmed = name
            output = folder / "timed_elbow_K1_sensitivity_candidate.csv"
            write_csv(output, columns, rows)
            checks[name]["trace"] = provenance(output)
            break
    save(folder / "timed_elbow_sensitivity.json", dict(
        generator=provenance(ROOT / "tools/study_s1_series_elbow_timed_sensitivity.py"),
        source_timing=provenance(timing_file), source_cycle=selected["trace"],
        source_protocol=provenance(ROOT / "tools/study_s1_series_elbow_release_timing.py"),
        model=baseline["model"], sample_step_s=.005, slider_limit_mm=float(slide_limit),
        K1_floor_deg=116., K1_target_deg=145., samples=checks, confirmed_candidate=confirmed,
        scope="Five-millisecond state samples. Initial screen only checks elbow pin versus rib and"
              " forearm versus extension stop. Nominal slider bounds are recorded separately."
              " Tested stiffness and force are provisional. Geometry of the wire and TPU is excluded."))


if __name__ == "__main__":
    main()
