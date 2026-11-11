"""Compare shoulder–elbow release timing under the same local force model."""

from itertools import combinations

import mujoco
import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from simulate_elbow_series_slider import SLIDE, SLIDE_DAMPING, SLIDE_STIFFNESS
from study_s1_local_force_geometry import load_reference


def run(model, initial, channels, preload, delta, shoulder_release, elbow_release, flexion_axis):
    data = mujoco.MjData(model)
    slide = model.joint(SLIDE)
    slide_q, slide_v = slide.qposadr.item(), slide.dofadr.item()
    gh_q = model.joint("GH").qposadr.item()
    elbow_q = model.joint("elbow").qposadr.item()
    elbow = model.actuator("FJ1486").id
    model.dof_frictionloss[slide_v] = .15
    mujoco.mj_setConst(model, data)
    data.qpos[gh_q:gh_q + 4] = initial
    initial_turn = Rotation.from_quat(np.roll(initial, -1))
    dt = model.opt.timestep
    final = round(4. / dt)
    cadence = round(.005 / dt)
    rows = []
    for step in range(final + 1):
        time = step * dt
        mujoco.mj_step1(model, data)
        shoulder_level = np.interp(time, [.2, .65, shoulder_release, shoulder_release + .4], [0., 1., 1., 0.])
        elbow_level = np.interp(time, [.2, .65, elbow_release, elbow_release + .4], [0., 1., 1., 0.])
        data.ctrl[channels] = preload + shoulder_level * delta
        data.ctrl[elbow] = 25. * elbow_level
        data.qfrc_applied[slide_v] = -max(0., .25 + SLIDE_STIFFNESS * data.qpos[slide_q]
                                       + SLIDE_DAMPING * data.qvel[slide_v])
        if step % cadence == 0 or step == final:
            current = Rotation.from_quat(np.roll(data.qpos[gh_q:gh_q + 4], -1))
            relative = current * initial_turn.inv()
            flexion = np.rad2deg(relative.as_rotvec()) @ flexion_axis
            rows.append([time, np.rad2deg(relative.magnitude()), flexion, np.rad2deg(data.qpos[elbow_q]),
                         data.qpos[slide_q] * 1000., *data.ctrl, *data.qpos])
        if step < final:
            mujoco.mj_step2(model, data)
    if not np.isfinite(rows).all() or any(warning.number for warning in data.warning):
        raise RuntimeError("Invalid dynamic response")
    return np.array(rows)


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    active_path = folder / "gravity_balanced_active_cycles.json"
    active = load_reference(provenance(active_path))
    start = load_reference(active["start"])
    baseline = load_reference(active["baseline"])
    geometry = load_reference(baseline["geometry_assembly"])
    shoulder = load_reference(baseline["shoulder_endpoints"])
    model_file = ROOT / baseline["model"]["file"]
    if provenance(model_file) != baseline["model"]:
        raise ValueError("Source MuJoCo model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    origins = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    owners = geometry["owners"]
    keys = list(geometry["sources"])
    solids = [load_component(geometry["sources"][key]) for key in keys]
    pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                      if owners[keys[i]] != owners[keys[j]]
                      and (owners[keys[i]] in MOBILE or owners[keys[j]] in MOBILE)])
    pin = keys.index("elbow:elbow_shoulder_pin")
    rib = keys.index("thorax_reference:FJ3330")
    first_pair = np.array([[0, 1]])
    pair_keys = [keys[pin], keys[rib]]
    channels = [model.actuator(name).id for name in active["shoulder_channels"]]
    axis = np.array(shoulder["axes_source"][0])
    preload = np.array(active["preload_N"])
    delta = .6 * np.array(active["full_flexion_increment_N"])
    initial = np.array(start["candidate"]["quaternion_wxyz"])
    names = ["time_s", "GH_delta_deg", "flexion_projection_deg", "elbow_deg", "slide_mm",
             *(f"{model.actuator(i).name}_input_N" for i in range(model.nu)),
             *(f"qpos_{i}" for i in range(model.nq))]

    def first_interference(rows, only_pin_rib):
        for row in rows:
            data.qpos[:] = row[-model.nq:]
            mujoco.mj_kinematics(model, data)
            positions = transforms(model, data, origins)
            if only_pin_rib:
                placed = [solids[pin].transform(positions[owners[keys[pin]]]), solids[rib]]
                hits = find_hits(placed, pair_keys, first_pair)
            else:
                placed = [part.transform(positions[owners[key]]) if owners[key] in positions else part
                          for key, part in zip(keys, solids, strict=True)]
                hits = find_hits(placed, keys, pairs)
            if hits:
                return dict(time_s=row[0], GH_delta_deg=row[1], flexion_projection_deg=row[2],
                            elbow_deg=row[3], hits=hits)
        return None

    results, trajectories = {}, {}
    options = [
        ("synchronous", 1.15, 1.15, 0.),
        ("shoulder_early", .85, 1.15, 0.),
        ("elbow_early", 1.15, .85, 0.),
        ("elbow_delay_025", 1.15, 1.40, 0.),
        ("elbow_delay_050", 1.15, 1.65, 0.),
        ("elbow_delay_075", 1.15, 1.90, 0.),
        ("elbow_delay_100", 1.15, 2.15, 0.),
        ("anterior_pull_2N", 1.15, 1.15, 2.),
    ]
    for name, shoulder_release, elbow_release, extra in options:
        input_delta = delta.copy()
        input_delta[active["shoulder_channels"].index("FJ1468")] += extra
        rows = run(model, initial, channels, preload, input_delta, shoulder_release, elbow_release, axis)
        first = first_interference(rows, only_pin_rib=True)
        trajectories[name] = rows
        results[name] = dict(shoulder_release_start_s=shoulder_release, elbow_release_start_s=elbow_release,
                             added_anterior_peak_N=extra, sampled_frames=len(rows), first_pin_rib_collision=first,
                             model_peak_flexion_projection_deg=rows[:, 2].max(), model_peak_elbow_deg=rows[:, 3].max(),
                             final_GH_offset_deg=rows[-1, 1], final_elbow_deg=rows[-1, 3])
        print(name, "target collision", None if first is None else first["time_s"], flush=True)
    if results["synchronous"]["first_pin_rib_collision"] is None:
        raise RuntimeError("Synchronous control did not reproduce prior interference")
    candidates = sorted((name for name, row in results.items() if row["first_pin_rib_collision"] is None),
                        key=lambda name: (results[name]["added_anterior_peak_N"] != 0.,
                                          abs(results[name]["elbow_release_start_s"]
                                              - results[name]["shoulder_release_start_s"])))
    selected = None
    for name in candidates:
        rows = trajectories[name]
        first = first_interference(rows, only_pin_rib=False)
        file = folder / f"gravity_balanced_timing_{name}.csv"
        write_csv(file, names, rows)
        results[name].update(first_collision=first, trace=provenance(file),
                             fully_audited=True,
                             clear_frames=len(rows) if first is None else int(np.sum(rows[:, 0] < first["time_s"])))
        print(name, "full audit", "clear" if first is None else first, flush=True)
        if first is None:
            selected = name
            break
    save(folder / "gravity_balanced_release_timing.json", dict(
        generator=provenance(ROOT / "tools/study_s1_series_elbow_release_timing.py"),
        input_cycles=provenance(active_path), baseline=active["baseline"], cases=results, selected=selected,
        period_s=.005, cross_body_pair_count=len(pairs),
        scope="One fixed starting pose and provisional cable force parameters. A target pair is screened"
              " first, then the shortlisted trajectories are checked for rigid cross-body interference."
              " Wire/TPU shape, contact dynamics outside the elbow stop, and physical repeatability are absent."))


if __name__ == "__main__":
    main()
