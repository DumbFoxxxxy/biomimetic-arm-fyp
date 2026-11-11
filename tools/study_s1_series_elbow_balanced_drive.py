"""Compare active shoulder and combined cycles from the gravity-balanced start."""

from itertools import combinations

import mujoco
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation

from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from simulate_elbow_series_slider import SLIDE, SLIDE_DAMPING, SLIDE_STIFFNESS
from study_s1_local_force_geometry import load_reference


def flexion_increment(model, data, indices, preload, axis):
    gh = model.joint("GH").dofadr.item()
    influence = []
    for actuator in indices:
        data.ctrl[:] = 0.
        data.ctrl[actuator] = 1.
        mujoco.mj_forward(model, data)
        influence.append(data.qfrc_actuator[gh:gh+3].copy())
    matrix = np.array(influence).T
    off_axis = null_space(axis.reshape(1, 3))
    fit = linprog(-(axis @ matrix), A_eq=off_axis.T @ matrix, b_eq=np.zeros(2),
                  bounds=list(zip(-preload, 25. - preload)), method="highs")
    if not fit.success or fit.fun >= 0.:
        raise RuntimeError("The existing channels cannot generate an isolated flexion increment")
    return fit.x, -fit.fun


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    start_path = folder / "gravity_balanced_start.json"
    start = load_reference(provenance(start_path))
    baseline = load_reference(start["baseline"])
    shoulder = load_reference(baseline["shoulder_endpoints"])
    geometry = load_reference(baseline["geometry_assembly"])
    xml_path = ROOT / baseline["model"]["file"]
    if provenance(xml_path) != baseline["model"]:
        raise ValueError("MuJoCo model changed")
    model = mujoco.MjModel.from_xml_string(xml_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    origins = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    keys = list(geometry["sources"])
    owners = geometry["owners"]
    solids = [load_component(geometry["sources"][key]) for key in keys]
    pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                      if owners[keys[i]] != owners[keys[j]]
                      and (owners[keys[i]] in MOBILE or owners[keys[j]] in MOBILE)])
    gh_q = model.joint("GH").qposadr.item()
    elbow_q = model.joint("elbow").qposadr.item()
    slider = model.joint(SLIDE)
    slide_q, slide_v = slider.qposadr.item(), slider.dofadr.item()
    channels = start["channels"]
    ids = [model.actuator(name).id for name in channels]
    elbow_id = model.actuator("FJ1486").id
    initial = np.array(start["candidate"]["quaternion_wxyz"])
    initial_turn = Rotation.from_quat(np.roll(initial, -1))
    initial_rotvec = initial_turn.as_rotvec()
    axis = np.array(shoulder["axes_source"][0])
    preload = np.array(start["candidate"]["shoulder_preload_N"])

    def collision():
        mujoco.mj_kinematics(model, data)
        poses = transforms(model, data, origins)
        placed = [part.transform(poses[owners[key]]) if owners[key] in poses else part
                  for key, part in zip(keys, solids, strict=True)]
        return find_hits(placed, keys, pairs)

    entry_hits = []
    entry_span = np.rad2deg(initial_turn.magnitude())
    stations = np.linspace(0., 1., int(np.ceil(entry_span / .1)) + 1)
    for fraction in stations:
        data.qpos[gh_q:gh_q + 4] = np.roll(Rotation.from_rotvec(fraction * initial_rotvec).as_quat(), 1)
        hits = collision()
        if hits:
            entry_hits.append(dict(rotation_from_original_deg=float(fraction * entry_span), hits=hits))
    print("Geometric entry:", len(entry_hits), "of", len(stations), "with intersections", flush=True)
    data.qpos[gh_q:gh_q+4] = initial
    increment, torque = flexion_increment(model, data, ids, preload, initial_turn.inv().apply(axis))
    dt = model.opt.timestep
    end_time = 4.
    cadence = round(.01 / dt)
    step_count = round(end_time / dt)
    trials = [
        ("shoulder_low", .25, 0., .15),
        ("shoulder_high", .60, 0., .15),
        ("combined_low", .25, 25., .15),
        ("combined_high", .60, 25., .15),
        ("combined_low_zero_rail_friction", .25, 25., 0.),
    ]
    outcomes = {}
    for name, scale, elbow_force, friction in trials:
        mujoco.mj_resetData(model, data)
        model.dof_frictionloss[slide_v] = friction
        mujoco.mj_setConst(model, data)
        data.qpos[gh_q:gh_q + 4] = initial
        trace = []
        first = None
        for step in range(step_count + 1):
            time = step * dt
            mujoco.mj_step1(model, data)
            command = np.interp(time, [0., .2, .65, 1.15, 1.55, end_time], [0., 0., 1., 1., 0., 0.])
            data.ctrl[ids] = preload + command * scale * increment
            data.ctrl[elbow_id] = command * elbow_force
            data.qfrc_applied[slide_v] = -max(0., .25 + SLIDE_STIFFNESS * data.qpos[slide_q]
                                           + SLIDE_DAMPING * data.qvel[slide_v])
            if step % cadence == 0 or step == step_count:
                current = Rotation.from_quat(np.roll(data.qpos[gh_q:gh_q + 4], -1))
                delta = current * initial_turn.inv()
                turn_deg = np.rad2deg(delta.magnitude())
                flex = np.rad2deg(delta.as_rotvec()) @ axis
                elbow = np.rad2deg(data.qpos[elbow_q])
                trace.append([time, turn_deg, flex, elbow, data.qpos[slide_q] * 1000., *data.ctrl, *data.qpos])
                if first is None:
                    hits = collision()
                    if hits:
                        first = dict(time_s=time, GH_delta_deg=turn_deg, flexion_projection_deg=flex,
                                     elbow_deg=elbow, hits=hits)
            if step < step_count:
                mujoco.mj_step2(model, data)
        values = np.array(trace)
        if not np.isfinite(values).all() or any(w.number for w in data.warning):
            raise RuntimeError(f"Unstable response: {name}")
        columns = ["time_s", "GH_delta_deg", "flexion_projection_deg", "elbow_deg", "slide_mm",
                   *(f"{model.actuator(i).name}_input_N" for i in range(model.nu)),
                   *(f"qpos_{i}" for i in range(model.nq))]
        saved = folder / f"gravity_balanced_{name}.csv"
        write_csv(saved, columns, trace)
        clear = values if first is None else values[values[:, 0] < first["time_s"]]
        outcomes[name] = dict(shoulder_increment_fraction=scale, elbow_input_peak_N=elbow_force,
                              rail_friction_N=friction, first_collision=first, trace=provenance(saved),
                              sampled_frames=len(values), clear_frames=len(clear),
                              clear_peak_flexion_projection_deg=float(clear[:, 2].max()) if len(clear) else None,
                              clear_peak_elbow_deg=float(clear[:, 3].max()) if len(clear) else None,
                              terminal_shoulder_offset_deg=float(values[-1, 1]),
                              terminal_elbow_deg=float(values[-1, 3]))
        print(name, outcomes[name], flush=True)
    save(folder / "gravity_balanced_active_cycles.json", dict(
        generator=provenance(ROOT / "tools/study_s1_series_elbow_balanced_drive.py"),
        start=provenance(start_path), baseline=start["baseline"], cases=outcomes, entry_sample_count=len(stations),
        entry_intersections=entry_hits, entry_step_upper_deg=.1, sample_period_s=.01,
        shoulder_channels=channels, preload_N=preload.tolist(), full_flexion_increment_N=increment.tolist(),
        full_increment_neutral_flexion_moment_Nm=float(torque),
        scope="Slerp entry is a geometric check with elbow held at neutral. After initialization the"
              " joint motion is solved under forces, gravity and the existing constraints."
              " The optimized flexion increment has zero off-axis torque only at the start."
              " Rigid audit is independent and only samples until first overlap."
              " Zero friction refers to the elbow slider; wire friction, wire geometry and TPU contact are absent."))


if __name__ == "__main__":
    main()
