"""Find gravity-balanced fixed-girdle starts and recheck coupled force trajectories."""

from itertools import combinations, product

import mujoco
import numpy as np
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation

from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from simulate_elbow_series_slider import SLIDE, SLIDE_DAMPING, SLIDE_STIFFNESS
from study_s1_local_force_geometry import load_reference


def shoulder_balance(model, data, channels, block, bound=25.):
    data.ctrl[:] = 0.
    mujoco.mj_forward(model, data)
    required = (data.qfrc_bias - data.qfrc_passive)[block].copy()
    vectors = []
    for index in channels:
        data.ctrl[:] = 0.
        data.ctrl[index] = 1.
        mujoco.mj_forward(model, data)
        vectors.append(data.qfrc_actuator[block].copy())
    matrix = np.array(vectors).T
    result = linprog(np.ones(len(channels)), A_eq=matrix, b_eq=required, bounds=(0., bound), method="highs")
    return result.x if result.success and np.linalg.norm(matrix @ result.x - required) < 1e-8 else None


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    baseline_path = folder / "shoulder_elbow_free_drive.json"
    baseline = read(baseline_path)
    model_file = ROOT / baseline["model"]["file"]
    if provenance(model_file) != baseline["model"]:
        raise ValueError("Base simulation changed")
    geometry = load_reference(baseline["geometry_assembly"])
    shoulder = load_reference(baseline["shoulder_endpoints"])
    xml = model_file.read_text(encoding="utf-8")
    model = mujoco.MjModel.from_xml_string(xml)
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    origins = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    owners = geometry["owners"]
    keys = list(geometry["sources"])
    solids = [load_component(geometry["sources"][key]) for key in keys]
    pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                      if owners[keys[i]] != owners[keys[j]]
                      and (owners[keys[i]] in MOBILE or owners[keys[j]] in MOBILE)])
    gh = model.joint("GH")
    q = gh.qposadr.item()
    v = gh.dofadr.item()
    block = slice(v, v + 3)
    frame = np.array(shoulder["axes_source"])
    channels = [model.actuator(name).id for name in shoulder["neutral"]]
    elbow = model.joint("elbow")
    slide = model.joint(SLIDE)
    elbow_q, elbow_v = elbow.qposadr.item(), elbow.dofadr.item()
    slide_q, slide_v = slide.qposadr.item(), slide.dofadr.item()
    model.dof_frictionloss[slide_v] = .15
    mujoco.mj_setConst(model, data)

    def collision():
        mujoco.mj_kinematics(model, data)
        local = transforms(model, data, origins)
        placed = [solid.transform(local[owners[key]]) if owners[key] in local else solid
                  for key, solid in zip(keys, solids, strict=True)]
        return find_hits(placed, keys, pairs)

    def stability(loads):
        initial = data.qpos.copy()
        data.ctrl[:] = 0.
        data.ctrl[channels] = loads
        directions = np.eye(model.nv)[v:v + 3] * 1e-5
        tangents = []
        for velocity in directions:
            moments = []
            for sign in (-1., 1.):
                data.qpos[:] = initial
                mujoco.mj_integratePos(model, data.qpos, velocity, sign)
                mujoco.mj_forward(model, data)
                moments.append((data.qfrc_bias - data.qfrc_passive - data.qfrc_actuator)[block].copy())
            tangents.append((moments[1] - moments[0]) / (2e-5))
        data.qpos[:] = initial
        jacobian = np.array(tangents).T
        return np.linalg.eigvalsh((jacobian + jacobian.T) / 2.).tolist()

    angular_grid = [list(range(-20, 51, 5)), list(range(-15, 46, 5)), list(range(-80, 81, 10))]
    samples = []
    for parameters in product(*angular_grid):
        angles = np.deg2rad(parameters)
        turn = Rotation.from_rotvec(frame[0] * angles[0]) * Rotation.from_rotvec(frame[1] * angles[1])
        turn *= Rotation.from_rotvec(frame[2] * angles[2])
        data.qpos[q:q+4] = np.roll(turn.as_quat(), 1)
        preload = shoulder_balance(model, data, channels, block)
        if preload is None:
            continue
        hits = collision()
        eigenvalues = stability(preload) if not hits else None
        mujoco.mj_forward(model, data)
        extension_moment = float((data.qfrc_passive - data.qfrc_bias)[elbow_v])
        samples.append(dict(parameters_deg=list(parameters), quaternion_wxyz=data.qpos[q:q+4].tolist(),
                            rotation_from_original_deg=float(np.rad2deg(turn.magnitude())),
                            shoulder_preload_N=preload.tolist(), rigid_hits=hits,
                            shoulder_stiffness_eigenvalues_Nm_per_rad=eigenvalues,
                            unconstrained_neutral_elbow_moment_Nm=extension_moment))
    candidates = [row for row in samples if not row["rigid_hits"]
                  and min(row["shoulder_stiffness_eigenvalues_Nm_per_rad"]) > 0.
                  and row["unconstrained_neutral_elbow_moment_Nm"] <= 0.]
    candidates.sort(key=lambda row: (row["rotation_from_original_deg"], max(row["shoulder_preload_N"])))
    print("Tested", np.prod([len(axis) for axis in angular_grid]), "balanced", len(samples),
          "collision-free", sum(not row["rigid_hits"] for row in samples), "stable", len(candidates), flush=True)
    cases = {}
    selected = candidates[0] if candidates else None
    if selected:
        print("Start candidate", selected, flush=True)
        for name, elbow_force in (("hold", 0.), ("elbow_load_release", 25.)):
            mujoco.mj_resetData(model, data)
            data.qpos[q:q+4] = selected["quaternion_wxyz"]
            rows = []
            first_collision = None
            dt = model.opt.timestep
            step_count = round(2. / dt)
            cadence = round(.01 / dt)
            for index in range(step_count + 1):
                time = index * dt
                mujoco.mj_step1(model, data)
                data.ctrl[channels] = selected["shoulder_preload_N"]
                data.ctrl[model.actuator("FJ1486").id] = elbow_force * np.interp(
                    time, [0., .2, .65, 1.15, 1.55, 2.], [0., 0., 1., 1., 0., 0.])
                data.qfrc_applied[slide_v] = -max(0., .25 + SLIDE_STIFFNESS * data.qpos[slide_q]
                                                + SLIDE_DAMPING * data.qvel[slide_v])
                if index % cadence == 0 or index == step_count:
                    rows.append([time, np.rad2deg(data.qpos[elbow_q]), data.qpos[slide_q] * 1000., *data.qpos])
                    if first_collision is None:
                        hits = collision()
                        if hits:
                            first_collision = dict(time_s=time, elbow_deg=rows[-1][1], hits=hits)
                if index < step_count:
                    mujoco.mj_step2(model, data)
            if any(warning.number for warning in data.warning) or not np.isfinite(rows).all():
                raise RuntimeError("Unstable model run")
            trace = folder / f"gravity_balanced_{name}.csv"
            write_csv(trace, ["time_s", "elbow_deg", "slide_mm", *(f"qpos_{i}" for i in range(model.nq))], rows)
            cases[name] = dict(peak_elbow_deg=float(max(row[1] for row in rows)),
                               final_elbow_deg=float(rows[-1][1]), first_collision=first_collision,
                               trace=provenance(trace))
            print(name, cases[name], flush=True)
    save(folder / "gravity_balanced_start.json", dict(
        generator=provenance(ROOT / "tools/study_s1_series_elbow_start_pose.py"),
        baseline=provenance(baseline_path), angular_grid_deg=angular_grid, channels=list(shoulder["neutral"]),
        balanced_samples=samples, candidate=selected, cases=cases,
        scope="GH position grid used only to choose the initial pose. The prescribed shoulder preloads remain"
              " constant and elbow input is tension-only. Subsequent motion is solved, not prescribed."
              " Geometry audits rigid cross-body overlaps, excluding wire and TPU. A positive GH stiffness"
              " test is local at fixed elbow, not a complete stability proof. Girdle is fixed, loads assumed."))


if __name__ == "__main__":
    main()
