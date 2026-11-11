"""Screen axial load and release at the existing elbow muscle interfaces."""

import argparse
import csv
import json

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read


def unit(vector):
    return vector / np.linalg.norm(vector)


def run_cycle(direction, stiffness, damping, mass, peak):
    axis = " ".join(str(value) for value in direction)
    model = mujoco.MjModel.from_xml_string(f"""
        <mujoco model="local muscle axial bench">
          <option timestep="0.0005" gravity="0 0 0" integrator="implicitfast"/>
          <worldbody>
            <body name="moving_cap">
              <joint type="slide" axis="{axis}" stiffness="{stiffness}" damping="{damping}"/>
              <inertial pos="0 0 0" mass="{mass}" diaginertia="1e-6 1e-6 1e-6"/>
            </body>
          </worldbody>
        </mujoco>
    """)
    data = mujoco.MjData(model)
    cap = model.body("moving_cap").id
    times = [0., 1., 2., 3., 4.]
    loads = [0., peak, peak, 0., 0.]
    trace = []
    for step in range(8001):
        time = step * model.opt.timestep
        force = float(np.interp(time, times, loads))
        data.xfrc_applied[cap, :3] = force * direction
        displacement = float(data.qpos[0])
        velocity = float(data.qvel[0])
        tpu_force = stiffness * displacement + damping * velocity
        if step % 20 == 0:
            trace.append((time, force, displacement * 1000., tpu_force, tpu_force - force))
        if step < 8000:
            mujoco.mj_step(model, data)
    values = np.asarray(trace)
    expected = peak / stiffness * 1000.
    error = abs(values[200, 2] - expected)
    if not np.isfinite(values).all() or any(w.number for w in data.warning):
        raise RuntimeError("Invalid axial bench simulation")
    return trace, dict(
        stiffness_N_per_m=stiffness, damping_Ns_per_m=damping, moving_mass_kg=mass,
        peak_force_N=peak, equilibrium_shortening_mm=expected,
        simulated_shortening_at_end_of_hold_mm=float(values[200, 2]),
        end_of_hold_error_mm=float(error), peak_shortening_mm=float(values[:, 2].max()),
        shortening_at_end_of_release_mm=float(values[300, 2]),
        shortening_at_end_of_recovery_mm=float(values[-1, 2]),
        max_abs_TPU_reaction_N=float(abs(values[:, 3]).max()),
        max_abs_combined_fixed_fixture_load_N=float(abs(values[:, 4]).max()),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-force", type=float, default=1., help="Assumed tensile load, N")
    args = parser.parse_args()
    if not np.isfinite(args.peak_force) or args.peak_force <= 0:
        parser.error("--peak-force must be finite and positive")

    passive_path = DEST / "elbow_passive_model/model_check.json"
    passive = read(passive_path)
    assembly_path = ROOT / passive["assembly"]["file"]
    if provenance(assembly_path) != passive["assembly"]:
        raise RuntimeError("Local interface assembly fingerprint changed")
    assembly = read(assembly_path)
    output = DEST.parent / "elbow_local_muscle_load"
    output.mkdir(exist_ok=True)
    results = {}
    for label, muscle in assembly["muscles"].items():
        head = muscle["anatomical_id"]
        proximal, distal = (
            np.asarray(next(port["source_position_mm"] for port in passive["ports"]
                            if port["head"] == head and port["role"] == role))
            for role in ("origin", "distal")
        )
        direction = unit(proximal - distal)
        terminal_tangent = unit(np.asarray(muscle["wire_path_source_mm"][-3]) - distal)
        difference = np.rad2deg(np.arccos(np.clip(direction @ terminal_tangent, -1., 1.)))
        results[head] = dict(
            role=label, proximal_fixed_bore_source_mm=proximal.tolist(),
            distal_moving_bore_source_mm=distal.tolist(), distal_pull_direction=direction.tolist(),
            proximal_pull_direction=(-direction).tolist(),
            neutral_bore_span_mm=float(np.linalg.norm(proximal - distal)),
            nominal_internal_tangent_angle_to_axial_pull_deg=float(difference), cases=[],
        )
        for stiffness in (100., 200., 400.):
            trace, metrics = run_cycle(direction, stiffness, 1., .025, args.peak_force)
            trace_path = output / f"{head}_k{stiffness:g}.csv"
            with trace_path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(("time_s", "input_N", "shortening_mm", "TPU_fixture_load_N",
                                 "combined_fixed_fixture_load_N"))
                writer.writerows(trace)
            results[head]["cases"].append(dict(**metrics, trace=provenance(trace_path)))
    report = dict(
        interface_report=provenance(passive_path), interface_assembly=provenance(assembly_path),
        solver=mujoco.__version__, timestep_s=.0005, trace_interval_s=.01,
        load_profile="1 s ramp, 1 s hold, 1 s release, 1 s recovery",
        assumptions="One cap fixed; other cap guided along bore chord in a virtual bench. "
        "Equal and opposite actuator forces; ideal linear axial spring and dashpot represent TPU. "
        "Stiffness, damping and moving mass are assumed. No resolved lattice, cable, contact or elbow motion.",
        muscles=results,
    )
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for head, result in results.items():
        print(f"{head}: span {result['neutral_bore_span_mm']:.2f} mm, "
              f"distal tangent offset {result['nominal_internal_tangent_angle_to_axial_pull_deg']:.2f} deg")
        for case in result["cases"]:
            print(f"  k={case['stiffness_N_per_m']:g} N/m: "
                  f"hold {case['simulated_shortening_at_end_of_hold_mm']:.4f} mm, "
                  f"residual {case['shortening_at_end_of_recovery_mm']:.6f} mm")


if __name__ == "__main__":
    main()
