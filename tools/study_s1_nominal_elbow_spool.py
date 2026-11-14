"""Replay the nominal spring with a position-commanded, unilateral flexor cable."""

import argparse
import csv

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_spring_integration import load_candidate
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from simulate_elbow_antagonistic_load import write_csv
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference
from study_s1_series_elbow_release_timing import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--friction", action="store_true", help="Compare directional sliding with and without compensation")
    mode.add_argument("--seating", action="store_true", help="Replay the compensated command with seating shortening")
    args = parser.parse_args()
    folder = DEST.parent / "integration_review/s1_series_elbow"
    source_file = folder / "return_spring_nominal_wire_travel.json"
    source = read(source_file)
    table = ROOT / source["table"]["file"]
    if provenance(table) != source["table"]:
        raise ValueError("Reference travel table changed")
    with table.open(encoding="utf8", newline="") as stream:
        reference = [row for row in csv.DictReader(stream) if row["muscle"] == "FJ1486"]
    times, lengths, force = np.array([[float(row[key]) for key in ("time_s", "model_tendon_mm", "active_tension_N")]
                                      for row in reference]).T
    _, extended, _, spring, sources, owners, revisions, _ = load_candidate(folder)
    if source["model"] != extended["model"]:
        raise ValueError("Model differs from reference trajectory")
    model_file = ROOT / source["model"]["file"]
    if provenance(model_file) != source["model"]:
        raise ValueError("Reference model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {key: data.body(key).xpos.copy() * 1000. for key in MOBILE}
    timing = load_reference(load_reference(extended["budget"])["timing"])
    active = load_reference(timing["input_cycles"])
    start = load_reference(active["start"])
    shoulder = load_reference(load_reference(timing["baseline"])["shoulder_endpoints"])
    schedule = timing["cases"][timing["selected"]]
    spring_size = load_reference(spring["size"])
    rate = spring["design"]["rate_N_per_mm"][1]
    elbow = model.actuator("FJ1486").id
    tendon = model.tendon("FJ1486").id
    gear = abs(model.actuator_gear[elbow, 0])
    rib = "thorax_reference:FJ3330"
    targets = bounded({rib: load_component(sources[rib])})
    parts = {key: load_component({**sources, **revisions}[key]) for key in
             ("bone:humerus", "elbow:elbow_washer_left", "elbow:elbow_shoulder_pin")}
    stiffness = 10000.  # N/m, assumed effective axial drive stiffness
    # Invert the previous force-controlled trajectory, adding a small
    # payout allowance as commanded tension falls to zero.
    available = lengths / 1000. - force / stiffness + .0005 * np.clip(1. - force / 2., 0., 1.)
    cases = []
    seating_ref = None
    if args.seating:
        stem = "nominal_spring_spool_mu_30_adjusted"
        seating_ref = provenance(folder / f"{stem}_flexor_wall_refined_segment.json")
        seating = load_reference(seating_ref)
        wall = load_reference(seating["wall_screen"])
        travel = load_reference(wall["travel"])
        with (ROOT / travel["trace"]["file"]).open(encoding="utf8", newline="") as stream:
            if provenance(ROOT / travel["trace"]["file"]) != travel["trace"] or travel["model"] != source["model"]:
                raise ValueError("Seating correction uses a different model or trace")
            saved = list(csv.DictReader(stream))
        np.testing.assert_allclose([float(row["time_s"]) for row in saved], times, rtol=0., atol=1e-10)
        available = np.array([float(row["available_length_mm"]) for row in saved]) / 1000.
        saved_transmission = np.array([float(row["transmission"]) for row in saved])
        saved_tension = np.array([float(row["flexor_tension_N"]) for row in saved])
        # The fixed-time loss reduces extension by the same amount as adding
        # it to the effective available length. The actual command is unchanged.
        plans = [("seating_reference", np.zeros_like(times), .3, False)]
        for rate_case in (5., 20., 50.):
            selected = sorted((row for row in seating["cases"] if row["subdivisions"] == 16
                               and row["nominal_node_stiffness_N_per_mm"] == rate_case), key=lambda row: row["time_s"])
            if not selected:
                raise ValueError(f"No refined seating samples at {rate_case:g} N/mm")
            # Preserve the reduced model only over its sampled time window.
            # Taper to zero at the surrounding baseline slack intervals.
            anchors = [row["time_s"] for row in selected]
            slack = times[saved_tension <= 1e-6]
            left = slack[slack < anchors[0]][-1]
            right = slack[slack > anchors[-1]][0]
            loss = np.interp(times, [left, *anchors, right],
                             [0., *(row["arc_shortening_mm"] / 1000. for row in selected), 0.])
            plans.append((f"seating_k{rate_case:g}", loss, .3, False))
    elif args.friction:
        turn = np.array([float(row["polyline_turn_rad"]) for row in reference])
        velocity = np.array([float(row["internal_arc_speed_mm_s"]) for row in reference])
        # Prescribed local sliding direction from the reference; transition
        # smoothly through zero speed to avoid an instantaneous tension jump.
        exponent = turn * np.tanh(velocity / 10.)
        plans = [("mu_00", 0., 0., False)]
        plans.extend((f"mu_{round(mu * 100):02d}_{mode}", 0., mu, mode == "adjusted")
                     for mu in (.1, .2, .3) for mode in ("original", "adjusted"))
    else:
        exponent = np.zeros_like(times)
        plans = [(label, offset, 0., False) for label, offset in
                 (("nominal", 0.), ("short_0p5mm", -.0005), ("long_0p5mm", .0005))]
    for label, offset, mu, adjusted in plans:
        transmission = saved_transmission if args.seating else np.exp(mu * exponent)
        command = available + offset
        if adjusted:
            command += force / stiffness * (1. - 1. / transmission)
        rows = run(model, np.array(start["candidate"]["quaternion_wxyz"]),
                   [model.actuator(name).id for name in active["shoulder_channels"]],
                   np.array(active["preload_N"]), .6 * np.array(active["full_flexion_increment_N"]),
                   schedule["shoulder_release_start_s"], schedule["elbow_release_start_s"],
                   np.array(shoulder["axes_source"][0]), return_preload_N=spring_size["installed_preload_N"],
                   return_stiffness_N_per_m=1000. * rate, flexor_spool=(times, command, stiffness),
                   flexor_transmission=transmission if args.friction or args.seating else None)
        np.testing.assert_allclose(rows[:, 0], times, atol=1e-10, rtol=0.)
        measured = []
        first = None
        for row, paid in zip(rows, command, strict=True):
            data.qpos[:] = row[-model.nq:]
            mujoco.mj_forward(model, data)
            extension = float(data.ten_length[tendon] - paid)
            measured.append(extension)
            if first is None:
                frames = transforms(model, data, neutral)
                hits = {key: hit for key, part in parts.items()
                        if (hit := overlaps(part.transform(frames[owners[key]]), targets))}
                if hits:
                    first = dict(time_s=float(row[0]), elbow_deg=float(row[3]), hits_mm3=hits)
        extension = np.asarray(measured)
        tension = rows[:, 5 + elbow] * gear
        inlet = np.maximum(stiffness * extension, 0.)
        np.testing.assert_allclose(tension, np.clip(inlet * transmission, 0., gear * model.actuator_ctrlrange[elbow, 1]),
                                   atol=1e-8, rtol=0.)
        output = folder / f"nominal_spring_spool_{label}.csv"
        write_csv(output, ["time_s", "elbow_deg", "slide_mm", "available_length_mm",
                           "extension_mm", "flexor_tension_N", "inlet_tension_N", "transmission",
                           *(f"qpos_{i}" for i in range(model.nq))],
                  np.column_stack((rows[:, [0, 3, 4]], command * 1000., extension * 1000.,
                                   tension, inlet, transmission, rows[:, -model.nq:])))
        result = dict(case=label, mu=mu, command_adjusted=adjusted or args.seating,
                      peak_inlet_tension_N=float(inlet.max()), peak_elbow_deg=float(rows[:, 3].max()),
                      end_elbow_deg=float(rows[-1, 3]), peak_tension_N=float(tension.max()),
                      command_travel_mm=float(np.ptp(command) * 1000.),
                      minimum_rail_reserve_mm=float(extended["rail_stroke_mm"] - rows[:, 4].max()),
                      slack_samples=int(np.count_nonzero(extension <= 0.)),
                      peak_slack_mm=float(max(0., -extension.min()) * 1000.),
                      end_slack_mm=float(max(0., -extension[-1]) * 1000.),
                      first_known_rib_hit=first, trace=provenance(output))
        if args.seating:
            result.update(maximum_prescribed_arc_loss_mm=float(np.max(offset) * 1000.),
                          physical_command_travel_mm=float(np.ptp(available) * 1000.))
            if label == "seating_reference":
                np.testing.assert_allclose(tension, saved_tension, atol=1e-7, rtol=0.)
        else:
            result["offset_mm"] = offset * 1000.
        print(result, flush=True)
        cases.append(result)
    filename = ("nominal_spring_spool_seating.json" if args.seating else
                "nominal_spring_spool_friction.json" if args.friction else "nominal_spring_spool_dynamics.json")
    save(folder / filename, dict(
        generator=provenance(ROOT / "tools/study_s1_nominal_elbow_spool.py"), source=provenance(source_file),
        drive_model=source["model"], stiffness_N_per_m=stiffness, cases=cases, seating_study=seating_ref,
        scope="Only FJ1486 uses an open-loop available-length command derived from a prior reference. "
              "Joint motion is solved by dynamics using the equivalent tendon path. Directional "
              "friction, when enabled, is a prescribed full-slip sensitivity using reference arc "
              "velocity and wrap, with tanh smoothing at 10 mm/s. It is not a stick-slip solver. "
              "Drive damping and external routes are omitted. Three known rib pairs checked at 5 ms."
              + (" Seating is a time-prescribed sensitivity with the recorded compensated command and"
                 " friction transmission held fixed. Available length in the case CSV is the effective"
                 " length including the interpolated loss, not additional physical payout. The tendon"
                 " Jacobian and local contact are not recomputed during integration." if args.seating else ""),
    ))


if __name__ == "__main__":
    main()
