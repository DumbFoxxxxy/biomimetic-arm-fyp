"""Measure internal cable travel, moment arms and bend loads on checked elbow traces."""

import argparse
import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_spring_integration import load_candidate
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from simulate_elbow_antagonistic_load import write_csv
from study_s1_local_force_geometry import load_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spool", choices=("mu_30_original", "mu_30_adjusted", "seating_k5"))
    args = parser.parse_args()
    folder = DEST.parent / "integration_review/s1_series_elbow"
    prior, extended, _, _, _, _, _, retention = load_candidate(folder)
    audit_stem = f"nominal_spring_spool_{args.spool}" if args.spool else "return_spring_nominal_rate"
    audit_file = folder / f"{audit_stem}_flexible_audit.json"
    audit_ref = provenance(audit_file)
    audit = load_reference(audit_ref)
    rigid_ref = provenance(folder / f"{audit_stem}_rigid_audit.json")
    rigid = load_reference(rigid_ref)
    if (audit["trace"] != rigid["trace"] or audit["retention"] != retention
            or rigid["retention"] != retention or rigid["findings"]
            or (args.spool and audit["spool_rigid_audit"] != rigid_ref)
            or (not args.spool and audit["findings"])):
        raise ValueError("The trajectory has missing or inconsistent geometry checks")
    trace = checked_trace(audit["trace"])
    if len(trace) != audit["nominal_rigid_and_wire_frames"] or len(trace) != rigid["sampled_frames"]:
        raise ValueError("Geometry reports do not cover the recorded trace")
    model_ref = extended["model"]
    model_path = ROOT / model_ref["file"]
    if audit["model"] != model_ref or provenance(model_path) != model_ref:
        raise ValueError("Model does not match the checked geometry")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    elbow_q = model.joint("elbow").qposadr.item()
    assembly = load_reference(load_reference(prior["passive"])["assembly"])
    channels = []
    for label, muscle in assembly["muscles"].items():
        path = np.array(muscle["wire_path_source_mm"])[1:-1]
        axis = np.array(load_reference(muscle["source_definition"])["basis_columns_display"])[:, 2]
        weight = weigh(path, path[0], axis / ((path[-1] - path[0]) @ axis))
        head = muscle["anatomical_id"]
        channels.append((label, head, path, weight, model.tendon(head).id))

    def measure():
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, neutral)
        humerus = frames["humerus"]
        rotation, translation = humerus[:, :3].T, humerus[:, 3]
        values = []
        bends = []
        for label, head, path, weights, tendon in channels:
            distal = frames[f"{head}_distal_pivot"]
            shift = (rotation @ (frames["FJ1477_carriage"][:, 3] - translation)
                     if label == "extensor" else np.zeros(3))
            posed = displacement_map(path, weights, rotation @ distal[:, :3],
                                     rotation @ (distal[:, 3] - translation), shift)
            segments = np.diff(posed, axis=0)
            lengths = np.linalg.norm(segments, axis=1)
            tangents = segments / lengths[:, None]
            turn = np.arctan2(np.linalg.norm(np.cross(tangents[:-1], tangents[1:]), axis=1),
                              np.sum(tangents[:-1] * tangents[1:], axis=1))
            values.append([lengths.sum(), data.ten_length[tendon] * 1000., turn.sum()])
            bends.append((posed[1:-1], turn, np.diff(tangents, axis=0)))
        return np.array(values), bends

    samples = []
    bend_samples = []
    epsilon = 1e-5
    for frame in trace:
        data.qpos[:] = [float(frame[f"qpos_{i}"]) for i in range(model.nq)]
        length, bends = measure()
        angle = data.qpos[elbow_q]
        data.qpos[elbow_q] = angle - epsilon
        minus, _ = measure()
        data.qpos[elbow_q] = angle + epsilon
        plus, _ = measure()
        # Positive arm pulls toward flexion. Hold slide and shoulder fixed.
        arms = (minus[:, :2] - plus[:, :2]) / (2. * epsilon)
        if args.spool:
            tensions = [float(frame["flexor_tension_N"]) if label == "flexor" else 0.
                        for label, *_ in channels]
        else:
            tensions = [float(frame[f"{head}_tension_N"]) for _, head, *_ in channels]
        samples.append(np.column_stack((length, arms, tensions)))
        if args.spool:
            bend_samples.append(bends)
    samples = np.array(samples)
    rows = []
    summary = {}
    times = np.array([float(frame["time_s"]) for frame in trace])
    for index, (label, head, *_rest) in enumerate(channels):
        arc, chord, turn, wire_arm, tendon_arm, force = samples[:, index].T
        delta = arc - arc[0]
        speed = np.gradient(arc, times)
        pull_power = force * -speed / 1000.
        moving_with_active_force = force > 1e-6
        pull_peak_index = int(np.argmax(pull_power))
        friction = []
        for mu in (0., .1, .2, .3):
            # Forward full-slip scenario only. Apply where line is being
            # taken up and the recorded active force is positive.
            take_up = (speed < 0.) & moving_with_active_force
            inlet = force[take_up] * np.exp(mu * turn[take_up])
            friction.append(dict(mu=mu, forward_samples=int(take_up.sum()),
                                 inlet_tension_peak_N=float(inlet.max()) if inlet.size else None,
                                 inlet_mechanical_power_peak_W=float(max(inlet * -speed[take_up] / 1000.))
                                 if inlet.size else None))
        summary[head] = dict(
            muscle=label, local_neutral_arc_mm=float(arc[0]), local_min_arc_mm=float(arc.min()),
            local_max_arc_mm=float(arc.max()), local_peak_to_peak_mm=float(np.ptp(arc)),
            local_peak_shortening_mm=float(max(0., -delta.min())),
            local_peak_extension_mm=float(max(0., delta.max())),
            local_return_delta_mm=float(delta[-1]), model_tendon_excursion_mm=float(np.ptp(chord)),
            local_arc_arm_mm=[float(wire_arm.min()), float(wire_arm.max())],
            model_tendon_arm_mm=[float(tendon_arm.min()), float(tendon_arm.max())],
            active_tension_peak_N=float(force.max()),
            peak_take_up_speed_mm_s=float(-speed.min()), peak_pay_out_speed_mm_s=float(speed.max()),
            peak_active_take_up_power_W=float(pull_power[pull_peak_index]),
            peak_active_take_up_power_time_s=float(times[pull_peak_index]) if moving_with_active_force.any() else None,
            local_polyline_turn_rad=[float(turn.min()), float(turn.max())],
            hypothetical_forward_full_slip=friction,
        )
        if args.spool:
            positions, angles, tangent_delta = (np.stack(items) for items in
                                                  zip(*(frame[index] for frame in bend_samples), strict=True))
            if not all(np.isfinite(field).all() for field in (positions, angles, tangent_delta)):
                raise ValueError(f"Degenerate cable segment in {head}")
            # Equal tension on both sides gives the force the cable applies
            # to a support at a polyline vertex. This does not locate contact.
            support = force[:, None, None] * tangent_delta
            load = np.linalg.norm(support, axis=2)
            i, j = np.unravel_index(np.argmax(load), load.shape)
            summary[head]["maximum_polyline_turn_deg"] = float(np.degrees(angles.max()))
            summary[head]["constant_distal_tension_bend_load"] = (dict(
                magnitude_N=float(load[i, j]), time_s=float(times[i]), internal_path_vertex=int(j + 1),
                turn_deg=float(np.degrees(angles[i, j])),
                position_in_humerus_source_axes_mm=positions[i, j].tolist(),
                vector_in_humerus_source_axes_N=support[i, j].tolist(),
                peak_sum_of_bend_load_magnitudes_N=float(load.sum(axis=1).max()),
            ) if moving_with_active_force.any() else None)
            if label == "flexor":
                applied = np.array([float(frame["transmission"]) for frame in trace])
                replay = np.exp(.3 * turn * np.tanh(speed / 10.))
                active = moving_with_active_force
                summary[head]["transmission_replay_sensitivity"] = dict(
                    maximum_active_absolute_change=float(np.max(np.abs(replay[active] - applied[active]))),
                    active_samples_opposite_slip_direction=int(np.count_nonzero(
                        active & (np.log(replay) * np.log(applied) < 0.) & (np.abs(speed) > 10.))),
                )
        for frame, values, change, velocity in zip(trace, samples[:, index], delta, speed, strict=True):
            rows.append([frame["time_s"], head, frame["elbow_deg"], *values, change, velocity])
    stem = f"{audit_stem}_wire_travel" if args.spool else "return_spring_nominal_wire_travel"
    output = folder / f"{stem}.csv"
    write_csv(output, ["time_s", "muscle", "elbow_deg", "internal_arc_mm", "model_tendon_mm", "polyline_turn_rad",
                       "internal_arc_arm_mm", "model_tendon_arm_mm", "active_tension_N",
                       "internal_arc_change_mm", "internal_arc_speed_mm_s"], rows)
    spool_scope = (
        " For spool traces, bend loads are a uniform distal-tension sensitivity at polyline vertices; "
        "wall contact and load sharing are unknown. Replay transmission uses postprocessed arc velocity "
        "and does not feed back into the saved dynamics. Passive extensor tension is excluded."
        if args.spool else ""
    )
    save(folder / f"{stem}.json", dict(
        generator=provenance(ROOT / "tools/summarize_s1_nominal_elbow_wire_travel.py"),
        flexible_audit=audit_ref, rigid_audit=rigid_ref, trace=audit["trace"], model=model_ref,
        table=provenance(output), sample_count=len(trace), channels=summary,
        scope="Internal steel path only, using the audited interpolation between its terminal supports. "
              "Moment arms are partial length derivatives holding the shoulder and slider fixed. "
              "Recorded tension is the active force command after actuator gearing. Passive TPU "
              "restoring force and external cable travel are excluded. The forward full-slip "
              "capstan scenario assumes every polyline bend participates in sliding; it does not "
              "predict contact, actual friction, or reverse sliding. Power uses simultaneous "
              "recorded active force and differentiated internal arc speed." + spool_scope,
        **({"flexible_findings": len(audit["findings"])} if args.spool else {}),
    ))
    print(summary)


if __name__ == "__main__":
    main()
