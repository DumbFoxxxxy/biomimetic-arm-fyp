"""Screen mounting and travel space for a proximal FJ1477 series slider."""

import csv
from pathlib import Path

import mujoco
import numpy as np

from audit_elbow_cap_travel import load_geometry
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from simulate_elbow_antagonistic_load import JOINTS


PREFIX = "elbow:FJ1477_"
ADAPTER = PREFIX + "bone_adapter"


def overlaps(part, solids, threshold=.001):
    bounds = np.asarray(part.bounding_box())
    result = {}
    for name, (target, other) in solids.items():
        if np.all(bounds[:3] <= other[3:]) and np.all(other[:3] <= bounds[3:]):
            volume = float(A.I.intersection_volume(part, target))
            if volume > threshold:
                result[name] = volume
    return result


def bounded(solids):
    return {name: (part, np.asarray(part.bounding_box())) for name, part in solids.items()}


def union(solids):
    mf = A.I._manifold()
    return mf.Manifold.batch_boolean(list(solids), mf.OpType.Add)


def pose_solids(data, solids, origins):
    transforms = {}
    for owner, origin in origins.items():
        body = data.body(owner)
        rotation = body.xmat.reshape(3, 3)
        transforms[owner] = np.c_[rotation, body.xpos * 1000. - rotation @ origin]
    return {name: part.transform(transforms[owner]) if owner in transforms else part
            for name, (owner, part) in solids.items()}


def maximum_rigid_travel(package, solids, travel, axis):
    for index, offset in enumerate(np.linspace(0., travel, max(1, int(np.ceil(travel))) + 1)):
        hit = overlaps(package.translate(offset * axis), solids)
        if hit:
            return dict(first_blocked_mm=float(offset), blocked_by_mm3=hit, preceding_samples=index)
    return None


def carriage_travel(delta, direction, relief):
    projection = delta @ direction
    target = np.linalg.norm(delta) - relief
    radicand = target**2 - (delta @ delta - projection**2)
    if radicand < 0. or projection <= 0.:
        return None
    return max(0., float(projection - np.sqrt(radicand)))


def main():
    output = DEST.parent / "elbow_local_muscle_load"
    response_file = output / "elbow_attachment_response.json"
    response = read(response_file)
    long_stroke_file = output / "elbow_extensor_long_stroke.json"
    long_stroke = read(long_stroke_file)
    model_file = ROOT / response["local_model"]["file"]
    trace_file = ROOT / long_stroke["trace"]["file"]
    assembly_file = DEST / "elbow_extension_stop/assembly.json"
    passive_file = DEST / "elbow_passive_model/model_check.json"
    for path, expected in ((model_file, response["local_model"]), (trace_file, long_stroke["trace"]),
                           (assembly_file, long_stroke["rigid_bones"])):
        if provenance(path) != expected:
            raise RuntimeError(f"Source changed: {path.name}")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    origin, distal = (data.site(f"FJ1477_{role}").xpos.copy() * 1000. for role in ("origin", "distal"))
    axis = distal - origin
    span = np.linalg.norm(axis)
    axis /= span
    side = np.cross(axis, [1., 0., 0.])
    side /= np.linalg.norm(side)
    radial = np.cross(axis, side)
    assembly = read(assembly_file)
    passive = read(passive_file)
    solids, origins = load_geometry(passive, assembly)
    moving_names = [name for name in solids if name.startswith(PREFIX + "proximal_")]
    package = union(solids[name][1] for name in moving_names)
    fixed = {name: solid for name, solid in solids.items() if name not in moving_names}
    fixed_adapter = fixed[ADAPTER][1]
    fixed_bounds = bounded(pose_solids(data, fixed, origins))
    direct = maximum_rigid_travel(package, fixed_bounds, 50., axis)
    # The original adapter may need a dedicated relief or replacement. Keep the
    # rest of the actual mounting and fasteners as solid obstacles in this case.
    without_adapter = maximum_rigid_travel(package,
                                          {name: row for name, row in fixed_bounds.items() if name != ADAPTER},
                                          50., axis)
    shaft_trials = []
    for offset in (14., 18., 22., 26.):
        for azimuth in range(0, 360, 15):
            angle = np.deg2rad(azimuth)
            perpendicular = offset * (np.cos(angle) * radial + np.sin(angle) * side)
            base = origin + perpendicular - 9. * axis
            shaft = rod(base, base + 69. * axis, 2.)
            socket = rod(base, base + 6. * axis, 3.5)
            socket_fusion = float(A.I.intersection_volume(socket, fixed_adapter))
            shaft_hits = overlaps(shaft, fixed_bounds)
            socket_hits = {key: value for key, value in overlaps(socket, fixed_bounds).items() if key != ADAPTER}
            # If a shaft occupies the adapter away from the socket, it needs
            # a separate mounting cut and is not an immediately usable trial.
            shaft_beyond_socket = shaft - socket
            outside_socket_hit = float(A.I.intersection_volume(shaft_beyond_socket, fixed_adapter))
            unseated_hits = {key: value for key, value in shaft_hits.items() if key != ADAPTER}
            moving_hit = maximum_rigid_travel(package, bounded({"shaft": shaft}), 50., axis)
            shaft_trials.append(dict(radial_offset_mm=offset, azimuth_deg=azimuth,
                                     base_source_mm=base.tolist(), socket_adapter_overlap_mm3=socket_fusion,
                                     shaft_adapter_overlap_outside_socket_mm3=outside_socket_hit,
                                     shaft_other_intersections_mm3=unseated_hits,
                                     socket_other_intersections_mm3=socket_hits, moving_cap_first_overlap=moving_hit))
    with trace_file.open(newline="", encoding="utf-8") as stream:
        trace = list(csv.DictReader(stream))
    peak_index = max(range(len(trace)), key=lambda idx: float(trace[idx]["elbow_deg"]))
    sample_indices = {0, peak_index, len(trace) - 1}
    for angle in (10., 30., 60., 90.):
        for window in (range(peak_index + 1), range(peak_index, len(trace))):
            sample_indices.add(min(window, key=lambda idx: abs(float(trace[idx]["elbow_deg"]) - angle)))
    samples = []
    sample_geometry = []
    for index in sorted(sample_indices):
        row = trace[index]
        angles = [float(row[f"{name}_deg"]) for name in JOINTS]
        data.qpos[addresses] = np.deg2rad(angles)
        mujoco.mj_kinematics(model, data)
        delta = data.site("FJ1477_distal").xpos * 1000. - origin
        length = np.linalg.norm(delta)
        rail_motion = np.clip((length - span) * 50. / long_stroke["series_design_budget_at_peak"]["total_travel_mm"],
                              0., 50.)
        residual = np.linalg.norm(delta - rail_motion * axis) - span
        posed = bounded(pose_solids(data, fixed, origins))
        sample_geometry.append((float(row["time_s"]), angles[0], delta, posed))
        hit = overlaps(package.translate(rail_motion * axis), posed)
        samples.append(dict(time_s=float(row["time_s"]), elbow_deg=angles[0],
                            input_rail_travel_mm=float(rail_motion), remaining_bore_growth_mm=float(residual),
                            rigid_package_overlaps_mm3=hit))
    peak_angle = float(trace[peak_index]["elbow_deg"])
    peak_delta = next(delta for _, angle, delta, _ in sample_geometry if angle == peak_angle)
    tilted_trials = []
    for tilt in (0, 10, 20, 30, 40):
        for azimuth in (range(0, 360, 30) if tilt else (0,)):
            phase = np.deg2rad(azimuth)
            tilt_rad = np.deg2rad(tilt)
            direction = np.cos(tilt_rad) * axis + np.sin(tilt_rad) * (np.cos(phase) * radial + np.sin(phase) * side)
            required = carriage_travel(peak_delta, direction, 50.)
            if required is None:
                continue
            fixed_no_adapter = {key: row for key, row in fixed_bounds.items() if key != ADAPTER}
            neutral_block = maximum_rigid_travel(package, fixed_no_adapter, required, direction)
            failures = []
            for time, elbow, delta, posed in sample_geometry:
                requested_relief = np.clip((np.linalg.norm(delta) - span) * 50.
                                           / long_stroke["series_design_budget_at_peak"]["total_travel_mm"], 0., 50.)
                travel = carriage_travel(delta, direction, requested_relief)
                if travel is None:
                    failures.append(dict(time_s=time, elbow_deg=elbow, unreachable_relief_mm=float(requested_relief)))
                    continue
                hit = overlaps(package.translate(travel * direction),
                               {key: row for key, row in posed.items() if key != ADAPTER})
                if hit:
                    failures.append(dict(time_s=time, elbow_deg=elbow, travel_mm=travel, overlaps_mm3=hit))
            tilted_trials.append(dict(tilt_deg=tilt, azimuth_deg=azimuth, direction_source=direction.tolist(),
                                     peak_translation_for_50mm_relief=required,
                                     neutral_travel_first_overlap_with_adapter_omitted=neutral_block,
                                     sampled_cycle_failures_with_adapter_omitted=failures))
    mf = A.I._manifold()
    fitting_keepout = mf.Manifold.cube([.5] * 3, True)
    fasteners = [fixed[PREFIX + f"bone_{index}_front_washer"][1] for index in (0, 1)]
    protected_seats = [fixed_adapter ^ washer.hull().minkowski_sum(mf.Manifold.cube([3.] * 3, True))
                       for washer in fasteners]
    relief_trials = []
    viable = [trial for trial in tilted_trials
              if trial["neutral_travel_first_overlap_with_adapter_omitted"] is None
              and not trial["sampled_cycle_failures_with_adapter_omitted"]]
    for trial in sorted(viable, key=lambda row: row["peak_translation_for_50mm_relief"]):
        direction = np.asarray(trial["direction_source"])
        travel = trial["peak_translation_for_50mm_relief"]
        delta = travel * direction
        sweep = union(mf.Manifold.batch_hull([solids[key][1], solids[key][1].translate(delta)])
                      for key in moving_names)
        cut = sweep.minkowski_sum(fitting_keepout)
        candidate = (fixed_adapter - cut).as_original().simplify(1e-6)
        components = [float(part.volume()) for part in candidate.decompose() if part.volume() > 1e-5]
        seat_loss = [float((seat - candidate).volume()) for seat in protected_seats]
        result = dict(tilt_deg=trial["tilt_deg"], azimuth_deg=trial["azimuth_deg"], travel_mm=travel,
                      removed_adapter_material_mm3=float((fixed_adapter - candidate).volume()),
                      remaining_adapter_component_volumes_mm3=components,
                      protected_washer_seat_material_removed_mm3=seat_loss)
        relief_trials.append(result)
        if len(components) != 1 or max(seat_loss) > .001:
            continue
        bounds = bounded({**{key: solid for key, (solid, _) in fixed_bounds.items() if key != ADAPTER},
                          ADAPTER: candidate})
        result["modified_adapter_travel_block"] = maximum_rigid_travel(package, bounds, travel, direction)
        if result["modified_adapter_travel_block"] is not None:
            continue
        try:
            _, source = export(candidate, output / "elbow_series_adapter_relief.obj", origin)
        except (AssertionError, ValueError) as error:
            result["export_error"] = str(error)
            continue
        result["prototype_relief_source"] = source
        reloaded = load_component(source)
        result["exported_adapter_added_mm3"] = float((reloaded - fixed_adapter).volume())
        result["exported_adapter_travel_block"] = maximum_rigid_travel(
            package, bounded({**{key: solid for key, (solid, _) in fixed_bounds.items() if key != ADAPTER},
                              ADAPTER: reloaded}), travel, direction)
        revised = {**fixed, ADAPTER: (fixed[ADAPTER][0], reloaded)}
        cycle_first_overlap = None
        maximum_slider_step = 0.
        maximum_joint_step = 0.
        prior_angles = None
        prior_motion = 0.
        for index, row in enumerate(trace):
            angles = np.asarray([float(row[f"{name}_deg"]) for name in JOINTS])
            data.qpos[addresses] = np.deg2rad(angles)
            mujoco.mj_kinematics(model, data)
            separation = data.site("FJ1477_distal").xpos * 1000. - origin
            requested_relief = np.clip((np.linalg.norm(separation) - span) * 50.
                                       / long_stroke["series_design_budget_at_peak"]["total_travel_mm"], 0., 50.)
            motion = carriage_travel(separation, direction, requested_relief)
            if motion is None or motion > travel + .001:
                cycle_first_overlap = dict(time_s=float(row["time_s"]), travel_mm=motion,
                                           reason="prescribed carriage exceeds checked travel")
                break
            if prior_angles is not None:
                maximum_joint_step = max(maximum_joint_step, float(np.max(abs(angles - prior_angles))))
                maximum_slider_step = max(maximum_slider_step, abs(motion - prior_motion))
            hit = overlaps(package.translate(motion * direction), bounded(pose_solids(data, revised, origins)))
            if hit:
                cycle_first_overlap = dict(time_s=float(row["time_s"]), elbow_deg=float(angles[0]),
                                           travel_mm=motion, overlaps_mm3=hit)
                break
            prior_angles, prior_motion = angles, motion
        result["exported_adapter_cycle"] = dict(poses_checked=index + 1,
                                                 first_failure=cycle_first_overlap,
                                                 max_joint_step_deg=maximum_joint_step,
                                                 max_slider_step_mm=maximum_slider_step)
        break
    save(output / "elbow_series_package.json", dict(
        generator=provenance(Path(__file__).resolve()), stop_assembly=provenance(assembly_file),
        pose_trace=provenance(trace_file), long_stroke=provenance(long_stroke_file),
        proximal_moving_parts=moving_names, neutral_slide_axis=axis.tolist(), bore_origin_mm=origin.tolist(),
        neutral_actual_assembly_travel=direct, neutral_adapter_replacement_assumption=without_adapter,
        shaft_radius_mm=2., shaft_length_mm=69., receiver_radius_mm=3.5, receiver_length_mm=6.,
        shaft_trials=shaft_trials, sampled_cycle=samples, tilted_trials=tilted_trials, relief_trials=relief_trials,
        protected_original_washer_seat_volumes_mm3=[float(seat.volume()) for seat in protected_seats],
        scope="50 mm prescribed sliding clearance feasibility. The moving part is the original rigid "
        "proximal cap and its fasteners. Single-shaft space and root socket overlap are checked in the "
        "neutral assembly, including existing fasteners. No complete rail pair, carriage, flexible "
        "muscle motion, spring, mounting connection, or powered response is inferred. Tilted directions "
        "initially screen with adapter omitted. The first viable exported relief is reloaded and "
        "checked with all listed rigid obstacles and original washer seat volumes."))
    print("direct", direct, "without adapter", without_adapter, flush=True)
    print("peak", next(row for row in samples if row["elbow_deg"] == float(trace[peak_index]["elbow_deg"])), flush=True)
    print("tilted directions clear with adapter omitted:", len(viable), "out of", len(tilted_trials), flush=True)
    print("adapter relief candidates retaining seats and exported:",
          [result for result in relief_trials if "prototype_relief_source" in result], flush=True)


if __name__ == "__main__":
    main()
