"""Recheck the exported series guide parts over the prior elbow response."""

import argparse
import csv
from itertools import combinations
from pathlib import Path

import mujoco
import numpy as np

from audit_elbow_cap_travel import load_geometry
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from simulate_elbow_antagonistic_load import JOINTS
from study_elbow_series_package import ADAPTER, PREFIX, bounded, carriage_travel, overlaps, pose_solids, union


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--free-slider-case", type=int, choices=range(5))
    args = parser.parse_args()
    folder = DEST.parent / "elbow_local_muscle_load"
    design_file = folder / "elbow_series_guide.json"
    design = read(design_file)
    chosen = design["selected"]
    package_file = ROOT / design["package"]["file"]
    assembly_file = ROOT / design["stop_assembly"]["file"]
    for path, expected in ((package_file, design["package"]), (assembly_file, design["stop_assembly"])):
        if provenance(path) != expected:
            raise RuntimeError(f"Input fingerprint changed: {path.name}")
    package = read(package_file)
    stroke = read(ROOT / package["long_stroke"]["file"])
    run_reference = None
    trace_reference = package["pose_trace"]
    if args.free_slider_case is not None:
        run_file = folder / "elbow_series_free_slider.json"
        free = read(run_file)
        if free["guide"] != provenance(design_file):
            raise RuntimeError("Free slider trace uses a different guide")
        run_reference = provenance(run_file)
        trace_reference = free["cases"][args.free_slider_case]["trace"]
    trace_file = ROOT / trace_reference["file"]
    if provenance(trace_file) != trace_reference:
        raise RuntimeError("Pose trace changed")
    response = read(ROOT / stroke["response"]["file"])
    model_file = ROOT / response["local_model"]["file"]
    if provenance(model_file) != response["local_model"]:
        raise RuntimeError("Local stop model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    origin = data.site("FJ1477_origin").xpos.copy() * 1000.
    initial_span = np.linalg.norm(data.site("FJ1477_distal").xpos * 1000. - origin)
    axis = np.asarray(design["axis_source"])
    check = read(DEST / "elbow_passive_model/model_check.json")
    originals, origins = load_geometry(check, read(assembly_file))
    guide = {name: load_component(source) for name, source in chosen["sources"].items()}
    original_adapter = originals[ADAPTER][1]
    cap_name = PREFIX + "proximal_rigid_cap"
    old_cap = originals[cap_name][1]
    old_fasteners = {name: originals[name][1] for name in package["proximal_moving_parts"] if name != cap_name}
    cap_addition = guide["moving_cap"] - old_cap
    cap_addition_hit = overlaps(cap_addition, bounded(old_fasteners))
    new_fixed = {name: part for name, part in guide.items() if name != "moving_cap"}
    neutral = {name: part for name, (owner, part) in originals.items()
               if name not in package["proximal_moving_parts"] and name != ADAPTER and owner not in origins}
    static_checks = {}
    for name, part in new_fixed.items():
        static_checks[name] = overlaps(part - original_adapter if name == "adapter" else part, bounded(neutral))
    fixed_part_pairs = {}
    for left, right in combinations(new_fixed, 2):
        amount = float(A.I.intersection_volume(new_fixed[left], new_fixed[right]))
        if amount > .001:
            fixed_part_pairs[f"{left}|{right}"] = amount
    mount_path = ROOT / check["assembly"]["file"]
    if provenance(mount_path) != check["assembly"]:
        raise RuntimeError("Original elbow mount geometry changed")
    original_wire = read(mount_path)["muscles"]["extensor"]["wire_path_source_mm"]
    local_wire = A.D.outer_envelope(np.asarray(original_wire), .3)
    neutral_wire_hits = overlaps(local_wire, bounded({
        **new_fixed, "adapter": new_fixed["adapter"] - original_adapter, "added_cap": cap_addition,
    }))
    moving_unit = union([guide["moving_cap"], *old_fasteners.values()])
    changing = {name: row for name, row in originals.items() if row[0] in origins}
    new_parts_bounds = bounded(new_fixed)
    with trace_file.open(newline="", encoding="utf-8") as stream:
        trace = list(csv.DictReader(stream))
    max_motion, max_step, max_joint_step, prior = 0., 0., 0., None
    first_failure = None
    slide_contact = None
    for index, row in enumerate(trace):
        angles = np.asarray([float(row[f"{name}_deg"]) for name in JOINTS])
        data.qpos[addresses] = np.deg2rad(angles)
        mujoco.mj_kinematics(model, data)
        if args.free_slider_case is None:
            delta = data.site("FJ1477_distal").xpos * 1000. - origin
            relief = np.clip((np.linalg.norm(delta) - initial_span) * 50.
                             / stroke["series_design_budget_at_peak"]["total_travel_mm"], 0., 50.)
            motion = carriage_travel(delta, axis, relief)
            if motion is None or motion > design["travel_mm"] + .001:
                first_failure = dict(index=index, motion_mm=motion, reason="travel outside assessed range")
                break
        else:
            motion = float(row["carriage_mm"])
        if prior is not None:
            max_step = max(max_step, abs(motion - prior[0]))
            max_joint_step = max(max_joint_step, float(np.max(abs(angles - prior[1]))))
        max_motion = max(max_motion, motion)
        prior = (motion, angles)
        moving_existing = bounded(pose_solids(data, changing, origins))
        hits = overlaps(moving_unit.translate(motion * axis),
                        {**bounded(neutral), **moving_existing, **new_parts_bounds})
        if args.free_slider_case is not None and motion > design["travel_mm"]:
            stop_hits = {name: hits.pop(name) for name in ("collar_a", "collar_b") if name in hits}
            if stop_hits and slide_contact is None:
                slide_contact = dict(index=index, time_s=float(row["time_s"]), carriage_mm=motion,
                                     intersections_mm3=stop_hits)
        fixed_hits = {name: found for name, part in new_fixed.items()
                      if (found := overlaps(part, moving_existing))}
        if hits or fixed_hits:
            first_failure = dict(index=index, time_s=float(row["time_s"]), elbow_deg=float(angles[0]),
                                 carriage_mm=motion, moving_unit_intersections=hits,
                                 fixed_guide_intersections=fixed_hits)
            break

    end_checks = {}
    for label, direction, travel in (("proximal", -1., 0.), ("distal", 1., design["travel_mm"])):
        end_checks[label] = []
        for beyond in (0., .05, .1, .25, .5):
            motion = travel + direction * beyond
            hit = overlaps(moving_unit.translate(motion * axis), new_parts_bounds)
            end_checks[label].append(dict(overtravel_mm=beyond, intersections_mm3=hit))
    target = ("elbow_series_guide_audit.json" if args.free_slider_case is None
              else f"elbow_series_free_slider_case{args.free_slider_case}_audit.json")
    save(folder / target, dict(
        generator=provenance(Path(__file__).resolve()), design=provenance(design_file),
        assembly=provenance(assembly_file), trace=provenance(trace_file), poses_checked=index + 1,
        threshold_mm3=.001, maximum_carriage_mm=max_motion, max_carriage_step_mm=max_step,
        max_joint_step_deg=max_joint_step, first_cycle_failure=first_failure,
        dynamics_source=run_reference, first_slider_collar_intersection=slide_contact,
        moving_cap_stock_vs_existing_fasteners=cap_addition_hit,
        fixed_parts_vs_existing_static=static_checks, fixed_parts_pair_overlaps=fixed_part_pairs,
        neutral_existing_extensor_wire_vs_new_material=neutral_wire_hits, end_contact_probes=end_checks,
        scope="Exported rigid geometry only. Series travel is read from the selected pose trace. "
        "In free-slider cases, end-collar overlap during software-limit overtravel is reported separately. "
        "New guide parts and the translated cap are checked, including actual rigid elbow parts. "
        "Stop probes indicate geometric contact. Collar and shaft retention, spring, TPU, wire, "
        "force and torque capacity and assembly access remain unverified."))
    print("checked", index + 1, "maximum carriage", max_motion, "first failure", first_failure, flush=True)
    print("static pairs", static_checks, fixed_part_pairs, "cap additions", cap_addition_hit, flush=True)
    print("neutral extensor wire vs new material", neutral_wire_hits, flush=True)
    print("end probes", end_checks, flush=True)


if __name__ == "__main__":
    main()
