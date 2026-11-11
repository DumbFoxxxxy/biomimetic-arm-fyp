"""Check returned spring and retention geometry over the local slider response."""

import csv
from pathlib import Path

import mujoco
import numpy as np

from audit_elbow_cap_travel import load_geometry
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_spring_mount import spring_envelope
from build_s1_teres_origin_mount import A, load_component, rod, save
from simulate_elbow_antagonistic_load import JOINTS
from study_elbow_series_package import ADAPTER, PREFIX, bounded, overlaps, pose_solids, union


def insertion_first_hit(part, axis, offsets, obstacles):
    for offset in offsets:
        hits = overlaps(part.translate(offset * axis), obstacles)
        if hits:
            return dict(offset_mm=float(offset), intersections_mm3=hits)
    return None


def main():
    folder = DEST.parent / "elbow_local_muscle_load"
    mount_file = folder / "elbow_series_spring_mount.json"
    guide_file = folder / "elbow_series_guide.json"
    response_file = folder / "elbow_series_free_slider.json"
    mount, guide, response = map(read, (mount_file, guide_file, response_file))
    if mount["guide"] != provenance(guide_file) or response["guide"] != provenance(guide_file):
        raise RuntimeError("Spring geometry and dynamics refer to different guide versions")
    assembly_file = ROOT / guide["stop_assembly"]["file"]
    package_file = ROOT / guide["package"]["file"]
    if guide["stop_assembly"] != provenance(assembly_file) or guide["package"] != provenance(package_file):
        raise RuntimeError("Stop or package source changed")
    originals, origins = load_geometry(read(DEST / "elbow_passive_model/model_check.json"), read(assembly_file))
    axis = np.asarray(guide["axis_source"])
    bases = np.asarray(guide["selected"]["rail_origins_source_mm"])
    retention = {side: {name: load_component(source) for name, source in row["sources"].items()}
                 for side, row in mount["retention"].items() if "sources" in row}
    old_cap = PREFIX + "proximal_rigid_cap"
    moving = read(package_file)["proximal_moving_parts"]
    parts = {name: load_component(source) for name, source in guide["selected"]["sources"].items()}
    fixed = {name: part for name, (_, part) in originals.items() if name not in moving and name != ADAPTER}
    moving_fasteners = {name: originals[name][1] for name in moving if name != old_cap}
    spring = mount["spring_mount"]
    installed = None
    if spring:
        installed = {name: load_component(source) for name, source in spring["sources"].items()}
        adapter, cap = installed["adapter"], installed["cap"]
    else:
        adapter, cap = parts["adapter"], parts["moving_cap"]

    neutral = {**fixed, **moving_fasteners, "adapter": adapter, "cap": cap}
    installed_pins = {name: installed[name] for name in ("fixed_pin", "moving_pin")} if installed else {}
    fitting = {}
    for side, base in zip(("a", "b"), bases, strict=True):
        if side not in retention:
            continue
        pair = retention[side]
        other_rails = {f"{other}_{name}": part for other, row in retention.items()
                       if other != side for name, part in row.items()}
        fitting[side] = dict(
            shaft_insertion_first_hit=insertion_first_hit(pair["shaft"], axis, np.linspace(18., 0., 37),
                                                          bounded({**neutral, **other_rails, **installed_pins})),
            screw_insertion_first_hit=insertion_first_hit(
                pair["retainer_screw"], axis, np.linspace(-6., 0., 13),
                bounded({**neutral, **other_rails, **installed_pins, "shaft": pair["shaft"]})),
            driver_envelope_intersections_mm3=overlaps(
                rod(base - 41. * axis, base - 29.01 * axis, 2.5),
                bounded({**neutral, **other_rails, **installed_pins, **pair})),
        )
        blocked_out = overlaps(pair["retainer_screw"].translate(.1 * axis), bounded({"adapter": adapter}))
        fitting[side]["retainer_bearing_probe_mm3"] = blocked_out
    checked = 0
    first_failure = None
    coil_peak_span = None
    if installed:
        model_file = ROOT / response["local_model"]["file"]
        trace_file = ROOT / response["cases"][2]["trace"]["file"]
        if (provenance(model_file) != response["local_model"]
                or provenance(trace_file) != response["cases"][2]["trace"]):
            raise RuntimeError("Free slide source changed")
        model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
        data = mujoco.MjData(model)
        addresses = [model.joint(name).qposadr.item() for name in JOINTS]
        fixed_new = {"adapter": adapter, "fixed_pin": installed["fixed_pin"]}
        fixed_new.update({f"{side}_{name}": solid for side, row in retention.items()
                          for name, solid in row.items()})
        changing = {name: row for name, row in originals.items() if row[0] in origins}
        statics = {name: part for name, part in fixed.items() if name not in changing}
        other_fixed = bounded({**statics, **fixed_new})
        moving_unit = union([cap, installed["moving_pin"], *moving_fasteners.values()])
        start = np.asarray(spring["fixed_eye_source_mm"])
        end = np.asarray(spring["moving_eye_neutral_source_mm"])
        coil_peak_span = 0.
        with trace_file.open(newline="", encoding="utf-8") as stream:
            for checked, row in enumerate(csv.DictReader(stream), start=1):
                data.qpos[addresses] = np.deg2rad([float(row[f"{name}_deg"]) for name in JOINTS])
                mujoco.mj_kinematics(model, data)
                offset = float(row["carriage_mm"]) * axis
                posed_end = end + offset
                coil_peak_span = max(coil_peak_span, float(np.linalg.norm(posed_end - start)))
                envelope = spring_envelope(start, posed_end, axis)
                obstacles = bounded(pose_solids(data, changing, origins))
                posed_cap = moving_unit.translate(offset)
                coil_hits = overlaps(envelope, {**other_fixed, **obstacles, **bounded({"moving_unit": posed_cap})})
                carriage_hits = overlaps(posed_cap, {**other_fixed, **obstacles})
                fixed_hits = {name: hit for name, solid in fixed_new.items() if (hit := overlaps(solid, obstacles))}
                if coil_hits or carriage_hits or fixed_hits:
                    first_failure = dict(time_s=float(row["time_s"]), elbow_deg=float(row["elbow_deg"]),
                                         carriage_mm=float(row["carriage_mm"]), coil_hits_mm3=coil_hits,
                                         carriage_hits_mm3=carriage_hits, fixed_hits_mm3=fixed_hits)
                    break
    save(folder / "elbow_series_spring_mount_audit.json", dict(
        generator=provenance(Path(__file__).resolve()), mount=provenance(mount_file),
        guide=provenance(guide_file), response=provenance(response_file),
        retention_fit=fitting, checked_dynamic_poses=checked,
        first_failure=first_failure, maximum_spring_eye_span_mm=coil_peak_span,
        scope="Geometry only. Neutral axial insertion and simple driver envelope, conditional former "
        "free slider trajectory with spring and pin packages added. Hook topology, actual spring, "
        "pin/retainer engagement, forces, and assembly order in the complete arm are unverified."))
    print("retention checks", fitting, flush=True)
    print("dynamic poses", checked, "first failure", first_failure, flush=True)


if __name__ == "__main__":
    main()
