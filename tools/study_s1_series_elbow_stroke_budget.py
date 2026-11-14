"""Estimate required slider extension and screen rebuilt retention shafts."""

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_spring_mount import spring_envelope
from build_s1_teres_origin_mount import load_component, rod, save
from simulate_elbow_series_slider import SLIDE
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference
from study_s1_series_elbow_release_timing import run


def rail(base, axis, travel):
    def cylinder(start, end, radius):
        return rod(base + start * axis, base + end * axis, radius)

    return (cylinder(-27., travel - 6., 2.) + cylinder(travel - 10., travel - 6., 4.)
            - cylinder(-27.1, -22., 1.23))


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
    mount = load_reference(assembly["mounts"])
    guide = load_reference(mount["guide"])
    model_file = ROOT / baseline["model"]["file"]
    if provenance(model_file) != baseline["model"]:
        raise ValueError("Current shoulder–elbow model changed")
    xml = model_file.read_text(encoding="utf-8")
    nominal = guide["travel_mm"]
    pose = np.array(start["candidate"]["quaternion_wxyz"])
    preload = np.array(active["preload_N"])
    delta = .6 * np.array(active["full_flexion_increment_N"])
    flexion = np.array(shoulder["axes_source"][0])
    checks = []
    trajectories = {}
    numerical_upper = 80.
    for force, stiffness in ((25., 200.), (30., 200.), (30., 160.), (35., 240.)):
        model = mujoco.MjModel.from_xml_string(xml)
        model.jnt_range[model.joint(SLIDE).id, 1] = numerical_upper / 1000.
        model.tendon_stiffness[model.tendon("FJ1486").id] = stiffness
        model.actuator_gear[model.actuator("FJ1486").id, 0] = -force / 25.
        ids = [model.actuator(name).id for name in active["shoulder_channels"]]
        rows = run(model, pose, ids, preload, delta, selected["shoulder_release_start_s"],
                   selected["elbow_release_start_s"], flexion)
        first = np.flatnonzero(rows[:, 3] >= 116.)
        last = first[0] if len(first) else len(rows)-1
        peak = float(rows[:, 4].max())
        checks.append(dict(flexor_force_N=force, flexor_stiffness_N_per_m=stiffness,
                           peak_elbow_deg=float(rows[:, 3].max()), peak_slide_mm=peak,
                           extra_for_full_cycle_mm=max(0., peak - nominal),
                           reached_116_sample_time_s=float(rows[last, 0]) if len(first) else None,
                           peak_slide_through_first_116_mm=float(rows[:last+1, 4].max()) if len(first) else None,
                           reached_numerical_upper=bool(peak >= numerical_upper)))
        trajectories[force, stiffness] = rows
        print("Stroke demand", checks[-1], flush=True)
    # Probe the original axial layout. The screw stays at the existing socket,
    # and both positive axial stop faces move by the same candidate extension.
    model = mujoco.MjModel.from_xml_string(xml)
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    owners = assembly["owners"]
    solids = {name: load_component(row) for name, row in assembly["sources"].items()}
    shaft_keys = [f"series:{name}_shaft" for name in ("a", "b")]
    axis = np.array(guide["axis_source"])
    bases = np.array(guide["selected"]["rail_origins_source_mm"])
    for base, name in zip(bases, shaft_keys, strict=True):
        original = solids[name]
        rebuilt = rail(base, axis, nominal)
        difference = (rebuilt - original).volume() + (original - rebuilt).volume()
        if difference > .001:
            raise ValueError(f"Original shaft construction did not reproduce for {name}: {difference}")
    spring = mount["spring_mount"]
    spring_a = np.array(spring["fixed_eye_source_mm"])
    spring_b = np.array(spring["moving_eye_neutral_source_mm"])
    benchmark = trajectories[30., 200.]
    geometry = []
    for extra in (2., 4., 6., 10., 15., 20.):
        travel = nominal + extra
        modified = {name: rail(base, axis, travel) for name, base in zip(shaft_keys, bases, strict=True)}
        addition = {name: part - solids[name] for name, part in modified.items()}
        first = None
        checked = 0
        for index in sorted(set(range(0, len(benchmark), 4)) | {int(benchmark[:, 4].argmax())}):
            row = benchmark[index]
            if row[4] > travel:
                continue
            checked += 1
            data.qpos[:] = row[-model.nq:]
            mujoco.mj_kinematics(model, data)
            frames = transforms(model, data, neutral)
            context = bounded({name: part.transform(frames[owners[name]]) if owners[name] in frames else part
                               for name, part in solids.items() if name not in shaft_keys})
            context.update(bounded({"series:return_spring_envelope": spring_envelope(
                spring_a, spring_b + row[4] * axis, axis).transform(frames["humerus"])}))
            hits = {name: overlaps(part.transform(frames["humerus"]), context) for name, part in addition.items()}
            hits = {name: collisions for name, collisions in hits.items() if collisions}
            if hits:
                first = dict(time_s=float(row[0]), elbow_deg=float(row[3]), slide_mm=float(row[4]), hits=hits)
                break
        geometry.append(dict(extra_travel_mm=extra, nominal_stroke_mm=travel,
                             checked_samples=checked, first_added_material_interference=first,
                             covers_benchmark_peak=bool(travel >= benchmark[:, 4].max())))
        print("Rail", geometry[-1], flush=True)
    save(folder / "elbow_series_stroke_budget.json", dict(
        generator=provenance(ROOT / "tools/study_s1_series_elbow_stroke_budget.py"),
        timing=provenance(timing_file), model=baseline["model"], assembly=baseline["geometry_assembly"],
        guide=mount["guide"], nominal_stroke_mm=nominal, virtual_upper_limit_mm=numerical_upper,
        dynamics=checks, rail_geometry=geometry,
        scope="Only the software upper limit is extended for the stroke-demand simulations."
              " The axial rail probes reproduce the old shafts and extend both together."
              " The 30 N, 200 N/m trajectory screens new shaft material at 20 ms plus peak travel."
              " The extension geometry, retention strength, longer spring selection, flexible clearances"
              " and physical capability remain unverified."))


if __name__ == "__main__":
    main()
