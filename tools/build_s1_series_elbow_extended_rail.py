"""Export and replay the 10 mm axial rail-extension candidate."""

from itertools import combinations
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import export, load_component, save
from simulate_elbow_antagonistic_load import write_csv
from simulate_elbow_series_slider import SLIDE
from study_s1_local_force_geometry import load_reference
from study_s1_series_elbow_release_timing import run
from study_s1_series_elbow_stroke_budget import rail


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    budget_file = folder / "elbow_series_stroke_budget.json"
    budget = read(budget_file)
    timing = load_reference(budget["timing"])
    base = load_reference(timing["baseline"])
    assembly = load_reference(budget["assembly"])
    guide = load_reference(budget["guide"])
    active = load_reference(timing["input_cycles"])
    start = load_reference(active["start"])
    shoulder = load_reference(base["shoulder_endpoints"])
    timing_case = timing["cases"][timing["selected"]]
    extension = next(row for row in budget["rail_geometry"]
                     if row["covers_benchmark_peak"] and row["first_added_material_interference"] is None)
    stroke = extension["nominal_stroke_mm"]
    axis = np.array(guide["axis_source"])
    metadata = {}
    for label, origin in zip(("a", "b"), guide["selected"]["rail_origins_source_mm"], strict=True):
        key = f"series:{label}_shaft"
        _, metadata[key] = export(rail(np.array(origin), axis, stroke),
                                  folder / f"extended_series_retention_{label}_shaft.obj", np.array(origin))
    sources = {**assembly["sources"], **metadata}
    source_model = ROOT / budget["model"]["file"]
    if provenance(source_model) != budget["model"]:
        raise ValueError("Current model changed")
    root = ET.parse(source_model).getroot()
    root.set("model", "S1 shoulder and elbow with 10 mm rail extension")
    root.find(f".//joint[@name='{SLIDE}']").set("range", f"0 {stroke/1000.}")
    root.find(".//motor[@name='FJ1486']").set("gear", "-1.2")
    xml = ET.tostring(root, encoding="unicode")
    model_file = folder / "elbow_extended_rail_free_drive.xml"
    model_file.write_text(xml, encoding="utf-8")
    model = mujoco.MjModel.from_xml_string(xml)
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    channels = [model.actuator(name).id for name in active["shoulder_channels"]]
    rows = run(model, np.array(start["candidate"]["quaternion_wxyz"]), channels,
               np.array(active["preload_N"]), .6 * np.array(active["full_flexion_increment_N"]),
               timing_case["shoulder_release_start_s"], timing_case["elbow_release_start_s"],
               np.array(shoulder["axes_source"][0]))
    rows[:, 5 + model.actuator("FJ1486").id] *= 1.2
    expected = next(row for row in budget["dynamics"]
                    if row["flexor_force_N"] == 30. and row["flexor_stiffness_N_per_m"] == 200.)
    np.testing.assert_allclose(rows[:, 3].max(), expected["peak_elbow_deg"], atol=1e-8, rtol=0.)
    np.testing.assert_allclose(rows[:, 4].max(), expected["peak_slide_mm"], atol=1e-8, rtol=0.)
    keys = list(sources)
    solids = [load_component(sources[key]) for key in keys]
    owners = assembly["owners"]
    pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                      if owners[keys[i]] != owners[keys[j]]
                      and (owners[keys[i]] in MOBILE or owners[keys[j]] in MOBILE)])
    first = None
    for i, row in enumerate(rows):
        data.qpos[:] = row[-model.nq:]
        mujoco.mj_kinematics(model, data)
        frames = transforms(model, data, neutral)
        placed = [part.transform(frames[owners[key]]) if owners[key] in frames else part
                  for key, part in zip(keys, solids, strict=True)]
        hits = find_hits(placed, keys, pairs)
        if hits:
            first = dict(time_s=float(row[0]), elbow_deg=float(row[3]), slide_mm=float(row[4]), hits=hits)
            break
        if i % 200 == 0:
            print("Rigid sample", i, flush=True)
    columns = ["time_s", "GH_delta_deg", "flexion_projection_deg", "elbow_deg", "slide_mm",
               *(f"{model.actuator(i).name}_input_N" for i in range(model.nu)),
               *(f"qpos_{i}" for i in range(model.nq))]
    csv = folder / "elbow_extended_rail_free_drive.csv"
    write_csv(csv, columns, rows)
    clear = rows if first is None else rows[rows[:, 0] < first["time_s"]]
    result = dict(generator=provenance(ROOT / "tools/build_s1_series_elbow_extended_rail.py"),
                  budget=provenance(budget_file), model=provenance(model_file), trajectory=provenance(csv),
                  shafts=metadata, rail_stroke_mm=stroke, flexor_force_N=30., checked_pair_count=len(pairs),
                  sampled_count=len(rows), clear_samples=len(clear), first_interference=first,
                  clear_peak_elbow_deg=float(clear[:, 3].max()) if len(clear) else None,
                  diagnostic_peak_elbow_deg=float(rows[:, 3].max()), peak_slide_mm=float(rows[:, 4].max()),
                  remaining_peak_travel_mm=float(stroke - rows[:, 4].max()), final_elbow_deg=float(rows[-1, 3]),
                  scope="Two longer fixed flanged shafts. Neutral sockets and moving cap unchanged."
                        " Rigid cross-body screening at 5 ms stops at the first failed sample."
                        " Does not validate spring selection, bending stiffness, retention strength, TPU, or wire.")
    save(folder / "elbow_extended_rail_free_drive.json", result)
    print(result, flush=True)


if __name__ == "__main__":
    main()
