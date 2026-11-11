"""Read back shared humeral supports and compare with the original assembly."""

import csv
from itertools import combinations
from pathlib import Path

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, assert_source, compose, find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_spring_mount import spring_envelope
from build_s1_teres_origin_mount import A, load_component, save
from refine_elbow_short_guide_relief import indexed_solid
from simulate_elbow_antagonistic_load import JOINTS
from study_elbow_series_package import bounded, overlaps


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    report_file = folder / "shared_anchor_bridges.json"
    report = read(report_file)
    parent_file = assert_source(report["baseline"])
    parent = read(parent_file)
    removed = set(report["archived_collar_keys"])
    replace = report["prototype_replacements"]
    if len([key for key in replace if key.startswith("elbow:")]) != 2:
        print("Both shared bridges are required before evaluating the combined replacement.", flush=True)
        return
    sources = {**{key: row for key, row in parent["sources"].items() if key not in removed}, **replace}
    owners = {key: owner for key, owner in parent["owners"].items() if key in sources}
    solids = {key: load_component(row) for key, row in sources.items()}
    keys = list(solids)
    neutral_pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                                if owners[keys[i]] != owners[keys[j]] or keys[i] in replace or keys[j] in replace])
    neutral = find_hits(list(solids.values()), keys, neutral_pairs)
    unchanged = [row for row in neutral if not {row["left"], row["right"]} & replace.keys()]
    introduced = [row for row in neutral if row not in unchanged]
    flexible, unreadable = {}, {}
    for key, row in parent["flexible_geometry_sources"].items():
        path = assert_source(row)
        try:
            flexible[key] = indexed_solid(row)[0] if path.suffix.lower() == ".obj" else load_component(row)
        except (AssertionError, ValueError, RuntimeError) as error:
            unreadable[key] = str(error) or type(error).__name__
    flexible_bounds = bounded(flexible)
    neutral_flexible_contacts = {}
    for key in replace:
        if not key.startswith("elbow:"):
            continue
        added = solids[key] - load_component(parent["sources"][key])
        neutral_flexible_contacts[key] = overlaps(added, flexible_bounds)
    response_path = DEST.parent / "elbow_local_muscle_load/elbow_series_free_slider.json"
    load_path = folder / "held_shoulder_elbow_load.json"
    response, held = map(read, (response_path, load_path))
    frame_file = assert_source(held["shoulder_frame"])
    frame = read(frame_file)
    gh = np.asarray(frame["joint_centre_source_mm"])
    model_file = assert_source(response["local_model"])
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    centres = {owner: data.body(owner).xpos.copy() * 1000. for owner in MOBILE}
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    slide = model.joint("FJ1477_series_travel")
    mount = read(assert_source(parent["mounts"]))["spring_mount"]
    start, end = (np.asarray(mount[key]) for key in ("fixed_eye_source_mm", "moving_eye_neutral_source_mm"))
    dynamic_keys = {key for key in keys if owners[key] in MOBILE and owners[key] != "humerus"}
    dynamic = {key: solids[key] for key in dynamic_keys}
    new_humeral = {key: solids[key] for key in replace}
    cases = []
    elevated = next(row for row in held["cases"] if row["axis"] == "flexion" and row["held_GH_deg"] == 30)
    for angle, source in ((0., response["cases"][2]["trace"]), (30., elevated["trace"])):
        trace_file = assert_source(source)
        with trace_file.open(newline="", encoding="utf-8") as stream:
            trace = list(csv.DictReader(stream))
        rotation = A.I.rotation_about(gh, frame["axes_source"]["flexion"], angle)
        rotated = {key: part.transform(rotation) for key, part in new_humeral.items()}
        fixed = bounded({key: part for key, part in solids.items() if owners[key] not in MOBILE})
        fixed_hits = {key: found for key, part in rotated.items() if (found := overlaps(part, fixed))}
        peak_angle = max(float(row["elbow_deg"]) for row in trace)
        first = None
        sampled = 0
        if not fixed_hits:
            for sampled, row in enumerate(trace, 1):
                q = np.deg2rad([float(row[f"{joint}_deg"]) for joint in JOINTS])
                data.qpos[addresses] = q
                travel = float(row["carriage_mm"])
                data.qpos[slide.qposadr.item()] = travel * .001
                mujoco.mj_kinematics(model, data)
                matrices = {name: compose(rotation, local) for name, local in transforms(model, data, centres).items()}
                posed = {key: part.transform(matrices[owners[key]]) for key, part in dynamic.items()}
                posed["series:spring_envelope"] = spring_envelope(start, end + slide.axis * travel,
                                                                 slide.axis).transform(rotation)
                obstacles = bounded(posed)
                collisions = {key: hit for key, part in rotated.items() if (hit := overlaps(part, obstacles))}
                if collisions:
                    first = dict(time_s=float(row["time_s"]), elbow_deg=float(row["elbow_deg"]),
                                 carriage_mm=travel, hits_mm3=collisions)
                    break
        cases.append(dict(held_GH_flexion_deg=angle, trace=provenance(trace_file), trace_peak_elbow_deg=peak_angle,
                          poses_checked=sampled, shoulder_fixed_hits_mm3=fixed_hits, first_moving_interference=first))
        print(cases[-1], flush=True)
    target = folder / "assembly_shared_anchors.json"
    save(target, dict(
        **{key: value for key, value in parent.items() if key not in ("sources", "owners", "scope")},
        sources=sources, owners=owners, shoulder_humerus_revision=provenance(report_file),
        scope="Combined neutral reference with two humeral adapter branches replacing the S1 collar. "
        "Shoulder/thorax contacts, updated tendon and flexible connections are still in review."))
    save(folder / "shared_anchor_bridges_audit.json", dict(
        generator=provenance(Path(__file__).resolve()), assembly=provenance(target), previous=provenance(parent_file),
        solid_count=len(solids), neutral_new_support_interferences=introduced,
        remaining_neutral_cross_body_interferences=unchanged, new_support_dynamic_cases=cases,
        new_support_material_neutral_flexible_overlaps_mm3=neutral_flexible_contacts,
        unchecked_flexible_sources=unreadable,
        threshold_mm3=.001, scope="Readback assessment. Dynamic solid checks target changed humeral supports "
        "against the elbow and spring at two held GH poses. They do not clear the overall arm, TPU, "
        "free terminal links, hardware installation or strength."))
    print("neutral new", len(introduced), "unchanged cross-body", len(unchanged), flush=True)


if __name__ == "__main__":
    main()
