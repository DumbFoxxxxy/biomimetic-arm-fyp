"""Assemble S1 and the current elbow on one humerus and screen GH poses."""

import csv
from itertools import combinations
from pathlib import Path

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_spring_mount import spring_envelope
from build_s1_teres_origin_mount import A, load_component, save
from simulate_elbow_antagonistic_load import JOINTS


MOBILE = {"humerus", "forearm", "FJ1486_distal_pivot", "FJ1477_distal_pivot", "FJ1477_carriage"}


def assert_source(row):
    path = ROOT / row["file"]
    if provenance(path)["sha256"] != row["sha256"]:
        raise RuntimeError(f"Input fingerprint changed: {path.name}")
    return path


def combined_sources(stop, passive, shoulder, spring, short_guides, mapped):
    sources = dict(stop["sources"])
    owners = {row["source_key"]: row["owner"] for row in passive["geom_owners"]}
    owners.update({key: "humerus" for key in sources if key.startswith("stop:")})
    skipped = {}
    for key, owner in shoulder["intended_body_owners"].items():
        if key in sources:
            continue
        if owner in ("humerus_ideal_attachment", "scapula_ideal_fastened"):
            owners[key] = "humerus" if owner == "humerus_ideal_attachment" else "scapula"
        elif key in ("cap:FJ1500/origin", "cap:FJ1504/origin", "cap:FJ1507/origin"):
            owners[key] = "scapula"
        else:
            skipped[key] = owner
            continue
        sources[key] = shoulder["sources"][key]
    for name, role in (("adapter", "humerus"), ("cap", "FJ1477_carriage")):
        key = "elbow:FJ1477_bone_adapter" if name == "adapter" else "elbow:FJ1477_proximal_rigid_cap"
        sources[key] = spring["spring_mount"]["sources"][name] if name == "adapter" else short_guides["trimmed_cap"]["source"]
        owners[key] = role
    for key in sources:
        if key.startswith("elbow:FJ1477_proximal_"):
            owners[key] = "FJ1477_carriage"
    for side, candidate in spring["retention"].items():
        for name, source in candidate["sources"].items():
            key = f"series:{side}_{name}"
            sources[key], owners[key] = source, "humerus"
    for name, owner in (("fixed_pin", "humerus"), ("moving_pin", "FJ1477_carriage")):
        key = f"series:{name}"
        sources[key] = spring["spring_mount"]["sources"][name]
        owners[key] = owner
    selected = next(row for row in mapped["candidates"] if row["allowance_mm"] == mapped["selected_allowance_mm"])
    return sources, owners, skipped, selected["source"]


def transforms(model, data, original_centres):
    matrices = {}
    for name, centre in original_centres.items():
        body = data.body(name)
        rotation = body.xmat.reshape(3, 3).copy()
        matrices[name] = np.c_[rotation, body.xpos * 1000. - rotation @ centre]
    return matrices


def compose(parent, child):
    return np.c_[parent[:, :3] @ child[:, :3], parent[:, :3] @ child[:, 3] + parent[:, 3]]


def find_hits(placed, keys, pairs):
    bounds = np.asarray([part.bounding_box() for part in placed])
    left, right = bounds[pairs[:, 0]], bounds[pairs[:, 1]]
    relevant = (np.all(left[:, :3] <= right[:, 3:] + 1e-8, axis=1)
                & np.all(right[:, :3] <= left[:, 3:] + 1e-8, axis=1))
    result = []
    for i, j in pairs[relevant]:
        volume = float(A.I.intersection_volume(placed[i], placed[j]))
        if volume > .001:
            result.append(dict(left=keys[i], right=keys[j], overlap_mm3=volume))
    return sorted(result, key=lambda row: -row["overlap_mm3"])


def main():
    folder = DEST.parent / "elbow_local_muscle_load"
    output = DEST.parent / "integration_review/s1_series_elbow"
    output.mkdir(exist_ok=True)
    selection_path = ROOT / "design/S1_CURRENT_CANDIDATE.json"
    selection = read(selection_path)
    shoulder_file = assert_source({"file": selection["assembly_report"], "sha256": selection["assembly_sha256"]})
    bare_file = assert_source({"file": selection["bare_motion_report"], "sha256": selection["bare_motion_report_sha256"]})
    shoulder, bare = read(shoulder_file), read(bare_file)
    passive_file = DEST / "elbow_passive_model/model_check.json"
    passive = read(passive_file)
    response_file = folder / "elbow_series_free_slider.json"
    mount_file = folder / "elbow_series_spring_mount.json"
    short_file = folder / "elbow_slider_short_guides.json"
    mapped_file = folder / "elbow_short_guide_mapped_relief.json"
    response, mount, short, mapped = map(read, (response_file, mount_file, short_file, mapped_file))
    if mapped["parent"] != provenance(short_file) or short["spring_mount"] != provenance(mount_file):
        raise RuntimeError("Relief versions do not match the installed series candidate")
    guide_file = assert_source(response["guide"])
    stop_file = assert_source(read(guide_file)["stop_assembly"])
    stop = read(stop_file)
    sources, owners, omitted, flex_source = combined_sources(stop, passive, shoulder, mount, short, mapped)
    assert_source(flex_source)
    solids = {key: load_component(source) for key, source in sources.items()}
    # Check shared GH material independently of the local shaft edits at the elbow.
    centre = np.asarray(bare["joint_centre_source_mm"])
    if not np.allclose(centre, np.asarray(passive["neutral_origins_m"]["humerus"]) * 1000., atol=1e-7):
        raise RuntimeError("GH frames differ")
    reference_humerus = load_component(shoulder["sources"]["bone:humerus"])
    actual_humerus = solids["bone:humerus"]
    mf = A.I._manifold()
    gh_zone = mf.Manifold.cube([160., 160., 200.], True).translate(centre + [0., 0., 60.])
    shared_check = dict(centre_source_mm=centre.tolist(), proximal_comparison_min_z_mm=float(centre[2] - 40.),
                        added_in_GH_region_mm3=float(((actual_humerus - reference_humerus) ^ gh_zone).volume()),
                        removed_in_GH_region_mm3=float(((reference_humerus - actual_humerus) ^ gh_zone).volume()))
    context_differences = [key for key in sources.keys() & shoulder["sources"].keys()
                           if sources[key]["sha256"] != shoulder["sources"][key]["sha256"]]
    model_file = assert_source(response["local_model"])
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    slider = model.joint("FJ1477_series_travel").qposadr.item()
    centres = {owner: data.body(owner).xpos.copy() * 1000. for owner in MOBILE}
    case = response["cases"][2]
    trace_file = assert_source(case["trace"])
    with trace_file.open(newline="", encoding="utf-8") as stream:
        trace = list(csv.DictReader(stream))
    selected_frames = {
        "neutral": 0, "elbow_flexed": max(range(len(trace)), key=lambda i: float(trace[i]["elbow_deg"])),
        "return_at_45": min(range(1000, len(trace)), key=lambda i: abs(float(trace[i]["elbow_deg"]) - 45.)),
    }
    new_elbow = {key for key in solids if key.startswith("series:")}
    new_elbow.update(key for key in solids if owners[key] == "FJ1477_carriage" or key == "elbow:FJ1477_bone_adapter")
    shoulder_added = set(sources) - set(stop["sources"]) - new_elbow
    # The neutral scene includes all declared rigid components and records same-body
    # overlaps separately from changing GH pairs.
    keys = list(solids)
    neutral_pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                                if keys[i] in new_elbow | shoulder_added or keys[j] in new_elbow | shoulder_added])
    neutral_hits = find_hits(list(solids.values()), keys, neutral_pairs)
    earlier = {}
    for row in neutral_hits:
        left, right = row["left"], row["right"]
        row["owners"] = [owners[left], owners[right]]
        elbow_key = next((key for key in (left, right) if key.startswith("elbow:") and key in stop["sources"]), None)
        if elbow_key and (left in shoulder_added or right in shoulder_added):
            if elbow_key not in earlier:
                earlier[elbow_key] = load_component(stop["sources"][elbow_key])
            shoulder_key = right if elbow_key == left else left
            row["previous_elbow_overlap_mm3"] = float(
                A.I.intersection_volume(earlier[elbow_key], solids[shoulder_key]))
    moving_pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                               if (owners[keys[i]] in MOBILE) != (owners[keys[j]] in MOBILE)])
    axes = {key: np.asarray(value) for key, value in bare["axes_source"].items()}
    spring_definition = mount["spring_mount"]
    spring_start = np.asarray(spring_definition["fixed_eye_source_mm"])
    spring_end = np.asarray(spring_definition["moving_eye_neutral_source_mm"])
    slide_axis = model.joint("FJ1477_series_travel").axis.copy()
    results = {}
    for name, index in selected_frames.items():
        pose = trace[index]
        data.qpos[addresses] = np.deg2rad([float(pose[f"{joint}_deg"]) for joint in JOINTS])
        travel = float(pose["carriage_mm"])
        data.qpos[slider] = travel * .001
        mujoco.mj_kinematics(model, data)
        local = transforms(model, data, centres)
        spring_solid = spring_envelope(spring_start, spring_end + travel * slide_axis, slide_axis)

        def hits(axis, angle):
            gh = A.I.rotation_about(centre, axis, angle)
            world = {owner: compose(gh, matrix) for owner, matrix in local.items()}
            posed = [solid.transform(world[owners[key]]) if owners[key] in world else solid
                     for key, solid in solids.items()]
            found = find_hits(posed, keys, moving_pairs)
            spring_pose = spring_solid.transform(gh)
            spring_keys = [key for key in keys if owners[key] not in MOBILE]
            spring_solids = [spring_pose, *(solids[key] for key in spring_keys)]
            spring_pairs = np.asarray([(0, i + 1) for i in range(len(spring_keys))])
            return found + find_hits(spring_solids, ["series:spring_envelope", *spring_keys], spring_pairs)

        zero = hits(axes["flexion"], 0.)
        baseline_pairs = {(row["left"], row["right"]) for row in zero}
        motions = []
        for spec in bare["results"]:
            label, sign = spec["axis"], spec["sign"]
            # Use the existing ideal S1 axis and earlier naked-joint clear range.
            limit = abs(float(spec["last_contiguous_sample_clear_deg"]))
            last_clear, first = 0., None
            for angle in range(1, int(limit) + 1):
                def new_hits(value):
                    return [row for row in hits(axes[label], sign * value)
                            if (row["left"], row["right"]) not in baseline_pairs]

                obstacles = new_hits(angle)
                if obstacles:
                    low, high = last_clear, float(angle)
                    while high - low > .1:
                        mid = (high + low) / 2.
                        intermediate = new_hits(mid)
                        if intermediate:
                            high, obstacles = mid, intermediate
                        else:
                            low = mid
                    last_clear = low
                    first = dict(angle_deg=sign * high, hits=obstacles)
                    break
                last_clear = float(angle)
            motions.append(dict(axis=label, direction=sign, screened_limit_deg=limit,
                                last_deg_before_new_pair=sign * last_clear, first_new_pair=first))
        results[name] = dict(trace_index=index, elbow_deg=float(pose["elbow_deg"]),
                             carriage_mm=travel, gh_neutral_hits=zero, motion=motions)
        print(name, "zero", len(zero), "motions",
              [(row["axis"], row["direction"], row["last_deg_before_new_pair"]) for row in motions], flush=True)
    save(output / "assembly.json", dict(
        generator=provenance(Path(__file__).resolve()), shoulder_selection=provenance(selection_path),
        shoulder=provenance(shoulder_file), elbow_base=provenance(stop_file),
        mounts=provenance(mount_file), cap_relief=provenance(short_file), tpu_relief=provenance(mapped_file),
        sources=sources, owners={key: owners[key] for key in sources}, omitted_shoulder_parts=omitted,
        flexible_geometry_sources={"tpu:elbow_extensor": flex_source,
                                   "tpu:elbow_flexor": stop["display_only_TPU_sources"]["tpu:elbow_flexor"],
                                   **{key: source for key, source in shoulder["sources"].items() if key.startswith("tpu:")}},
        source_frame="BP3D neutral mm",
        scope="Integration candidate and neutral rigid inventory. Selected S1 ideal mounts with declared "
        "owners only. Historic bench, unresolved shoulder terminals, and flexible contacts are outside "
        "the rigid screening. Current pointers not replaced."))
    save(output / "gh_clearance.json", dict(
        assembly=provenance(output / "assembly.json"), bare_s1_motion=provenance(bare_file),
        trace=provenance(trace_file), shared_humerus=shared_check, source_differences=sorted(context_differences),
        rigid_count=len(solids), pairs_neutral_with_new_components=len(neutral_pairs),
        gh_cross_body_pairs=len(moving_pairs), neutral_new_part_overlaps=neutral_hits,
        intersection_threshold_mm3=.001, angle_step_deg=1., boundary_resolution_deg=.1, samples=results,
        scope="Prescribed single GH axes at three elbow states, with SC and AC held neutral. Initial "
        "overlaps are listed and never treated as a clearance pass; additional pairs are tracked "
        "separately. No K acceptance, shoulder actuation, or flexible collision inference."))
    print("neutral new-part overlaps", neutral_hits[:10], flush=True)


if __name__ == "__main__":
    main()
