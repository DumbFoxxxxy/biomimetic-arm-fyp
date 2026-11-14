"""Check the two existing elbow muscle channels along the coupled replay."""

import argparse
import csv

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_spring_mount import spring_envelope
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_extensor_long_stroke import polyline_lengths, turn_radius
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def checked_trace(source):
    path = ROOT / source["file"]
    if provenance(path) != source:
        raise ValueError("Trajectory changed after rigid audit")
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    variants = parser.add_mutually_exclusive_group()
    variants.add_argument("--extended", action="store_true")
    variants.add_argument("--interface-relief", action="store_true")
    args = parser.parse_args()
    extended = args.extended or args.interface_relief
    folder = DEST.parent / "integration_review/s1_series_elbow"
    timing_path = folder / "gravity_balanced_release_timing.json"
    timing = read(timing_path)
    case = timing["cases"][timing["selected"]]
    if not case["fully_audited"] or case["first_collision"] is not None:
        raise ValueError("Need the rigid-screened timed cycle")
    baseline = load_reference(timing["baseline"])
    geometry = load_reference(baseline["geometry_assembly"])
    rigid_sources = dict(geometry["sources"])
    flexible_sources = dict(geometry["flexible_geometry_sources"])
    model_source = baseline["model"]
    trace_source = case["trace"]
    evidence = {}
    if extended:
        extension_path = folder / "elbow_extended_rail_free_drive.json"
        extension = read(extension_path)
        budget = load_reference(extension["budget"])
        if (extension["first_interference"] is not None or extension["clear_samples"] != extension["sampled_count"]
                or budget["assembly"] != baseline["geometry_assembly"]):
            raise ValueError("Extended guide did not pass the selected assembly rigid screen")
        relief_path = folder / "FJ1486_return_relief.json"
        relief = read(relief_path)
        if relief["original"] != flexible_sources["tpu:elbow_flexor"]:
            raise ValueError("Flexor relief belongs to a different assembly")
        flexible_sources["tpu:elbow_flexor"] = relief["selected"]["mesh"]
        rigid_sources.update(extension["shafts"])
        model_source = extension["model"]
        trace_source = extension["trajectory"]
        evidence = dict(extended_drive=provenance(extension_path), flexor_relief=provenance(relief_path))
        mount = load_reference(geometry["mounts"])
        guide = load_reference(mount["guide"])
        spring_axis = np.array(guide["axis_source"])
        spring_start = np.array(mount["spring_mount"]["fixed_eye_source_mm"])
        spring_end = np.array(mount["spring_mount"]["moving_eye_neutral_source_mm"])
        if args.interface_relief:
            original_audit = provenance(folder / "extended_elbow_muscle_geometry.json")
            for label, name in (("flexor", "FJ1486_interface_relief_TPU.json"),
                                ("extensor", "FJ1477_extended_relief_TPU.json")):
                study = read(folder / name)
                key = f"tpu:elbow_{label}"
                if study["input_audit"] != original_audit or study["original"] != flexible_sources[key] or study["selected"] is None:
                    raise ValueError(f"Missing compatible relief for {label}")
                flexible_sources[key] = study["selected"]["mesh"]
                evidence[f"{label}_matched_relief"] = provenance(folder / name)
    passive_file = DEST / "elbow_passive_model/model_check.json"
    passive = read(passive_file)
    elbow = load_reference(passive["assembly"])
    model_file = ROOT / model_source["file"]
    if provenance(model_file) != model_source:
        raise ValueError("Dynamics model hash changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    owners = geometry["owners"]
    fixed = bounded({key: load_component(source) for key, source in rigid_sources.items()
                     if owners[key] not in MOBILE})
    mobile = {key: load_component(source) for key, source in rigid_sources.items()
              if owners[key] in MOBILE}
    muscles = {}
    for label, specification in elbow["muscles"].items():
        head = specification["anatomical_id"]
        path = np.array(specification["wire_path_source_mm"])[1:-1]
        definition = specification["source_definition"]
        if isinstance(definition, dict):
            definition = load_reference(definition)
        else:
            definition = read(ROOT / definition)
        axis = np.array(definition["basis_columns_display"])[:, 2]
        gradient = axis / ((path[-1] - path[0]) @ axis)
        source = flexible_sources[f"tpu:elbow_{label}"]
        mesh, vertices, faces = indexed_solid(source)
        _, boundary, nonmanifold = mesh_geom.edge_census(faces)
        if boundary or nonmanifold or mesh.status() != A.I._manifold().Error.NoError:
            raise ValueError(f"Invalid native TPU mesh: {label}")
        muscles[label] = dict(head=head, path=path, path_weights=weigh(path, path[0], gradient),
                              vertices=vertices, weights=weigh(vertices, path[0], gradient), gradient=gradient,
                              faces=faces.astype(np.uint64), rest_path_mm=polyline_lengths(path).sum(),
                              source=source)

    trace = checked_trace(trace_source)
    elbow_angles = np.array([float(row["elbow_deg"]) for row in trace])
    peak = int(elbow_angles.argmax())
    samples = set(range(0, len(trace), 20)) | {peak, len(trace) - 1,
             int(np.argmax([float(row["slide_mm"]) for row in trace]))}
    for angle in range(5, int(elbow_angles[peak]), 5):
        samples.add(int(abs(elbow_angles[:peak+1] - angle).argmin()))
        samples.add(peak + int(abs(elbow_angles[peak:] - angle).argmin()))
    if args.interface_relief:
        previous = load_reference(original_audit)
        seeds = [index for index, row in enumerate(previous["records"])
                 if any(muscle.get("rigid_TPU_contacts") for muscle in row["muscles"].values())]
        samples.update(i for seed in seeds for i in range(max(0, seed - 7), min(len(trace), seed + 8)))
    manifold = A.I._manifold()
    records = []
    for index, row in enumerate(trace):
        data.qpos[:] = [float(row[f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        positions = transforms(model, data, neutral)
        body_matrix = positions["humerus"]
        humerus_rotation, humerus_translation = body_matrix[:, :3], body_matrix[:, 3]
        def local_pose(name):
            matrix = positions[name]
            return humerus_rotation.T @ matrix[:, :3], humerus_rotation.T @ (matrix[:, 3] - humerus_translation)

        rigid = {**fixed, **bounded({key: part.transform(positions[owners[key]]) for key, part in mobile.items()})}
        spring_report = None
        if extended:
            moving_eye = spring_end + float(row["slide_mm"]) * spring_axis
            spring = spring_envelope(spring_start, moving_eye, spring_axis).transform(body_matrix)
            spring_report = dict(eye_distance_mm=float(np.linalg.norm(moving_eye - spring_start)),
                                 rigid_hits=overlaps(spring, rigid))
            rigid.update(bounded({"series:return_spring_envelope": spring}))
        wires, tpu, muscle_results = {}, {}, {}
        for label, muscle in muscles.items():
            head = muscle["head"]
            distal_rotation, distal_shift = local_pose(f"{head}_distal_pivot")
            proximal_shift = (np.zeros(3) if label == "flexor" else local_pose("FJ1477_carriage")[1])
            path = displacement_map(muscle["path"], muscle["path_weights"], distal_rotation, distal_shift, proximal_shift)
            world_path = path @ humerus_rotation.T + humerus_translation
            ends = np.array([data.site(f"{head}_{name}").xpos * 1000. for name in ("origin", "distal")])
            endpoint_error = float(np.linalg.norm(world_path[[0, -1]] - ends, axis=1).max())
            if endpoint_error > 1e-7:
                raise ValueError(f"Internal path disconnected at {head}, sample {index}: {endpoint_error}")
            wire = A.D.outer_envelope(world_path, .3)
            wires[label] = wire
            muscle_results[label] = dict(path_change_mm=float(polyline_lengths(world_path).sum() - muscle["rest_path_mm"]),
                                         min_three_point_radius_mm=turn_radius(world_path),
                                         endpoint_max_error_mm=endpoint_error, rigid_hits=overlaps(wire, rigid))
            if index in samples:
                posed = displacement_map(muscle["vertices"], muscle["weights"], distal_rotation, distal_shift, proximal_shift)
                interior = (muscle["weights"] > 0.) & (muscle["weights"] < 1.)
                delta = distal_rotation - np.eye(3)
                gradients = (np.eye(3) + muscle["weights"][interior, None, None] * delta
                             + (muscle["vertices"][interior] @ delta.T + distal_shift - proximal_shift)[:, :, None]
                             * muscle["gradient"][None, None, :])
                determinant = float(np.linalg.det(gradients).min())
                if determinant <= 0.:
                    raise ValueError(f"Non-invertible TPU interpolation: {head}, {index}")
                world = posed @ humerus_rotation.T + humerus_translation
                solid = manifold.Manifold(manifold.Mesh64(np.ascontiguousarray(world), muscle["faces"]))
                if solid.status() != manifold.Error.NoError:
                    raise ValueError(f"Invalid TPU interpolation: {head}, {index}")
                tpu[label] = solid
                muscle_results[label].update(min_interpolation_determinant=determinant,
                                             rigid_TPU_contacts=overlaps(solid, rigid))
        if tpu:
            flexible = bounded(tpu)
            for label, wire in wires.items():
                muscle_results[label]["wire_TPU_hits"] = overlaps(wire, flexible)
        record = dict(time_s=float(row["time_s"]), elbow_deg=elbow_angles[index],
                      carriage_travel_mm=float(row["slide_mm"]), TPU_sample=index in samples,
                      muscles=muscle_results)
        if spring_report is not None:
            record["spring"] = spring_report
        if args.interface_relief and tpu:
            record["elbow_TPU_pair_overlap_mm3"] = overlaps(tpu["flexor"], bounded({"extensor": tpu["extensor"]})).get("extensor", 0.)
        records.append(record)
        if index % 200 == 0:
            print("Checked frame", index, flush=True)
    summaries = {}
    for label in muscles:
        rows = [row["muscles"][label] for row in records]
        hits = [r for r in records if r["muscles"][label]["rigid_hits"]]
        soft = [r for r in records if r["muscles"][label].get("wire_TPU_hits")]
        summaries[label] = dict(
            path_change_range_mm=[min(row["path_change_mm"] for row in rows), max(row["path_change_mm"] for row in rows)],
            rigid_hit_frames=len(hits), first_rigid_hit=hits[0] if hits else None,
            TPU_hit_frames=len(soft), first_TPU_hit=soft[0] if soft else None)
    spring_summary = {}
    if extended:
        spring_hits = [row for row in records if row["spring"]["rigid_hits"]]
        spring_summary = dict(spring=dict(
            eye_distance_range_mm=[min(row["spring"]["eye_distance_mm"] for row in records),
                                   max(row["spring"]["eye_distance_mm"] for row in records)],
            rigid_hit_frames=len(spring_hits), first_rigid_hit=spring_hits[0] if spring_hits else None))
    filename = "timed_elbow_muscle_geometry.json"
    if extended:
        filename = "interface_relief_elbow_muscle_geometry.json" if args.interface_relief else "extended_elbow_muscle_geometry.json"
    save(folder / filename, dict(
        generator=provenance(ROOT / "tools/audit_s1_series_elbow_timed_muscle_geometry.py"),
        timing=provenance(timing_path), passive=provenance(passive_file), trace=trace_source, **evidence,
        mesh_sources={label: row["source"] for label, row in muscles.items()},
        rigid_wire_samples=len(records), interpolated_TPU_samples=len(samples), wire_diameter_mm=.6,
        scope="The local wire uses the source internal centreline and follows imposed cap-based displacement."
              " Local bending, elasticity, friction and gravity are not solved for the wire."
              " TPU meshes are displaced kinematically; overlap volumes describe candidate contacts,"
              " not measured compression, failure or compatibility. External routing is excluded.",
        results=summaries, records=records, **spring_summary))
    print(summaries, spring_summary, flush=True)


if __name__ == "__main__":
    main()
