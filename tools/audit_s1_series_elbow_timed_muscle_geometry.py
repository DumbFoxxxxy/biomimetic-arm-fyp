"""Check the two existing elbow muscle channels along the coupled replay."""

import csv

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
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
    folder = DEST.parent / "integration_review/s1_series_elbow"
    timing_path = folder / "gravity_balanced_release_timing.json"
    timing = read(timing_path)
    case = timing["cases"][timing["selected"]]
    if not case["fully_audited"] or case["first_collision"] is not None:
        raise ValueError("Need the rigid-screened timed cycle")
    baseline = load_reference(timing["baseline"])
    geometry = load_reference(baseline["geometry_assembly"])
    passive_file = DEST / "elbow_passive_model/model_check.json"
    passive = read(passive_file)
    elbow = load_reference(passive["assembly"])
    model_file = ROOT / baseline["model"]["file"]
    if provenance(model_file) != baseline["model"]:
        raise ValueError("Dynamics model hash changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    owners = geometry["owners"]
    fixed = bounded({key: load_component(source) for key, source in geometry["sources"].items()
                     if owners[key] not in MOBILE})
    mobile = {key: load_component(source) for key, source in geometry["sources"].items()
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
        source = geometry["flexible_geometry_sources"][f"tpu:elbow_{label}"]
        mesh, vertices, faces = indexed_solid(source)
        _, boundary, nonmanifold = mesh_geom.edge_census(faces)
        if boundary or nonmanifold or mesh.status() != A.I._manifold().Error.NoError:
            raise ValueError(f"Invalid native TPU mesh: {label}")
        muscles[label] = dict(head=head, path=path, path_weights=weigh(path, path[0], gradient),
                              vertices=vertices, weights=weigh(vertices, path[0], gradient), gradient=gradient,
                              faces=faces.astype(np.uint64), rest_path_mm=polyline_lengths(path).sum(),
                              source=source)

    trace = checked_trace(case["trace"])
    elbow_angles = np.array([float(row["elbow_deg"]) for row in trace])
    peak = int(elbow_angles.argmax())
    samples = set(range(0, len(trace), 20)) | {peak, len(trace) - 1,
             int(np.argmax([float(row["slide_mm"]) for row in trace]))}
    for angle in range(5, 100, 5):
        samples.add(int(abs(elbow_angles[:peak+1] - angle).argmin()))
        samples.add(peak + int(abs(elbow_angles[peak:] - angle).argmin()))
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
        records.append(dict(time_s=float(row["time_s"]), elbow_deg=elbow_angles[index],
                            carriage_travel_mm=float(row["slide_mm"]), TPU_sample=index in samples,
                            muscles=muscle_results))
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
    save(folder / "timed_elbow_muscle_geometry.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_series_elbow_timed_muscle_geometry.py"),
        timing=provenance(timing_path), passive=provenance(passive_file), trace=case["trace"],
        mesh_sources={label: row["source"] for label, row in muscles.items()},
        rigid_wire_samples=len(records), interpolated_TPU_samples=len(samples), wire_diameter_mm=.6,
        scope="The local wire uses the source internal centreline and follows imposed cap-based displacement."
              " Local bending, elasticity, friction and gravity are not solved for the wire."
              " TPU meshes are displaced kinematically; overlap volumes describe candidate contacts,"
              " not measured compression, failure or compatibility. External routing is excluded.",
        results=summaries, records=records))
    print(summaries, flush=True)


if __name__ == "__main__":
    main()
