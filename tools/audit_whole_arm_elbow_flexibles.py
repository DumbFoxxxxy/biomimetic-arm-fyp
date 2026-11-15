"""Check the existing elbow TPU and retained return spring in the whole-arm replay."""

from itertools import combinations

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import find_hits, transforms
from audit_s1_series_elbow_spring_integration import load_candidate
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_series_elbow_return_spring import spring_path, sweep
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from project_runtime import load_mujoco_model
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    replay_ref = provenance(folder / "whole_arm_passive_elbow_audit.json")
    replay = load_reference(replay_ref)
    scene = load_reference(replay["scene"])
    transfer = load_reference(scene["transfer"])
    base = load_reference(replay["base"])
    elbow_folder = DEST.parent / "integration_review/s1_series_elbow"
    prior, extended, mount, spring, _, revised_owners, revisions, retention = load_candidate(elbow_folder)
    if extended != load_reference(replay["extended_rails"]):
        raise ValueError("Return spring and elbow replay require matching extended rails")
    reference = load_reference(load_reference(prior["passive"])["assembly"])
    mf = A.I._manifold()
    muscles = {}
    for label, row in reference["muscles"].items():
        definition = load_reference(row["source_definition"])
        points = np.asarray(row["wire_path_source_mm"])[1:-1]
        axis = np.asarray(definition["basis_columns_display"])[:, 2]
        gradient = axis / ((points[-1] - points[0]) @ axis)
        source = prior["mesh_sources"][label]
        solid, vertices, faces = indexed_solid(source)
        _, boundary, nonmanifold = mesh_geom.edge_census(faces)
        if boundary or nonmanifold or solid.status() != mf.Error.NoError:
            raise ValueError(f"Invalid indexed TPU mesh for {label}")
        muscles[label] = (row["anatomical_id"], points, weigh(points, points[0], gradient),
                          vertices, weigh(vertices, points[0], gradient), faces)
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
               **extended["shafts"], **revisions, "base:T12": base["source"]}
    owners = {**transfer["owner_by_source"], **{key: revised_owners[key] for key in revisions},
              "base:T12": "thorax"}
    parts = {key: load_component(meta) for key, meta in sources.items()}
    keys = list(parts)
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                        if keys[i] in revisions or keys[j] in revisions])
    model = load_mujoco_model(ROOT / scene["model"]["file"])
    data = mujoco.MjData(model)
    origins = {key: np.asarray(scene["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    start = np.asarray(mount["fixed_eye_source_mm"])
    end = np.asarray(mount["moving_eye_neutral_source_mm"])
    length = np.linalg.norm(end - start)
    axis = (end - start) / length
    pin = np.asarray(mount["pin_axis_source"])
    spring_frame = np.c_[pin, np.cross(axis, pin), axis, start]
    spec = spring["design"]
    radius = spec["wire_diameter_mm"] / 2.
    mean = (spec["coil_OD_mm"] - spec["wire_diameter_mm"]) / 2.
    turns = spec["active_turns"] + 2
    accepted = replay["samples"][:replay["contiguous_samples"]]
    selected = {row["trace_frame"]: row for row in
                [*accepted[::8], accepted[-1], max(accepted, key=lambda row: row["input_deg"]),
                 max(accepted, key=lambda row: row["carriage_mm"])]}
    def inspect(row):
        data.qpos[:] = row["qpos"]
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        posed = {key: part.transform(frames[owners[key]]) for key, part in parts.items()}
        hard = bounded(posed)
        hardware_hits = find_hits(list(posed.values()), keys, pairs)
        hum = frames["humerus"]
        rotation, offset = hum[:, :3], hum[:, 3]
        coil = sweep(spring_path(length + max(0., row["carriage_mm"]), turns, radius, mean), radius)
        coil = coil.transform(spring_frame).transform(hum)
        spring_hits = overlaps(coil, hard)
        targets = {**hard, **bounded({"series:return_spring": coil})}
        bodies, wires, findings = {}, {}, {}
        for label, (head, points, weights, vertices, vertex_weights, faces) in muscles.items():
            distal = frames[f"{head}_distal_pivot"]
            shift = (rotation.T @ (frames["FJ1477_carriage"][:, 3] - offset)
                     if label == "extensor" else np.zeros(3))
            local = (rotation.T @ distal[:, :3], rotation.T @ (distal[:, 3] - offset), shift)
            path = displacement_map(points, weights, *local) @ rotation.T + offset
            body_points = displacement_map(vertices, vertex_weights, *local) @ rotation.T + offset
            solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(body_points), faces))
            if solid.status() != mf.Error.NoError:
                raise ValueError(f"Invalid interpolated {label} shape at frame {row['trace_frame']}")
            bodies[label] = solid
            wires[label] = A.D.outer_envelope(path, .30)
            findings[label] = dict(body_intersections_mm3=overlaps(solid, targets),
                                   wire_intersections_mm3=overlaps(wires[label], targets))
        for label, wire in wires.items():
            findings[label]["wire_TPU_intersections_mm3"] = overlaps(wire, bounded(bodies))
        print(row["trace_frame"], len(hardware_hits), len(spring_hits), flush=True)
        return dict(trace_frame=row["trace_frame"], input_deg=row["input_deg"],
                    hardware_hits=hardware_hits, spring_hits_mm3=spring_hits,
                    muscles=findings, muscle_pair_overlap_mm3=A.I.intersection_volume(*bodies.values()))

    results = {index: inspect(row) for index, row in selected.items()}
    contact_seeds = [index for index, row in results.items()
                     if any(data["body_intersections_mm3"] for data in row["muscles"].values())]
    for row in accepted:
        index = row["trace_frame"]
        if index not in results and any(abs(index - seed) <= 32 for seed in contact_seeds):
            results[index] = inspect(row)
    save(folder / "whole_arm_elbow_flexible_audit.json", dict(
        replay=replay_ref, retention=retention, source_meshes=prior["mesh_sources"], revisions=revisions,
        hardware_checked_pair_count=len(pairs), samples=[results[index] for index in sorted(results)],
        generator=provenance(ROOT / "tools/audit_whole_arm_elbow_flexibles.py"),
        scope="Prescribed existing elbow replay, every eighth previously audited pose plus peaks, endpoint and contact neighborhoods. "
              "Indexed original lattices with existing local relief revisions; prescribed cap interpolation, "
              "not compliant equilibrium. Shared rigid arm and base used; shoulder TPU remains unposed."))


if __name__ == "__main__":
    main()
