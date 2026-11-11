"""Pull the timed flexor–humerus contact into the source TPU and test local relief."""

import numpy as np
import mujoco

from audit_elbow_slider_short_guides import displacement_map, export_lattice, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from refine_elbow_short_guide_relief import indexed_solid, inverse_map
from study_elbow_series_package import union
from study_s1_local_force_geometry import load_reference


def mesh_data(part, first, gradient):
    points, faces = mesh_geom.manifold_mesh_to_arrays(part.to_mesh64())
    return points, weigh(points, first, gradient), faces.astype(np.uint64)


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    audit_file = folder / "timed_elbow_muscle_geometry.json"
    audit = read(audit_file)
    timing = load_reference(audit["timing"])
    baseline = load_reference(timing["baseline"])
    assembly = load_reference(baseline["geometry_assembly"])
    passive = load_reference(audit["passive"])
    flexor = load_reference(passive["assembly"])["muscles"]["flexor"]
    definition = load_reference(flexor["source_definition"])
    source = audit["mesh_sources"]["flexor"]
    original, vertices, faces = indexed_solid(source)
    humerus = load_component(assembly["sources"]["bone:humerus"])
    path = np.array(flexor["wire_path_source_mm"])[1:-1]
    axis = np.array(definition["basis_columns_display"])[:, 2]
    length = float((path[-1] - path[0]) @ axis)
    gradient = axis / length
    weights = weigh(vertices, path[0], gradient)
    path_weights = weigh(path, path[0], gradient)
    model_file = ROOT / baseline["model"]["file"]
    if provenance(model_file) != baseline["model"]:
        raise ValueError("Source model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    trace = checked_trace(audit["trace"])
    manifold = A.I._manifold()

    def pose(index):
        row = trace[index]
        data.qpos[:] = [float(row[f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        matrix = transforms(model, data, neutral)
        bone, end = matrix["humerus"], matrix["FJ1486_distal_pivot"]
        return bone[:, :3].T @ end[:, :3], bone[:, :3].T @ (end[:, 3] - bone[:, 3])

    def deform(mesh, rotation, translation):
        points, stations, triangles = mesh
        posed = displacement_map(points, stations, rotation, translation, np.zeros(3))
        part = manifold.Manifold(manifold.Mesh64(np.ascontiguousarray(posed), triangles))
        if part.status() != manifold.Error.NoError:
            raise ValueError("Invalid displaced flexor mesh")
        return part

    contacts = {index: row["muscles"]["flexor"].get("rigid_TPU_contacts", {}).get("bone:humerus", 0.)
                for index, row in enumerate(audit["records"])}
    training = {index for index, overlap in contacts.items() if overlap > .001}
    validation = {index for index, row in enumerate(audit["records"])
                  if row["TPU_sample"] or 1.20 <= row["time_s"] <= 1.50}
    poses = {index: pose(index) for index in validation}
    patches, baseline_hits = [], {}
    for index, (rotation, translation) in poses.items():
        overlap = deform((vertices, weights, faces), rotation, translation) ^ humerus
        volume = float(overlap.volume())
        if volume > .001:
            baseline_hits[index] = volume
        if index not in training:
            continue
        if abs(volume - contacts[index]) > 1e-5:
            raise ValueError("Mapped contact is inconsistent with the prior audit")
        for piece in overlap.decompose():
            if piece.volume() <= 1e-6:
                continue
            points, _ = mesh_geom.manifold_mesh_to_arrays(piece.to_mesh64())
            mapped = inverse_map(points, rotation, translation, np.zeros(3), path[0], gradient)
            patches.append(manifold.Manifold.hull_points(np.ascontiguousarray(mapped)))
    if not patches:
        raise ValueError("No source contact patches to guide a relief")
    tube = A.D.outer_envelope(path, definition["tube_outer_radius_mm"] + .05)
    end_stock = union([original.trim_by_plane(-axis, -float(path[0] @ axis + 8.)),
                       original.trim_by_plane(axis, float(path[-1] @ axis - 8.))])
    protected = original ^ union([tube, end_stock])
    trials, selected = [], None
    for margin in (.10, .20, .35, .50):
        pad = manifold.Manifold.cube([2. * margin] * 3, True)
        cutter = union(part.minkowski_sum(pad) for part in patches) - protected
        candidate = (original - cutter).as_original().simplify(1e-6)
        protected_loss = float((protected - candidate).volume())
        pieces = [float(piece.volume()) for piece in candidate.decompose() if piece.volume() > 1e-5]
        trial = dict(margin_mm=margin, removed_mm3=float((original - candidate).volume()),
                     protected_loss_mm3=protected_loss, positive_pieces_mm3=pieces)
        trials.append(trial)
        if protected_loss > .001 or len(pieces) != 1:
            continue
        points = mesh_data(candidate, path[0], gradient)
        overlaps = {}
        for index, (rotation, translation) in poses.items():
            volume = float((deform(points, rotation, translation) ^ humerus).volume())
            if volume > .001:
                overlaps[index] = volume
        trial["humerus_overlap_frames"] = len(overlaps)
        trial["max_remaining_overlap_mm3"] = max(overlaps.values(), default=0.)
        print(margin, trial["removed_mm3"], len(overlaps), trial["max_remaining_overlap_mm3"], flush=True)
        if overlaps:
            continue
        actual, output = export_lattice(candidate, folder / "FJ1486_return_relief_TPU.obj",
                                        np.array(source["translation_source_mm"]))
        actual_mesh = mesh_data(actual, path[0], gradient)
        for index, (rotation, translation) in poses.items():
            posed = deform(actual_mesh, rotation, translation)
            wire = A.D.outer_envelope(displacement_map(path, path_weights, rotation, translation, np.zeros(3)), .3)
            if max(float((posed ^ humerus).volume()), float((posed ^ wire).volume())) > .001:
                raise ValueError(f"Readback failed at frame {index}")
        selected = dict(**trial, mesh=output)
        break
    save(folder / "FJ1486_return_relief.json", dict(
        generator=provenance(ROOT / "tools/build_s1_series_elbow_flexor_relief.py"),
        reference_audit=provenance(audit_file), original=source,
        fitting_frame_count=len(training), checked_frame_count=len(poses),
        additional_frames_outside_fit=len(set(poses) - training),
        baseline_contact_count=len(baseline_hits), baseline_peak_overlap_mm3=max(baseline_hits.values()),
        protected_channel_radius_mm=definition["tube_outer_radius_mm"] + .05,
        protected_end_span_mm=8., trials=trials, selected=selected,
        scope="The contact is mapped in the humerus frame using the prior geometric interpolation."
              " This is a subtractive TPU candidate. The same recorded pose sequence is rechecked."
              " The local tube wall and end stock are guarded. Stiffness and contact mechanics are not inferred."))


if __name__ == "__main__":
    main()
