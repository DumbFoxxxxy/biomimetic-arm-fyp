"""Pull sampled shaft contacts back into the undeformed short-guide mesh."""

import csv
from pathlib import Path

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, export_lattice, weigh
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, rod, save
from simulate_elbow_antagonistic_load import JOINTS
from study_elbow_series_package import union


def indexed_solid(source):
    path = ROOT / source["file"]
    if provenance(path)["sha256"] != source["sha256"]:
        raise RuntimeError(f"Mesh changed: {path.name}")
    vertices, faces = mesh_geom.load_obj(str(path))
    vertices += source["translation_source_mm"]
    faces = np.ascontiguousarray(faces, dtype=np.uint64)
    return A.I._manifold().Manifold(A.I._manifold().Mesh64(np.ascontiguousarray(vertices), faces)), vertices, faces


def inverse_map(target, rotation, translation, shift, origin, gradient):
    points = target - shift
    delta = rotation - np.eye(3)
    for _ in range(16):
        station = weigh(points, origin, gradient)
        error = displacement_map(points, station, rotation, translation, shift) - target
        if np.linalg.norm(error, axis=1).max() < 1e-7:
            return points
        active = (station > 0.) & (station < 1.)
        difference = points @ delta.T + translation - shift
        jacobian = np.eye(3) + station[:, None, None] * delta
        jacobian[active] += difference[active, :, None] * gradient[None, None, :]
        points -= np.linalg.solve(jacobian, error[:, :, None])[:, :, 0]
    raise RuntimeError("Could not pull shaft overlap back to the source TPU frame")


def main():
    folder = DEST.parent / "elbow_local_muscle_load"
    parent_file = folder / "elbow_slider_short_guides.json"
    parent = read(parent_file)
    mount = read(ROOT / parent["spring_mount"]["file"])
    response = read(ROOT / parent["response"]["file"])
    definition_file = ROOT / parent["short_guide_design"]["file"]
    if provenance(definition_file) != parent["short_guide_design"]:
        raise RuntimeError("Short guide definition changed")
    definition = read(definition_file)
    definition_offset = np.array([-140., -56., 1100.])
    centres = np.asarray([row["centre_display_mm"] for row in definition["actual_short_tube_guides"]]) + definition_offset
    muscle_axis = np.asarray(definition["basis_columns_display"])[:, 2]
    origin = centres[0]
    gradient = muscle_axis / ((centres[-1] - origin) @ muscle_axis)
    model_path = ROOT / response["local_model"]["file"]
    if provenance(model_path) != response["local_model"]:
        raise RuntimeError("Free slider model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    slider = model.joint("FJ1477_series_travel").qposadr.item()
    distal = data.body("FJ1477_distal_pivot")
    neutral_distal = distal.xpos.copy() * 1000.
    neutral_slide = data.body("FJ1477_carriage").xpos.copy() * 1000.
    shaft = load_component(mount["retention"]["a"]["sources"]["shaft"])
    base_row = next(row for row in parent["shaft_relief_options"] if row["cutter_radius_mm"] == 5.5)
    base, vertices, faces = indexed_solid(base_row["candidate_source"])
    weights = weigh(vertices, origin, gradient)
    mf = A.I._manifold()
    protected = union(rod(centre - np.asarray(row["axis_display"]) * row["length_mm"] / 2.,
                          centre + np.asarray(row["axis_display"]) * row["length_mm"] / 2.,
                          row["outer_radius_mm"] + .25)
                      for centre, row in zip(centres, definition["actual_short_tube_guides"], strict=True)) ^ base
    pose_samples, intermediate, heldout, local_cuts = [], [], [], []
    for label, case in parent["cases"].items():
        trace_file = ROOT / case["trace"]["file"]
        if provenance(trace_file) != case["trace"]:
            raise RuntimeError("Trace changed after original collision check")
        with trace_file.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        samples = {row["source_pose_index"]: row for row in case["mesh_samples"]}
        angles = np.asarray([float(row["elbow_deg"]) for row in rows])
        peak = int(angles.argmax())
        validation = set()
        for angle in (2.5, 7.5, 12.5, 17.5, 22.5, 27.5):
            validation.add(int(abs(angles[:peak + 1] - angle).argmin()))
            validation.add(peak + int(abs(angles[peak:] - angle).argmin()))
        independent = set()
        for angle in np.arange(1.25, 30., 2.5):
            independent.add(int(abs(angles[:peak + 1] - angle).argmin()))
            independent.add(peak + int(abs(angles[peak:] - angle).argmin()))
        independent -= set(samples) | validation
        for frame in sorted(set(samples) | validation | independent):
            sample = samples.get(frame)
            row = rows[frame]
            data.qpos[addresses] = np.deg2rad([float(row[f"{name}_deg"]) for name in JOINTS])
            data.qpos[slider] = float(row["carriage_mm"]) * .001
            mujoco.mj_kinematics(model, data)
            rotation = distal.xmat.reshape(3, 3).copy()
            translation = distal.xpos * 1000. - rotation @ neutral_distal
            shift = data.body("FJ1477_carriage").xpos * 1000. - neutral_slide
            if frame in independent:
                heldout.append((label, frame, rotation, translation, shift))
                continue
            if sample is None:
                intermediate.append((label, frame, rotation, translation, shift))
            else:
                pose_samples.append((label, frame, rotation, translation, shift))
                if "rail_a_shaft" not in sample["relief_probes"]["5.5"]["tpu_rigid_overlap_mm3"]:
                    continue
            posed = displacement_map(vertices, weights, rotation, translation, shift)
            intersection = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), faces)) ^ shaft
            for part in intersection.decompose():
                if part.volume() <= .001:
                    continue
                points, _ = mesh_geom.manifold_mesh_to_arrays(part.to_mesh64())
                local_points = inverse_map(points, rotation, translation, shift, origin, gradient)
                local_cuts.append(mf.Manifold.hull_points(np.ascontiguousarray(local_points)))
    if not local_cuts:
        raise RuntimeError("No residual contacts found to guide local correction")
    options = []
    selected = None
    for allowance in (.2, .4, .7, 1.):
        pad = mf.Manifold.cube([2. * allowance] * 3, True)
        trial = (base - union(part.minkowski_sum(pad) for part in local_cuts)).as_original().simplify(1e-6)
        parts = sorted((float(part.volume()) for part in trial.decompose() if part.volume() > 1e-5), reverse=True)
        row = dict(allowance_mm=allowance, added_cut_mm3=float((base - trial).volume()),
                   protected_guide_removed_mm3=float((protected - trial).volume()), component_volumes_mm3=parts)
        options.append(row)
        if len(parts) != 1 or row["protected_guide_removed_mm3"] > .001:
            continue
        try:
            exported, source = export_lattice(trial, folder / f"elbow_short_guides_mapped_relief_{allowance:g}.obj", origin)
        except (AssertionError, ValueError) as error:
            row["export_error"] = str(error) or type(error).__name__
            continue
        points, triangles = mesh_geom.manifold_mesh_to_arrays(exported.to_mesh64())
        new_weights = weigh(points, origin, gradient)
        triangles = np.ascontiguousarray(triangles, dtype=np.uint64)
        contacts = []
        for label, frame, rotation, translation, shift in pose_samples + intermediate + heldout:
            posed = displacement_map(points, new_weights, rotation, translation, shift)
            body = mf.Manifold(mf.Mesh64(np.ascontiguousarray(posed), triangles))
            amount = float((body ^ shaft).volume())
            if amount > .001 or body.status() != mf.Error.NoError:
                contacts.append(dict(case=label, frame=frame, shaft_overlap_mm3=amount,
                                     geometry_status=str(body.status())))
        row.update(source=source, rechecked_poses=len(pose_samples),
                   additional_cut_fit_poses=len(intermediate), held_out_poses=len(heldout),
                   remaining_shaft_contacts=contacts)
        if not contacts:
            selected = allowance
            break
    save(folder / "elbow_short_guide_mapped_relief.json", dict(
        generator=provenance(Path(__file__).resolve()), parent=provenance(parent_file),
        source_lattice=base_row["candidate_source"], back_mapped_local_cut_count=len(local_cuts),
        candidates=options, selected_allowance_mm=selected,
        scope="Extensor rail A only. Pullback uses an imposed map, and sampled contact changes do not "
        "establish a physical TPU deformation, strength, or material capability. Rigid cap trim "
        "and all other previous components retain their separate verification scope."))
    for row in options:
        print(row, flush=True)


if __name__ == "__main__":
    main()
