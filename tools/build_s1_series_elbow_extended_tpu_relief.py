"""Protect the channels while relieving local TPU contacts on the extended elbow cycle."""

import argparse

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, export_lattice, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, rod, save
from refine_elbow_short_guide_relief import indexed_solid, inverse_map
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("muscle", choices=("flexor", "extensor"))
    parser.add_argument("--interface-guard", action="store_true")
    args = parser.parse_args()
    label = args.muscle
    if args.interface_guard and label != "flexor":
        parser.error("The interface guard applies to the flexor only")
    folder = DEST.parent / "integration_review/s1_series_elbow"
    audit_path = folder / "extended_elbow_muscle_geometry.json"
    audit = read(audit_path)
    extended = load_reference(audit["extended_drive"])
    budget = load_reference(extended["budget"])
    assembly = load_reference(budget["assembly"])
    elbow = load_reference(load_reference(audit["passive"])["assembly"])["muscles"][label]
    definition = load_reference(elbow["source_definition"])
    source = audit["mesh_sources"][label]
    original, vertices, triangles = indexed_solid(source)
    path = np.array(elbow["wire_path_source_mm"])[1:-1]
    axis = np.array(definition["basis_columns_display"])[:, 2]
    gradient = axis / ((path[-1] - path[0]) @ axis)
    head = elbow["anatomical_id"]
    first = path[0]
    trace = checked_trace(audit["trace"])
    model_path = ROOT / extended["model"]["file"]
    if provenance(model_path) != extended["model"]:
        raise ValueError("Extended drive model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    manifold = A.I._manifold()
    hits = {index: row["muscles"][label]["rigid_TPU_contacts"] for index, row in enumerate(audit["records"])
            if row["muscles"][label].get("rigid_TPU_contacts")}
    sources = {**assembly["sources"], **extended["shafts"]}
    blockers = {key: load_component(sources[key]) for key in sorted({k for row in hits.values() for k in row})}
    samples = {index for index, row in enumerate(audit["records"]) if row["TPU_sample"]}
    samples.update(index for seed in hits for index in range(max(0, seed - 7), min(len(trace), seed + 8)))

    def mesh_data(part):
        points, faces = mesh_geom.manifold_mesh_to_arrays(part.to_mesh64())
        return points, weigh(points, first, gradient), faces.astype(np.uint64)

    def solid_at(mesh, rotation, translation, shift):
        vertices, weights, faces = mesh
        points = displacement_map(vertices, weights, rotation, translation, shift)
        solid = manifold.Manifold(manifold.Mesh64(np.ascontiguousarray(points), faces))
        if solid.status() != manifold.Error.NoError:
            raise ValueError("Invalid displaced lattice")
        return solid

    def sample(index):
        data.qpos[:] = [float(trace[index][f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        transforms_world = transforms(model, data, neutral)
        humerus = transforms_world["humerus"]
        world_rotation, world_translation = humerus[:, :3], humerus[:, 3]

        def local(owner):
            world = transforms_world.get(owner, np.c_[np.eye(3), np.zeros(3)])
            return np.c_[world_rotation.T @ world[:, :3], world_rotation.T @ (world[:, 3] - world_translation)]

        distal = local(f"{head}_distal_pivot")
        shift = local("FJ1477_carriage")[:, 3] if label == "extensor" else np.zeros(3)
        obstacles = {key: part if assembly["owners"][key] == "humerus"
                     else part.transform(local(assembly["owners"][key])) for key, part in blockers.items()}
        return (distal[:, :3], distal[:, 3], shift), bounded(obstacles)

    states = {index: sample(index) for index in sorted(samples)}
    original_mesh = vertices, weigh(vertices, first, gradient), triangles
    fitting_hits = dict(hits)
    if label == "extensor" or args.interface_guard:
        for index in sorted(samples - hits.keys()):
            if index % 2 or not any(abs(index - seed) <= 7 for seed in hits):
                continue
            affine, obstacles = states[index]
            contact = overlaps(solid_at(original_mesh, *affine), obstacles)
            if contact:
                fitting_hits[index] = contact
    patches = []
    for index, row in fitting_hits.items():
        affine, obstacles = states[index]
        posed = solid_at(original_mesh, *affine)
        for key, volume in row.items():
            intersection = posed ^ obstacles[key][0]
            if abs(intersection.volume() - volume) > 1e-5:
                raise ValueError(f"Baseline contact did not reproduce: {key} at {index}")
            for component in intersection.decompose():
                if component.volume() <= 1e-6:
                    continue
                points, _ = mesh_geom.manifold_mesh_to_arrays(component.to_mesh64())
                source_points = inverse_map(points, *affine, first, gradient)
                patches.append(manifold.Manifold.hull_points(np.ascontiguousarray(source_points)))
    evidence = {}
    suffix = "interface_relief_TPU" if args.interface_guard else "extended_relief_TPU"
    patch_union = union(patches).as_original().simplify(1e-6)
    if label == "flexor":
        channel = A.D.outer_envelope(path, definition["tube_outer_radius_mm"] + .05)
    else:
        short = read(DEST.parent / "elbow_local_muscle_load/elbow_slider_short_guides.json")
        guides = load_reference(short["short_guide_design"])["actual_short_tube_guides"]
        if len(guides) != len(path):
            raise ValueError("Extensor guide count does not match the internal path")
        source_centres = np.array([row["centre_display_mm"] for row in guides])
        np.testing.assert_allclose(source_centres - source_centres[0], path - first, atol=1e-7, rtol=0.)
        half = np.array([row["axis_display"] for row in guides]) * np.array([row["length_mm"] for row in guides])[:, None] / 2.
        channel = union(rod(p - h, p + h, row["outer_radius_mm"] + .05)
                        for p, h, row in zip(path, half, guides, strict=True))
    proximal_stock = original.trim_by_plane(-axis, -float(path[0] @ axis + 8.))
    distal_stock = original.trim_by_plane(axis, float(path[-1] @ axis - 8.))
    if args.interface_guard:
        # Protect material at all ten physical cap, crimp and clamp solids,
        # including a guard extending beyond their actual contact surfaces.
        keys = sorted(key for key, owner in assembly["owners"].items() if owner == f"{head}_distal_pivot")
        if f"elbow:{head}_cap_rotor" not in keys or f"elbow:{head}_active_wire_crimp_barrel" not in keys:
            raise ValueError("Incomplete distal interface")
        pad = manifold.Manifold.cube([1.5] * 3, True)
        interface = union(load_component(sources[key]).minkowski_sum(pad) for key in keys)
        distal_stock = original ^ interface
        evidence["distal_interface_guard"] = dict(source_parts=keys, cube_dilation_half_width_mm=.75,
                                                   material_volume_mm3=float(distal_stock.volume()))
    end_stock = union([proximal_stock, distal_stock])
    channel_stock = original ^ channel
    protected = union([channel_stock, end_stock])
    necessary_guard_overlap = float((patch_union ^ protected).volume())
    trials = []
    selected = None
    path_weights = weigh(path, first, gradient)
    margins = (.10, .25, .45) if args.interface_guard else ((.10,) if label == "flexor" else (.45, .70, 1.20))
    for margin in margins:
        pad = manifold.Manifold.cube([2. * margin] * 3, True)
        removed = union(part.minkowski_sum(pad) for part in patches) - protected
        candidate = (original - removed).as_original().simplify(1e-6)
        components = sorted(candidate.decompose(), key=lambda part: part.volume(), reverse=True)
        pieces = [float(part.volume()) for part in components if part.volume() > 0.]
        detached_volume = sum(pieces[1:])
        if 0. < detached_volume <= .001 and len(components) == len(pieces):
            candidate = components[0]
        lost = float((protected - candidate).volume())
        row = dict(margin_mm=margin, removed_volume_mm3=float((original - candidate).volume()),
                   protected_loss_mm3=lost, positive_component_volumes_mm3=pieces,
                   detached_volume_mm3=detached_volume)
        trials.append(row)
        if not pieces or (len(pieces) > 1 and candidate is not components[0]) or lost > .001:
            continue
        mesh = mesh_data(candidate)
        residual = {}
        for index, (affine, obstacles) in states.items():
            intersection = overlaps(solid_at(mesh, *affine), obstacles)
            if intersection:
                residual[index] = intersection
        row.update(overlap_frame_count=len(residual), residual_overlaps_mm3=residual)
        print(label, margin, "residual frames", len(residual), flush=True)
        if residual:
            continue
        actual, meta = export_lattice(candidate, folder / f"{head}_{suffix}.obj",
                                      np.array(source["translation_source_mm"]))
        mesh = mesh_data(actual)
        for index, (affine, obstacles) in states.items():
            posed = solid_at(mesh, *affine)
            wire = A.D.outer_envelope(displacement_map(path, path_weights, *affine), .3)
            if overlaps(posed, obstacles) or (posed ^ wire).volume() > .001:
                raise ValueError(f"Exported relief failed readback at frame {index}")
        selected = dict(**row, mesh=meta)
        break
    diagnosis = {}
    if selected is None and "residual_overlaps_mm3" in trials[-1]:
        guard_meshes = {"channel_or_short_guides": mesh_data(channel_stock), "end_stock": mesh_data(end_stock)}
        old_shaft = bounded({"series:a_shaft": load_component(assembly["sources"]["series:a_shaft"])}) if label == "extensor" else None
        for index, remaining in trials[-1]["residual_overlaps_mm3"].items():
            affine, obstacles = states[index]
            targets = {key: obstacles[key] for key in remaining}
            posed_original = solid_at(original_mesh, *affine)
            diagnosis[index] = dict(
                original=overlaps(posed_original, targets),
                protected_stock={key: overlaps(solid_at(mesh, *affine), targets) for key, mesh in guard_meshes.items()})
            if old_shaft is not None:
                diagnosis[index]["original_with_short_shaft"] = overlaps(posed_original, old_shaft)
    save(folder / f"{head}_{suffix}.json", dict(
        generator=provenance(ROOT / "tools/build_s1_series_elbow_extended_tpu_relief.py"),
        input_audit=provenance(audit_path), muscle=label, original=source, **evidence,
        original_contact_frames=len(hits), training_frames=len(fitting_hits),
        additional_training_frame_indices=sorted(fitting_hits.keys() - hits.keys()),
        checked_frames=len(states), held_out_frames=len(states)-len(fitting_hits),
        targeted_obstacles=list(blockers), guard_overlap_of_inverse_mapped_hulls_mm3=necessary_guard_overlap,
        trials=trials, selected=selected, residual_diagnosis=diagnosis,
        scope="Projected overlap hulls cut only outside a nominal channel/guide and end-interface guard."
              " A detached fragment may be discarded only within the 0.001 mm3 overlap tolerance"
              " and after rechecking guard loss. Clearance uses kinematic TPU interpolation and nominal wire."))
    print(label, "selected", selected is not None, "guard overlap", necessary_guard_overlap, flush=True)


if __name__ == "__main__":
    main()
