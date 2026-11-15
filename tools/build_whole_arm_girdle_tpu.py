"""Build channel-centered girdle lattices and screen their nominal motion envelope."""

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, mesh_geom, rod, save
from project_runtime import load_mujoco_model
from study_elbow_series_package import bounded, overlaps


def wire_frame(start, end, reference_side):
    axis = end - start
    length = np.linalg.norm(axis)
    axis /= length
    side = reference_side - axis * (axis @ reference_side)
    side /= np.linalg.norm(side)
    return np.column_stack((side, np.cross(axis, side), axis)), length


def lattice(length, widths, count):
    """Use staggered rings to make a tapered crossed-strut lattice around a sleeve."""
    mf = A.I._manifold()
    end = length - 2.
    heights = np.linspace(2., end, count + 1)
    rings = []
    pitch = 2 * np.pi / 8
    for i, z in enumerate(heights):
        t = i / count
        radius = 1.85 + np.asarray(widths) * np.sin(np.pi * t)
        angles = np.arange(8) * pitch + (i % 2) * pitch / 2
        rings.append(np.column_stack((radius[0] * np.cos(angles),
                                      radius[1] * np.sin(angles), np.full(8, z))))
    pieces = [rod([0., 0., 2.], [0., 0., end], 1.35, 24)]
    for z in (2., end - 1.1):
        pieces.append(rod([0., 0., z], [0., 0., z + 1.1], 2.15, 24))
    for i, ring in enumerate(rings):
        for j, point in enumerate(ring):
            # Radial links attach the taper to the line guide throughout its length.
            pieces.append(rod([0., 0., point[2]], point, .45, 12))
            if i < count:
                for k in (j, (j + (-1 if i % 2 == 0 else 1)) % len(ring)):
                    pieces.append(rod(point, rings[i + 1][k], .45, 12))
            pieces.append(mf.Manifold.sphere(.48, 12).translate(point))
    material = mf.Manifold.batch_boolean(pieces, mf.OpType.Add)
    bore = rod([0., 0., 0.], [0., 0., length], .8, 24)
    return (material - bore).as_original().simplify(1e-6)


def relieve_neutral(body, obstacles, hits):
    # Keep the original lattice and sleeve around the wire wherever room exists.
    # Only the opposing rigid surfaces are used as cutters.
    if not hits:
        return body, []
    trimmed = body
    for key in hits:
        trimmed -= obstacles[key][0]
    trimmed = trimmed.as_original().simplify(1e-6)
    positive = [part for part in trimmed.decompose() if part.volume() > 1e-8]
    if len(positive) != 1 or trimmed.volume() < .90 * body.volume():
        return body, []
    return trimmed, sorted(hits)


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    scene_file = folder / "whole_arm_shared_scene_check.json"
    scene = read(scene_file)
    motion_file = folder / "whole_arm_shared_motion.json"
    motion = read(motion_file)
    if motion["model"] != scene["model"]:
        raise ValueError("The motion and assembly refer to different scenes")
    transfer = read(ROOT / scene["transfer"]["file"])
    seat_file = folder / "girdle_seat_pair.json"
    seats = read(seat_file)["seats"]
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"]}
    rigid = {key: load_component(source) for key, source in sources.items()}
    owners = transfer["owner_by_source"]
    origins = {key: np.asarray(scene["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    model = load_mujoco_model(ROOT / scene["model"]["file"])
    data = mujoco.MjData(model)
    names = ("scapular_SC1_pos", "K7_SC1_pos")
    paths = {}
    neutral = bounded(rigid)
    original_file = DEST.parent / "integration_review/s1_series_elbow/assembly_FJ1504_local_interface_bidirectional.json"
    available_tpu, unavailable_tpu = {}, {}
    for key, source in read(original_file)["flexible_geometry_sources"].items():
        try:
            available_tpu[key] = load_component(source)
        except ValueError as error:
            unavailable_tpu[key] = str(error)
    original_tpu = bounded(available_tpu)
    definitions = {}
    for head, key in (("FJ1459", "FJ1459_g3"), ("FJ1536", "FJ1536_g2")):
        ports = [seats[f"{key}_{end}"] for end in ("origin", "insertion")]
        start, end = map(np.asarray, (port["port_source_mm"] for port in ports))
        anatomy = DEST.parents[1] / "输入/muscle" / f"{head}.obj"
        vertices, _ = mesh_geom.load_obj(str(anatomy))
        axis = (end - start) / np.linalg.norm(end - start)
        station = (vertices - start) @ axis
        section = vertices[(station > 0.) & (station < np.linalg.norm(end - start))]
        if len(section) < 8:
            raise ValueError(f"Not enough reference anatomy at {key}")
        lateral = section - start - np.outer((section - start) @ axis, axis)
        _, _, vectors = np.linalg.svd(lateral - lateral.mean(axis=0), full_matrices=False)
        frame, length = wire_frame(start, end, vectors[0])
        widths = np.clip(np.percentile(abs(lateral @ frame[:, :2]), 60, axis=0) * .22, 1., 5.)
        count = int(np.clip(np.ceil(length / 7.), 4, 16))
        trials = []
        candidates = []
        for scale in (1., .65, .35):
            material = lattice(length, widths * scale, count)
            world = material.transform(np.column_stack((frame, start)))
            hits = overlaps(world, neutral)
            trials.append(dict(scale=scale, rigid_overlap_mm3=hits))
            candidates.append((sum(hits.values()), -scale, world))
            if not hits:
                break
        _, chosen, world = min(candidates, key=lambda item: item[:2])
        original_volume = world.volume()
        world, trimmed_against = relieve_neutral(world, neutral, overlaps(world, neutral))
        body, meta = export(world, folder / f"{key}_TPU_lattice.obj", start)
        definitions[key] = (body, ports, frame, start, length, dict(
            source=meta, anatomy=provenance(anatomy), ring_count=count + 1,
            chosen_width_scale=-chosen, taper_increments_mm=(widths * -chosen).tolist(),
            channel_diameter_mm=1.6, terminal_free_length_mm=2., neutral_port_distance_mm=length.item(),
            neutral_trials=trials, neutral_relief_against=trimmed_against,
            neutral_removed_volume_mm3=float(original_volume - world.volume()),
            neutral_reloaded_overlap_mm3=overlaps(body, neutral),
            neutral_existing_tpu_overlap_mm3=overlaps(body, original_tpu)))
        print(key, "neutral rigid overlap trials", trials, flush=True)

    for name in names:
        samples = []
        for row in motion["paths"][name]["accepted"]:
            data.qpos[:] = row["qpos"]
            mujoco.mj_forward(model, data)
            matrices = transforms(model, data, origins)
            obstacles = bounded({key: part.transform(matrices[owners[key]]) for key, part in rigid.items()})
            placed = {}
            wires = {}
            results = {}
            for key, (solid, ports, rest_frame, rest_start, length, _) in definitions.items():
                points = [matrices[p["owner"]] @ np.r_[p["port_source_mm"], 1.] for p in ports]
                frame, span = wire_frame(*points, rest_frame[:, 0])
                stretch = (span - 4.) / (length - 4.)
                mapping = frame @ np.diag([1., 1., stretch]) @ rest_frame.T
                offset = points[0] + frame[:, 2] * (2. - 2. * stretch) - mapping @ rest_start
                placed[key] = solid.transform(np.column_stack((mapping, offset)))
                wires[key] = rod(*points, .30)
                hit = overlaps(placed[key], obstacles)
                results[key] = dict(port_distance_mm=span.item(), axial_scale=stretch.item(),
                                    rigid_overlaps_mm3=hit,
                                    own_wire_overlap_mm3=A.I.intersection_volume(placed[key], wires[key]))
            for key, part in placed.items():
                results[key]["neighbor_wire_overlap_mm3"] = sum(
                    A.I.intersection_volume(part, wire) for other, wire in wires.items() if other != key)
            between = A.I.intersection_volume(*placed.values())
            samples.append(dict(input_deg=row["input_deg"], muscles=results, pair_overlap_mm3=between))
        paths[name] = samples
        print(name, len(samples), "nominal TPU poses", flush=True)
    save(folder / "girdle_tpu_geometry.json", dict(
        scene=provenance(scene_file), motion=provenance(motion_file), seats=provenance(seat_file),
        existing_tpu_reference=provenance(original_file),
        existing_tpu_solids_checked=list(original_tpu), existing_tpu_solids_unavailable=unavailable_tpu,
        generator=provenance(ROOT / "tools/build_whole_arm_girdle_tpu.py"),
        muscles={key: entry[-1] for key, entry in definitions.items()}, paths=paths,
        scope="Neutral connected lattice geometry plus axial-only, fixed-cross-section affine preview. "
              "Overlaps flag where compliant contact and the actual channel support need testing. "
              "One sampling pass covers geometry shared by outward and return poses. No material "
              "fit, force balance, compressive instability or physical performance verified."))


if __name__ == "__main__":
    main()
