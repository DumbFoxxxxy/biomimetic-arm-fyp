"""Screen smaller flexor keeper nuts at the existing through-bolt and washer."""

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, save
from refine_elbow_short_guide_relief import indexed_solid
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    parent_path = folder / "FJ1486_extended_relief_TPU.json"
    parent = read(parent_path)
    audit = load_reference(parent["input_audit"])
    extended = load_reference(audit["extended_drive"])
    assembly = load_reference(load_reference(extended["budget"])["assembly"])
    elbow = load_reference(load_reference(audit["passive"])["assembly"])["muscles"]["flexor"]
    definition = load_reference(elbow["source_definition"])
    nut_source = assembly["sources"]["elbow:FJ1486_keeper_nut"]
    design_path = (ROOT / nut_source["file"]).with_name("design.json")
    feature = next(row for row in read(design_path)["features"] if row["muscle"] == "FJ1486")
    axis = np.array(feature["keeper_bolt_axis_display"])
    # The OBJ has its own local origin. Recover the source-plane centre from
    # its two section centroids, and verify the historical 1.6 mm axial span.
    nut = load_component(nut_source)
    mf = A.I._manifold()
    points = np.asarray(nut.to_mesh64().vert_properties)[:, :3]
    z = points @ axis
    rear, seat = float(z.min()), float(z.max())
    np.testing.assert_allclose(seat - rear, 1.6, atol=.002, rtol=0.)
    origin = np.unique(points[np.isclose(z, seat, atol=.002, rtol=0.)], axis=0).mean(axis=0)
    ring_radii = np.linalg.norm(points - origin - ((points - origin) @ axis)[:, None] * axis, axis=1)
    if abs(ring_radii.max() - 4. / np.sqrt(3.)) > .005:
        raise ValueError("Nut section centre does not reproduce the original hex")
    u = A.D.unit(np.cross(axis, [1., 0., 0.]))
    frame = np.c_[u, np.cross(axis, u), axis, origin]
    washer_source = assembly["sources"]["elbow:FJ1486_keeper_back_washer"]
    washer = load_component(washer_source)
    seat_witness = washer.translate(-.02 * axis)
    original_seat = (nut ^ seat_witness).volume()
    if original_seat < .05:
        raise ValueError("Back washer does not meet the nut seating face")

    source = parent["original"]
    _, vertices, faces = indexed_solid(source)
    path = np.array(elbow["wire_path_source_mm"])[1:-1]
    muscle_axis = np.array(definition["basis_columns_display"])[:, 2]
    weights = weigh(vertices, path[0], muscle_axis / ((path[-1] - path[0]) @ muscle_axis))
    model_path = ROOT / extended["model"]["file"]
    if provenance(model_path) != extended["model"]:
        raise ValueError("Drive model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    trace = checked_trace(audit["trace"])
    fitted = {int(index) for index in parent["residual_diagnosis"]}
    samples = {index for index, row in enumerate(audit["records"]) if row["TPU_sample"]}
    samples |= {index for seed in fitted for index in range(max(0, seed - 7), min(len(trace), seed + 8))}
    candidates = [(1.6, 2.1), (1.6, 2.0), (1.6, 1.9), (1.6, 1.8), (1.4, 2.0), (1.4, 1.9), (1.2, 1.9)]
    shapes = [nut ^ mf.Manifold.cylinder(height, radius, radius, 96).translate([0., 0., -height]).transform(frame)
              for height, radius in candidates]
    rows = [dict(height_mm=height, max_radius_mm=radius, removed_mm3=float((nut - part).volume()),
                 relative_seat_area=float((part ^ seat_witness).volume() / original_seat), overlaps={})
            for (height, radius), part in zip(candidates, shapes, strict=True)]
    original_hits = {}
    contacts = {}
    for index in sorted(samples):
        data.qpos[:] = [float(trace[index][f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        pose = transforms(model, data, neutral)
        reference = pose["forearm"]
        rotation, translation = reference[:, :3], reference[:, 3]
        h = pose["humerus"]
        d = pose["FJ1486_distal_pivot"]
        points = displacement_map(vertices, weights, h[:, :3].T @ d[:, :3], h[:, :3].T @ (d[:, 3] - h[:, 3]), np.zeros(3))
        points = (points @ h[:, :3].T + h[:, 3] - translation) @ rotation
        muscle = mf.Manifold(mf.Mesh64(np.ascontiguousarray(points), faces))
        if muscle.status() != mf.Error.NoError:
            raise ValueError(f"Invalid flexor interpolation at {index}")
        contact = (muscle ^ nut).as_original().simplify(1e-6)
        amount = float(contact.volume())
        if amount > .001:
            original_hits[index] = amount
            contacts[index] = contact
            for row, part in zip(rows, shapes, strict=True):
                volume = float((contact ^ part).volume())
                if volume > .001:
                    row["overlaps"][index] = volume
    core = nut ^ mf.Manifold.cylinder(1.6, 1.50, 1.50, 96).translate([0., 0., -1.6]).transform(frame)
    contact_hulls = [part.hull() for index, contact in contacts.items() if index % 2 == 0
                     for part in contact.decompose() if part.volume() > 1e-6]
    for allowance in (.06, .12, .20):
        pad = mf.Manifold.cube([2. * allowance] * 3, True)
        cutter = mf.Manifold.batch_boolean([patch.minkowski_sum(pad) for patch in contact_hulls], mf.OpType.Add)
        part = (nut - (cutter - core)).as_original().simplify(1e-6)
        if sum(piece.volume() > 1e-6 for piece in part.decompose()) != 1 or (core - part).volume() > .001:
            continue
        overlaps = {index: float((contact ^ part).volume()) for index, contact in contacts.items()}
        shapes.append(part)
        rows.append(dict(height_mm=1.6, local_allowance_mm=allowance, core_radius_mm=1.50,
                         removed_mm3=float((nut - part).volume()),
                         relative_seat_area=float((part ^ seat_witness).volume() / original_seat),
                         overlaps={index: volume for index, volume in overlaps.items() if volume > .001}))
    selected = None
    clear = [i for i, row in enumerate(rows) if not row["overlaps"] and row["relative_seat_area"] >= .70]
    if clear:
        index = max(clear, key=lambda i: (rows[i]["height_mm"], rows[i]["relative_seat_area"]))
        actual, mesh = export(shapes[index], folder / "FJ1486_low_profile_keeper_nut.obj", origin)
        if ((actual - nut).volume() > .001 or (actual ^ seat_witness).volume() < .70 * original_seat
                or any((actual ^ contact).volume() > .001 for contact in contacts.values())):
            raise ValueError("Readback added material or lost the required bearing seat")
        selected = dict(**rows[index], mesh=mesh)
    save(folder / "FJ1486_low_profile_keeper.json", dict(
        generator=provenance(ROOT / "tools/build_s1_series_elbow_flexor_keeper.py"),
        parent=provenance(parent_path), historical_design=provenance(design_path), original_nut=nut_source,
        original_washer=washer_source, original_flexor=source, original_overlap_frames=original_hits,
        tested_frames=len(samples), local_fit_frame_indices=[i for i in contacts if i % 2 == 0],
        trials=rows, selected=selected,
        scope="Subset of original nut material on existing bore and washer. Seat area measured with a 0.02 mm translated washer."
              " 70% area screen is provisional; threads, torque, preload, tools and retention load are unverified."))
    print("Nut candidate", selected, "baseline hits", len(original_hits),
          [(len(row["overlaps"]), row["relative_seat_area"]) for row in rows], flush=True)


if __name__ == "__main__":
    main()
