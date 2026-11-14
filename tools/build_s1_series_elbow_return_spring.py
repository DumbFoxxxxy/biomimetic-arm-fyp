"""Generate a continuous coil and captive end eyes on the existing spring pins."""

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def bezier(start, end, start_tangent, end_tangent):
    controls = np.array([start, start + 1.1 * start_tangent, end - 1.1 * end_tangent, end])
    t = np.linspace(0., 1., 28)[:, None]
    return (1.-t)**3 * controls[0] + 3.*t*(1.-t)**2 * controls[1] + 3.*t**2*(1.-t) * controls[2] + t**3 * controls[3]


def spring_path(span, turns, wire_radius, mean_radius):
    eye_radius = .9 + .15 + wire_radius
    angle = np.linspace(-np.pi/6., -11.*np.pi/6., 61)
    eye = np.c_[np.zeros(len(angle)), eye_radius*np.sin(angle), eye_radius*np.cos(angle)]
    eye_tangent = np.array([0., -np.sqrt(3.)/2., .5])
    # Two transition turns are included in the coil model. Space left at
    # each end accommodates the pin eye and the curved transition.
    coil_start = 4.45
    coil_span = span - 2.*coil_start
    phase = np.linspace(0., 2.*np.pi*turns, turns*64+1)
    coil = np.c_[-mean_radius*np.sin(phase), mean_radius*np.cos(phase), coil_start + coil_span*phase/phase[-1]]
    tangent = A.D.unit(np.array([-mean_radius, 0., coil_span/phase[-1]]))
    lead = bezier(eye[-1], coil[0], eye_tangent, tangent)
    tail = bezier(eye[-1], np.array([0., mean_radius, coil_start]), eye_tangent, tangent * [-1., 1., 1.])
    far_end = np.vstack([eye, tail[1:-1]])
    far_end[:, 2] = span - far_end[:, 2]
    return np.vstack([eye, lead[1:-1], coil, far_end[::-1]])


def sweep(points, radius):
    tangents = np.gradient(points, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    # Carry the previous transverse direction along the curved centreline.
    reference = np.eye(3)[np.argmin(abs(tangents[0]))]
    normals = []
    for tangent in tangents:
        reference = A.D.unit(reference - (reference @ tangent)*tangent)
        normals.append(reference)
    normal = np.array(normals)
    binormal = np.cross(tangents, normal)
    sides = 16
    theta = np.arange(sides)*2.*np.pi/sides
    rings = points[:, None] + radius*(np.cos(theta)[None, :, None]*normal[:, None]
                                    + np.sin(theta)[None, :, None]*binormal[:, None])
    vertices = np.vstack([rings.reshape(-1, 3), points[[0, -1]]])
    a = np.arange((len(points)-1)*sides).reshape(-1, sides)
    b = np.roll(a, -1, axis=1)
    side_a = np.stack([a, b, b+sides], axis=-1).reshape(-1, 3)
    side_b = np.stack([a, b+sides, a+sides], axis=-1).reshape(-1, 3)
    circle = np.arange(sides)
    cap_a = np.column_stack([np.full(sides, len(vertices)-2), np.roll(circle, -1), circle])
    end = circle + (len(points)-1)*sides
    cap_b = np.column_stack([np.full(sides, len(vertices)-1), end, np.roll(end, -1)])
    mf = A.I._manifold()
    solid = mf.Manifold(mf.Mesh64(vertices, np.vstack([side_a, side_b, cap_a, cap_b]).astype(np.uint64)))
    if solid.status() != mf.Error.NoError or solid.volume() <= 0.:
        raise ValueError("Invalid swept spring")
    return solid


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    size_path = folder / "return_spring_size_screen.json"
    size = read(size_path)
    audit = load_reference(size["geometry"])
    extended = load_reference(audit["extended_drive"])
    assembly = load_reference(load_reference(extended["budget"])["assembly"])
    mount = load_reference(size["spring_mount"])["spring_mount"]
    design = size["rates_within_5_percent_at_central_modulus"][0]
    radius = .5 * design["wire_diameter_mm"]
    mean = .5 * (design["coil_OD_mm"] - design["wire_diameter_mm"])
    turns = design["active_turns"] + 2
    start = np.array(mount["fixed_eye_source_mm"])
    end = np.array(mount["moving_eye_neutral_source_mm"])
    axis = A.D.unit(end-start)
    pin_axis = np.array(mount["pin_axis_source"])
    basis = np.c_[pin_axis, np.cross(axis, pin_axis), axis, start]
    rest = float(np.linalg.norm(end-start))
    stroke = size["full_rail_stroke_mm"]
    sources = {**assembly["sources"], **extended["shafts"]}
    parts = {key: load_component(value) for key, value in sources.items()}
    owners = assembly["owners"]
    mf = A.I._manifold()

    def geometry(travel):
        path = spring_path(rest+travel, turns, radius, mean)
        return sweep(path, radius).transform(basis)

    def envelope(travel):
        # The wire around each end is bounded by the first/last transitions.
        # Include the hook behind the pin as well as the continuous coil.
        path = spring_path(rest+travel, turns, radius, mean)
        first, last = 87, len(path)-87
        tail_points = (path[:first] @ basis[:, :3].T + start, path[last:] @ basis[:, :3].T + start)
        block = lambda pts: mf.Manifold.cube((np.ptp(pts, axis=0)+2.*radius).tolist()).translate((pts.min(axis=0)-radius).tolist())
        coil = rod(start + (4.45-radius)*axis, end + (travel-4.45+radius)*axis, mean+radius+.015, 64)
        return union([coil, *(block(pts) for pts in tail_points)])

    rigid_mount = {key: value for key, value in parts.items() if owners[key] in ("humerus", "FJ1477_carriage")}
    export_rows = {}
    for label, travel in (("neutral", 0.), ("full_stroke", stroke)):
        part = geometry(travel)
        actual, mesh = export(part, folder / f"return_spring_{label}.obj", start)
        if (actual - envelope(travel)).volume() > .001:
            raise ValueError("Conservative envelope does not contain spring mesh")
        targets = bounded({key: solid.translate(travel*axis) if owners[key] == "FJ1477_carriage" else solid
                           for key, solid in rigid_mount.items()})
        export_rows[label] = dict(mesh=mesh, mount_hits_mm3=overlaps(actual, targets))

    model_path = ROOT / extended["model"]["file"]
    if provenance(model_path) != extended["model"]:
        raise ValueError("Drive model changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    trace = checked_trace(audit["trace"])
    records = []
    for index, row in enumerate(trace):
        data.qpos[:] = [float(row[f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        poses = transforms(model, data, neutral)
        travel = max(0., float(row["slide_mm"]))
        h = poses["humerus"]
        spring_bound = envelope(travel).transform(h)
        bounds = np.array(spring_bound.bounding_box())
        candidates = {}
        for key, part in parts.items():
            current = part.transform(poses[owners[key]]) if owners[key] in poses else part
            other = np.array(current.bounding_box())
            if np.all(bounds[:3] <= other[3:]) and np.all(other[:3] <= bounds[3:]):
                candidates[key] = current
        # Only run the full swept wire on parts the envelope can meet. Pin
        # envelopes overlap by construction, so their seating is checked
        # exactly across stroke in the separate mounting screen below.
        pins = ("series:fixed_pin", "series:moving_pin")
        candidates = {k: v for k, v in candidates.items() if k not in pins}
        envelope_hits = overlaps(spring_bound, bounded(candidates))
        if envelope_hits:
            exact = geometry(travel).transform(h)
            hits = overlaps(exact, bounded({key: candidates[key] for key in envelope_hits}))
            records.append(dict(frame=index, time_s=float(row["time_s"]), slide_mm=travel, hits_mm3=hits))
        if index % 200 == 0:
            print("Spring path frame", index, flush=True)

    pin_results = []
    for travel in np.linspace(0., stroke, 21):
        targets = bounded({"fixed_pin": parts["series:fixed_pin"],
                           "moving_pin": parts["series:moving_pin"].translate(travel*axis)})
        pin_results.append(dict(travel_mm=float(travel), hits_mm3=overlaps(geometry(travel), targets)))
    spring = geometry(0.)
    open_mount = bounded({key: part for key, part in rigid_mount.items() if key not in pins})
    paths = []
    for sign in (-1., 1.):
        failures = []
        for distance in np.linspace(8., 0., 33):
            hit = overlaps(spring.translate(sign*distance*pin_axis), open_mount)
            if hit:
                failures.append(dict(distance_mm=float(distance), hits_mm3=hit))
        paths.append(dict(pin_axis_sign=sign, steps=33, rigid_hits=failures))
    save(folder / "return_spring_geometry.json", dict(
        generator=provenance(ROOT / "tools/build_s1_series_elbow_return_spring.py"), size=provenance(size_path),
        design=design, exports=export_rows, trajectory_frames=len(trace), exact_trajectory_candidates=records,
        pin_overlap_screen=pin_results, spring_approach_with_pins_removed=paths,
        nominal_hook_eye_ID_mm=2.1, nominal_hook_throat_mm=.9 + .15 - radius,
        scope="Continuous swept wire with analytic pitch and end eyes at two pin centres. Rigid surroundings only."
              " Installed spring geometry is imposed; transition and hook stresses, pin axial retention,"
              " tool access and manufacturing remain separate checks."))
    print("Exact candidate frames", len(records), "with collision", sum(bool(r["hits_mm3"]) for r in records),
          "insertion failures", [len(r["rigid_hits"]) for r in paths], flush=True)


if __name__ == "__main__":
    main()
