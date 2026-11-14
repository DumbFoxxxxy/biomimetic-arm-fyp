"""Compare washer relief and right AC fastener access at the seated joint."""

import numpy as np

from audit_s1_infraspinatus_split_installation import insertion
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_AC_closure_screen.json"
    audit = read(source)
    scene = {key: load_component(meta) for key, meta in load_reference(audit["source"])["sources"].items()}
    hardware = {key.removeprefix("AC:AC_"): load_component(meta)
                for key, meta in audit["hardware_sources"].items()}
    axis = np.asarray(audit["axis_source"])
    washer_key = "bottom_washer_left"
    washer = hardware[washer_key]
    scapula = scene["bone:scapula"]
    clearance_cube = A.I._manifold().Manifold.cube([.12] * 3, True)
    swept_contacts = []
    for offset in np.linspace(0., 1., 41):
        contact = (washer.translate(-offset * axis) ^ scapula).translate(offset * axis)
        if contact.volume() > 1e-8:
            swept_contacts.append(contact.minkowski_sum(clearance_cube))
    candidate = (washer - union(swept_contacts)).as_original().simplify(1e-6)
    fixed = bounded({**scene, **{key: part for key, part in hardware.items()
                                if key not in (washer_key, "nut_left")}})
    trace = insertion(candidate, fixed, -axis, 12.)
    relief = dict(removed_mm3=float(washer.volume() - candidate.volume()),
                  original_mm3=float(washer.volume()), trace=trace)
    if not trace["hits"]:
        center = np.asarray(washer.bounding_box()).reshape(2, 3).mean(axis=0)
        actual, mesh = export(candidate, folder / "AC_left_rear_washer_relief.obj", center)
        relief["mesh"] = mesh
        relief["readback"] = insertion(actual, fixed, -axis, 12.)
        if relief["readback"]["hits"]:
            raise ValueError("Relieved washer path changed after readback")

    right = []
    oblique = []
    fixed_tools = bounded({**scene, **{key: part for key, part in hardware.items()
                                      if key != "bolt_right"}})
    for name, direction, radii in (("bolt_right", axis, (1., 1.5, 2.)),
                                  ("nut_right", -axis, (3.4, 4., 4.5))):
        part = hardware[name]
        vertices = np.asarray(part.to_mesh64().vert_properties)[:, :3]
        coordinates = vertices @ direction
        face = vertices[coordinates >= coordinates.max() - .001]
        start = face.mean(axis=0) + .05 * direction
        if name == "bolt_right":
            basis = A.D.unit(np.cross(direction, np.eye(3)[np.argmin(abs(direction))]))
            other = np.cross(direction, basis)
            for tilt_deg in (5., 10., 15., 20., 25., 30.):
                tilt = np.deg2rad(tilt_deg)
                for azimuth in np.linspace(0., 2. * np.pi, 24, endpoint=False):
                    vector = (np.cos(tilt) * direction
                              + np.sin(tilt) * (np.cos(azimuth) * basis + np.sin(azimuth) * other))
                    short_tool = rod(start, start + 20. * vector, 1.5, 48)
                    oblique.append(dict(tilt_deg=tilt_deg, direction_source=vector.tolist(),
                                        hits_mm3=overlaps(short_tool, fixed_tools)))
        # Test coaxial tools from both faces, including the reverse fastener
        # arrangement's approach side. These are access envelopes only.
        for sign in (1., -1.):
            local_direction = sign * direction
            if sign > 0:
                point = start
            else:
                point = start + (coordinates.min() - start @ direction - .05) * direction
            for radius in radii:
                body = rod(point, point + 20. * local_direction, radius, 64)
                if name.startswith("nut_"):
                    body -= rod(point - local_direction, point + 21. * local_direction, 2.1, 64)
                hits = overlaps(body, bounded(scene))
                right.append(dict(part=name, side=sign, radius_mm=radius,
                                  start_source_mm=point.tolist(), hits_mm3=hits))
    save(folder / "FJ1504_AC_fastener_access.json", dict(
        generator=provenance(ROOT / "tools/study_s1_ac_fastener_access.py"),
        source=provenance(source), left_washer=relief, right_axial_tool_screen=right,
        right_oblique_driver_screen=oblique,
        scope="Nominal axis of the existing fasteners. Left washer relief changes only the "
              "washer. Right probes test bare access from each face, without assuming new "
              "thread engagement, captive hardware or a changed clamp geometry.",
    ))
    print("left washer", relief["removed_mm3"], len(trace["hits"]), "exported", "mesh" in relief)
    print("right", [(row["part"], row["side"], row["radius_mm"], list(row["hits_mm3"])) for row in right])
    print("oblique clear by tilt", [(angle, sum(row["tilt_deg"] == angle and not row["hits_mm3"]
                                               for row in oblique))
                                    for angle in sorted({row["tilt_deg"] for row in oblique})])


if __name__ == "__main__":
    main()
