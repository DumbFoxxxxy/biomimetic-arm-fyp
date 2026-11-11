"""Export the smallest screened FJ1500 lateral rest-shape candidate."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, mesh_geom, save
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_local_force_geometry import load_reference


TPU = "tpu:infraspinatus"


def warp(points, origin, direction, length, offset):
    points = np.asarray(points)
    fraction = np.clip(((points - origin) @ direction / length - .30) / .48, 0., 1.)
    weights = fraction**2 * (3. - 2. * fraction)
    return points + weights[:, None] * offset


def main():
    folder = DEST.parent / "s1_local_muscle_load"
    output = folder / "FJ1500_cap_relief"
    screen_file = output / "tpu_rest_offset_screen.json"
    screen = read(screen_file)
    equilibrium = load_reference(screen["reference_equilibrium"])
    geometry = load_reference(equilibrium["local_geometry"])
    assembly = load_reference(equilibrium["shoulder_assembly"])
    route = next(row for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    sleeve = load_reference(screen["compact_sleeve"])["selected"]["source"]
    stable_rows = [row for row in screen["results"] if row["GH_flexion_deg"] in (0., 30., 60.)]
    offsets = sorted({row["rest_offset_mm"] for row in stable_rows})
    distance = next(value for value in offsets if all(not row["rigid_conflict_states"] for row in stable_rows
                                                       if row["rest_offset_mm"] == value))
    offset = distance * np.array(screen["trial_direction_source"])
    origin = np.array(route["ports"]["origin"]["inner_face_source_mm"])
    end = np.array(route["ports"]["distal"]["outer_face_source_mm"])
    axis = (end - origin) / np.linalg.norm(end - origin)
    length = np.linalg.norm(end - origin)
    if abs(axis @ offset) > 1e-10:
        raise ValueError("The offset must be lateral to keep the warp invertible")
    original_source = assembly["sources"][TPU]
    original = load_component(original_source)
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(original.to_mesh64())
    moved = warp(vertices, origin, axis, length, offset)
    if not np.allclose(warp(moved, origin, axis, length, -offset), vertices, atol=1e-10, rtol=0.):
        raise RuntimeError("Lateral warp did not invert")
    manifold = A.I._manifold()
    warped = manifold.Manifold(manifold.Mesh64(moved, faces.astype(np.uint32)))
    actual, meta = export(warped, output / "FJ1500_lateral_rest_TPU.obj",
                          np.array(original_source["translation_source_mm"]))
    route_points = np.array(route["full_wire_path_source_mm"])
    warped_route = warp(route_points, origin, axis, length, offset)
    wire = A.D.outer_envelope(warped_route, .3)
    collisions = {}
    for key, source in assembly["sources"].items():
        if not (key.startswith("bone:") or key.startswith("GH_cap_")
                or key.startswith((f"cap:{HEAD}/", f"free_terminal:{HEAD}", f"stop:{HEAD}/",
                                   f"origin:{HEAD}/", f"link:{HEAD}/muscle_"))):
            continue
        part = load_component(source)
        positioned = load_component(sleeve) if key == f"link:{HEAD}/muscle_uncrimped_sleeve" else part
        if assembly["intended_body_owners"][key] == f"free_terminal:{HEAD}" or key.startswith(f"link:{HEAD}/muscle_"):
            positioned = positioned.translate(offset)
        baseline = float(A.I.intersection_volume(original, part))
        current = float(A.I.intersection_volume(actual, positioned))
        if max(baseline, current) > .001:
            collisions[key] = dict(original_mm3=baseline, warped_rest_mm3=current)
    save(output / "lateral_rest_tpu_geometry.json", dict(
        generator=provenance(ROOT / "tools/build_s1_infraspinatus_offset_tpu.py"),
        bias_screen=provenance(screen_file), original=original_source, candidate=meta,
        rest_offset_mm=distance, rest_offset_source_mm=offset.tolist(), input_vertices=len(vertices),
        original_volume_mm3=float(original.volume()), warped_volume_mm3=float(actual.volume()),
        original_mesh_deflection_fraction_limits=[.30, .78], reference_wire_radius_mm=.3,
        rest_wire_intersection_mm3=float(A.I.intersection_volume(wire, actual)),
        rest_path_source_mm=warped_route.tolist(), flexible_rest_overlap_with_rigid_mm3=collisions,
        scope="Rest-shape geometry only. Original proximal region is fixed and the terminal region translates."
              " The two exterior ends of the local actuation route are warped with the muscle."
              " The original CAD strain, axial/lateral stiffness and actual deformed geometry are not inferred."))
    print("Rest offset", distance, "mm; original / warped volumes", original.volume(), actual.volume(),
          "wire overlap", A.I.intersection_volume(wire, actual), "rest overlaps", collisions, flush=True)


if __name__ == "__main__":
    main()
