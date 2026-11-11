"""Replace the competing shaft collar with supported S1 eyes on elbow adapters."""

from itertools import product
from pathlib import Path

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_guide import unit
from build_s1_teres_origin_mount import A, export, load_component, mesh_geom, rod, save
from study_elbow_series_package import bounded, overlaps, union


def one_component(solid):
    return sum(piece.volume() > 1e-5 for piece in solid.decompose()) == 1


def guard_seats(key, base, solids):
    prefix = key.removesuffix("adapter")
    pad = A.I._manifold().Manifold.cube([3.] * 3, True)
    return [base ^ solid.hull().minkowski_sum(pad) for name, solid in solids.items()
            if name.startswith(prefix) and "washer" in name and "bone_" in name]


def bridge_candidates(base, eye, point, pin_axis, obstacles, hardware):
    mf = A.I._manifold()
    if overlaps(eye - base, obstacles):
        return None, dict(tested=0, original_eye_blocked=True)
    vertices, _ = mesh_geom.manifold_mesh_to_arrays(base.to_mesh64())
    nearest = np.argsort(np.linalg.norm(vertices - point, axis=1))
    first = unit(np.cross(pin_axis, np.eye(3)[np.argmin(abs(pin_axis))]))
    second = np.cross(pin_axis, first)
    endpoints = []
    for angle in np.linspace(0., 2. * np.pi, 16, endpoint=False):
        direction = np.cos(angle) * first + np.sin(angle) * second
        endpoints.append(point + 2.5 * direction)
    tested = 0
    for radius in (2.6, 2.2, 1.8):
        for i, tip in product(nearest[:24], endpoints):
            surface = vertices[i]
            reach = np.linalg.norm(surface - tip)
            if reach < .01:
                continue
            direction = (surface - tip) / reach
            root = surface + 1.7 * direction
            midpoint = (tip + root) / 2.
            across = unit(np.cross(direction, np.eye(3)[np.argmin(abs(direction))]))
            other = np.cross(direction, across)
            middle_options = [None, *(midpoint + side * distance
                                      for distance in (3., 6., 9.) for side in (across, -across, other, -other))]
            for mid in middle_options:
                tested += 1
                if mid is None:
                    link = rod(tip, root, radius)
                else:
                    link = union([rod(tip, mid, radius), rod(mid, root, radius),
                                  mf.Manifold.sphere(radius, 16).translate(mid)])
                link -= hardware
                if (A.I.intersection_volume(link, base) <= 1.
                        or A.I.intersection_volume(link, eye) <= 1.
                        or overlaps((link + eye) - base, obstacles)):
                    continue
                candidate = union([base, link, eye]).as_original().simplify(1e-6)
                if one_component(candidate):
                    return candidate, dict(radius_mm=radius, tip_mm=tip.tolist(), root_mm=root.tolist(),
                                           midpoint_mm=None if mid is None else mid.tolist(), tested=tested)
    return None, dict(tested=tested)


def prepare_placement(base, point, eye, hardware, guard_neutral, other_bounds, protected, offset):
    posed_hardware = {name: solid.translate(offset) for name, solid in hardware.items()}
    impacts = {name: found for name, solid in posed_hardware.items()
               if (found := overlaps(solid, other_bounds))}
    if impacts:
        return None, dict(offset_mm=offset.tolist(), translated_hardware_hits_mm3=impacts)
    moved = point + offset
    guard = guard_neutral.translate(offset)
    reduced = (base - guard).as_original().simplify(1e-6)
    loss = [float((seat - reduced).volume()) for seat in protected]
    row = dict(offset_mm=offset.tolist(), protected_seat_loss_mm3=loss)
    if max(loss, default=0.) > .001 or not one_component(reduced):
        return None, row
    return (reduced, moved, eye.translate(offset), guard, posed_hardware), row


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    parent = folder / "assembly.json"
    assembly = read(parent)
    sources = assembly["sources"]
    solids = {key: load_component(row) for key, row in sources.items()}
    collar_file = A.DESIGN / "s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json"
    collar = read(collar_file)
    context = read(ROOT / assembly["shoulder"]["file"])
    mf = A.I._manifold()
    pad = mf.Manifold.cube([.3] * 3, True)
    removed_collar = {key for key in solids if key.startswith("collar:")}
    working = {key: solid for key, solid in solids.items() if key not in removed_collar}
    reports = []
    modified = {}
    for head, label, target in (("FJ1500", "upper", "FJ1477"), ("FJ1504", "lower", "FJ1486")):
        key = f"elbow:{target}_bone_adapter"
        base = working[key]
        anchor = collar["anchors"][head]
        point, axis = np.asarray(anchor["eye_centre_source_mm"]), unit(np.asarray(anchor["bore_axis_source"]))
        original_eye = solids[f"collar:jaw_{label}"] ^ rod(point - 1.5 * axis, point + 1.5 * axis, 3.23)
        hardware_keys = [name for name in working if name.startswith(f"link:{head}/bone_")]
        groove_key = f"link:{head}/bone_groove_ring"
        hardware = union([*(working[name] for name in hardware_keys),
                          load_component(context["sources"][groove_key])]).minkowski_sum(pad)
        bore = rod(point - 12. * axis, point + 12. * axis, 1.1)
        hardware += bore
        conflicts = [hardware]
        if target == "FJ1486":
            # The original deltoid cap stays in place. Relieve the elbow
            # adapter only where it intersected this existing cap.
            conflicts.append(working["cap:FJ1468/distal"].minkowski_sum(pad))
        reduced = (base - union(conflicts)).as_original().simplify(1e-6)
        protected = guard_seats(key, base, solids)
        seat_losses = [float((seat - reduced).volume()) for seat in protected]
        result = dict(head=head, target=key, eye_centre_source_mm=point.tolist(), pin_axis_source=axis.tolist(),
                      adapter_relief_mm3=float((base - reduced).volume()),
                      protected_seat_loss_mm3=seat_losses)
        reports.append(result)
        original_usable = one_component(reduced) and max(seat_losses, default=0.) <= .001
        result["original_trim_invalid"] = not original_usable
        others = {name: solid for name, solid in working.items() if name != key and name not in hardware_keys}
        offset_options = []
        if original_usable:
            offset_options.append((reduced, point, original_eye, hardware, {name: working[name] for name in hardware_keys}))
        if head == "FJ1500":
            hardware_parts = {name: working[name] for name in hardware_keys}
            hardware_parts[groove_key] = load_component(context["sources"][groove_key])
            other_bounds = bounded(others)
            tested = []
            directions = [unit(np.array([i, j, k], dtype=float)) for i, j, k in product((-1, 0, 1), repeat=3)
                          if i or j or k]
            for radius in (.5, 1., 1.5, 2., 2.5, 3., 4.):
                for direction in directions:
                    placement, row = prepare_placement(base, point, original_eye, hardware_parts,
                                                       hardware, other_bounds, protected, radius * direction)
                    tested.append(row)
                    if placement is not None:
                        offset_options.append(placement)
                if len(offset_options) >= 8:
                    break
            result["offset_placement_trials"] = tested
        result["bridge_trials"] = []
        for prepared, position, eye, keepout, own_parts in offset_options:
            obstacles = bounded({**others, **own_parts})
            candidate, details = bridge_candidates(prepared, eye, position, axis, obstacles, keepout)
            result["bridge_trials"].append(dict(offset_mm=(position - point).tolist(), **details))
            if candidate is None:
                continue
            try:
                reloaded, meta = export(candidate, folder / f"shared_anchor_{head}_adapter.obj", position)
            except (AssertionError, ValueError) as error:
                result["export_error"] = str(error) or type(error).__name__
                continue
            roundtrip = overlaps(reloaded, obstacles)
            eye_loss = float((eye - reloaded).volume())
            if roundtrip or eye_loss > .001:
                continue
            result.update(candidate_source=meta, rigid_neutral_hits_mm3=roundtrip, bridge=details,
                          final_eye_centre_source_mm=position.tolist(), eye_offset_mm=(position - point).tolist(),
                          exact_eye_material_removed_mm3=eye_loss,
                          final_adapter_relief_mm3=float((base - prepared).volume()),
                          final_protected_seat_loss_mm3=[float((seat - reloaded).volume()) for seat in protected])
            working[key] = reloaded
            modified[key] = meta
            for name in hardware_keys:
                shifted = dict(sources[name])
                shifted["translation_source_mm"] = (np.asarray(shifted["translation_source_mm"]) + position - point).tolist()
                modified[name] = shifted
                working[name] = own_parts[name]
            break
    save(folder / "shared_anchor_bridges.json", dict(
        generator=provenance(Path(__file__).resolve()), baseline=provenance(parent),
        original_eye_reference=provenance(collar_file), archived_collar_keys=sorted(removed_collar),
        reports=reports, prototype_replacements=modified,
        scope="Alternative to the full humeral S1 collar. Uses existing elbow bone mounts, holds "
        "the original Y pin axes and retains or locally shifts the pin centre for rigid clearance. "
        "No shoulder tension, flexible contact, tendon dynamics, assembly access or strength acceptance."))
    print(reports, flush=True)


if __name__ == "__main__":
    main()
