"""Build and recheck seat prototypes for the shared scapular antagonist pair."""

from itertools import combinations

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, rod, save
from screen_whole_arm_girdle_routes import SharedAssembly
from study_elbow_series_package import overlaps


PAIR = ("FJ1459_g3", "FJ1536_g2")


def build_seat(base, port, inward, bone, radius):
    axis = port - base
    axis /= np.linalg.norm(axis)
    # Extend to the source surface for the local fixture footprint.
    post = rod(base - .7 * axis, port, radius)
    # A short eye leaves room for the wire to change angle along the motion.
    sleeve = rod(port - .6 * inward, port + .6 * inward, 1.05)
    hole = rod(port - 1.8 * inward, port + 1.8 * inward, .65)
    raw = (post + sleeve) - hole
    return raw - bone, float((raw ^ bone).volume())


def inspect_state(seats, routes, pose):
    frames, solids, _ = pose
    ports = {}
    placed = {}
    hits = []
    for name, (part, info) in seats.items():
        rotation, offset = frames[info["owner"]]
        placed[name] = part.transform(np.column_stack((rotation, offset)))
        ports[name] = rotation @ info["port_source_mm"] + offset
        for key, volume in overlaps(placed[name], solids).items():
            hits.append(dict(left=name, right=key, overlap_mm3=volume))
    cables = {}
    for route in routes:
        a, b = (ports[f"{route['id']}_{role}"] for role in ("origin", "insertion"))
        cables[route["id"]] = rod(a, b, .30)
        for key, volume in overlaps(cables[route["id"]], solids).items():
            hits.append(dict(left=route["id"], right=key, overlap_mm3=volume))
    for (left, a), (right, b) in combinations({**placed, **cables}.items(), 2):
        volume = A.I.intersection_volume(a, b)
        if volume > .001:
            hits.append(dict(left=left, right=right, overlap_mm3=volume))
    return hits


def main():
    output = DEST.parent / "whole_arm_first_pass"
    screen_file = output / "girdle_route_screen.json"
    screen = read(screen_file)
    shared = SharedAssembly(DEST / "unified_SC_AC")
    if screen["assembly"] != provenance(shared.assembly_file):
        raise ValueError("The anchor screen uses a different assembly")
    routes = [next(row for row in screen["routes"] if row["id"] == key) for key in PAIR]
    neutral = shared.pose()
    seats = {}
    for route in routes:
        selected = route["neutral_selection"]
        base, ports = map(np.asarray, (selected["bone_attach_reference_mm"], selected["source_points_mm"]))
        vector = ports[1] - ports[0]
        vector /= np.linalg.norm(vector)
        for index, role in enumerate(("origin", "insertion")):
            name = f"{route['id']}_{role}"
            owner = route["owners"][index]
            footprint = selected["seat_axis_probe_hits_mm3"][index]
            if not footprint:
                raise ValueError(f"No indicated support at {name}")
            bone_key = max(footprint, key=footprint.get)
            if shared.solids[bone_key][0] != owner:
                raise ValueError(f"Wrong bone owner at {name}")
            bone = shared.solids[bone_key][1]
            port = ports[index]
            direction = vector if index == 0 else -vector
            trials = []
            for radius in (1., .85, .75):
                candidate, fit = build_seat(base[index], port, direction, bone, radius)
                parts = [part for part in candidate.decompose() if part.volume() > 1e-8]
                checks = overlaps(candidate, neutral[1])
                trials.append(dict(radius_mm=radius, components=len(parts), bone_trim_mm3=fit,
                                   neutral_overlaps_mm3=checks))
                if len(parts) == 1 and fit > .001 and not checks:
                    break
            else:
                raise ValueError(f"No seated interface candidate for {name}: {trials}")
            rebuilt, source = export(candidate, output / f"{name}_eye.obj", port)
            seats[name] = (rebuilt, dict(owner=owner, bone=bone_key, source=source,
                                        attachment_source_mm=base[index].tolist(), port_source_mm=port.tolist(),
                                        nominal_eye_diameter_mm=1.3, trials=trials))
    branch = [next(row for row in route["motion"] if row["axis"] == 0 and row["direction"] == -1)
              for route in routes]
    paths = [[(row["angle_deg"], row["ac_rotvec_deg"]) for row in item["outgoing"]] for item in branch]
    common = paths[0][:min(map(len, paths))]
    if any(path[:len(common)] != common for path in paths):
        raise ValueError("The two routes use incompatible coordination samples")
    accepted = []
    first_limit = None
    for angle, ac in common:
        hits = inspect_state(seats, routes, shared.pose(0, angle, ac))
        if hits:
            first_limit = dict(angle_deg=angle, ac_rotvec_deg=ac, hits=hits)
            break
        accepted.append(dict(angle_deg=angle, ac_rotvec_deg=ac))
    return_samples = 0
    for row in reversed(accepted):
        hits = inspect_state(seats, routes, shared.pose(0, row["angle_deg"], row["ac_rotvec_deg"]))
        if hits:
            raise ValueError(f"Return interface check failed at {row['angle_deg']}: {hits}")
        return_samples += 1
    save(output / "girdle_seat_pair.json", dict(
        assembly=screen["assembly"], model=screen["model"], route_screen=provenance(screen_file),
        generator=provenance(ROOT / "tools/build_whole_arm_girdle_seats.py"),
        seats={key: info for key, (_, info) in seats.items()}, checked_joint_axis=0,
        outward=accepted, return_samples=return_samples, first_seat_or_wire_limit=first_limit,
        reference_limits=[item["first_limit"] for item in branch],
        scope="Shared-assembly solid check of four conformal seat candidates and two 0.60 mm chords. "
              "Bone-conforming footprints allow a future bonded or keyed interface; bond, detachable "
              "fastening, crimp/stop detail and strength have not been designed or validated. "
              "No deformation, friction or actuation model."))
    print(f"Four seats exported; {len(accepted)} outbound and {return_samples} return samples", flush=True)
    if first_limit:
        print(first_limit, flush=True)


if __name__ == "__main__":
    main()
