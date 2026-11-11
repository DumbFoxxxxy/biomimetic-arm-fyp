"""Bound stroke relief from existing cap swivels and check rigid travel."""

from itertools import combinations
from pathlib import Path

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from simulate_elbow_antagonistic_load import HEADS, build_local_model


def unit(vector):
    return vector / np.linalg.norm(vector)


def swivel_orbit(origin, bore, centre, axis):
    radius = bore - centre
    plane_centre = centre + axis * (axis @ radius)
    radial = bore - plane_centre
    tangent = np.cross(axis, radial)
    projected = origin - plane_centre
    angle = np.arctan2(projected @ tangent, projected @ radial)
    rho = np.linalg.norm(projected - axis * (projected @ axis))
    h = projected @ axis
    r = np.linalg.norm(radial)
    return plane_centre, radial, tangent, angle, np.hypot(h, rho - r), np.hypot(h, rho + r)


def geometry_bounds(model, data, head, elbow_deg):
    data.qpos[:] = 0.
    data.joint("elbow").qpos[:] = np.deg2rad(elbow_deg)
    mujoco.mj_kinematics(model, data)
    origin = data.site(f"{head}_origin").xpos.copy() * 1000.
    bore = data.site(f"{head}_distal").xpos.copy() * 1000.
    pivot = data.joint(f"{head}_swivel")
    centre, radial, tangent, angle, minimum, maximum = swivel_orbit(
        origin, bore, pivot.xanchor * 1000., unit(pivot.xaxis))
    elbow = data.joint("elbow")
    swivel_samples = np.deg2rad(np.arange(-180., 180., .1))
    positions = centre + np.outer(np.cos(swivel_samples), radial) + np.outer(np.sin(swivel_samples), tangent)
    lengths = np.linalg.norm(positions - origin, axis=1)
    arms = (np.cross(positions - elbow.xanchor * 1000., (origin - positions) / lengths[:, None])
            @ elbow.xaxis)
    pivot.qpos[:] = angle
    mujoco.mj_kinematics(model, data)
    measured = np.linalg.norm(data.site(f"{head}_distal").xpos * 1000. - origin)
    np.testing.assert_allclose(measured, minimum, atol=1e-7)
    pivot.qpos[:] = angle + np.pi
    mujoco.mj_kinematics(model, data)
    measured = np.linalg.norm(data.site(f"{head}_distal").xpos * 1000. - origin)
    np.testing.assert_allclose(measured, maximum, atol=1e-7)
    return dict(elbow_deg=elbow_deg, neutral_cap_length_mm=float(np.linalg.norm(bore - origin)),
                min_length_mm=float(minimum), max_length_mm=float(maximum),
                bore_orbit_radius_mm=float(np.linalg.norm(radial)),
                best_length_swivel_deg=float(np.rad2deg(angle)),
                sampled_min_signed_moment_arm_mm=float(arms.min()),
                sampled_max_signed_moment_arm_mm=float(arms.max()),
                sampled_max_arm_swivel_deg=float(np.rad2deg(swivel_samples[arms.argmax()])))


def rigid_hits(data, owner_solids, origins, pairs):
    keys = list(owner_solids)
    transforms = {}
    for name in origins:
        body = data.body(name)
        rotation = body.xmat.reshape(3, 3)
        transforms[name] = np.c_[rotation, body.xpos * 1000. - rotation @ origins[name]]
    placed = [solid.transform(transforms[name]) if name in transforms else solid
              for name, solid in owner_solids.values()]
    bounds = np.asarray([part.bounding_box() for part in placed])
    left, right = bounds[pairs[:, 0]], bounds[pairs[:, 1]]
    admitted = np.all(left[:, :3] <= right[:, 3:] + 1e-8, axis=1)
    admitted &= np.all(right[:, :3] <= left[:, 3:] + 1e-8, axis=1)
    result = []
    for i, j in pairs[admitted]:
        volume = float(A.I.intersection_volume(placed[i], placed[j]))
        if volume > .001:
            result.append(dict(parts=[keys[i], keys[j]], overlap_mm3=volume))
    return result


def rigid_travel(model, data, owner_solids, origins, head, elbow_deg, target_deg):
    owner = f"{head}_distal_pivot"
    owners = [name for name, _ in owner_solids.values()]
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(owners)), 2)
                        if (owners[i] == owner) != (owners[j] == owner)])
    data.qpos[:] = 0.
    data.joint("elbow").qpos[:] = np.deg2rad(elbow_deg)

    def hits(swivel):
        data.joint(f"{head}_swivel").qpos[:] = np.deg2rad(swivel)
        mujoco.mj_kinematics(model, data)
        return rigid_hits(data, owner_solids, origins, pairs)

    first = hits(0.)
    if first:
        return dict(last_clear_deg=None, first_blocked_deg=0., intersections=first)
    last = 0.
    for swivel in np.linspace(0., target_deg, max(1, int(np.ceil(abs(target_deg)))) + 1)[1:]:
        collisions = hits(swivel)
        if collisions:
            blocked = swivel
            while abs(blocked - last) > .1:
                mid = (last + blocked) / 2.
                middle_hits = hits(mid)
                if middle_hits:
                    blocked, collisions = mid, middle_hits
                else:
                    last = mid
            return dict(last_clear_deg=float(last), first_blocked_deg=float(blocked), intersections=collisions)
        last = swivel
    return dict(last_clear_deg=float(last), first_blocked_deg=None, intersections=[])


def load_geometry(check, stop):
    owner_map = {row["source_key"]: row["owner"] for row in check["geom_owners"]}
    owner_map.update({key: "humerus" for key in stop["sources"] if key.startswith("stop:")})
    origins = {name: np.asarray(check["neutral_origins_m"][name]) * 1000.
               for name in ("forearm", *(f"{head}_distal_pivot" for head in HEADS))}
    owner_solids = {key: (owner_map[key], load_component(source)) for key, source in stop["sources"].items()}
    return owner_solids, origins


def main():
    check_path = DEST / "elbow_passive_model/model_check.json"
    check = read(check_path)
    model_path = ROOT / check["model"]["file"]
    stop_path = DEST / "elbow_extension_stop/assembly.json"
    stop = read(stop_path)
    if provenance(model_path) != check["model"]:
        raise RuntimeError("Passive elbow model fingerprint changed")
    model = mujoco.MjModel.from_xml_string(build_local_model(model_path, check))
    data = mujoco.MjData(model)
    owner_solids, origins = load_geometry(check, stop)
    selected = (0, 10, 20, 30, 60, 90, 116, 145)
    results = {}
    for head in HEADS:
        rows = [geometry_bounds(model, data, head, angle) for angle in range(151)]
        baseline = rows[0]["neutral_cap_length_mm"]
        for angle in selected:
            row = rows[angle]
            row["min_length_change_from_neutral_mm"] = row["min_length_mm"] - baseline
            print(head, angle, "min length change", round(row["min_length_change_from_neutral_mm"], 3),
                  "max signed arm", round(row["sampled_max_signed_moment_arm_mm"], 3), flush=True)
            if angle > 30:
                continue
            target = row["best_length_swivel_deg"]
            travel = rigid_travel(model, data, owner_solids, origins, head, angle, target)
            if travel["last_clear_deg"] is not None:
                data.joint(f"{head}_swivel").qpos[:] = np.deg2rad(travel["last_clear_deg"])
                mujoco.mj_kinematics(model, data)
                length = np.linalg.norm(data.site(f"{head}_distal").xpos - data.site(f"{head}_origin").xpos)
                travel["last_sampled_length_change_mm"] = float(length * 1000. - baseline)
            row["local_rigid_screen_toward_min_length"] = travel
        results[head] = rows
    save(DEST.parent / "elbow_local_muscle_load/elbow_cap_travel.json", dict(
        generator=provenance(Path(__file__).resolve()), interfaces=provenance(check_path),
        geometry=provenance(stop_path), poses=results, rigid_parts=len(owner_solids),
        overlap_threshold_mm3=.001, swivel_step_max_deg=1., blocked_bracket_max_deg=.1,
        moment_arm_swivel_step_deg=.1, elbow_geometry_step_deg=1.,
        scope="Length bounds allow a full cap turn as a mathematical envelope; moment arms sample it. "
        "Rigid screening varies one swivel at a time from neutral along the shortest arc toward minimum "
        "length, with elbow fixed at 0/10/20/30 deg and the other swivel neutral. Only pairs affected "
        "by that swivel are checked, including the existing stop geometry. No forces or TPU contact."))


if __name__ == "__main__":
    main()
