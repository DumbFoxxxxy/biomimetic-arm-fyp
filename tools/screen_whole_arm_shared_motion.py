"""Measure first rigid limits of prescribed motions in the assembled arm."""

from itertools import combinations
import argparse

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import find_hits, transforms
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from project_runtime import load_mujoco_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alternatives-only", action="store_true")
    args = parser.parse_args()
    folder = DEST.parent / "whole_arm_first_pass"
    check_file = folder / "whole_arm_shared_scene_check.json"
    check = read(check_file)
    transfer = read(ROOT / check["transfer"]["file"])
    scene = ROOT / check["model"]["file"]
    if provenance(scene) != check["model"]:
        raise ValueError("The whole-arm scene changed since the load check")
    model = load_mujoco_model(scene)
    data = mujoco.MjData(model)
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"]}
    keys = list(sources)
    owners = transfer["owner_by_source"]
    parts = [load_component(sources[key]) for key in keys]
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                        if owners[keys[i]] != owners[keys[j]]])
    origins = {key: np.asarray(check["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    gh = model.joint("GH").dofadr.item()
    elbow = model.joint("elbow").dofadr.item()
    sc = model.joint("SC").dofadr.item()
    ac = model.joint("AC").dofadr.item()
    pointer = read(ROOT / "design/S1_CURRENT_CANDIDATE.json")
    axes_file = ROOT / pointer["bare_motion_report"]
    axes = {name: np.asarray(value) for name, value in read(axes_file)["axes_source"].items()}
    neutral_humerus = origins["forearm"] - origins["humerus"]
    neutral_humerus /= np.linalg.norm(neutral_humerus)
    seats_file = folder / "girdle_seat_pair.json"
    seat_report = read(seats_file)
    seat_routes = {
        key: [seat_report["seats"][f"{key}_{role}"] for role in ("origin", "insertion")]
        for key in ("FJ1459_g3", "FJ1536_g2")
    }
    neutral_scapular_line = origins["humerus"] - origins["scapula"]

    def sample(delta, check_girdle_wires=False):
        mujoco.mj_resetData(model, data)
        mujoco.mj_integratePos(model, data.qpos, delta, 1.)
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        posed = [part.transform(frames[owners[key]]) for key, part in zip(keys, parts)]
        hits = find_hits(posed, keys, pairs)
        direction = data.body("forearm").xpos - data.body("humerus").xpos
        direction /= np.linalg.norm(direction)
        # Absolute direction change is retained for combined paths.
        angle = np.rad2deg(np.arccos(np.clip(direction @ neutral_humerus, -1., 1.)))
        # Signed change of the AC-to-GH reference line projected into the source XZ plane.
        line = frames["scapula"][:, :3] @ neutral_scapular_line
        inclination = np.arctan2(line[2], -line[0]) - np.arctan2(
            neutral_scapular_line[2], -neutral_scapular_line[0])
        result = dict(qpos=data.qpos.tolist(), hits=hits, humerus_direction_change_deg=float(angle),
                      coronal_scapular_line_change_deg=float(np.rad2deg(np.arctan2(
                          np.sin(inclination), np.cos(inclination)))))
        if check_girdle_wires and not hits:
            placed_wires = []
            for route, seats in seat_routes.items():
                points = [frames[seat["owner"]] @ np.r_[seat["port_source_mm"], 1.] for seat in seats]
                placed_wires.append(rod(*points, .30))
            all_solids = [*posed, *placed_wires]
            wire_keys = [*keys, *seat_routes]
            wire_pairs = np.asarray([(i, j) for j in range(len(keys), len(wire_keys)) for i in range(j)])
            result["hits"] = find_hits(all_solids, wire_keys, wire_pairs)
            result["girdle_wire_min_gap_capped_3mm"] = min(
                A.I.gap(all_solids[i], all_solids[j], 3.) for i, j in wire_pairs)
        return result

    neutral = sample(np.zeros(model.nv))
    if neutral["hits"]:
        raise ValueError("Verified neutral scene could not be reproduced")
    specifications = []
    for name, axis, sign, maximum in (
        ("K2_GH_flexion", "flexion", 1, 120),
        ("K3_GH_abduction", "abduction", 1, 120),
        ("GH_axial_positive", "longitudinal", 1, 90),
        ("GH_axial_negative", "longitudinal", -1, 90),
    ):
        delta = np.zeros(model.nv)
        delta[gh:gh + 3] = axes[axis] * sign
        specifications.append((name, delta, maximum))
    delta = np.zeros(model.nv)
    delta[elbow] = 1.
    specifications.insert(0, ("K1_elbow_fixed_passive_interfaces", delta, 145))
    results = {}
    for name, delta, maximum in ([] if args.alternatives_only else specifications):
        accepted = [dict(input_deg=0., **neutral)]
        first = None
        for angle in range(1, maximum + 1):
            result = sample(np.deg2rad(angle) * delta)
            if result["hits"]:
                for tenth in range(1, 10):
                    fine = angle - 1 + tenth / 10.
                    refined = sample(np.deg2rad(fine) * delta)
                    if refined["hits"]:
                        first = dict(input_deg=fine, **refined)
                        break
                    accepted.append(dict(input_deg=fine, **refined))
                first = first or dict(input_deg=angle, **result)
                break
            accepted.append(dict(input_deg=float(angle), **result))
        for row in reversed(accepted):
            result = sample(np.deg2rad(row["input_deg"]) * delta)
            if result["hits"]:
                raise ValueError(f"Return path differs for {name}")
        results[name] = dict(accepted=accepted, first_limit=first, checked_return=len(accepted),
                             checked_upper_bound_deg=maximum)
        print(f"{name}: {accepted[-1]['input_deg']:.1f} deg, "
              f"first hit {first['hits'][0] if first else None}", flush=True)

    # The saved SC/AC trajectory was originally verified on the bare chain.
    # Check where the transferred arm interfaces first restrict that same path.
    seat_path = seat_report["outward"]
    coordinated = []
    for row in seat_path:
        delta = np.zeros(model.nv)
        delta[sc] = np.deg2rad(row["angle_deg"])
        delta[ac:ac + 3] = np.deg2rad(row["ac_rotvec_deg"])
        coordinated.append((row["angle_deg"], delta))
    prior_branches = (
        ("SC_AC_coordinated", 0., 0.), ("SC_AC_with_elbow", .8, 0.),
        ("K7_coordinated_with_GH_elevation", .8, 1.),
    )
    for name, elbow_fraction, shoulder_fraction in (() if args.alternatives_only else prior_branches):
        accepted = []
        first = None
        samples = []
        for angle, delta in coordinated:
            delta = delta.copy()
            delta[elbow] = np.deg2rad(abs(angle) * elbow_fraction)
            delta[gh:gh + 3] = np.deg2rad(abs(angle) * shoulder_fraction) * axes["abduction"]
            result = dict(input_deg=angle, **sample(delta))
            if result["hits"]:
                first = result
                break
            samples.append(delta)
            accepted.append(result)
        for delta in reversed(samples):
            if sample(delta)["hits"]:
                raise ValueError(f"Return path differs for {name}")
        results[name] = dict(accepted=accepted, first_limit=first, checked_return=len(accepted),
                             checked_upper_bound_deg=abs(coordinated[-1][0]))
        print(f"{name}: {accepted[-1]['input_deg']:.3f} deg, "
              f"first hit {first['hits'][0] if first else None}", flush=True)

    if args.alternatives_only:
        existing = read(folder / "whole_arm_shared_motion.json")
        if existing["model"] != check["model"]:
            raise ValueError("Existing action ranges use a different scene")
        results = existing["paths"]

    for axis, ratio, sign in ((0, 0., 1), (1, 0., 1), (2, .5, 1), (1, -2., -1)):
        for combined in (False, True):
            name = f"{'K7' if combined else 'scapular'}_SC{axis}_{'pos' if sign > 0 else 'neg'}"
            delta = np.zeros(model.nv)
            delta[sc + axis] = sign
            delta[ac + axis] = sign * ratio
            if combined:
                delta[elbow] = 2.
                delta[gh:gh + 3] = 3. * axes["abduction"]
            samples = []
            first = None
            steps_per_degree = max(1, int(np.ceil(np.max(np.abs(delta)))))
            ticks = 10 * steps_per_degree
            for tick in range(ticks + 1):
                amount = tick / steps_per_degree
                state = sample(np.deg2rad(amount) * delta, check_girdle_wires=True)
                if state["hits"]:
                    if tick:
                        for fraction in range(1, 10):
                            fine_amount = (tick - 1 + fraction / 10.) / steps_per_degree
                            fine = sample(np.deg2rad(fine_amount) * delta, check_girdle_wires=True)
                            if fine["hits"]:
                                state, amount = fine, fine_amount
                                break
                            samples.append(dict(input_deg=sign * fine_amount, **fine))
                    first = dict(input_deg=sign * amount, **state)
                    break
                samples.append(dict(input_deg=sign * amount, **state))
            if not samples:
                raise ValueError(f"Two-wire neutral assembly is blocked: {first}")
            for row in reversed(samples):
                if sample(np.deg2rad(abs(row["input_deg"])) * delta, check_girdle_wires=True)["hits"]:
                    raise ValueError(f"Return path differs for {name}")
            results[name] = dict(accepted=samples, first_limit=first, checked_return=len(samples),
                                 checked_upper_bound_deg=10., AC_per_SC=ratio,
                                 GH_abduction_per_abs_SC=3. if combined else 0.,
                                 elbow_flexion_per_abs_SC=2. if combined else 0.)
            end = samples[-1]
            print(name, end["input_deg"], end["coronal_scapular_line_change_deg"],
                  end["humerus_direction_change_deg"], first["hits"][0] if first else None, flush=True)

    save(folder / "whole_arm_shared_motion.json", dict(
        model=check["model"], scene_check=provenance(check_file), anatomical_axes=provenance(axes_file),
        coordinated_input=provenance(seats_file), cross_owner_pair_count=len(pairs), paths=results,
        scope="Prescribed poses, 1 degree steps and 0.1 degree refinement for single-axis probes. "
              "SC/AC reference paths follow supplied samples. Alternate paths also check both girdle wires "
              "against all rigid geometry and each other; the reported AC-to-GH coronal angle is a scapular orientation proxy. "
              "Elbow swivel angles and slider fixed at zero; "
              "this deliberately screens that interface condition only. Fixed TPU references and cable "
              "actuation are outside this rigid motion audit. Axial signs are source convention."))


if __name__ == "__main__":
    main()
