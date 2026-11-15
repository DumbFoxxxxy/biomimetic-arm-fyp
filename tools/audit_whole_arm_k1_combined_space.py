"""Screen the eight TPU bodies and elbow return spring in the relieved shared arm."""

import argparse
from itertools import combinations

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import export_lattice, weigh
from audit_s1_series_elbow_integration import transforms
from audit_s1_series_elbow_spring_integration import load_candidate
from audit_whole_arm_k1_shared import CAP, section
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_series_elbow_return_spring import spring_path, sweep
from build_s1_teres_origin_mount import A, export, load_component, save
from project_runtime import load_mujoco_model
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def affine(points, frame):
    return points @ frame[:, :3].T + frame[:, 3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--spring-relief", action="store_true", help="Try local spring clearance in the existing parts")
    options.add_argument("--refine-spring-relief", action="store_true", help="Clear remaining spring contacts from the first trial")
    args = parser.parse_args()
    folder = DEST.parent / "whole_arm_first_pass"
    path_ref = provenance(folder / "whole_arm_k1_cap_relief.json")
    path = load_reference(path_ref)
    if not path["cap_relief"]["geometry_accepted"] or any(row["hits"] for row in path["K7_relieved_cap"]):
        raise ValueError("Requires the cap already screened on both trajectories")
    combined = load_reference(path["combined"])
    tpu_ref = provenance(folder / "whole_arm_k7_tpu_moving_cuff_audit.json")
    tpu = load_reference(tpu_ref)
    if tpu["combined"] != path["combined"]:
        raise ValueError("Muscle meshes use another shared arm")
    motion = load_reference(combined["motion"])
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    routes = load_reference(combined["routed"])["routes"]
    support = load_reference(combined["support"])
    revision = load_reference(load_reference(support["audit"])["humeral_revision"])["selected"]["meshes"]
    base = load_reference(path["base"])
    rigid_sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
                     "elbow:FJ1477_bone_adapter": support["source"], "base:T12": base["source"],
                     CAP: path["cap_relief"]["selected_source"]}
    rigid_sources.update({key: source for key, source in revision.items()
                          if key.startswith("link:") and key in rigid_sources})
    owners = {**transfer["owner_by_source"], "base:T12": "thorax"}
    rigid = {key: load_component(source) for key, source in rigid_sources.items()}
    statics = {key: load_component(source) for key, source in tpu["tpu_sources"].items()
               if not key.startswith("tpu:elbow_")}
    elbow_folder = DEST.parent / "integration_review/s1_series_elbow"
    prior, _, mount, spring, _, _, retention_revision, retention = load_candidate(elbow_folder)
    spring_ref = provenance(elbow_folder / "return_spring_geometry.json")
    if spring != load_reference(spring_ref):
        raise ValueError("Spring definition differs from its swept geometry record")
    mf = A.I._manifold()
    active = {}
    for role, channel in (("flexor", "FJ1486"), ("extensor", "FJ1477")):
        name = f"tpu:elbow_{role}"
        source = tpu["tpu_sources"][name]
        if source != prior["mesh_sources"][role]:
            raise ValueError(f"Muscle mesh no longer matches the elbow reference: {role}")
        neutral, vertices, faces = indexed_solid(source)
        if neutral.status() != mf.Error.NoError:
            raise ValueError(f"Invalid neutral {role} geometry")
        route = routes[channel]
        first, last = np.asarray(route["points_source_mm"])[[0, -1]]
        vector = last - first
        weights = weigh(vertices, first, vector / (vector @ vector))[:, None]
        active[name] = (vertices, faces, weights, route["owners"])

    start = np.asarray(mount["fixed_eye_source_mm"])
    neutral_end = np.asarray(mount["moving_eye_neutral_source_mm"])
    pin_axis = np.asarray(mount["pin_axis_source"])
    for role, point in (("fixed", start), ("moving", neutral_end)):
        pin = rigid[f"series:{role}_pin"]
        marker = mf.Manifold.sphere(.1, 12).translate(point)
        if (marker - pin).volume() > 1e-5:
            raise ValueError(f"Spring {role} eye is not centred on the shared pin")
    spec = spring["design"]
    radius = spec["wire_diameter_mm"] / 2.
    coil_radius = spec["coil_OD_mm"] / 2. - radius
    turns = spec["active_turns"] + 2
    model_path = ROOT / path["assembly"]["file"]
    if provenance(model_path) != path["assembly"]:
        raise ValueError("The checked assembly changed")
    model = load_mujoco_model(model_path)
    data = mujoco.MjData(model)
    origins = {key: np.asarray(value) for key, value in scene["neutral_body_origins_mm"].items()}
    trace = path["accepted"]
    non_elbow_joints = [model.joint(name).qposadr.item() for name in ("SC", "AC", "GH")]
    parked = [index for start_index in non_elbow_joints for index in range(start_index, start_index + 4)]
    for row in trace:
        np.testing.assert_allclose(np.asarray(row["qpos"])[parked], model.qpos0[parked], atol=1e-10, rtol=0.)
    angles = np.asarray([row["elbow_deg"] for row in trace])
    tpu_indices = {int(np.argmin(abs(angles - goal))) for goal in range(0, 31, 5)} | {len(trace) - 1}
    fixed_mount = "elbow:FJ1477_bone_adapter"
    deltoid_name = "tpu:deltoid"

    def posed_spring(frames):
        moving_end = affine(neutral_end, frames["FJ1477_carriage"])
        vector = moving_end - start
        length = np.linalg.norm(vector)
        axis = vector / length
        orthogonal_pin = pin_axis - axis * (pin_axis @ axis)
        orthogonal_pin /= np.linalg.norm(orthogonal_pin)
        basis = np.c_[orthogonal_pin, np.cross(axis, orthogonal_pin), axis, start]
        coil = sweep(spring_path(length, turns, radius, coil_radius), radius).transform(basis)
        if coil.status() != mf.Error.NoError:
            raise ValueError("Invalid spring mesh")
        return coil, float(length)

    trial = {}
    parent_trial = None
    if args.spring_relief or args.refine_spring_relief:
        # Cut only local contact patches. The moving end patch is mapped back
        # into carriage coordinates before any union or stock subtraction.
        parts = {fixed_mount: rigid[fixed_mount], CAP: rigid[CAP], deltoid_name: statics[deltoid_name]}
        training = sorted(tpu_indices)
        roles = {fixed_mount: "fixed_seat", CAP: "moving_seat", deltoid_name: "deltoid"}
        if args.refine_spring_relief:
            parent_trial = provenance(folder / "whole_arm_k1_spring_relief.json")
            previous = load_reference(parent_trial)
            if previous["path"] != path_ref or previous["muscle_meshes"] != tpu_ref:
                raise ValueError("First spring relief trial is stale")
            if len(previous["samples"]) != len(trace):
                raise ValueError("First relief trial has incomplete pose coverage")
            trial = {key: previous["trial"][key] for key in (fixed_mount, CAP)}
            for key, result in trial.items():
                rigid[key] = load_component(result["source"])
            parts = {deltoid_name: load_component(previous["trial"][deltoid_name]["source"])}
            names = {row["trace_frame"] for row in previous["samples"] if row["spring_to_static_TPU_mm3"]}
            training = [i for i, row in enumerate(trace) if row["trace_frame"] in names]
            roles = {deltoid_name: "deltoid_refined"}
        cuts = {name: [] for name in parts}
        pad = mf.Manifold.cube([.4] * 3, True)
        for index in training:
            data.qpos[:] = trace[index]["qpos"]
            mujoco.mj_forward(model, data)
            frames = transforms(model, data, origins)
            coil, _ = posed_spring(frames)
            for key, part in parts.items():
                posed = coil
                if key == CAP:
                    carriage = frames[owners[CAP]]
                    rotation, offset = carriage[:, :3], carriage[:, 3]
                    posed = coil.transform(np.c_[rotation.T, -rotation.T @ offset])
                for patch in (part ^ posed).decompose():
                    if patch.volume() > .001:
                        cuts[key].append(patch.hull().minkowski_sum(pad))
        cable_guards = union(A.D.outer_envelope(np.asarray(routes[key]["points_source_mm"]), .72)
                             for key in ("FJ1468", "FJ1467", "FJ1513"))
        pin_guard = mf.Manifold.cube([.12] * 3, True)
        guards = {fixed_mount: rigid["series:fixed_pin"].minkowski_sum(pin_guard),
                  CAP: rigid["series:moving_pin"].minkowski_sum(pin_guard), deltoid_name: cable_guards}
        for key, label in roles.items():
            original = parts[key]
            raw_cut = union(cuts[key])
            clash = A.I.intersection_volume(original ^ raw_cut, guards[key])
            candidate = (original - (raw_cut - guards[key])).as_original().simplify(1e-7)
            pieces = sorted((p for p in candidate.decompose() if p.volume() > 1e-5), key=lambda p: -p.volume())
            discarded = 0.
            if key == deltoid_name and len(pieces) > 1:
                small = pieces[1:]
                discarded = sum(p.volume() for p in small)
                if discarded <= .002 * original.volume() and all(
                        A.I.intersection_volume(p, guards[key]) <= .001 for p in small):
                    candidate = (candidate ^ pieces[0]).as_original().simplify(1e-7)
                else:
                    discarded = 0.
            removed = original - candidate
            protected_loss = A.I.intersection_volume(removed, guards[key])
            components = [piece.volume() for piece in candidate.decompose() if piece.volume() > 1e-5]
            selected = None
            if len(components) == 1 and protected_loss <= .001:
                meta = tpu["tpu_sources"][key] if key == deltoid_name else rigid_sources[key]
                writer = export_lattice if key == deltoid_name else export
                _, selected = writer(candidate, folder / f"K1_spring_relief_{label}.obj",
                                     np.asarray(meta["translation_source_mm"]))
                checked = load_component(selected)
                added = max(0., (checked - original).volume())
                if added > .001:
                    raise ValueError(f"Relief added material after export: {key}")
                if key == deltoid_name:
                    statics[key] = checked
                else:
                    rigid[key] = checked
            trial[key] = dict(source=selected, components_mm3=components, cut_count=len(cuts[key]),
                              raw_cut_in_guard_mm3=clash, removed_mm3=removed.volume(),
                              removed_in_guard_mm3=protected_loss, discarded_fragments_mm3=discarded,
                              net_stock_volume_change_mm3=original.volume() - candidate.volume())
        print("K1 spring relief candidates:", [(key, row["source"] is not None, row["removed_mm3"])
                                                for key, row in trial.items()], flush=True)
        mujoco.mj_resetData(model, data)
    # Limit the costly lattice booleans to the sampled angles, while checking
    # the spring and the wires at every already approved eight-line pose.
    samples = []
    static_tpu = bounded(statics)
    moving = {"forearm", "FJ1486_distal_pivot", "FJ1477_distal_pivot", "FJ1477_carriage"}
    fixed_rigid = bounded({key: part for key, part in rigid.items() if owners[key] not in moving})
    dynamic_rigid = {key: part for key, part in rigid.items() if owners[key] in moving}
    mujoco.mj_forward(model, data)
    neutral_frames = transforms(model, data, origins)
    fixed_wires = bounded({key: section(route, neutral_frames)[0] for key, route in routes.items()
                           if not any(owner in moving for owner in route["owners"])})
    dynamic_wires = {key: route for key, route in routes.items() if key not in fixed_wires}
    for index, row in enumerate(trace):
        data.qpos[:] = row["qpos"]
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        posed = {key: part.transform(frames[owners[key]]) for key, part in dynamic_rigid.items()}
        hard = {**fixed_rigid, **bounded(posed)}
        coil, length = posed_spring(frames)
        wires = {**fixed_wires, **bounded({key: section(route, frames)[0] for key, route in dynamic_wires.items()})}
        sample = dict(trace_frame=row["trace_frame"], elbow_deg=row["elbow_deg"], separation_mm=float(length),
                      spring_to_rigid_mm3=overlaps(coil, hard), spring_to_wire_mm3=overlaps(coil, wires),
                      spring_to_static_TPU_mm3=overlaps(coil, static_tpu))
        if index in tpu_indices:
            active_shapes = {}
            for key, (vertices, faces, weight, bodies) in active.items():
                a, b = (frames[body] for body in bodies)
                moved = (1. - weight) * affine(vertices, a) + weight * affine(vertices, b)
                solid = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), np.ascontiguousarray(faces, dtype=np.uint64)))
                if solid.status() != mf.Error.NoError or solid.volume() <= 0.:
                    raise ValueError(f"Invalid posed muscle {key} in frame {index}")
                active_shapes[key] = solid
            muscles = {**statics, **active_shapes}
            target_tpu = bounded(muscles)
            sample["moving_TPU_to_rigid_mm3"] = {key: overlaps(part, hard) for key, part in active_shapes.items()}
            sample["wire_to_TPU_mm3"] = {key: overlaps(wire, target_tpu) for key, (wire, _) in wires.items()}
            sample["spring_to_moving_TPU_mm3"] = overlaps(coil, bounded(active_shapes))
            collisions = {}
            for a, b in combinations(muscles, 2):
                if a in active_shapes or b in active_shapes:
                    contact = overlaps(muscles[a], {b: target_tpu[b]})
                    if contact:
                        collisions[f"{a} / {b}"] = contact[b]
            sample["moving_TPU_to_other_TPU_mm3"] = collisions
            print(f"K1 flexible geometry {row['elbow_deg']:.2f} degrees", flush=True)
        samples.append(sample)
    target = ("whole_arm_k1_spring_relief_refined.json" if args.refine_spring_relief else
              "whole_arm_k1_spring_relief.json" if args.spring_relief else "whole_arm_k1_combined_space.json")
    save(folder / target, dict(
        path=path_ref, muscle_meshes=tpu_ref, spring_definition=spring_ref, prior_retention=retention,
        prior_retention_installed={key: rigid_sources.get(key) == source for key, source in retention_revision.items()},
        neutral_spring_eye_source_mm=[start.tolist(), neutral_end.tolist()], samples=samples, trial=trial,
        parent_spring_relief=parent_trial,
        generator=provenance(ROOT / "tools/audit_whole_arm_k1_combined_space.py"),
        scope="Eight-wire relieved-cap K1 trajectory with a geometrically prescribed return spring and eight TPU meshes. "
              "Spring checked at every outward sample; interpolated elbow TPU checked near every 5 degrees and at the end. "
              "Shoulder and girdle joints are parked. Shoulder cuff TPU uses unloaded neutral geometry. The original return "
              "spring pins are retained and the donor dog-screw revision is absent. Contacts are recorded at 0.001 mm3. "
              "The return visits the exact same prescribed pose set. No physical elastic equilibrium or spring retention validation."))
    print(f"K1: {len(samples)} spring poses and {len(tpu_indices)} lattice poses", flush=True)


if __name__ == "__main__":
    main()
