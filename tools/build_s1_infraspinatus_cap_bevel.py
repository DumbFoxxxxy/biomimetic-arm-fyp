"""Find a local FJ1500 cap bevel that preserves functional material."""

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import export, load_component, mesh_geom, save
from study_s1_infraspinatus_cap_relief import CAP, HEAD, guard_solids, inverse, union
from study_s1_local_force_geometry import load_reference


def vertices(part):
    return mesh_geom.manifold_mesh_to_arrays(part.to_mesh64())[0]


def summarize(rows, cap):
    limits = {}
    for load in (0., 5., 10.):
        last, first = None, None
        for row, patches, other in rows:
            if row["input_N"] != load:
                continue
            hits = list(other)
            for obstruction, patch in patches:
                volume = float((cap ^ patch).volume())
                if volume > .001:
                    hits.append(dict(moving=CAP, obstruction=obstruction, overlap_mm3=volume))
            if hits:
                first = dict(angle_deg=row["angle_deg"], hits=sorted(hits, key=lambda hit: -hit["overlap_mm3"]))
                break
            last = row["angle_deg"]
        limits[str(int(load))] = dict(last_clear_deg=last, first_collision=first)
    return limits


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1500_flexion_path.json"
    path = read(path_file)
    rest = load_reference(path["rest_geometry"])
    screen = load_reference(rest["bias_screen"])
    equilibrium = load_reference(screen["reference_equilibrium"])
    assembly = load_reference(equilibrium["shoulder_assembly"])
    geometry = load_reference(equilibrium["local_geometry"])
    model = next(row for row in load_reference(equilibrium["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    axis = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    centre = np.array(geometry["gh_centre_source_mm"])
    pin = np.array(model["neutral_cap_pin_source_mm"])
    original = load_component(assembly["sources"][CAP])
    guards = {name: original ^ part for name, part in guard_solids(model, route).items()}
    guarded = vertices(union(guards.values())) - pin
    obstacles, rows = {}, []
    for row in path["outward"]:
        patches, other = [], []
        for hit in row["hits"]:
            if CAP not in (hit["moving"], hit["obstruction"]):
                other.append(hit)
                continue
            obstruction = hit["obstruction"] if hit["moving"] == CAP else hit["moving"]
            if obstruction not in obstacles:
                obstacles[obstruction] = load_component(assembly["sources"][obstruction])
            body = obstacles[obstruction]
            owner = assembly["intended_body_owners"][obstruction]
            if owner == "humerus_ideal_attachment":
                gh = Rotation.from_rotvec(axis * np.deg2rad(row["angle_deg"]))
                body = body.transform(rigid_transform(gh, centre, centre))
            elif owner not in ("fixed_context", "scapula_ideal_fastened"):
                raise ValueError(f"Unmapped cap obstruction: {obstruction} ({owner})")
            cap = Rotation.from_rotvec(np.deg2rad(row["cap_rotation_vector_deg"]))
            transform = rigid_transform(cap, pin, np.array(row["cap_pin_source_mm"]))
            patch = original ^ body.transform(inverse(transform))
            if abs(patch.volume() - hit["overlap_mm3"]) > .0001:
                raise ValueError(f"Mapped interference does not match {obstruction}")
            patches.append((obstruction, patch))
        rows.append((row, patches, other))

    # Fibonacci directions give deterministic coverage of possible bevel normals.
    stations = np.arange(768)
    z = 1. - 2. * (stations + .5) / len(stations)
    phase = stations * np.pi * (3. - np.sqrt(5.))
    radial = np.sqrt(1. - z**2)
    directions = np.c_[radial * np.cos(phase), radial * np.sin(phase), z]
    floors = (guarded @ directions.T).max(axis=0) + .03
    target_patches = [vertices(patch) - pin for row, patches, _ in rows
                      if row["input_N"] == 10. and row["angle_deg"] <= 72.
                      for _, patch in patches]
    limits = [(patch @ directions.T).min(axis=0) - .05 for patch in target_patches]
    candidates = []
    for idx, normal in enumerate(directions):
        cut = floors[idx]
        covered = 0
        for limit in limits:
            if limit[idx] < floors[idx]:
                break
            cut = limit[idx] if covered == 0 else min(cut, limit[idx])
            covered += 1
        candidates.append((covered, idx, cut))
    max_covered = max(row[0] for row in candidates)
    if max_covered == 0:
        bounds = []
        for angle in (64.4, 65., 66., 67., 68., 70.):
            required = union(patch for row, patches, _ in rows if row["angle_deg"] <= angle
                             for _, patch in patches)
            bounds.append(dict(angle_deg=angle, required_volume_mm3=float(required.volume()),
                               protected_volume_mm3={key: float((required ^ guard).volume())
                                                     for key, guard in guards.items()}))
        print("No plane clears first patch; protected material bounds", bounds, flush=True)
        save(folder / "FJ1500_guarded_cap_bevel.json", dict(
            generator=provenance(ROOT / "tools/build_s1_infraspinatus_cap_bevel.py"),
            source_path=provenance(path_file), original=assembly["sources"][CAP],
            trial_normals=len(directions), candidate=None, protected_bounds=bounds,
            scope="Finite plane search plus inverse-mapped cap interference volumes. No bevel selected."))
        return
    trials = []
    for covered, idx, cut in candidates:
        if covered != max_covered:
            continue
        normal = directions[idx]
        cap = original.trim_by_plane(-normal, -float(cut + pin @ normal))
        if len(cap.decompose()) != 1:
            continue
        removed = float(original.volume() - cap.volume())
        trials.append((removed, idx, cut, cap))
    if not trials:
        raise RuntimeError("No connected bevel retained")
    trials.sort(key=lambda row: row[0])
    ranked = [(summarize(rows, cap), removed, idx, cut, cap) for removed, idx, cut, cap in trials[:12]]
    limits_after, removed, idx, cut, cap = max(ranked, key=lambda row: (row[0]["10"]["last_clear_deg"],
                                                                       row[0]["5"]["last_clear_deg"], -row[1]))
    removed_material = original - cap
    lost = {key: float((removed_material ^ guard).volume()) for key, guard in guards.items()}
    if any(value > 1e-6 for value in lost.values()):
        raise RuntimeError(f"Cut entered protected stock: {lost}")
    actual, source = export(cap, folder / "FJ1500_guarded_cap_bevel.obj",
                            np.array(assembly["sources"][CAP]["translation_source_mm"]))
    verified = summarize(rows, actual)
    for load in limits_after:
        if verified[load]["last_clear_deg"] != limits_after[load]["last_clear_deg"]:
            raise RuntimeError("Bevel export changed the diagnosed range")
    report = dict(
        generator=provenance(ROOT / "tools/build_s1_infraspinatus_cap_bevel.py"),
        source_path=provenance(path_file), original=assembly["sources"][CAP], candidate=source,
        removed_mm3=removed, removed_fraction=removed / original.volume(), protected_removal_mm3=lost,
        plane_normal_source=directions[idx].tolist(), plane_origin_source_mm=(pin + cut * directions[idx]).tolist(),
        kept_plane_side="normal projection less than plane position",
        outward_limits_before={key: row["last_clear"]["angle_deg"] for key, row in path["limits"].items()},
        outward_limits_after=verified, trial_normals=len(directions),
        candidate_planes_with_max_consecutive_patch_coverage=len(trials),
        checked_candidate_count=len(ranked),
        scope="Existing sampled outward equilibrium states are replayed with a subtractive cap bevel."
              " Every recorded prior collision involving the cap is re-intersected."
              " Other contacts are carried forward unchanged. New continuous equilibrium and strength are not verified."
              " A plane is selected from a finite normal set. The stock guards are geometric design references.")
    save(folder / "FJ1500_guarded_cap_bevel.json", report)
    print("Removed", removed, "mm³; protected", lost, "limits", verified, flush=True)


if __name__ == "__main__":
    main()
