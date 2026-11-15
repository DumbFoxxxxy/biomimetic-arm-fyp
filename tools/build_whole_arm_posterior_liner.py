"""Rebuild the posterior channel in the existing deltoid lattice."""

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import export_lattice
from audit_s1_infraspinatus_distributed_deformation import blend
from audit_s1_series_elbow_integration import transforms
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, save
from project_runtime import load_mujoco_model
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    audit_ref = provenance(folder / "whole_arm_posterior_detour_audit.json")
    audit = load_reference(audit_ref)
    scene = load_reference(audit["scene"])
    motion = load_reference(audit["motion"])
    transfer = load_reference(scene["transfer"])
    local = DEST.parent / "s1_local_muscle_load"
    historical_ref = provenance(local / "FJ1467_distal_guide_TPU.json")
    historical = load_reference(historical_ref)
    parent, _, _ = indexed_solid(historical["parent_TPU"])
    route = audit["routes"]["FJ1467"]
    path = np.asarray(route["points_source_mm"])
    first, last = audit["posterior_detour"]["edited_node_indices"]
    section = slice(first - 2, last + 3)
    centerline = path[section]
    envelope = A.D.outer_envelope
    shell = envelope(centerline, .95)
    lumen = envelope(centerline, .52)
    # Preserve the existing shared belly. Rebuild only the edited liner and bore.
    candidate = ((parent + shell) - lumen).as_original().simplify(1e-7)
    if (parent ^ (shell - lumen)).volume() <= .001:
        raise ValueError("The new liner does not bond to the original lattice")
    exported, mesh = export_lattice(candidate, folder / "deltoid_posterior_detour_TPU.obj", path[first])
    lost = parent - exported
    neighbors = {}
    for key in ("FJ1468", "FJ1513"):
        points = np.asarray(audit["routes"][key]["points_source_mm"])
        guard = envelope(points, 1.10)
        neighbors[key] = dict(removed_in_guard_mm3=A.I.intersection_volume(lost, guard),
                              new_liner_wire_overlap_mm3=A.I.intersection_volume(exported - parent, envelope(points, .30)))

    section_zone = exported ^ envelope(centerline, 1.14)
    mf = A.I._manifold()
    zone_mesh = section_zone.to_mesh64()
    vertices = np.asarray(zone_mesh.vert_properties)[:, :3]
    faces = np.array(zone_mesh.tri_verts, dtype=np.uint64, order="C", copy=True)
    front = np.asarray(audit["routes"]["FJ1468"]["points_source_mm"])
    axis = front[-1] - front[0]
    axis /= np.linalg.norm(axis)
    window = load_reference(route["guide_motion_reference"])["deformation_window_mm"]
    vertex_weights = blend((vertices - front[0]) @ axis, *window, .1)[0]
    # Near the edited span, the taper at the anchored proximal terminal has
    # ended. Interpolate the line using the exact saved per-node weights.
    wire_weights = np.asarray(route["distal_weights"])[section]
    parts = {key: load_component(source) for key, source in
             {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"]}.items()}
    owners = transfer["owner_by_source"]
    origins = {k: np.asarray(scene["neutral_body_origins_mm"][k]) for k in set(owners.values())}
    model = load_mujoco_model(ROOT / scene["model"]["file"])
    data = mujoco.MjData(model)
    samples = []
    wire_audit = audit["paths"]["K7_SC1_pos"]
    accepted_count = min(row["initial_clear_samples"] for row in wire_audit["summary"].values())
    for row in motion["paths"]["K7_SC1_pos"]["accepted"][:accepted_count]:
        data.qpos[:] = row["qpos"]
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        a, b = (frames[owner] for owner in route["owners"])
        def deform(points, weights):
            return ((1. - weights[:, None]) * (points @ a[:, :3].T + a[:, 3])
                    + weights[:, None] * (points @ b[:, :3].T + b[:, 3]))
        posed = mf.Manifold(mf.Mesh64(np.ascontiguousarray(deform(vertices, vertex_weights)), faces))
        if posed.status() != mf.Error.NoError:
            raise ValueError(f"Invalid posed liner at SC {row['input_deg']}")
        obstacles = bounded({key: part.transform(frames[owners[key]]) for key, part in parts.items()})
        wire = envelope(deform(centerline, wire_weights), .30)
        samples.append(dict(input_deg=row["input_deg"], local_TPU_rigid_hits_mm3=overlaps(posed, obstacles),
                            local_TPU_wire_overlap_mm3=A.I.intersection_volume(posed, wire)))
    save(folder / "whole_arm_posterior_liner.json", dict(
        audit=audit_ref, historical_guide=historical_ref, source=mesh,
        parent=historical["parent_TPU"], edited_with_transitions_node_indices=[section.start, section.stop - 1],
        outer_radius_mm=.95, lumen_radius_mm=.52, existing_neighbor_checks=neighbors, samples=samples,
        generator=provenance(ROOT / "tools/build_whole_arm_posterior_liner.py"),
        scope="Connected neutral lattice and local geometric interpolation at the posterior channel. "
              "Includes sampled poses before the first wire limit. No load, fatigue, or full-belly deformation validation."))
    print(f"Exported connected deltoid mesh; {len(samples)} local poses", flush=True)


if __name__ == "__main__":
    main()
