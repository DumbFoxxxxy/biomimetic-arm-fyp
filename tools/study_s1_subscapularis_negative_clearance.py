"""Measure local relief for the first negative scapular motion."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_AC_blind_coupled_motion.json"
    coupled = read(source)
    closure = load_reference(coupled["closure_audit"])
    design = load_reference(closure["source"])
    access = load_reference(load_reference(design["source"])["source"])
    entry = load_reference(load_reference(access["source"])["source"])
    cap = load_reference(entry["source"])
    blind = load_reference(cap["source"])
    split = load_reference(load_reference(load_reference(blind["parent"])["split_audit"])["split"])
    removable = load_reference(load_reference(load_reference(split["parent"])["parent"])["source"])
    assembly = load_reference(load_reference(removable["integrated_sleeve"])["mount"])
    route = load_reference(assembly["FJ1504_cap_revision"])
    cap_key = "cap:FJ1504/origin"
    bolt_key = "origin:FJ1504/bolt_0"
    bone_key = "bone:scapula"
    moving = {key: load_component(entry["sources"][key]) for key in (cap_key, bolt_key)}
    moving[bone_key] = load_component(design["selected"]["meshes"][bone_key])
    ribs = bounded({key: load_component(meta) for key, meta in entry["sources"].items()
                    if key.startswith("thorax_reference:")})
    pivot = np.asarray(closure["pivot_source_mm"])
    axis = np.array([0., 1., 0.])
    wire = A.D.outer_envelope(np.array(route["selected"]["local_wire_path_source_mm"]),
                              route["bore_radius_mm"] + .55)
    washer = load_component(entry["sources"]["origin:FJ1504/front_washer_0"])
    verts = np.asarray(moving[bolt_key].to_mesh64().vert_properties)[:, :3]
    normal = np.linalg.svd(verts - verts.mean(axis=0), full_matrices=False)[2][0]
    margin = A.I._manifold().Manifold.cube([.1] * 3, True)
    result = []
    for target in (-.25, -.5, -1., -2.):
        offsets = np.linspace(0., target, round(abs(target) / .05) + 1)
        contacts = {key: [] for key in moving}
        sample_hits = []
        for angle in offsets:
            pose = A.I.rotation_about(pivot, axis, angle)
            undo = A.I.rotation_about(pivot, axis, -angle)
            row = {}
            for key, solid in moving.items():
                posed = solid.transform(pose)
                hits = overlaps(posed, ribs, threshold=1e-6)
                if hits:
                    row[key] = hits
                    for rib in hits:
                        contact = (posed ^ ribs[rib][0]).transform(undo)
                        contacts[key].append(contact.minkowski_sum(margin))
            sample_hits.append(dict(angle_deg=float(angle), hits_mm3=row))
        reliefs = {}
        revised = {}
        for key, parts in contacts.items():
            if not parts:
                continue
            cut = moving[key] ^ union(parts)
            remaining = moving[key] - cut
            revised[key] = remaining.as_original().simplify(1e-6)
            reliefs[key] = dict(
                removed_mm3=float(cut.volume()), bounds_source_mm=list(cut.bounding_box()),
                wire_wall_guard_loss_mm3=float((cut ^ wire).volume()),
                washer_proximity_loss_mm3=float((cut ^ union(
                    washer.translate(normal * shift) for shift in (-.25, 0., .25))).volume()),
                surviving_components=sum(p.volume() > 1e-6 for p in remaining.decompose()),
                remaining_rib_hits=[
                    dict(angle_deg=float(angle), hits_mm3=found)
                    for angle in offsets
                    if (found := overlaps(remaining.transform(A.I.rotation_about(pivot, axis, angle)), ribs))
                ],
            )
        case = dict(target_deg=target, samples=sample_hits, relief=reliefs)
        if target == -.25 and revised and all(
                row["wire_wall_guard_loss_mm3"] < .001 and row["washer_proximity_loss_mm3"] < .001
                and row["surviving_components"] == 1 and not row["remaining_rib_hits"]
                for row in reliefs.values()):
            meshes = {}
            for key, solid in revised.items():
                name = key.split(":")[-1].replace("/", "_")
                actual, meshes[key] = export(solid, folder / f"FJ1504_negative_relief_{name}.obj", pivot)
                for angle in offsets:
                    posed = actual.transform(A.I.rotation_about(pivot, axis, angle))
                    if overlaps(posed, ribs):
                        raise ValueError("Local relief changed after export")
            case["local_relief_meshes"] = meshes
        result.append(case)
        print(target, {key: (round(row["removed_mm3"], 4), round(row["wire_wall_guard_loss_mm3"], 4),
                             round(row["washer_proximity_loss_mm3"], 4), len(row["remaining_rib_hits"]))
                       for key, row in reliefs.items()}, flush=True)
    save(folder / "FJ1504_negative_rib_clearance.json", dict(
        generator=provenance(ROOT / "tools/study_s1_subscapularis_negative_clearance.py"),
        source=provenance(source), axis_source=axis.tolist(), pivot_source_mm=pivot.tolist(), cases=result,
        scope="Pure AC rotation with clavicle fixed. Three selected solids against all thorax "
        "reference meshes, 0.05 degree samples. Contact cuts are numerical witnesses, "
        "with a 0.10 mm box added before intersection. Only the first interval is exported "
        "when the local protection and readback checks pass.",
    ))


if __name__ == "__main__":
    main()
