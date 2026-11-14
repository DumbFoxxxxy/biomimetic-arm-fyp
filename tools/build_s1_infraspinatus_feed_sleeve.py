"""Fit a supported FJ1500 feed sleeve on the same-side scapular mount."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    route_file = folder / "FJ1500_proximal_rebore_screen.json"
    route = read(route_file)["exported"]
    mount_file = DEST.parent / "integration_review/s1_series_elbow/assembly_FJ1500_front_mount.json"
    mount = read(mount_file)
    path = np.array(route["trial"]["new_path_source_mm"])
    entry, exit = path[1:3]
    length = float(np.linalg.norm(exit-entry))
    axis = (exit-entry)/length
    cap_key = "cap:FJ1500/origin"
    cap = load_component(route["parts"]["cap"])
    tpu = load_component(route["parts"]["TPU"])
    original_cap = load_component(mount["sources"][cap_key])
    obstacles = bounded({k: load_component(src) for k, src in mount["sources"].items() if k != cap_key})
    wire = A.D.outer_envelope(path,.30)
    bore = rod(entry-2.*axis,exit+.02*axis,.86,64)
    mf = A.I._manifold()
    fit_pad = mf.Manifold.cube([.20]*3, True)
    fastener_guard = union(part.minkowski_sum(mf.Manifold.cube([.6]*3,True))
                           for key,(part,_) in obstacles.items() if key.startswith("origin:FJ1500/")) ^ tpu
    baseline_hits = overlaps(wire,obstacles)
    trials = []
    selected = None
    for radius, front, collar in product((1.41,1.61,1.86), (0., .75, 1.5), (0., .55, 1.)):
        stock = rod(entry-front*axis,exit-.10*axis,radius,64)
        if collar:
            stock = union([stock,rod(entry-front*axis,entry+(1.5-front)*axis,radius+collar,64)])
        sleeve = stock-bore
        bond = float((sleeve^cap).volume())
        actual = ((cap+sleeve)-bore).as_original().simplify(1e-6)
        added = actual-original_cap
        hits = overlaps(added,obstacles)
        tpu_overlap = float((added^tpu).volume())
        wire_overlap = float((wire^actual).volume())
        pieces = sum(s.volume()>1e-5 for s in actual.decompose())
        trial = dict(outer_radius_mm=radius, proximal_projection_mm=front, flange_extra_radius_mm=collar,
                     fusion_volume_mm3=bond, added_mm3=float(added.volume()), rigid_hits_mm3=hits,
                     added_stock_TPU_overlap_mm3=tpu_overlap, nominal_wire_overlap_mm3=wire_overlap, components=pieces)
        trials.append(trial)
        if hits or wire_overlap>.001 or bond<1. or pieces!=1 or selected is not None or baseline_hits:
            continue
        recess = added.minkowski_sum(fit_pad)
        recessed = (tpu-recess).as_original().simplify(1e-6)
        guard_loss = float((fastener_guard-recessed).volume())
        connected = sum(s.volume()>1e-5 for s in recessed.decompose()) == 1
        trial.update(TPU_recess_removed_mm3=float((tpu-recessed).volume()),
                     TPU_guard_loss_mm3=guard_loss, recessed_TPU_connected=connected)
        if connected and guard_loss<=.001 and (recessed^actual).volume() <= (tpu^cap).volume()+.001:
            reread, meta = export(actual, folder/"FJ1500_supported_feed_cap.obj", entry)
            flexible, flex_meta = export(recessed, folder/"FJ1500_supported_feed_TPU.obj", entry)
            if (overlaps(reread-original_cap,obstacles) or (reread^wire).volume()>.001
                    or ((reread-original_cap)^flexible).volume()>.001 or (wire^flexible).volume()>.001):
                raise ValueError("Candidate fit changed after readback")
            selected = dict(**trial, mesh=meta, recessed_TPU=flex_meta)
    installation = []
    if selected is not None:
        cap = load_component(selected["mesh"])
        tpu = load_component(selected["recessed_TPU"])
        present = {k:v for k,v in obstacles.items() if not k.startswith("origin:FJ1500/")}
        directions = {"feed_axis":-axis, "reverse_feed_axis":axis}
        for label, direction in directions.items():
            hits = []
            for distance in np.linspace(8.,0.,41):
                current = cap.translate(distance*direction)
                rigid_hits = overlaps(current,present)
                flexible_overlap = float((current^tpu).volume())
                if rigid_hits or flexible_overlap>.001:
                    hits.append(dict(distance_mm=float(distance),rigid_hits_mm3=rigid_hits,TPU_overlap_mm3=flexible_overlap))
            installation.append(dict(direction=label,vector=direction.tolist(),hits=hits))
    save(folder/"FJ1500_feed_sleeve_screen.json",dict(
        generator=provenance(ROOT/"tools/build_s1_infraspinatus_feed_sleeve.py"), route=provenance(route_file),
        mount=provenance(mount_file), bore_inner_radius_mm=.86, start_source_mm=entry.tolist(),
        bore_axis_source=axis.tolist(), span_mm=length, nominal_wire_hits_same_side_mount_mm3=baseline_hits,
        trials=trials, selected=selected, local_approach_trials=installation, scope="Integral annular sleeve fused to original cap."
        " Added stock screened at rest. Local rigid approach with the TPU present and fasteners removed. No load qualification."))
    print("Wire on same-side mount",baseline_hits,"selected",selected is not None,flush=True)


if __name__ == "__main__":
    main()
