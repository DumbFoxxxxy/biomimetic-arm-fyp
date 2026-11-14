"""Screen a rotated and translated FJ1500 proximal feed bore and TPU transition."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_infraspinatus_smooth_tendon import centreline
from study_s1_infraspinatus_wire_travel import route_metrics
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    parent_file = folder / "FJ1500_port_transition_screen.json"
    parent = read(parent_file)
    wire = load_reference(parent["wire_baseline"])
    rest = load_reference(load_reference(wire["deformation"])["rest_geometry"])
    prior = load_reference(load_reference(rest["bias_screen"])["reference_equilibrium"])
    assembly = load_reference(prior["shoulder_assembly"])
    path = np.array(rest["rest_path_source_mm"])
    tpu = load_component(rest["candidate"])
    key = "cap:FJ1500/origin"
    cap = load_component(assembly["sources"][key])
    parts = {k: load_component(src) for k, src in assembly["sources"].items()
             if k.startswith(("bone:", "GH_cap_", "origin:FJ1500/"))}
    obstacles = bounded(parts)
    bolts = [v for k, v in parts.items() if any(s in k for s in ("bolt", "washer", "nut"))]
    mf = A.I._manifold()
    entry = path[1]
    normal = A.D.unit(path[2]-entry)
    inward = A.D.unit(path[23]-path[2])
    lateral = A.D.unit(inward-normal*(inward@normal))
    cross = np.cross(normal,lateral)
    checked = []
    successes = []
    for lateral_mm, cross_mm, tilt_deg, join in product((-6.,-3.,0.,3.,6.),(-4.,0.,4.),(20.,40.,60.,75.),(18,22,26)):
        angle = np.deg2rad(tilt_deg)
        axis = np.cos(angle)*normal+np.sin(angle)*lateral
        start = entry+lateral_mm*lateral+cross_mm*cross
        inner = start+3.25/np.cos(angle)*axis
        end = path[join]
        chord = np.linalg.norm(end-inner)
        forward = A.D.unit(path[join+1]-end)
        curve = centreline([inner,inner+chord*.36*axis,end-chord*.36*forward,end],49)
        transition = np.vstack([start-5.*axis,start,curve,path[join+1]])
        _, turns, radii = route_metrics(transition)
        radius = float(radii.min())
        wire_solid = A.D.outer_envelope(transition,.30)
        hit = overlaps(wire_solid,obstacles)
        bore = rod(start-1.4*axis,inner+1.4*axis,.86)
        cap_cut = cap^bore
        bearing = [float((bore^bolt).volume()) for bolt in bolts]
        row = dict(lateral_mm=lateral_mm,cross_mm=cross_mm,tilt_deg=tilt_deg,join_index=join,
                   min_radius_mm=radius,turn_deg=float(np.rad2deg(turns.sum())), rigid_wire_hits_mm3=hit,
                   cap_bore_removed_mm3=float(cap_cut.volume()), protected_fastener_overlap_mm3=max(bearing,default=0.))
        checked.append(row)
        if not hit and row["cap_bore_removed_mm3"]>1. and row["protected_fastener_overlap_mm3"]<=.001:
            successes.append((row, transition, bore))
    verified = []
    meshes = None
    for row, transition, bore in sorted(successes,key=lambda q:-q[0]["min_radius_mm"])[:10]:
        path_new = np.vstack([transition[:-1],path[row["join_index"]+1:]])
        # Recut the lead as well as the curve so the old and new channel
        # cannot leave a lip at the pad seam.
        lumen = A.D.outer_envelope(transition,.86)
        modified_cap = (cap-bore-lumen).as_original().simplify(1e-6)
        modified_tpu = (tpu-lumen).as_original().simplify(1e-6)
        guard = union(v.minkowski_sum(mf.Manifold.cube([.6]*3,True)) for v in bolts)^tpu
        wire_solid = A.D.outer_envelope(path_new,.30)
        components = [sum(q.volume()>1e-5 for q in part.decompose()) for part in (modified_cap,modified_tpu)]
        volumes = [float((part^wire_solid).volume()) for part in (modified_cap,modified_tpu)]
        loss = float((guard-modified_tpu).volume())
        trial = dict(**row, cap_and_TPU_components=components, cap_and_TPU_wire_overlap_mm3=volumes,
                     TPU_fastener_guard_loss_mm3=loss, new_path_source_mm=path_new.tolist())
        bore_axis = A.D.unit(path_new[2]-path_new[1])
        station = np.linspace(0., np.linalg.norm(path_new[2]-path_new[1]), 31)
        rings = []
        for distance in station:
            p = path_new[1]+distance*bore_axis
            ring = rod(p-.05*bore_axis,p+.05*bore_axis,1.41)-rod(p-.1*bore_axis,p+.1*bore_axis,.86)
            rings.append(float((ring^modified_cap).volume()/ring.volume()))
        trial["cap_annular_stock_fraction"] = rings
        verified.append(trial)
        if meshes is None and components==[1,1] and max(volumes+[loss])<=.001:
            output = {}
            for name, part in (("cap", modified_cap),("TPU",modified_tpu)):
                actual, metadata = export(part,folder/f"FJ1500_proximal_rebore_{name}.obj",entry)
                if (actual^wire_solid).volume()>.001:
                    raise ValueError("Wire clearance did not reproduce after export")
                output[name] = metadata
            meshes = dict(trial=trial, parts=output)
    save(folder/"FJ1500_proximal_rebore_screen.json",dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_proximal_rebore.py"), parent=provenance(parent_file),
        trials=checked, initial_candidates=len(successes), recut_trials=verified, exported=meshes,
        scope="Static, original scapular reference. Tilted rigid bore and tangent-continuous path before entering existing channel."
              " Checks do not cover loaded deformation, complete TPU walls, modified front-mount version or retention load."))
    print(len(checked),len(successes),len(verified),bool(meshes),flush=True)


if __name__ == "__main__":
    main()
