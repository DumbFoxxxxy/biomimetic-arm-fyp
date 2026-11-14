"""Distinguish seated contact from obstructions to inserting the FJ1500 pair."""

from itertools import product

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    input_file = folder / "FJ1500_feed_installation_screen.json"
    previous = read(input_file)
    feed = load_reference(previous["source"])
    assembly = load_reference(feed["mount"])
    route = load_reference(feed["route"])["exported"]["trial"]["new_path_source_mm"]
    cap = load_component(feed["selected"]["mesh"])
    tpu = load_component(feed["selected"]["recessed_TPU"])
    original_targets = {key: load_component(src) for key, src in assembly["sources"].items()
                        if key != "cap:FJ1500/origin" and not key.startswith("origin:FJ1500/")}
    baseline = overlaps(tpu, bounded(original_targets), threshold=1e-6)
    scapula = original_targets["bone:scapula"]
    rest_contact = tpu ^ scapula
    protected = tpu ^ A.D.outer_envelope(np.array(route), 1.31)
    mf = A.I._manifold()
    relief_trials = []
    relieved = tpu
    for pad in (0., .005, .02):
        stock = rest_contact if pad == 0. else rest_contact.minkowski_sum(mf.Manifold.cube([2.*pad]*3,True))
        candidate = (tpu-(stock-protected)).as_original().simplify(1e-6)
        components = sorted((float(piece.volume()) for piece in candidate.decompose() if piece.volume()>1e-8),reverse=True)
        relief_trials.append(dict(padding_mm=pad,components_mm3=components,remaining_bone_overlap_mm3=float((candidate^scapula).volume())))
        if len(components)==1 and relief_trials[-1]["remaining_bone_overlap_mm3"]<=.001:
            relieved = candidate
            break
    # Defer the neighbouring fasteners only for this insertion screen.
    targets = bounded({key: value for key,value in original_targets.items() if not key.startswith("origin:FJ1504/")})
    centre = np.array(feed["start_source_mm"])
    normal = -np.array(feed["bore_axis_source"])
    candidates = [np.array(q["direction"]) for q in previous["direction_search"] if q["FJ1504_origin_hardware_later"]]

    def check(direction, tilt, rotation_axis, samples):
        hits = []
        for distance in samples:
            turn = Rotation.from_rotvec(rotation_axis*np.deg2rad(tilt*distance/8.))
            pose = rigid_transform(turn,centre,centre+distance*direction)
            cap_hits = overlaps(cap.transform(pose),targets)
            tpu_hits = overlaps(relieved.transform(pose),targets)
            if distance == 0. and tpu_hits.get("bone:scapula",0.) <= baseline.get("bone:scapula",0.) + 1e-8:
                tpu_hits.pop("bone:scapula",None)
            if cap_hits or tpu_hits:
                hits.append(dict(distance_mm=float(distance),cap_hits=cap_hits,TPU_hits=tpu_hits))
        return hits

    straight = []
    for direction in candidates:
        hit = check(direction,0.,normal,(.4,2.,4.,8.))
        cost = max((max([*q["cap_hits"].values(),*q["TPU_hits"].values()]) for q in hit), default=0.)
        straight.append((cost,direction))
    straight.sort(key=lambda item:item[0])
    trials = []
    for (_, direction), tilt, basis in product(straight[:10],(0.,-10.,10.,-20.,20.),np.eye(3)):
        if tilt == 0. and basis[0] != 1.:
            continue
        hits = check(direction,tilt,basis,(.2,.4,1.,2.,4.,6.,8.))
        score = max((max([*q["cap_hits"].values(),*q["TPU_hits"].values()]) for q in hits), default=0.)
        trials.append(dict(direction=direction.tolist(), tilt_deg=tilt, rotation_axis=basis.tolist(),
                           max_intersection_mm3=score, coarse_hits=hits))
    for row in sorted(trials,key=lambda r:r["max_intersection_mm3"])[:5]:
        row["dense_hits"] = check(np.array(row["direction"]),row["tilt_deg"],np.array(row["rotation_axis"]),np.linspace(0.,8.,81))
    points, _ = mesh_geom.manifold_mesh_to_arrays(rest_contact.to_mesh64())
    save(folder/"FJ1500_tilted_installation_screen.json",dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_tilted_installation.py"), parent=provenance(input_file),
        seated_original_TPU_contacts_mm3=baseline,
        seated_scapula_contact_bounds_mm=[points.min(axis=0).tolist(),points.max(axis=0).tolist()] if len(points) else None,
        seated_relief_removed_mm3=float((tpu-relieved).volume()), protected_loss_mm3=float((protected-relieved).volume()), relief_trials=relief_trials,
        seated_relief_TPU_contacts_mm3=overlaps(relieved,bounded(original_targets),threshold=1e-6),
        straight_directions_checked=len(straight), trials=trials,
        scope="Virtual seated micro-relief and prescribed tilt decreasing linearly toward the seat."
              " No material exported. Neighbouring FJ1504 fastener installation and TPU deformation are excluded."))
    print('Seated before',baseline,'after',overlaps(relieved,targets,threshold=1e-6),
          'best',[(r['max_intersection_mm3'],len(r['dense_hits'])) for r in trials if 'dense_hits' in r],flush=True)


if __name__ == "__main__":
    main()
