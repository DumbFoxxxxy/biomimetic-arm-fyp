"""Map the best FJ1500 side approaches back to the proximal TPU pad."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_side_entry_screen.json"
    screen = read(source)
    sleeve = load_reference(screen["source"])
    feed = load_reference(sleeve["integrated_sleeve"])
    assembly = load_reference(feed["mount"])
    route = np.array(load_reference(feed["route"])["exported"]["trial"]["new_path_source_mm"])
    tpu = load_component(feed["selected"]["recessed_TPU"])
    solids = {k:load_component(v) for k,v in assembly["sources"].items()
              if k != "cap:FJ1500/origin" and not k.startswith("origin:FJ1504/")}
    solids["cap:FJ1500/origin"] = load_component(sleeve["selected"]["cap"])
    targets = bounded(solids)
    axis = A.D.unit(route[-1]-route[2])
    channel_guard = A.D.outer_envelope(route,1.31)^tpu
    mf = A.I._manifold()
    attachment = union(v.minkowski_sum(mf.Manifold.cube([.6]*3,True)) for k,v in solids.items()
                       if k.startswith("origin:FJ1500/"))^tpu
    protected = channel_guard+attachment
    rows = []
    paths = sorted(screen["straight_traces"],key=lambda t:sum(sum(s["hits_mm3"].values()) for s in t["samples"]))[:3]
    for trace in paths:
        direction = np.array(trace["direction_source"])
        cuts = []
        bounds = []
        for d in np.linspace(.1,2.,20):
            for key in overlaps(tpu.translate(d*direction),targets):
                contact = tpu^solids[key].translate(-d*direction)
                if contact.volume()>.001:
                    cuts.append(contact)
                    points,_ = mesh_geom.manifold_mesh_to_arrays(contact.to_mesh64())
                    bounds.extend(((points-route[2])@axis).tolist())
        swept = union(cuts)
        candidate = (tpu-(swept-protected)).as_original().simplify(1e-6)
        residual = []
        for d in np.linspace(.05,2.,40):
            hits = overlaps(candidate.translate(d*direction),targets)
            if hits:
                residual.append(dict(offset_mm=float(d),hits_mm3=hits))
        rows.append(dict(direction_source=direction.tolist(),contact_station_range_mm=[min(bounds),max(bounds)],
                         channel_guard_overlap_mm3=float((swept^channel_guard).volume()),
                         fastener_guard_overlap_mm3=float((swept^attachment).volume()),
                         removed_mm3=float((tpu-candidate).volume()),
                         connected_parts=sum(q.volume()>1e-5 for q in candidate.decompose()),residual=residual))
    save(folder/"FJ1500_side_entry_relief_screen.json",dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_side_entry_relief.py"), parent=provenance(source),
        checks=rows, scope="Virtual relief on final 2 mm of three straight approaches. Protects wire and fastener"
        " envelopes. No modified geometry exported; full approach and TPU mechanics are unverified."))
    for row in rows:
        print(row["contact_station_range_mm"],row["channel_guard_overlap_mm3"],row["fastener_guard_overlap_mm3"],
              row["removed_mm3"],row["connected_parts"],len(row["residual"]),flush=True)


if __name__ == "__main__":
    main()
