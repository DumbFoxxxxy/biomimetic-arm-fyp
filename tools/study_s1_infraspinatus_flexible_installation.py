"""Screen cap-first FJ1500 installation with a temporarily bent TPU belly."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_feed_sleeve_screen.json"
    feed = read(source)
    assembly = load_reference(feed["mount"])
    route = load_reference(feed["route"])["exported"]["trial"]
    path = np.array(route["new_path_source_mm"])
    cap = load_component(feed["selected"]["mesh"])
    tpu = load_component(feed["selected"]["recessed_TPU"])
    points, triangles = mesh_geom.manifold_mesh_to_arrays(tpu.to_mesh64())
    triangles = np.ascontiguousarray(triangles, dtype=np.uint64)
    mf = A.I._manifold()
    fixed = bounded({key: load_component(src) for key, src in assembly["sources"].items()
                     if key != "cap:FJ1500/origin" and not key.startswith("origin:FJ1504/")})
    fixed.update(bounded({"cap:FJ1500/origin": cap}))
    axis = np.array(feed["bore_axis_source"])
    origin = path[2]
    span_axis = A.D.unit(path[-1]-origin)
    tangent = A.D.unit(np.cross(axis,span_axis))
    outward = A.D.unit(np.cross(span_axis,tangent))
    nominal = overlaps(tpu, fixed, threshold=1e-6)

    def weights(vertices, reach, start=0.):
        s = np.clip(((vertices-origin)@span_axis-start)/reach,0.,1.)
        return s*s*(3.-2.*s), 6.*s*(1.-s)/reach

    def check(direction, lift, reach, fractions, with_wire=False):
        w,dw = weights(points,reach,8.)
        wire_weight,_ = weights(path,reach,8.)
        rows = []
        for s in fractions:
            # Move the entry pad along the sleeve. The distal end remains
            # unconnected during insertion and can be held away by hand.
            translation = 8.*s*axis
            bend = lift*np.sin(s*np.pi/2.)*direction
            jacobian = 1.+dw*(bend@span_axis)
            placed = points+translation+w[:,None]*bend
            body = mf.Manifold(mf.Mesh64(np.ascontiguousarray(placed),triangles))
            if body.status()!=mf.Error.NoError:
                raise ValueError("Invalid assembly pose mesh")
            hits = overlaps(body,fixed)
            if s == 0.:
                hits = {k:v for k,v in hits.items() if v>nominal.get(k,0.)+1e-8}
            result = dict(entry_offset_mm=float(8.*s),hits_mm3=hits,
                          minimum_sampled_jacobian=float(jacobian.min()))
            if with_wire:
                # This line moves with the TPU only. During assembly the
                # wire is unthreaded; this checks the channel geometry.
                moved_path = path+translation+wire_weight[:,None]*bend
                result["nominal_wire_TPU_overlap_mm3"] = float((body^A.D.outer_envelope(moved_path,.30)).volume())
            rows.append(result)
        return rows

    trials = []
    for angle,lift,reach in product(np.arange(0.,360.,45.),(4.,10.,18.),(35.,70.)):
        direction = np.cos(np.deg2rad(angle))*outward+np.sin(np.deg2rad(angle))*tangent
        rows = check(direction,lift,reach,(.05,.25,.50,.75,1.))
        score = max((max(q["hits_mm3"].values(),default=0.) for q in rows))
        trials.append(dict(angle_deg=float(angle),direction_source=direction.tolist(),tip_lift_mm=lift,
                           blend_length_mm=reach,max_overlap_mm3=score,coarse=rows))
    for trial in sorted(trials,key=lambda q:q["max_overlap_mm3"])[:4]:
        trial["fine"] = check(np.array(trial["direction_source"]),trial["tip_lift_mm"],
                              trial["blend_length_mm"],np.linspace(0.,1.,41),with_wire=True)
    save(folder/"FJ1500_flexible_installation_screen.json",dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_flexible_installation.py"), feed=provenance(source),
        baseline_TPU_contacts_mm3=nominal, trials=trials,
        scope="Temporary TPU bend prescribed during insertion, with FJ1500 cap and fasteners fixed."
              " FJ1504 origin hardware deferred, distal tendon disconnected, wire unthreaded."
              " Nominal centreline checked inside deformed TPU only. No force, damage or measured deflection claim."))
    print('Best',[(r['angle_deg'],r['tip_lift_mm'],r['max_overlap_mm3']) for r in sorted(trials,key=lambda q:q['max_overlap_mm3'])[:4]],flush=True)


if __name__ == "__main__":
    main()
