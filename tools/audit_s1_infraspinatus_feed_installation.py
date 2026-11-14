"""Compare local installation sequences for the supported FJ1500 feed cap."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def probe(moving, fixed, direction, step):
    hits = []
    for distance in np.linspace(8., 0., round(8./step)+1):
        for name, part in moving.items():
            contact = overlaps(part.translate(distance*direction), fixed)
            if contact:
                hits.append(dict(distance_mm=float(distance), part=name, overlap_mm3=contact))
    return hits


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_feed_sleeve_screen.json"
    report = read(source)
    mount = load_reference(report["mount"])
    route = load_reference(report["route"])
    wire = load_reference(load_reference(route["parent"])["wire_baseline"])
    deformation = load_reference(wire["deformation"])
    combined = load_reference(deformation["combined"])
    baseline = load_reference(load_reference(combined["exit_revision"])["baseline"])
    port = next(row["ports"]["origin"] for row in load_reference(baseline["input_route"])["routes"] if row["head"] == "FJ1500")
    cap_key = "cap:FJ1500/origin"
    fasteners = {key: load_component(value) for key, value in mount["sources"].items() if key.startswith("origin:FJ1500/")}
    fixed = bounded({key: load_component(value) for key, value in mount["sources"].items()
                     if key != cap_key and key not in fasteners})
    cap = load_component(report["selected"]["mesh"])
    tpu = load_component(report["selected"]["recessed_TPU"])
    feed = -np.array(report["bore_axis_source"])
    normal = np.array(port["normal_source"])
    directions = {"feed_axis":feed,"cap_normal":normal,"reverse_cap_normal":-normal,
                  "halfway":A.D.unit(feed+normal)}
    results = []
    for name, direction in directions.items():
        for sequence, moving, obstacles in (
            ("bench_cap_into_TPU", {"cap":cap}, bounded({"TPU":tpu})),
            ("mount_preassembled_pair", {"cap":cap,"TPU":tpu}, fixed),
            ("mount_TPU_after_cap", {"TPU":tpu}, {**fixed,**bounded({"cap":cap})}),
        ):
            hits = probe(moving, obstacles, direction, .4)
            results.append(dict(sequence=sequence,direction=name,vector_source=direction.tolist(),step_mm=.4,hits=hits))
    qualified = []
    for sequence in ("bench_cap_into_TPU", "mount_preassembled_pair"):
        candidates = [r for r in results if r["sequence"] == sequence and not r["hits"]]
        for candidate in candidates:
            bench = sequence == "bench_cap_into_TPU"
            hits = probe({"cap":cap} if bench else {"cap":cap,"TPU":tpu}, bounded({"TPU":tpu}) if bench else fixed,
                         np.array(candidate["vector_source"]), .1)
            qualified.append(dict(sequence=sequence,direction=candidate["direction"],step_mm=.1,hits=hits))
    # Search a larger set of approach directions. The quick witnesses reject
    # directions early; every surviving direction is checked at 0.1 mm.
    search = []
    for index in range(180):
        z = 1.-2.*(index+.5)/180.
        azimuth = index*np.pi*(3.-np.sqrt(5.))
        direction = np.array([np.sqrt(1.-z*z)*np.cos(azimuth),np.sqrt(1.-z*z)*np.sin(azimuth),z])
        for later in (False, True):
            targets = {key:part for key,part in fixed.items() if not later or not key.startswith("origin:FJ1504/")}
            moving = {"cap":cap,"TPU":tpu}
            first = next((dict(distance_mm=distance,part=name,overlap_mm3=contact)
                          for distance in (.4,1.,2.,4.,8.) for name,part in moving.items()
                          if (contact:=overlaps(part.translate(distance*direction),targets))),None)
            dense = None if first else probe(moving,targets,direction,.1)
            search.append(dict(direction=direction.tolist(),FJ1504_origin_hardware_later=later,
                               witness=first,dense_hits=dense))
    seated = {**fixed, **bounded({"cap":cap,"TPU":tpu})}
    fastener_checks = []
    for index in (0,1):
        names = [f"origin:FJ1500/front_washer_{index}",f"origin:FJ1500/bolt_{index}"]
        for name in names:
            hits = probe({name:fasteners[name]},seated,normal,.2)
            fastener_checks.append(dict(part=name,step_mm=.2,nominal_overlap_mm3=overlaps(fasteners[name],seated),hits=hits))
            seated.update(bounded({name:fasteners[name]}))
    # Cut only the material swept by the original outside approach. Protect
    # the fastener neighbourhood and the complete channel wall near the seam.
    manifold = A.I._manifold()
    route_points = np.array(route["exported"]["trial"]["new_path_source_mm"])
    guard = union([A.D.outer_envelope(route_points,1.31),
                   *(p.minkowski_sum(manifold.Manifold.cube([.6]*3,True)) for p in fasteners.values())]) ^ tpu
    pad = manifold.Manifold.cube([.16]*3,True)
    bench_relief = []
    for direction_name in ("feed_axis","halfway","cap_normal"):
        direction = directions[direction_name]
        contact_stock = []
        for distance in np.linspace(0.,8.,41):
            contact = tpu^cap.translate(distance*direction)
            contact_stock.extend(q.hull().minkowski_sum(pad) for q in contact.decompose() if q.volume()>1e-5)
        if not contact_stock:
            continue
        cutter = union(contact_stock)-guard
        actual = (tpu-cutter).as_original().simplify(1e-6)
        components = sum(q.volume()>1e-5 for q in actual.decompose())
        loss = float((guard-actual).volume())
        hits = probe({"cap":cap},bounded({"TPU":actual}),direction,.1)
        row = dict(direction=direction_name, removed_mm3=float((tpu-actual).volume()),components=components,
                   protected_loss_mm3=loss, remaining_hits=hits)
        if components==1 and loss<=.001 and not hits:
            _, row["mesh"] = export(actual,folder/f"FJ1500_feed_bench_relief_{direction_name}.obj",np.array(report["start_source_mm"]))
            if probe({"cap":cap},bounded({"TPU":load_component(row["mesh"])}),direction,.1):
                raise ValueError("Bench path did not reproduce after readback")
        bench_relief.append(row)
    save(folder/"FJ1500_feed_installation_screen.json",dict(
        generator=provenance(ROOT/"tools/audit_s1_infraspinatus_feed_installation.py"), source=provenance(source),
        mount=report["mount"], approach_mm=8., results=results, finer_checks=qualified, fastener_checks=fastener_checks,
        direction_search=search, bench_relief=bench_relief,
        scope="Rest-shaped TPU and actual revised cap solids, nominal neutral rigid assembly. Straight local approaches."
              " Fastener checks use nominal axial translation and include any seating/thread overlap."
              " Does not simulate flexible assembly, helical threading, tools or whole-shoulder installation."))
    for result in results:
        print(result["sequence"],result["direction"],len(result["hits"]),flush=True)
    print('finer',[(r['sequence'],r['direction'],len(r['hits'])) for r in qualified],flush=True)
    print('clear directions by staging',[(later,sum(q['dense_hits']==[] and q['FJ1504_origin_hardware_later']==later for q in search))
                                         for later in (False,True)],flush=True)


if __name__ == "__main__":
    main()
