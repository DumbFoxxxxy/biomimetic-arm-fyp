"""Find a side-entry approach with the FJ1500 sleeve removed."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_removable_feed_screen.json"
    removable = read(source)
    feed = load_reference(removable["integrated_sleeve"])
    assembly = load_reference(feed["mount"])
    tpu = load_component(feed["selected"]["recessed_TPU"])
    cap = load_component(removable["selected"]["cap"])
    fixed = {key: load_component(value) for key, value in assembly["sources"].items()
             if key != "cap:FJ1500/origin"}
    fixed["cap:FJ1500/origin"] = cap
    targets = bounded({key: part for key, part in fixed.items() if not key.startswith("origin:FJ1504/")})
    baseline = overlaps(tpu, targets, threshold=1e-6)

    def check(displacement):
        return overlaps(tpu.translate(displacement), targets)

    def score(hits):
        return sum(hits.values())

    # Select directions by their actual first separation from the fitted
    # position. The previous screen ranked directions by all 8 mm at once.
    directions = []
    for i in range(144):
        z = 1. - 2.*(i+.5)/144.
        azimuth = i*np.pi*(3.-np.sqrt(5.))
        r = np.sqrt(1.-z*z)
        direction = np.array([r*np.cos(azimuth),r*np.sin(azimuth),z])
        hits = check(.3*direction)
        directions.append(dict(direction_source=direction.tolist(),first_offset_mm=.3,hits_mm3=hits))
    finalists = sorted(directions,key=lambda row:score(row["hits_mm3"]))[:12]
    traces = []
    for row in finalists:
        direction = np.array(row["direction_source"])
        hits = [dict(offset_mm=float(d),hits_mm3=check(d*direction)) for d in np.arange(.1,8.01,.1)]
        traces.append(dict(direction_source=direction.tolist(), samples=hits,
                           furthest_clear_local_step_mm=max((q["offset_mm"] for q in hits if not q["hits_mm3"]),default=0.)))
    # Try a dog-leg path from the best separation direction, holding the
    # orientation fixed and moving to the side before pulling farther out.
    best = min(traces,key=lambda row:(sum(bool(q["hits_mm3"]) for q in row["samples"][:10]),
                                     sum(score(q["hits_mm3"]) for q in row["samples"])))
    direction = np.array(best["direction_source"])
    basis_a = A.D.unit(np.cross(direction,np.eye(3)[np.argmin(abs(direction))]))
    basis_b = np.cross(direction,basis_a)
    bent = []
    for azimuth, travel in product(np.arange(0.,360.,45.),(2.,4.,8.)):
        angle = np.deg2rad(azimuth)
        side = np.cos(angle)*basis_a+np.sin(angle)*basis_b
        samples = []
        for d in np.arange(.25,8.01,.25):
            displacement = d*direction+min(d,2.)/2.*travel*side
            samples.append(dict(offset_source_mm=displacement.tolist(),hits_mm3=check(displacement)))
        bent.append(dict(azimuth_deg=float(azimuth),lateral_mm=travel,
                         hit_count=sum(bool(q["hits_mm3"]) for q in samples),samples=samples))
    save(folder/"FJ1500_side_entry_screen.json",dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_side_entry.py"), source=provenance(source),
        seated_contacts_mm3=baseline, first_step_directions=directions, straight_traces=traces, two_stage_traces=bent,
        scope="Nominal rigid placement of resting TPU with sleeve removed and FJ1504 origin hardware deferred."
              " Seated contacts recorded separately. No flexible contact force or later hardware installation validation."))
    print("Clear first steps",sum(not q["hits_mm3"] for q in directions),
          "straight",[(sum(bool(s["hits_mm3"]) for s in q["samples"]),q["furthest_clear_local_step_mm"]) for q in traces],
          "best dogleg",min(q["hit_count"] for q in bent),flush=True)


if __name__ == "__main__":
    main()
