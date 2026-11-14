"""Screen insertion of the scapula with the preassembled FJ1504 cap."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def screen(moving, fixed, directions):
    targets = bounded(fixed)

    def probe(direction, distance):
        rows = {}
        for key, solid in moving.items():
            hits = overlaps(solid.translate(direction * distance), targets)
            if hits:
                rows[key] = hits
        return rows

    def trace(direction, offsets):
        hits = []
        for distance in offsets:
            found = probe(direction, distance)
            if found:
                hits.append(dict(offset_mm=float(distance), moving_hits_mm3=found))
        return dict(outward_direction_source=direction.tolist(), sample_count=len(offsets), hits=hits)

    def rank(result):
        return (len(result["hits"]),
                sum(v for hit in result["hits"] for part in hit["moving_hits_mm3"].values()
                    for v in part.values()))

    coarse = [trace(direction, (.5, 2., 5., 10., 16., 24.)) for direction in directions]
    best = sorted(coarse, key=rank)[:3]
    dense = [trace(np.array(row["outward_direction_source"]), np.linspace(24., 0., 241)) for row in best]
    return dict(moving=list(moving), fixed=list(fixed), coarse=coarse, fine=dense)


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_cap_with_bolts_screen.json"
    cap_screen = read(source)
    blind = load_reference(cap_screen["source"])
    split_audit = load_reference(load_reference(blind["parent"])["split_audit"])
    split = load_reference(split_audit["split"])
    side = load_reference(load_reference(split["parent"])["parent"])
    removable = load_reference(side["source"])
    assembly = load_reference(load_reference(removable["integrated_sleeve"])["mount"])
    sources = dict(assembly["sources"])
    sources.update(blind["selected"]["meshes"])
    sources["cap:FJ1504/origin"] = cap_screen["local_cap_relief"]["mesh"]
    sources["cap:FJ1500/origin"] = removable["selected"]["cap"]
    sources["FJ1504/TPU"] = assembly["flexible_geometry_sources"]["tpu:subscapularis"]
    replaced = set(blind["retired_hardware"]) - set(blind["selected"]["meshes"])
    pending = [key for key, owner in assembly["owners"].items()
               if owner not in {"scapula", "clavicle", "thorax"} or key.startswith(("AC:", "GH"))]
    sources = {key: row for key, row in sources.items() if key not in replaced and key not in pending}
    solids = {key: load_component(row) for key, row in sources.items()}
    moving = {key: part for key, part in solids.items() if assembly["owners"].get(key) == "scapula"
              or key == "FJ1504/TPU"}
    fixed = {key: part for key, part in solids.items() if assembly["owners"].get(key) == "thorax"}
    clavicle = {"bone:clavicle": solids["bone:clavicle"]}
    sc_parts = {key: part for key, part in fixed.items() if key.startswith("SC:")}
    fixed_open = {key: part for key, part in fixed.items() if key not in sc_parts}
    directions = [axis * sign for axis in np.eye(3) for sign in (-1., 1.)]
    for i in range(48):
        z = 1. - 2. * (i + .5) / 48.
        angle = i * np.pi * (3. - np.sqrt(5.))
        r = np.sqrt(1. - z * z)
        directions.append(np.array([r * np.cos(angle), r * np.sin(angle), z]))
    setups = [
        ("clavicle_seated_AC_hardware_deferred", moving, {**fixed, **clavicle}),
        ("clavicle_carried_SC_and_AC_hardware_deferred", {**moving, **clavicle}, fixed_open),
    ]
    cases = {}
    for name, mobile, stationary in setups:
        cases[name] = screen(mobile, stationary, directions)
        print(name, [len(row["hits"]) for row in cases[name]["fine"]], flush=True)
    save(folder / "FJ1504_preassembled_scapula_entry.json", dict(
        generator=provenance(ROOT / "tools/study_s1_scapula_preassembled_entry.py"),
        source=provenance(source), sources=sources, cases=cases,
        deferred_keys=pending, retired_keys=sorted(replaced), threshold_mm3=.001,
        scope="24 mm fixed-orientation translations near final fit. FJ1504 neutral TPU moves "
              "with scapula, with the distal connection open. FJ1500 belly and sleeve are deferred. "
              "AC and GH hardware are deferred in both cases; SC hardware is deferred only "
              "when the clavicle is carried. No final joint closure or tool access validation.",
    ))


if __name__ == "__main__":
    main()
