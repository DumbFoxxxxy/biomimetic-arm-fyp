"""Recheck the modified AC closure and local rigid motion with exported parts."""

import numpy as np

from audit_s1_infraspinatus_split_installation import insertion
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_AC_right_blind_fastener.json"
    design = read(source)
    pose = load_reference(design["source"])
    access = load_reference(pose["source"])
    closure = load_reference(access["source"])
    entry = load_reference(closure["source"])
    prior = entry["cases"]["clavicle_seated_AC_hardware_deferred"]
    selected = design["selected"]
    scene = {key: load_component(meta) for key, meta in entry["sources"].items()}
    scene["bone:scapula"] = load_component(selected["meshes"]["bone:scapula"])
    scene["cap_right"] = load_component(selected["meshes"]["cap_right"])
    hardware = {key.removeprefix("AC:AC_"): load_component(meta)
                for key, meta in closure["hardware_sources"].items() if key.endswith("_left")}
    hardware["bottom_washer_left"] = load_component(access["left_washer"]["mesh"])
    hardware["screw_right"] = load_component(selected["meshes"]["screw_right"])
    old = {key.removeprefix("AC:AC_"): load_component(meta)
           for key, meta in closure["hardware_sources"].items()}
    axis = np.asarray(closure["axis_source"])
    stage_config = (
        ("cap_left", None, 12.), ("top_washer_left", axis, 12.), ("bolt_left", axis, 22.),
        ("bottom_washer_left", -axis, 12.), ("nut_left", -axis, 12.), ("screw_right", axis, 16.),
    )
    installed = dict(scene)
    stages = []
    for name, direction, travel in stage_config:
        if direction is None:
            directions = [np.asarray(trace["offset_direction_source"])
                          for trace in closure["cases"][0]["stages"][0]["traces"] if not trace["hits"]]
        else:
            directions = [direction]
        traces = [insertion(hardware[name], bounded(installed), vector, travel) for vector in directions]
        stages.append(dict(part=name, traces=traces, clear=any(not trace["hits"] for trace in traces)))
        installed[name] = hardware[name]
    print("closure", [(row["part"], row["clear"]) for row in stages], flush=True)

    pivot = np.asarray(pose["pivot_source_mm"])
    targets = bounded({key: scene[key] for key in prior["fixed"]})
    baseline = {key: load_component(entry["sources"][key]) for key in prior["moving"]}
    comparison = {"original_AC": {**baseline, **old}, "short_right_screw": {
        **{key: scene[key] for key in prior["moving"]}, **hardware, "cap_right": scene["cap_right"]}}
    motion = {}
    for name, mobile in comparison.items():
        def hits(angle, rotation_axis):
            transform = A.I.rotation_about(pivot, rotation_axis, angle)
            return {key: found for key, part in mobile.items() if key != "FJ1504/TPU"
                    and (found := overlaps(part.transform(transform), targets))}

        neutral = hits(0., np.eye(3)[0])
        rows = []
        for rotation_axis in np.eye(3):
            for sign in (-1., 1.):
                last = 0.
                first = None
                samples = np.linspace(.5, 20., 40) if not neutral else ()
                for magnitude in samples:
                    contact = hits(sign * magnitude, rotation_axis)
                    if contact:
                        lo, hi = last, magnitude
                        while hi - lo > .1:
                            mid = (lo + hi) / 2.
                            refined = hits(sign * mid, rotation_axis)
                            if refined:
                                hi, contact = mid, refined
                            else:
                                lo = mid
                        last = lo
                        first = dict(angle_deg=sign * hi, hits_mm3=contact)
                        break
                    last = magnitude
                rows.append(dict(axis_source=rotation_axis.tolist(), last_clear_deg=sign * last,
                                 first_interference=first, neutral_blocked=bool(neutral)))
        motion[name] = dict(neutral_hits_mm3=neutral, axes=rows)
        print(name, [(round(row["last_clear_deg"], 3),
                      list((row["first_interference"] or {}).get("hits_mm3", {}))) for row in rows], flush=True)
    save(folder / "FJ1504_AC_blind_closure_audit.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_ac_blind_closure.py"), source=provenance(source),
        closure_stages=stages, motion=motion, pivot_source_mm=pivot.tolist(),
        scope="Nominal AC closure with right cap held at its registered pose, including the "
              "neutral subscapularis in insertion checks. Motion rotates the scapula and AC "
              "hardware against seated clavicle and thorax; flexible contacts and later GH/arm "
              "hardware are omitted. Threshold 0.001 mm3. Thread and temporary retention unverified.",
    ))


if __name__ == "__main__":
    main()
