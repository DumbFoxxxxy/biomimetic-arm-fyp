"""Read back the FJ1504 curved cap and revised proximal TPU as one interface."""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np

from audit_s1_series_elbow_integration import MOBILE, assert_source, find_hits
from build_channel_first_inventory import DEST, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from study_elbow_series_package import bounded, overlaps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bidirectional", action="store_true")
    args = parser.parse_args()
    suffix = "_bidirectional" if args.bidirectional else ""
    folder = DEST.parent / "integration_review/s1_series_elbow"
    exit_file = folder / f"FJ1504_curved_local_exit{suffix}.json"
    pad_file = folder / "FJ1504_recessed_TPU_origin.json"
    guide, pad = read(exit_file), read(pad_file)
    source_file = assert_source(guide["assembly"])
    parent = read(source_file)
    if guide["assembly"] != pad["assembly"]:
        if not args.bidirectional:
            raise RuntimeError("The TPU and curved cap refer to different assemblies")
        pad_parent = read(assert_source(pad["assembly"]))
        for key in pad_parent["sources"]:
            if key.startswith("origin:FJ1504/") and parent["sources"].get(key) != pad_parent["sources"][key]:
                raise RuntimeError("FJ1504 rigid pad support changed")
    tpu_source, cap_source = pad["selected"]["source"], guide["selected"]["source"]
    cap_key = "cap:FJ1504/origin"
    sources = {**parent["sources"], cap_key: cap_source}
    rigid = {key: load_component(source) for key, source in sources.items()}
    tpu = load_component(tpu_source)
    neck = np.asarray(pad["selected"]["path_source_mm"])
    exterior = np.asarray(guide["selected"]["local_wire_path_source_mm"])
    np.testing.assert_allclose(neck[0], guide["selected"]["cap_inner_source_mm"], atol=1e-9)
    local_path = np.vstack([exterior[::-1], neck])
    nominal = A.D.outer_envelope(local_path, .3)
    tpu_hits = overlaps(tpu, bounded(rigid))
    wire_hits = overlaps(nominal, bounded({**rigid, "revised_TPU": tpu}))
    cap_hits = overlaps(rigid[cap_key], bounded({key: part for key, part in rigid.items() if key != cap_key}))
    stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(neck, axis=0), axis=1))]
    sections = []
    for distance in np.arange(1., 25.1, 1.):
        index = np.searchsorted(stations, distance) - 1
        axis = (neck[index + 1] - neck[index]) / (stations[index + 1] - stations[index])
        point = neck[index] + axis * (distance - stations[index])
        outer = rod(point - .1 * axis, point + .1 * axis, 1.5)
        inner = rod(point - .2 * axis, point + .2 * axis, .9)
        ring = outer - inner
        sections.append(dict(station_mm=float(distance), TPU_annulus_fraction=float((ring ^ tpu).volume() / ring.volume())))
    whole_probes = []
    if args.bidirectional:
        centres = read(DEST / "elbow_passive_model/model_check.json")["neutral_origins_m"]
        for joint, anchor, moving in (("SC", "clavicle", MOBILE | {"clavicle", "scapula"}),
                                      ("AC", "scapula", MOBILE | {"scapula"})):
            keys = list(rigid)
            active = [parent["owners"][key] in moving for key in keys]
            pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2) if active[i] != active[j]])
            for probe in guide["selected"]["whole_cap_SC_AC_probes"]:
                if probe["joint"] != joint:
                    continue
                matrix = A.I.rotation_about(np.asarray(centres[anchor]) * 1000.,
                                            np.eye(3)[probe["source_axis_index"]], probe["angle_deg"])
                posed = [solid.transform(matrix) if move else solid
                         for solid, move in zip(rigid.values(), active, strict=True)]
                whole_probes.append(dict(**probe, combined_intersections=find_hits(posed, keys, pairs)))
    target = folder / f"assembly_FJ1504_local_interface{suffix}.json"
    save(target, dict(
        **{key: value for key, value in parent.items() if key not in (
            "sources", "flexible_geometry_sources", "scope", "FJ1504_cap_revision", "FJ1504_TPU_revision")},
        sources=sources, flexible_geometry_sources={**parent["flexible_geometry_sources"], "tpu:subscapularis": tpu_source},
        FJ1504_cap_revision=provenance(exit_file), FJ1504_TPU_revision=provenance(pad_file),
        scope="Neutral combined candidate including the recessed TPU and cap-integrated FJ1504 outlet."))
    save(folder / f"FJ1504_local_interface{suffix}_audit.json", dict(
        generator=provenance(Path(__file__).resolve()), assembly=provenance(target),
        revised_TPU_vs_rigid_mm3=tpu_hits, curved_cap_vs_rigid_mm3=cap_hits,
        full_local_nominal_wire_vs_material_mm3=wire_hits, neutral_nominal_wire_diameter_mm=.6,
        combined_bidirectional_pose_probes=whole_probes,
        proximal_TPU_wall_samples=sections, local_path_source_mm=local_path.tolist(),
        scope="Neutral rigid and TPU mesh checks on the original internal FJ1504 route with the revised "
        "proximal portion and integral cap outlet. Annular occupancy diagnoses physical material near "
        "the proximal bore. No external supply routing, cable endurance, motion or drive acceptance."))
    print("TPU", tpu_hits, "cap", cap_hits, "wire", wire_hits, flush=True)
    if whole_probes:
        print("combined pose checks", len(whole_probes), "clear",
              sum(not pose["combined_intersections"] for pose in whole_probes), flush=True)
    print("annular material occupancy", [(row["station_mm"], round(row["TPU_annulus_fraction"], 3))
                                          for row in sections], flush=True)


if __name__ == "__main__":
    main()
