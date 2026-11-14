"""Screen paired blind fixing axes on the current FJ1467 scapula."""

from itertools import product

import numpy as np

from audit_s1_ac_blind_coupled_motion import load_assembly
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, rod, save, shoulder_pads
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    local = DEST.parent / "s1_local_muscle_load"
    folder = local / "FJ1500_cap_relief"
    route_ref = provenance(local / "FJ1467_posterior_rigid_locations.json")
    route = load_reference(route_ref)
    clearance = load_reference(load_reference(route["neighbor_replay"])["rigid_audit"])
    girdle = load_reference(clearance["girdle_path"])
    closure_ref = provenance(folder / "FJ1504_AC_blind_closure_audit.json")
    if girdle["closure_audit"] != closure_ref:
        raise ValueError("FJ1467 report and current scapula use different AC assemblies")
    _, _, sources = load_assembly(folder)
    bone = load_component(sources["bone:scapula"])
    ac_keys = [key for key in sources if key.startswith("AC:")]
    stock_targets = bounded({key: load_component(sources[key]) for key in ac_keys})
    port = route["source_ports"]["origin"]
    centre = np.asarray(port["inner_face_source_mm"])
    normal = np.asarray(port["normal_source"])
    basis = shoulder_pads.pad_basis(normal)
    line = A.D.outer_envelope(np.asarray(route["neutral_source_route_mm"])[:10], 1.10)
    targets = {**stock_targets, **bounded({"FJ1467_wire_guard": line})}
    rows = []
    for shift, pitch, azimuth in product((0., 3., -3.), (10., 14., 18.), range(0, 180, 10)):
        support_centre = centre + shift * basis[:, 0]
        angle = np.deg2rad(azimuth)
        lateral = basis[:, 0] * np.cos(angle) + basis[:, 1] * np.sin(angle)
        holes = []
        for sign in (-1, 1):
            point = support_centre + sign * pitch / 2. * lateral
            crossings = bone.ray_cast(point + normal * 14., point - normal * 45.)
            if len(crossings) < 2 or np.dot(crossings[0].normal, normal) <= 0:
                holes.append(dict(reference_mm=point.tolist(), fit=False, reason="No first entering bone span"))
                continue
            top, bottom = (float((np.asarray(hit.position) - point) @ normal) for hit in crossings[:2])
            # The first 1 mm is the entry transition. The next 4 mm is
            # insert length, followed by a 1 mm axial reserve.
            start = point + (top - 1.) * normal
            end = point + (top - 6.) * normal
            guard = rod(end, start, 3.)
            approach = rod(end, point + normal * 8., 2.)
            missing = max(0., float((guard - bone).volume()))
            hits = overlaps(approach, targets)
            holes.append(dict(reference_mm=point.tolist(), bone_entry_mm=(point + top * normal).tolist(),
                              insert_start_mm=start.tolist(), end_of_reserve_mm=end.tolist(),
                              stock_mm=top - bottom, reserve_after_4mm_insert_mm=top - bottom - 5.,
                              guard_missing_mm3=missing, approach_hits_mm3=hits,
                              fit=top - bottom >= 6. and missing <= .001 and not hits))
        rows.append(dict(shift_mm=shift, pitch_mm=pitch, azimuth_deg=azimuth, holes=holes,
                         fit=all(hole["fit"] for hole in holes)))
    eligible = [row for row in rows if row["fit"]]
    # A conforming support does not require symmetric blind axes. Screen
    # independently sampled parallel axes on a wider patch when needed.
    independent = []
    compact = []
    if not eligible:
        for u, v in product(np.arange(-14., 14.1, 2.), repeat=2):
            offset = np.array([u, v])
            point = centre + basis[:, :2] @ offset
            crossings = bone.ray_cast(point + normal * 14., point - normal * 45.)
            if len(crossings) < 2 or np.dot(crossings[0].normal, normal) <= 0:
                continue
            top, bottom = (float((np.asarray(hit.position) - point) @ normal) for hit in crossings[:2])
            if top > 7.:
                continue
            for accepted, depth, radius, approach_radius in ((independent, 6., 3., 2.), (compact, 4., 2., 1.5)):
                if top - bottom < depth:
                    continue
                start, end = (point + (top - distance) * normal for distance in (1., depth))
                if (rod(end, start, radius) - bone).volume() > .001:
                    continue
                if overlaps(rod(end, point + normal * 8., approach_radius), targets):
                    continue
                accepted.append(dict(offset_uv_mm=offset.tolist(), axis_source_mm=point.tolist(),
                                     bone_entry_mm=(point + top * normal).tolist(), stock_mm=top - bottom))
    independent_pairs = pair_axes(independent)
    compact_pairs = pair_axes(compact)
    if eligible:
        eligible.sort(key=lambda row: (abs(row["shift_mm"]),
                                       -min(hole["reserve_after_4mm_insert_mm"] for hole in row["holes"])))
    out = local / "FJ1467_current_scapula_blind_axes.json"
    save(out, dict(generator=provenance(ROOT / "tools/study_s1_posterior_mount_stock.py"),
                   route=route_ref, closure=closure_ref, scapula=sources["bone:scapula"],
                   normal_source=normal.tolist(), checked_AC_parts=ac_keys, cases=rows,
                   eligible_count=len(eligible), selected=eligible[0] if eligible else None,
                   independent_axes=independent, independent_pair_count=len(independent_pairs),
                   independent_pair_nearest_centre=independent_pairs[0] if independent_pairs else None,
                   compact_axes=compact, compact_pair_count=len(compact_pairs),
                   compact_pair_nearest_centre=compact_pairs[0] if compact_pairs else None,
                   scope="Neutral axes only. Standard insert is 4 mm deep after a 1 mm transition, "
                         "with 3 mm radial stock radius, 1 mm axial reserve and 2 mm approach radius. "
                         "Compact hypotheses use 2 mm insert depth, 2 mm stock radius, 1 mm reserve "
                         "and 1.5 mm approach radius. No material is removed. Support conformity, "
                         "actual hardware retention and dynamic clearance remain unverified."))
    print(f"Screened {len(rows)} symmetric pairs; {len(eligible)} passed. "
          f"Asymmetric grid: {len(independent_pairs)} standard and {len(compact_pairs)} compact pairs", flush=True)


def pair_axes(independent):
    independent_pairs = []
    for i, first in enumerate(independent):
        for second in independent[i + 1:]:
            p, q = (np.asarray(row["offset_uv_mm"]) for row in (first, second))
            spacing = np.linalg.norm(p - q)
            if 8. <= spacing <= 22.:
                independent_pairs.append(dict(first=first, second=second, spacing_mm=float(spacing),
                                              centre_offset_mm=float(np.linalg.norm((p + q) / 2.))))
    independent_pairs.sort(key=lambda row: row["centre_offset_mm"])
    return independent_pairs


if __name__ == "__main__":
    main()
