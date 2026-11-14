"""Screen load-side TPU stock along the elbow flexor's internal steel path."""

import argparse

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import MOBILE, transforms
from audit_s1_series_elbow_spring_integration import load_candidate
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, rod, save
from refine_elbow_short_guide_relief import indexed_solid
from study_s1_local_force_geometry import load_reference


def radial_stock(mesh, point, direction, reach):
    hits = mesh.ray_cast(point, point + reach * direction)
    if not hits:
        return dict(first_material_mm=None, first_stock_mm=None, stock_truncated=False)
    entry = float((np.asarray(hits[0].position) - point) @ direction)
    if np.dot(hits[0].normal, direction) > 0.:
        return dict(first_material_mm=0., first_stock_mm=entry, stock_truncated=False)
    exit_distance = float((np.asarray(hits[1].position) - point) @ direction) if len(hits) > 1 else reach
    return dict(first_material_mm=entry, first_stock_mm=exit_distance - entry, stock_truncated=len(hits) == 1)


def refine_segment(mesh, points, rays, segment, divisions):
    """Add radial springs along one segment with a conserved stiffness sum."""
    indexed = {ray["vertex"]: ray for ray in rays}
    first, last = indexed[segment], indexed[segment + 1]
    span = points[segment + 1] - points[segment]
    tangent = span / np.linalg.norm(span)
    fractions = np.arange(1, divisions) / divisions
    inserted = points[segment] + fractions[:, None] * span
    refined = np.vstack((points[:segment + 1], inserted, points[segment + 1:]))
    supports = []
    for vertex, ray in indexed.items():
        weight = .5 + .5 / divisions if vertex in (segment, segment + 1) else 1.
        supports.append(dict(ray, vertex=vertex + (divisions - 1 if vertex > segment else 0),
                             stiffness_weight=weight))
    for vertex, (fraction, point) in enumerate(zip(fractions, inserted, strict=True), start=segment + 1):
        radial = (1. - fraction) * np.array(first["direction_source"]) + fraction * np.array(last["direction_source"])
        radial -= (radial @ tangent) * tangent
        radial /= np.linalg.norm(radial)
        stock = radial_stock(mesh, point, radial, 5.)
        if stock["first_material_mm"] is None:
            raise ValueError("Refined radial support did not reach TPU")
        supports.append(dict(vertex=vertex, direction_source=radial.tolist(), stiffness_weight=1. / divisions,
                             **stock))
    return refined, sorted(supports, key=lambda ray: ray["vertex"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spool", choices=("mu_30_adjusted", "seating_k5"), default="mu_30_adjusted",
                        help="Select the reference or the rerun with prescribed seating loss")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-seating", action="store_true",
                        help="Recheck the solved wire envelopes against the interpolated wall")
    mode.add_argument("--refine-segment", action="store_true", help="Refine radial supports between nodes 8 and 9")
    args = parser.parse_args()
    folder = DEST.parent / "integration_review/s1_series_elbow"
    stem = f"nominal_spring_spool_{args.spool}"
    travel_ref = provenance(folder / f"{stem}_wire_travel.json")
    travel = load_reference(travel_ref)
    audit = load_reference(travel["flexible_audit"])
    prior, extended, _, _, _, _, _, retention = load_candidate(folder)
    if travel["model"] != extended["model"] or audit["retention"] != retention or audit["trace"] != travel["trace"]:
        raise ValueError("Wall screen sources do not match the audited assembly")
    trace = checked_trace(travel["trace"])
    model_path = ROOT / extended["model"]["file"]
    if provenance(model_path) != extended["model"] or len(trace) != travel["sample_count"]:
        raise ValueError("Wall screen model or trajectory changed")
    model = mujoco.MjModel.from_xml_string(model_path.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    assembly = load_reference(load_reference(prior["passive"])["assembly"])
    muscle = assembly["muscles"]["flexor"]
    path = np.asarray(muscle["wire_path_source_mm"])[1:-1]
    axis = np.asarray(load_reference(muscle["source_definition"])["basis_columns_display"])[:, 2]
    gradient = axis / ((path[-1] - path[0]) @ axis)
    source = prior["mesh_sources"]["flexor"]
    _, vertices, faces = indexed_solid(source)
    path_weights = weigh(path, path[0], gradient)
    vertex_weights = weigh(vertices, path[0], gradient)
    time = np.array([float(row["time_s"]) for row in trace])
    force = np.array([float(row["flexor_tension_N"]) for row in trace])
    peak = travel["channels"]["FJ1486"]["constant_distal_tension_bend_load"]
    indices = set(np.flatnonzero(force > 1e-6)[::10])
    indices.update((int(np.argmax(force)), int(np.argmin(np.abs(time - peak["time_s"])))))
    mf = A.I._manifold()
    results = []
    paths = []
    seating = []
    refinement = []
    if args.refine_segment:
        from study_s1_elbow_flexor_wall_compliance import solve

        wall_ref = provenance(folder / f"{stem}_flexor_wall.json")
        wall = load_reference(wall_ref)
        if wall["travel"] != travel_ref or wall["mesh"] != source:
            raise ValueError("Radial probes do not match the active travel and wall mesh")
    if args.check_seating:
        stem_wall = f"{stem}_flexor_wall"
        compliance_ref = provenance(folder / f"{stem_wall}_compliance.json")
        compliance = load_reference(compliance_ref)
        wall = load_reference(compliance["wall_screen"])
        if wall["travel"] != travel_ref or wall["mesh"] != source:
            raise ValueError("Seating study does not match the active travel and wall mesh")
        solved = {}
        for case in compliance["cases"]:
            if case["stations"]:
                solved.setdefault(case["time_s"], []).append(case)
    for index in sorted(indices):
        row = trace[index]
        data.qpos[:] = [float(row[f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        frames = transforms(model, data, neutral)
        humerus = frames["humerus"]
        distal = frames[f"{muscle['anatomical_id']}_distal_pivot"]
        rotation = humerus[:, :3].T
        affine = (rotation @ distal[:, :3], rotation @ (distal[:, 3] - humerus[:, 3]), np.zeros(3))
        posed = displacement_map(path, path_weights, *affine)
        paths.append(dict(time_s=float(time[index]), tension_N=float(force[index]), points_mm=posed.tolist()))
        moved = displacement_map(vertices, vertex_weights, *affine)
        mesh = mf.Manifold(mf.Mesh64(np.ascontiguousarray(moved), faces))
        if mesh.status() != mf.Error.NoError:
            raise ValueError(f"Invalid interpolated TPU at frame {index}")
        if args.refine_segment:
            rays = [ray for ray in wall["rays"] if ray["time_s"] == float(time[index]) and ray["offset_deg"] == 0]
            if not {8, 9}.issubset({ray["vertex"] for ray in rays}):
                continue
            for divisions in (1, 2, 4, 8, 16):
                refined, supports = refine_segment(mesh, posed, rays, 8, divisions)
                for stiffness in (5., 20., 50.):
                    response = solve(refined, supports, float(force[index]), stiffness)
                    local = [row for row in response["stations"] if 8 <= row["vertex"] <= 8 + divisions]
                    refinement.append(dict(time_s=float(time[index]), subdivisions=divisions,
                                           nominal_spacing_mm=float(np.linalg.norm(posed[9] - posed[8]) / divisions),
                                           nominal_node_stiffness_N_per_mm=stiffness,
                                           region_reaction_sum_N=sum(row["reaction_N"] for row in local),
                                           region_wall_shift_peak_mm=max(row["wall_displacement_mm"] for row in local),
                                           **response))
            continue
        if args.check_seating:
            cases = solved.get(float(time[index]), [])
            if not cases:
                continue
            directions = {ray["vertex"]: np.array(ray["direction_source"]) for ray in wall["rays"]
                          if ray["time_s"] == float(time[index]) and ray["offset_deg"] == 0}
            baseline = float((A.D.outer_envelope(posed, .30) ^ mesh).volume())
            for case in cases:
                shifted = posed.copy()
                loaded = set()
                for station in case["stations"]:
                    vertex = station["vertex"]
                    shifted[vertex] += station["centre_shift_mm"] * directions[vertex]
                    if station["reaction_N"] > 1e-4:
                        loaded.add(vertex)
                overlap = float((A.D.outer_envelope(shifted, .30) ^ mesh).volume())
                segment_hits = []
                if overlap > .001:
                    for segment in range(len(shifted) - 1):
                        volume = float((A.D.outer_envelope(shifted[segment:segment + 2], .30) ^ mesh).volume())
                        if volume > .001:
                            start, end = shifted[segment:segment + 2]
                            span = end - start
                            middle = rod(start + .25 * span, start + .75 * span, .30, 48)
                            segment_hits.append(dict(segment=segment, overlap_mm3=volume,
                                                     middle_half_overlap_mm3=float((middle ^ mesh).volume()),
                                                     touches_loaded_node=bool({segment, segment + 1} & loaded)))
                seating.append(dict(time_s=case["time_s"],
                                    stiffness_per_station_N_per_mm=case["stiffness_per_station_N_per_mm"],
                                    baseline_overlap_mm3=baseline, seated_overlap_mm3=overlap,
                                    segment_hits=segment_hits))
            continue
        tangents = np.diff(posed, axis=0)
        tangents /= np.linalg.norm(tangents, axis=1)[:, None]
        change = np.diff(tangents, axis=0)
        loads = force[index] * np.linalg.norm(change, axis=1)
        for bend in np.flatnonzero(loads >= 1.):
            vertex = int(bend + 1)
            radial = change[bend] / np.linalg.norm(change[bend])
            tangent = tangents[bend] + tangents[bend + 1]
            side = np.cross(tangent, radial)
            side /= np.linalg.norm(side)
            for degrees in (-20, 0, 20):
                angle = np.deg2rad(degrees)
                direction = radial * np.cos(angle) + side * np.sin(angle)
                results.append(dict(
                    time_s=float(time[index]), vertex=vertex, offset_deg=degrees,
                    estimated_load_N=float(loads[bend]), centre_source_mm=posed[vertex].tolist(),
                    direction_source=direction.tolist(), **radial_stock(mesh, posed[vertex], direction, 5.),
                ))
    if args.refine_segment:
        if not refinement:
            raise ValueError("No frames had both original radial supports")
        save(folder / f"{stem}_flexor_wall_refined_segment.json", dict(
            generator=provenance(ROOT / "tools/audit_s1_elbow_flexor_channel_wall.py"),
            solver=provenance(ROOT / "tools/study_s1_elbow_flexor_wall_compliance.py"),
            wall_screen=wall_ref, cases=refinement,
            scope="Subdivide segment 8–9 only, in poses with both original supports. Endpoints retain "
                  "their exterior half-weights and the interior weights use trapezoidal quadrature. "
                  "The total nominal spring rate for this region is twice the listed node rate at "
                  "all subdivisions. The summed reaction is a scalar sum, not a vector resultant. "
                  "The radial springs and guide geometry are a reduced frictionless model.",
        ))
        print(f"Solved {len(refinement)} segment refinement cases")
        return
    if args.check_seating:
        if len(seating) != sum(map(len, solved.values())) or not seating:
            raise ValueError("Missing solved seating cases")
        save(folder / f"{stem_wall}_seating_audit.json", dict(
            generator=provenance(ROOT / "tools/audit_s1_elbow_flexor_channel_wall.py"),
            compliance=compliance_ref, cases=seating,
            scope="Boolean overlap against the pose-interpolated wall before radial spring displacement. "
                  "Overlap at loaded nodes is expected in this comparison. Segment results indicate "
                  "where to refine the model and their volumes are not additive. A flat-ended cylinder "
                  "checks the middle half of each intersecting segment independently. An unloaded-node segment "
                  "can still intersect the adjacent moving wall. No physical pinch or load capacity is inferred.",
        ))
        print(f"Audited {len(seating)} seated cable envelopes")
        return
    if not results:
        raise ValueError("No loaded channel stations found")
    peak_ray = max(results, key=lambda ray: ray["estimated_load_N"])
    if not np.isclose(peak_ray["estimated_load_N"], peak["magnitude_N"], rtol=0., atol=1e-8):
        raise ValueError("Channel replay did not reproduce the recorded peak bend load")
    sampled_stock = [ray for ray in results if ray["first_material_mm"] is not None]
    save(folder / f"{stem}_flexor_wall.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_elbow_flexor_channel_wall.py"),
        travel=travel_ref, mesh=source, sampled_frames=len(indices), rays=results, paths=paths,
        missing_material_rays=len(results) - len(sampled_stock),
        truncated_stock_rays=sum(ray["stock_truncated"] for ray in sampled_stock),
        minimum_stock_lower_bound_mm=min((ray["first_stock_mm"] for ray in sampled_stock), default=None),
        scope="Every tenth loaded sample plus peak tension and peak bend-load frames. At bends above 1 N, "
              "probe the load direction and two directions offset by 20 degrees up to 5 mm from the steel "
              "centreline. First stock is the radial thickness of the first contiguous material interval "
              "and is a lower bound if that interval extends past the ray. "
              "Interpolated TPU geometry and uniform distal tension; this is not a contact or strength solution.",
    ))
    print(f"Screened {len(results)} rays in {len(indices)} loaded frames")


if __name__ == "__main__":
    main()
