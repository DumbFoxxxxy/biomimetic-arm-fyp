"""Quantify wire travel and friction sensitivity on the FJ1500 channel centreline."""

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_infraspinatus_distributed_deformation import blend
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from study_s1_local_force_geometry import load_reference


def route_metrics(points):
    segments = np.diff(points, axis=0)
    length = np.linalg.norm(segments, axis=1)
    if np.any(length < 1e-8):
        raise ValueError("Collapsed centreline segment")
    tangent = segments / length[:, None]
    cross = np.linalg.norm(np.cross(tangent[:-1], tangent[1:]), axis=1)
    dots = np.sum(tangent[:-1]*tangent[1:], axis=1)
    turns = np.arctan2(cross, dots)
    radii = np.linalg.norm(points[2:]-points[:-2], axis=1)/(2.*np.maximum(cross, 1e-15))
    return float(length.sum()), turns, radii


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_distributed_deformation_screen.json"
    audit = read(source)
    combined = load_reference(audit["combined"])
    rest = load_reference(audit["rest_geometry"])
    revision = load_reference(combined["exit_revision"])
    base = load_reference(revision["baseline"])
    route = next(r for r in load_reference(base["input_route"])["routes"] if r["head"] == "FJ1500")
    prior = load_reference(load_reference(rest["bias_screen"])["reference_equilibrium"])
    spec = next(r for r in load_reference(prior["kinematics"])["models"] if r["head"] == "FJ1500")
    origin = np.array(route["ports"]["origin"]["inner_face_source_mm"])
    axis = np.array(route["ports"]["distal"]["outer_face_source_mm"])-origin
    axis /= np.linalg.norm(axis)
    pin = np.array(spec["neutral_cap_pin_source_mm"])+rest["rest_offset_source_mm"]
    path = np.array(rest["rest_path_source_mm"])
    rest_span, rest_turn, rest_radii = route_metrics(path)
    lattice = load_component(rest["candidate"])
    smoothing = []
    candidate = path.copy()
    # Keep the port and adjacent lead stations fixed, then test the smoother
    # routes against the existing TPU material at rest before considering a recut.
    for iteration in range(1, 17):
        candidate[2:-2] = .25*candidate[1:-3]+.5*candidate[2:-2]+.25*candidate[3:-1]
        if iteration in (1, 2, 4, 8, 16):
            length, turn, radius = route_metrics(candidate)
            smoothing.append(dict(passes=iteration, length_mm=length, min_radius_mm=float(radius.min()),
                                  total_turn_deg=float(np.rad2deg(turn.sum())),
                                  max_station_displacement_mm=float(np.linalg.norm(candidate-path,axis=1).max()),
                                  nominal_wire_TPU_overlap_mm3=float((A.D.outer_envelope(candidate,.30)^lattice).volume())))
    states = next(q["equilibrium"] for q in combined["scenarios"] if q["label"] == "offset_16mm")
    results = []
    friction = (0., .1, .2, .3)
    for geometry in audit["results"]:
        if not geometry["name"].startswith("extended_ramps"):
            continue
        weight, _ = blend((path-origin)@axis, *geometry["blend_range_mm"], geometry["ramp_fraction"])
        checked = {(q["frame"],q["input_N"]):q for q in geometry["cases"]}
        records = []
        for row in states:
            load = row["commanded_muscle_tension_N"]
            pair = row["frame"],load
            if checked[pair]["nominal_internal_wire_overlap_mm3"]>.001:
                raise ValueError("Nominal route intersects TPU")
            rotation = Rotation.from_rotvec(np.deg2rad(row["cap_rotation_vector_deg"]))
            posed = path+weight[:,None]*(rotation.apply(path-pin)-(path-pin)+row["cap_translation_source_mm"])
            length, turning, radius = route_metrics(posed)
            tight = int(radius.argmin())
            wrap = float(turning.sum())
            rows = [dict(mu=mu, distal_over_proximal=float(np.exp(-mu*wrap)),
                         input_for_specified_distal_load_N=float(load*np.exp(mu*wrap))) for mu in friction]
            reaction = float((2.*load*np.sin(turning/2.)).max())
            records.append(dict(frame=row["frame"], time_s=row["time_s"], specified_distal_load_N=load,
                                wire_path_mm=length, take_up_from_rest_mm=rest_span-length,
                                straight_proxy_change_mm=row["muscle_proxy_chord_change_mm"],
                                total_polyline_turn_rad=wrap, tightest_three_point_radius_mm=float(radius[tight]),
                                tightest_rest_station_index=tight+1, tightest_point_source_mm=posed[tight+1].tolist(),
                                largest_sampled_turn_reaction_at_equal_tension_N=reaction, friction_sensitivity=rows))
        results.append(dict(distribution=geometry["name"], records=records))
    if len(results)!=2:
        raise ValueError("Both checked distributed routes are required")
    save(folder/"FJ1500_wire_travel_and_friction.json", dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_wire_travel.py"), deformation=provenance(source),
        rest_path_mm=rest_span, wire_radius_mm=.30, friction_probes=list(friction), scenarios=results,
        rest_route_turns=[dict(index=i, turn_deg=float(np.rad2deg(t)), radius_mm=float(r), point_source_mm=path[i].tolist())
                          for i,(t,r) in enumerate(zip(rest_turn,rest_radii,strict=True),1)], rest_smoothing_trials=smoothing,
        scope="Centreline lengths under imposed TPU deformation. Positive take-up means cable path shortens."
              " Three-point radii use stored stations and are local curvature indicators."
              " Capstan sensitivity assumes sliding contact over every sampled turn and is applied after equilibrium."
              " It does not predict static friction, contact distribution, return hysteresis or a new motion."))
    for scenario in results:
        high = [q for q in scenario["records"] if q["specified_distal_load_N"]==10.]
        print(scenario["distribution"], len(scenario["records"]), '10 N travel',
              min(q["take_up_from_rest_mm"] for q in high),max(q["take_up_from_rest_mm"] for q in high),
              'r_min',min(q["tightest_three_point_radius_mm"] for q in high),
              '10 N proximal at mu .3',max(q["friction_sensitivity"][-1]["input_for_specified_distal_load_N"] for q in high))


if __name__ == "__main__":
    main()
