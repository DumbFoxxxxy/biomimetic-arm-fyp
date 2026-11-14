"""Check tangent-continuous FJ1500 short-tendon routes at the solved cap poses."""

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_elbow_series_package import bounded, overlaps
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def centreline(points, count):
    a, b, c, d = points
    q = np.linspace(0., 1., count)[:, None]
    return (1.-q)**3*a + 3.*(1.-q)**2*q*b + 3.*(1.-q)*q**2*c + q**3*d


def fit_route(exits, tips, target_radius, available):
    leads = tips - exits
    directions = leads / np.linalg.norm(leads, axis=1)[:, None]
    chord = tips[1]-tips[0]
    span = np.linalg.norm(chord)
    if min(directions[0] @ chord, -directions[1] @ chord) <= 0.:
        return None
    gauss, weights = np.polynomial.legendre.leggauss(24)
    ts = np.linspace(0., 1., 129)[:, None]

    def controls(lengths):
        return np.array([tips[0], tips[0]+lengths[0]*directions[0], tips[1]+lengths[1]*directions[1], tips[1]])

    def first_second(points, t):
        a, b, c = 3. * np.diff(points, axis=0)
        return (1.-t)**2*a+2.*t*(1.-t)*b+t**2*c, 2.*((1.-t)*(b-a)+t*(c-b))

    def arc(lengths):
        velocity, _ = first_second(controls(lengths), (gauss[:, None]+1.)/2.)
        return float(.5 * weights @ np.linalg.norm(velocity, axis=1))

    def constraints(lengths):
        velocity, acceleration = first_second(controls(lengths), ts)
        speed = np.linalg.norm(velocity, axis=1)
        curvature = np.linalg.norm(np.cross(velocity, acceleration), axis=1) / np.maximum(speed**3, 1e-12)
        return np.r_[1./target_radius-curvature, (velocity @ chord) / span**2-.005]

    limits = [(1e-3, span), (1e-3, span)]
    def optimize(equal_length=None):
        conditions = [dict(type="ineq", fun=constraints)]
        if equal_length is not None:
            conditions.append(dict(type="eq", fun=lambda lengths: arc(lengths)-equal_length))
        candidates = []
        for ratio in (.20, .33, .55):
            solution = minimize(arc if equal_length is None else lambda x: np.sum((x-span/3.)**2),
                                np.full(2, span*ratio), method="SLSQP", bounds=limits, constraints=conditions,
                                options={"maxiter":100, "ftol":1e-10})
            if constraints(solution.x).min() >= -1e-7 and (equal_length is None or abs(arc(solution.x)-equal_length) < 1e-6):
                candidates.append(solution.x)
        return min(candidates, key=arc) if candidates else None

    minimum = optimize()
    if minimum is None:
        return None
    shortest = arc(minimum)
    fitted = optimize(available) if available >= shortest-1e-6 else None
    params = fitted if fitted is not None else minimum
    q = np.linspace(0., 1., 1025)[:, None]
    v, acc = first_second(controls(params), q)
    radii = np.linalg.norm(v, axis=1)**3 / np.maximum(np.linalg.norm(np.cross(v, acc), axis=1), 1e-15)
    points = centreline(controls(params), 65)
    return dict(minimum_cubic_arc_mm=shortest, added_length_over_chord_mm=shortest-span,
                shortage_mm=shortest-available, fitted_available_length=fitted is not None,
                verified_minimum_radius_mm=float(radii.min()), points=points)


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    screen_file = folder / "FJ1500_tpu_coupling_screen.json"
    screen = read(screen_file)
    path = load_reference(screen["prior_path"])
    rest = load_reference(path["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    prior = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(prior["local_geometry"])
    assembly = load_reference(prior["shoulder_assembly"])
    model = next(row for row in load_reference(prior["kinematics"])["models"] if row["head"] == HEAD)
    ports = next(row["ports"] for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    sleeve = load_reference(path["sleeve"])
    if sleeve["original"] != assembly["sources"][key]:
        raise ValueError("S1 geometry changed")
    sources = {**assembly["sources"], key: sleeve["selected"]["source"]}
    solids = {key: load_component(source) for key, source in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    owners = assembly["intended_body_owners"]
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    offset = np.array(rest["rest_offset_source_mm"])
    adjusted, adjusted_ports = rest_geometry(model, ports, offset)
    pins = [np.array(model[f"neutral_{label}_pin_source_mm"]) for label in ("cap", "bone")]
    axes = [np.array(model[f"{label}_pin_axis_source"]) for label in ("cap", "bone")]
    gain = screen["candidate"]["gain"]
    direction = np.array(screen["coupling_direction_source"])

    def wire_targets(state, cap_pin, gh):
        cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
        phi_a, phi_b = np.deg2rad(state["swivel_angles_deg"])
        transforms = dict(cap=rigid_transform(cap, pins[0], cap_pin),
                          arm=rigid_transform(gh, centre, centre),
                          eye_a=rigid_transform(cap*Rotation.from_rotvec(axes[0]*phi_a), pins[0], cap_pin),
                          eye_b=rigid_transform(gh*Rotation.from_rotvec(axes[1]*phi_b), pins[1], centre+gh.apply(pins[1]-centre)))
        targets = {}
        for name, part in solids.items():
            frame = None
            if name in (f"link:{HEAD}/muscle_groove_ring", f"link:{HEAD}/muscle_uncrimped_sleeve"):
                frame = transforms["eye_a"]
            elif name in (f"link:{HEAD}/bone_groove_ring", f"link:{HEAD}/bone_uncrimped_sleeve"):
                frame = transforms["eye_b"]
            elif owners[name] == f"free_terminal:{HEAD}":
                frame = transforms["cap"]
            elif owners[name] == "humerus_ideal_attachment":
                frame = transforms["arm"]
            targets[name] = part.transform(frame) if frame is not None else part
        return bounded(targets)

    trials = []
    smooth_states = []
    for angle in (60., 68., 75., 85., 96.):
        gh = Rotation.from_rotvec(flexion * np.deg2rad(angle))
        for added_length in (0., .25, .5):
            states = solve_branch(adjusted, adjusted_ports, centre, gh, bias["stiffness"], [0., 5., 10.],
                                  retain_leads=True, installed_length_offset=added_length,
                                  lateral_coupling=(direction, gain))
            for state in states:
                cap_pin, hits = collision_state(solids, owners, state, offset, gh, centre, *pins, *axes)
                exits, tips = (np.array(state[key]) for key in ("tendon_exits_source_mm", "lead_tips_source_mm"))
                chord = float(np.linalg.norm(tips[1]-tips[0]))
                available = chord + state["chord_gap_mm"]
                routes = []
                targets = wire_targets(state, cap_pin, gh)
                for limit in (2., 4., 6.):
                    route = fit_route(exits, tips, limit, available)
                    if route is None:
                        routes.append(dict(radius_probe_mm=limit, solved=False))
                        continue
                    points = route.pop("points")
                    tube = A.D.outer_envelope(points, .30)
                    routes.append(dict(radius_probe_mm=limit, solved=True, **route,
                                       radius_verified=route["verified_minimum_radius_mm"] >= limit-1e-3,
                                       rigid_wire_hits_mm3=overlaps(tube, targets), centerline_source_mm=points.tolist()))
                trials.append(dict(GH_deg=angle, input_N=state["commanded_muscle_tension_N"],
                                   assumed_installed_extra_mm=added_length, tension_N=state["short_tendon_tension_N"],
                                   chord_mm=chord, available_central_span_mm=available,
                                   rigid_terminal_hits=hits, routes=routes))
        for bending in (0., 10., 50.):
            states = solve_branch(adjusted, adjusted_ports, centre, gh, bias["stiffness"], [0., 5., 10., 5., 0.],
                                  retain_leads=True, lateral_coupling=(direction, gain), strand_bending_Nmm2=bending)
            for state in states:
                cap_pin, hits = collision_state(solids, owners, state, offset, gh, centre, *pins, *axes)
                exits, tips = (np.array(state[key]) for key in ("tendon_exits_source_mm", "lead_tips_source_mm"))
                lengths = np.linalg.norm(tips - exits, axis=1)
                tangents = (tips-exits) / lengths[:, None]
                h = np.linalg.norm(tips[1]-tips[0])/3.
                points = centreline([tips[0], tips[0]+h*tangents[0], tips[1]+h*tangents[1], tips[1]], 65)
                tube = A.D.outer_envelope(points, .30)
                smooth_states.append(dict(GH_deg=angle, assumed_bending_Nmm2=bending, cap_pin_source_mm=cap_pin.tolist(),
                                          rigid_terminal_hits=hits, central_curve_source_mm=points.tolist(),
                                          central_curve_rigid_hits_mm3=overlaps(tube, wire_targets(state, cap_pin, gh)), **state))
        print("Checked FJ1500 smooth routes at", angle, flush=True)
    def boundary_pose(angle):
        gh = Rotation.from_rotvec(flexion * np.deg2rad(angle))
        state = solve_branch(adjusted, adjusted_ports, centre, gh, bias["stiffness"], [0., 5., 10.],
                             retain_leads=True, lateral_coupling=(direction, gain), strand_bending_Nmm2=10.)[-1]
        cap, hits = collision_state(solids, owners, state, offset, gh, centre, *pins, *axes)
        return dict(GH_deg=angle, cap_pin_source_mm=cap.tolist(), rigid_terminal_hits=hits)

    coarse = [boundary_pose(float(angle)) for angle in range(68, 76)]
    first = next(i for i, row in enumerate(coarse) if row["rigid_terminal_hits"])
    if first == 0:
        raise ValueError("The lower boundary of the local screen is already obstructed")
    fine = [boundary_pose(round(coarse[first-1]["GH_deg"]+.1*i, 1)) for i in range(1, 10)]
    save(folder / "FJ1500_smooth_tendon_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_smooth_tendon.py"),
        solver=provenance(ROOT / "tools/study_s1_y_terminal_equilibrium.py"),
        coupling=provenance(screen_file), radius_probes_mm=[2.,4.,6.], added_installed_lengths_mm=[0.,.25,.5],
        cases=trials, arc_equilibrium=smooth_states, ten_newton_local_rigid_boundary=sorted(coarse+fine, key=lambda row: row["GH_deg"]),
        scope="Baseline cap equilibria use a chord inequality. At each fixed resulting pose,"
        " independently fit a monotone cubic span tangent to both retained crimp leads."
        " Minimum is a numerical result within that family. Assumed radii and installed lengths require calibration."
        " Separate arc equilibrium uses one cubic with handle lengths fixed to one third of the chord, an arc"
        " length inequality and assumed bending energy; this is a reduced cable model."
        " TPU contact, wire axial elasticity and full rod equilibrium are excluded."))


if __name__ == "__main__":
    main()
