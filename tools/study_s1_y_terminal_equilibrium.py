"""Solve the S1-Y2 free caps and eye swivels under local cable loading."""

import argparse
import json

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from study_s1_local_force_geometry import load_reference


ANGLE_SCALE = 10.  # mm-equivalent numerical coordinate per radian
GRADIENT_STEP = 1e-4


def derivative(function, state):
    directions = np.eye(len(state)) * GRADIENT_STEP
    return np.array([(function(state + step) - function(state - step)) / (2. * GRADIENT_STEP)
                     for step in directions])


def solve_branch(model, ports, gh_centre, humerus_rotation, stiffness, drive_levels,
                 retain_leads=False, installed_length_offset=0., lateral_coupling=None,
                 translation_energy=None):
    pin_a = np.array(model["neutral_cap_pin_source_mm"])
    pin_b = np.array(model["neutral_bone_pin_source_mm"])
    axis_a = np.array(model["cap_pin_axis_source"])
    axis_b = np.array(model["bone_pin_axis_source"])
    offset_a = np.array(model["neutral_crimp_outer_exit_A_source_mm"]) - pin_a
    offset_b = np.array(model["neutral_crimp_outer_exit_B_source_mm"]) - pin_b
    origin = np.array(ports["origin"]["inner_face_source_mm"])
    stop_offset = np.array(ports["distal"]["outer_face_source_mm"]) - pin_a
    bone_pin = gh_centre + humerus_rotation.apply(pin_b - gh_centre)
    strand_length = model["nominal_free_strand_material_length_mm"] + installed_length_offset
    path = np.array(model["active_strand_centreline_source_mm"])
    lead_a = path[1] - path[0] if retain_leads else np.zeros(3)
    lead_b = path[-2] - path[-1] if retain_leads else np.zeros(3)
    lead_length = np.linalg.norm(lead_a) + np.linalg.norm(lead_b)
    active_axis = origin - pin_a - stop_offset
    active_axis /= np.linalg.norm(active_axis)
    axial, transverse, cap_rotation, eye_rotation = stiffness
    coupling_axis = np.zeros(3)
    coupling_gain = 0.
    if lateral_coupling is not None:
        direction, coupling_gain = lateral_coupling
        coupling_axis = np.array(direction, dtype=float)
        coupling_axis -= (coupling_axis @ active_axis) * active_axis
        coupling_axis /= np.linalg.norm(coupling_axis)

    def geometry(state):
        angle = state[3:] / ANGLE_SCALE
        cap = Rotation.from_rotvec(angle[:3])
        free_pin = pin_a + state[:3]
        stop = free_pin + cap.apply(stop_offset)
        eye_a = Rotation.from_rotvec(axis_a * angle[3])
        eye_b = Rotation.from_rotvec(axis_b * angle[4])
        point_a = free_pin + cap.apply(eye_a.apply(offset_a))
        point_b = bone_pin + humerus_rotation.apply(eye_b.apply(offset_b))
        tip_a = point_a + cap.apply(eye_a.apply(lead_a))
        tip_b = point_b + humerus_rotation.apply(eye_b.apply(lead_b))
        return stop, point_a, point_b, tip_a, tip_b

    def gap(state):
        _, _, _, tip_a, tip_b = geometry(state)
        return strand_length - lead_length - np.linalg.norm(tip_a - tip_b)

    initial_length = np.linalg.norm(pin_a + stop_offset - origin)

    def potential(state, pull):
        displacement = state[:3]
        axial_motion = displacement @ active_axis
        lateral = displacement - axial_motion * active_axis
        lateral -= coupling_gain * axial_motion * coupling_axis
        angle = state[3:] / ANGLE_SCALE
        stop, *_ = geometry(state)
        elastic = (.5 * (axial * axial_motion**2 + transverse * (lateral @ lateral))
                   if translation_energy is None else translation_energy(displacement))
        return elastic + .5 * (cap_rotation * (angle[:3] @ angle[:3])
                     + eye_rotation * (angle[3:] @ angle[3:])) + pull * (np.linalg.norm(stop - origin) - initial_length)

    state = np.zeros(8)
    rows = []
    for pull in drive_levels:
        def objective(position):
            return potential(position, pull)

        result = minimize(objective, state, jac=lambda q: derivative(objective, q), method="SLSQP",
                          constraints=[dict(type="ineq", fun=gap, jac=lambda q: derivative(gap, q))],
                          options=dict(ftol=1e-11, maxiter=250))
        state = result.x
        clearance = gap(state)
        gradient = derivative(objective, state)
        gap_gradient = derivative(gap, state)
        tension = max(0., float(gradient @ gap_gradient / (gap_gradient @ gap_gradient))) if clearance < 1e-6 else 0.
        residual = np.linalg.norm(gradient - tension * gap_gradient)
        def lagrangian(position):
            return objective(position) - tension * gap(position)

        tangent = null_space(gap_gradient.reshape(1, -1)) if tension > 1e-7 else np.eye(len(state))
        steps = tangent.T * .001
        curvature = tangent.T @ np.column_stack([
            (derivative(lagrangian, state + step) - derivative(lagrangian, state - step)) / .002
            for step in steps])
        minimum_curvature = float(np.linalg.eigvalsh((curvature + curvature.T) / 2.)[0])
        if clearance < -1e-6 or residual > 5e-5 or minimum_curvature < -1e-6:
            raise RuntimeError(f"{model['head']} at {pull} N: {result.message}, gap {clearance}, residual {residual}")
        stop, point_a, point_b, tip_a, tip_b = geometry(state)
        humeral_force = tension * (tip_a - tip_b) / np.linalg.norm(tip_a - tip_b)
        active_force = pull * (origin - stop) / np.linalg.norm(origin - stop)
        torque = np.cross(tip_b - gh_centre, humeral_force) / 1000.
        rows.append(dict(
            commanded_muscle_tension_N=pull, short_tendon_tension_N=tension,
            chord_gap_mm=clearance, stationarity_residual_N=float(residual),
            minimum_tangent_curvature_N_per_mm=minimum_curvature, solver_status=result.message,
            complementarity_residual_Nmm=abs(clearance * tension),
            muscle_proxy_chord_change_mm=float(np.linalg.norm(stop - origin) - initial_length),
            cap_translation_source_mm=state[:3].tolist(),
            cap_rotation_vector_deg=np.rad2deg(state[3:6] / ANGLE_SCALE).tolist(),
            swivel_angles_deg=np.rad2deg(state[6:] / ANGLE_SCALE).tolist(),
            cap_stop_source_mm=stop.tolist(), tendon_exits_source_mm=[point_a.tolist(), point_b.tolist()],
            lead_tips_source_mm=[tip_a.tolist(), tip_b.tolist()],
            proximal_actuation_reaction_N=(-active_force).tolist(),
            scapula_net_reaction_N=(-humeral_force).tolist(),
            humerus_reaction_N=humeral_force.tolist(), GH_torque_source_Nm=torque.tolist()))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retain-leads", action="store_true")
    args = parser.parse_args()
    output = DEST.parent / "s1_local_muscle_load"
    local_file = output / "gh_endpoint_geometry.json"
    local = read(local_file)
    assembly_path = ROOT / local["shoulder_assembly"]["file"]
    load_reference(local["shoulder_assembly"])
    kinematics_path = assembly_path.parent / "free_terminal_kinematics/kinematics_audit.json"
    kinematics = read(kinematics_path)
    if kinematics["assembly_sha256"] != local["shoulder_assembly"]["sha256"]:
        raise ValueError("Free-terminal geometry and S1-Y2 assembly do not match")
    port_map = {row["head"]: row["ports"] for row in load_reference(local["port_geometry"])["routes"]}
    gh_centre = np.array(local["gh_centre_source_mm"])
    axes = np.array(local["axes_source"])
    probes = [("neutral", 0., axes[0]), ("flexion", 30., axes[0]), ("flexion", 60., axes[0]),
              ("flexion", 96., axes[0]), ("abduction", 30., axes[1]),
              ("longitudinal", -30., axes[2]), ("longitudinal", 30., axes[2])]
    stiffness_cases = [(axial, 2. * axial, cap, eye) for axial, cap, eye in
                       [(.2, 50., 1.), (.1, 50., 1.), (.4, 50., 1.), (.2, 25., .1), (.2, 100., 5.)]]
    length_offsets = [0., 1.] if args.retain_leads else [0.]
    if args.retain_leads:
        stiffness_cases = stiffness_cases[:1]
    loads = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 8., 6., 4., 2., 0.]
    cases = []
    for model in kinematics["models"]:
        head = model["head"]
        for motion, angle, axis in probes:
            humerus_rotation = Rotation.from_rotvec(axis * np.deg2rad(angle))
            for stiffness in stiffness_cases:
                for offset in length_offsets:
                    rows = solve_branch(model, port_map[head], gh_centre, humerus_rotation, stiffness, loads,
                                        args.retain_leads, offset)
                    for row in rows:
                        row["GH_torque_projection_Nm"] = (axes @ np.array(row["GH_torque_source_Nm"])).tolist()
                    cases.append(dict(head=head, prescribed_axis=motion, prescribed_GH_deg=angle,
                                      stiffness=list(stiffness), assumed_length_offset_mm=offset, rows=rows))
        print(f"Solved {head} cases", flush=True)
    strand_scope = ("The two source exit leads remain straight. The middle span is the chord between lead tips,"
                    " a length lower bound. It permits idealized point bends at the tips."
                    if args.retain_leads else "A straight chord ignores crimp exit tangent constraints and obstacles.")
    report = dict(
        generator=provenance(ROOT / "tools/study_s1_y_terminal_equilibrium.py"),
        local_geometry=provenance(local_file), kinematics=provenance(kinematics_path),
        shoulder_assembly=local["shoulder_assembly"], GH_projection_axes=local["axes"],
        stiffness_columns=["axial_N_per_mm", "transverse_N_per_mm", "cap_Nmm_per_rad", "eye_Nmm_per_rad"],
        stiffness_uncalibrated=True, nominal_strand_length_used_as_installed_length_hypothesis=True,
        retains_source_exit_leads=args.retain_leads,
        scope="Local frictionless static equilibrium with prescribed GH pose."
              " Cap translation and rotation and both swivel angles are solved."
              f" {strand_scope}"
              " No cap gravity, contact, assembly-wide collision, dynamic return, or stiffness calibration.",
        cases=cases)
    filename = "y_terminal_retained_lead_equilibrium.json" if args.retain_leads else "y_terminal_equilibrium.json"
    (output / filename).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Solved {sum(len(case['rows']) for case in cases)} quasi-static loading/unloading states")


if __name__ == "__main__":
    main()
