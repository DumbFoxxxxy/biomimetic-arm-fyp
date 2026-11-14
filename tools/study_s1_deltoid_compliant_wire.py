"""Solve a reduced, frictionless FJ1468 cable/TPU contact model."""

import argparse

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, minimize

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from channel_first_bone_contact import BoneContactSurface
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nonlinear", action="store_true", help="Resolve signed-distance contact during optimization")
    parser.add_argument("--distributed", action="store_true", help="Split the TPU patch into two coupled regions")
    args = parser.parse_args()
    nonlinear = args.nonlinear
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1468_contact_load_screen.json"
    loads = read(source)
    taut = load_reference(loads["source"])
    motion = load_reference(taut["source"])
    transition = load_reference(load_reference(motion["wall_audit"])["source"])
    mounts = load_reference(load_reference(transition["source"])["proximal_mount"])["sources"]
    reference = max(taut["cases"], key=lambda row: row["node_count"])
    full_reference = np.array(reference["route_source_mm"])
    end_index = 3 + reference["node_count"]
    points = full_reference[3:end_index]
    patch_force = np.array(loads["by_owner"]["TPU"]["resultant_N_at_10N"])
    mode = A.D.unit(patch_force)
    solids = {"TPU": load_component(motion["merged_TPU"]),
              "cap": load_component(mounts["cap:FJ1468/origin"])}
    lo, hi = points.min(axis=0) - 4., points.max(axis=0) + 4.
    box = A.I._manifold().Manifold.cube(hi - lo).translate(lo)
    surfaces = {name: BoneContactSurface(solid ^ box, band_mm=2., cell_mm=1.)
                for name, solid in solids.items()}

    tangent = points[2:] - points[:-2]
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    transverse = np.cross(tangent, np.eye(3)[np.argmin(abs(tangent), axis=1)])
    transverse /= np.linalg.norm(transverse, axis=1)[:, None]
    wall_count = 2 if args.distributed else 1
    count = 2 * len(tangent) + wall_count
    wall_slice = slice(count - wall_count, count)
    node_jac = np.zeros((len(points), 3, count))
    for i, (a, b) in enumerate(zip(transverse, np.cross(tangent, transverse)), start=1):
        node_jac[i, :, 2 * (i - 1):2 * i] = np.stack([a, b], axis=1)
    subdivisions = 8
    index = np.r_[np.repeat(np.arange(len(points) - 1), subdivisions), len(points) - 2]
    fraction = np.r_[np.tile(np.arange(subdivisions) / subdivisions, len(points) - 1), 1.]
    sample_points = (1. - fraction[:, None]) * points[index] + fraction[:, None] * points[index + 1]
    sample_jac = ((1. - fraction[:, None, None]) * node_jac[index]
                  + fraction[:, None, None] * node_jac[index + 1])
    # One effective translation mode for the near-port patch. Its amplitude
    # tapers to zero at the two fixed cable endpoints.
    stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    at = (1. - fraction) * stations[index] + fraction * stations[index + 1]
    weights = np.clip(np.minimum(at, stations[-1] - at) / 2., 0., 1.)
    if args.distributed:
        fraction_distal = np.clip((at / stations[-1] - .25) / .5, 0., 1.)
        shapes = weights[:, None] * np.column_stack((1. - fraction_distal, fraction_distal))
        # Equal displacements have the same total spring energy as the
        # one-mode model. A coupling spring resists differential displacement.
        stiffness_ratios = .5 * np.eye(2) + .25 * np.array([[1., -1.], [-1., 1.]])
    else:
        shapes = weights[:, None]
        stiffness_ratios = np.eye(1)
    matrices, gaps = [], []
    relative_jacobians = {}
    for name, surface in surfaces.items():
        dist, normal = surface.query(sample_points)
        jacobian = sample_jac.copy()
        if name == "TPU":
            jacobian[:, :, wall_slice] = -shapes[:, None, :] * mode[None, :, None]
        relative_jacobians[name] = jacobian
        matrices.append(np.einsum("ij,ijk->ik", normal, jacobian))
        gaps.append(dist - .30)
    constraint = LinearConstraint(np.vstack(matrices), -np.concatenate(gaps), np.inf)
    cached_state = None
    cached_contact = None

    def contact(state):
        nonlocal cached_state, cached_contact
        if cached_state is None or not np.array_equal(state, cached_state):
            distances, gradients = [], []
            for name, surface in surfaces.items():
                jacobian = relative_jacobians[name]
                distance, normal = surface.query(sample_points + jacobian @ state)
                distances.append(distance - .30)
                gradients.append(np.einsum("ij,ijk->ik", normal, jacobian))
            cached_state = state.copy()
            cached_contact = np.concatenate(distances), np.vstack(gradients)
        return cached_contact

    nonlinear_constraint = NonlinearConstraint(lambda x: contact(x)[0], 0., np.inf,
                                               jac=lambda x: contact(x)[1])
    cable_jac = np.diff(node_jac, axis=0)
    diff = np.diff(points, axis=0)
    neutral_length = float(np.linalg.norm(diff, axis=1).sum())
    cases = []
    load_cases = ((10., 20.), (10., 10.), (20., 20.)) if nonlinear or args.distributed else (
        (0., 20.), (5., 20.), (10., 40.), (10., 20.), (10., 10.), (20., 20.))
    for tension, stiffness in load_cases:
        def energy(x):
            segments = diff + cable_jac @ x
            lengths = np.linalg.norm(segments, axis=1)
            gradient = tension * np.einsum("ij,ijk->k", segments / lengths[:, None], cable_jac)
            wall = x[wall_slice]
            restoring_force = stiffness * (stiffness_ratios @ wall)
            gradient[wall_slice] += restoring_force
            return tension * (lengths.sum() - neutral_length) + .5 * (wall @ restoring_force), gradient

        bound = Bounds(np.r_[np.full(count - wall_count, -.5), np.zeros(wall_count)], np.full(count, .5))
        result = minimize(energy, np.zeros(count), jac=True, method="SLSQP", bounds=bound,
                          constraints=[constraint], options={"ftol": 1e-10, "maxiter": 150})
        if nonlinear:
            bound = Bounds(np.r_[np.full(count - wall_count, -1.25), np.zeros(wall_count)], np.full(count, 1.25))
            result = minimize(energy, result.x, jac=True, method="SLSQP", bounds=bound,
                              constraints=[nonlinear_constraint], options={"ftol": 1e-8, "maxiter": 60})
        exact = {}
        for name, surface in surfaces.items():
            dist, _ = surface.query(sample_points + relative_jacobians[name] @ result.x)
            exact[name] = float(dist.min() - .3)
        linear_residual = float((constraint.A @ result.x - constraint.lb).min())
        posed = points + node_jac @ result.x
        full = np.vstack((full_reference[:3], posed, full_reference[end_index:]))
        cap_hits = overlaps(A.D.outer_envelope(full, .30), bounded({"cap": solids["cap"]}))
        wall = result.x[wall_slice]
        row = dict(tension_N=tension, stiffness_N_per_mm=stiffness, converged=bool(result.success),
                   solver_message=result.message, wall_travel_mm=float(np.max(shapes @ wall)),
                   wall_mode_amplitudes_mm=wall.tolist(),
                   wall_spring_force_N=float(np.sum(stiffness * (stiffness_ratios @ wall))),
                   predicted_cable_change_mm=float(np.linalg.norm(diff + cable_jac @ result.x, axis=1).sum()
                                                   - neutral_length),
                   linear_min_clearance_mm=linear_residual, exact_sample_min_clearance_mm=exact,
                   full_cable_cap_overlap_mm3=cap_hits,
                   displacement_bound_hit=bool(np.any(abs(result.x) >= bound.ub - .0001)),
                   cable_points_source_mm=posed.tolist())
        if nonlinear:
            values, jacobian = contact(result.x)
            multipliers = getattr(result, "multipliers", None)
            if multipliers is not None:
                sample_count = len(sample_points)
                row["contact_reaction_scalar_N"] = {
                    name: float(multipliers[i * sample_count:(i + 1) * sample_count].sum())
                    for i, name in enumerate(surfaces)
                }
                row["stationarity_residual_N"] = float(np.max(abs(energy(result.x)[1] - jacobian.T @ multipliers)))
            row["minimum_contact_clearance_mm"] = float(values.min())
        cases.append(row)
        print({k: v for k, v in row.items() if k != "cable_points_source_mm"}, flush=True)
    if args.distributed:
        output = "FJ1468_compliant_wire_distributed.json" if nonlinear else "FJ1468_compliant_wire_distributed_linear.json"
    else:
        output = "FJ1468_compliant_wire_nonlinear.json" if nonlinear else "FJ1468_compliant_wire_screen.json"
    save(folder / output, dict(
        generator=provenance(ROOT / "tools/study_s1_deltoid_compliant_wire.py"), source=provenance(source),
        cases=cases, wall_direction_source=mode.tolist(), contact_subdivisions_per_segment=subdivisions,
        wall_modes=wall_count, wall_stiffness_matrix_per_unit_stiffness=stiffness_ratios.tolist(),
        contact_method="iterated_signed_distance" if nonlinear else "linearized_normals",
        scope="Neutral geometry with effective TPU displacement modes, fixed end positions "
              "and frictionless constant cable tension. Nonlinear mode reevaluates the signed "
              "distances throughout optimization. Bounds are numerical trust limits. "
              "No physical stroke, friction, or actual TPU stiffness is inferred.",
    ))


if __name__ == "__main__":
    main()
