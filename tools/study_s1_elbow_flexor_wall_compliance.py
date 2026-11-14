"""Estimate elbow flexor wire seating with unilateral radial TPU supports."""

import numpy as np
from scipy.optimize import minimize

from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from study_s1_local_force_geometry import load_reference


def solve(points, rays, tension, stiffness):
    indices = np.array([ray["vertex"] for ray in rays], dtype=int)
    directions = np.array([ray["direction_source"] for ray in rays])
    gaps = np.array([ray["first_material_mm"] for ray in rays]) - .30
    weights = np.array([ray.get("stiffness_weight", 1.) for ray in rays])
    spring_rates = stiffness * weights
    if np.any(gaps < 0.):
        raise ValueError("Nominal steel already overlaps sampled wall")
    if not len(indices):
        return dict(arc_shortening_mm=None, maximum_equilibrium_residual_N=None,
                    optimizer_success=None, stations=[])

    def evaluate(displacement):
        posed = points.copy()
        posed[indices] += displacement[:, None] * directions
        segments = np.diff(posed, axis=0)
        lengths = np.linalg.norm(segments, axis=1)
        tangents = segments / lengths[:, None]
        cable_gradient = -tension * np.einsum("ij,ij->i", np.diff(tangents, axis=0)[indices - 1], directions)
        compression = np.maximum(displacement - gaps, 0.)
        reaction = spring_rates * compression
        return lengths.sum(), compression, reaction, cable_gradient + reaction

    def energy(displacement):
        length, compression, _, gradient = evaluate(displacement)
        return tension * (length - initial_length) + .5 * (compression @ (spring_rates * compression)), gradient

    initial_length = np.linalg.norm(np.diff(points, axis=0), axis=1).sum()
    result = minimize(energy, np.zeros(len(indices)), jac=True, method="L-BFGS-B",
                      bounds=[(0., None)] * len(indices),
                      options=dict(ftol=1e-14, gtol=1e-7, maxiter=500))
    length, compression, reaction, gradient = evaluate(result.x)
    residual = np.where(result.x <= 1e-8, np.minimum(gradient, 0.), gradient)
    error = float(np.max(np.abs(residual)))
    # The energy is convex in these radial displacements. The projected
    # gradient checks equilibrium even when line search stalls near roundoff.
    if not np.isfinite(error) or error > 1e-4:
        raise ValueError(f"Radial contact failed to equilibrate: {result.message}, residual {error:g} N")
    return dict(
        arc_shortening_mm=float(initial_length - length),
        maximum_equilibrium_residual_N=error,
        optimizer_success=bool(result.success),
        stations=[dict(vertex=int(vertex), centre_shift_mm=float(shift),
                       wall_displacement_mm=float(movement), reaction_N=float(load),
                       initial_radial_wire_gap_mm=float(gap))
                  for vertex, shift, movement, load, gap
                  in zip(indices, result.x, compression, reaction, gaps, strict=True)],
    )


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    stem = "nominal_spring_spool_mu_30_adjusted_flexor_wall"
    source = provenance(folder / f"{stem}.json")
    report = load_reference(source)
    by_time = {}
    for ray in report["rays"]:
        if ray["offset_deg"] == 0:
            by_time.setdefault(ray["time_s"], []).append(ray)
    cases = []
    for stiffness in (5., 10., 20., 50.):
        for frame in report["paths"]:
            rays = by_time.get(frame["time_s"], [])
            points = np.array(frame["points_mm"])
            if rays and not np.allclose(points[[ray["vertex"] for ray in rays]],
                               [ray["centre_source_mm"] for ray in rays], atol=1e-8, rtol=0.):
                raise ValueError("Radial support positions do not match the cable path")
            response = solve(points, rays, frame["tension_N"], stiffness)
            cases.append(dict(time_s=frame["time_s"], tension_N=frame["tension_N"],
                              stiffness_per_station_N_per_mm=stiffness, **response))
    output = folder / f"{stem}_compliance.json"
    save(output, dict(
        generator=provenance(ROOT / "tools/study_s1_elbow_flexor_wall_compliance.py"),
        wall_screen=source, cases=cases,
        scope="Quasi-static energy minimum with fixed recorded tension and frictionless cable. "
              "Selected nodes slide only along the sampled load rays; endpoints and unprobed nodes "
              "stay fixed. Wall supports are independent unilateral radial springs at the selected "
              "discretization. Stiffness values are sensitivity inputs, not identified TPU properties. "
              "Segment-wise contact, deformation outside the support direction, pressure, strength "
              "and feedback into arm motion are not evaluated.",
    ))
    print(f"Solved {len(cases)} radial contact cases")


if __name__ == "__main__":
    main()
