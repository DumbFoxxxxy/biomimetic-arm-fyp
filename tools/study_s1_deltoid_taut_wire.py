"""Solve a frictionless tension-only wire in the fixed neutral FJ1468 guide."""

import numpy as np
from scipy.optimize import minimize

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from channel_first_bone_contact import BoneContactSurface
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def solve(points, surface, steps):
    """Minimize length with each interior node free within its transverse plane."""
    directions = points[2:] - points[:-2]
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    transverse = np.cross(directions, np.eye(3)[np.argmin(abs(directions), axis=1)])
    transverse /= np.linalg.norm(transverse, axis=1)[:, None]
    basis = np.stack([transverse, np.cross(directions, transverse)], axis=2)
    count = len(points) - 2
    node_jac = np.zeros((len(points), 3, 2 * count))
    for index in range(count):
        node_jac[index + 1, :, 2 * index:2 * index + 2] = basis[index]
    i = np.repeat(np.arange(len(points) - 1), steps)
    f = np.tile(np.arange(steps) / steps, len(points) - 1)
    sampling_jac = (1. - f[:, None, None]) * node_jac[i] + f[:, None, None] * node_jac[i + 1]
    ref_samples = (1. - f[:, None]) * points[i] + f[:, None] * points[i + 1]

    def route(x):
        return points + node_jac @ x

    def length(x):
        diff = np.diff(route(x), axis=0)
        lengths = np.linalg.norm(diff, axis=1)
        return lengths.sum(), np.einsum("ij,ijk->k", diff / lengths[:, None], np.diff(node_jac, axis=0))

    previous = None
    cached = None

    def contact(x):
        nonlocal previous, cached
        if previous is None or not np.array_equal(x, previous):
            previous = x.copy()
            dist, normals = surface.query(ref_samples + sampling_jac @ x)
            cached = dist - .30, np.einsum("ij,ijk->ik", normals, sampling_jac)
        return cached

    start = np.zeros(2 * count)
    initial_clearance = float(contact(start)[0].min())
    result = minimize(length, start, jac=True, method="SLSQP", bounds=[(-1., 1.)] * len(start),
                      constraints=[dict(type="ineq", fun=lambda x: contact(x)[0],
                                        jac=lambda x: contact(x)[1])],
                      options=dict(maxiter=60, ftol=1e-8))
    posed = route(result.x)
    clearance = contact(result.x)[0]
    reaction = []
    multipliers = getattr(result, "multipliers", None)
    if multipliers is not None and result.success and clearance.min() >= -1e-5:
        # KKT multipliers of a unit-tension length objective scale with tension.
        samples = ref_samples + sampling_jac @ result.x
        reaction = [dict(point_source_mm=samples[k].tolist(), normal_reaction_N_at_10N=float(10. * value))
                    for k, value in enumerate(multipliers) if value > 1e-5]
    return posed, dict(solver_success=bool(result.success), solver_message=result.message,
                       iterations=int(result.nit), initial_clearance_mm=initial_clearance,
                       minimum_constraint_clearance_mm=float(clearance.min()),
                       active_sample_count=int(np.count_nonzero(clearance < .005)),
                       arc_mm=float(result.fun), reaction_samples=reaction,
                       displacement_bound_hit=bool(np.max(abs(result.x)) > .999))


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1468_revised_channel_motion.json"
    motion = read(source)
    audit = load_reference(motion["wall_audit"])
    transition = load_reference(audit["source"])
    path = np.asarray(next(r["centreline_source_mm"] for r in transition["recuts"] if "mesh" in r))
    proximal = load_reference(load_reference(transition["source"])["proximal_mount"])["sources"]
    solids = {"TPU": load_component(motion["merged_TPU"]),
              "cap": load_component(proximal["cap:FJ1468/origin"])}
    local = path[3:52]
    lo, hi = local.min(axis=0) - 4., local.max(axis=0) + 4.
    region = A.I._manifold().Manifold.cube(hi - lo).translate(lo)
    # Clip for signed-distance queries only. Verify against full meshes.
    surface = BoneContactSurface(union(solids.values()) ^ region, band_mm=2., cell_mm=1.)
    cases = []
    for stride, stations in ((4, 4), (3, 5)):
        nodes = local[np.r_[np.arange(0, len(local) - 1, stride), len(local) - 1]]
        solved, report = solve(nodes, surface, stations)
        full = np.vstack([path[:3], solved, path[52:]])
        report.update(node_count=len(nodes), samples_per_segment=stations, route_source_mm=full.tolist(),
                      full_wire_hits_mm3=overlaps(A.D.outer_envelope(full, .3), bounded(solids)),
                      full_arc_mm=float(np.linalg.norm(np.diff(full, axis=0), axis=1).sum()))
        cases.append(report)
        print({key: value for key, value in report.items() if key not in ("route_source_mm", "reaction_samples")},
              flush=True)
    save(folder / "FJ1468_taut_wire_neutral.json", dict(
        generator=provenance(ROOT / "tools/study_s1_deltoid_taut_wire.py"), source=provenance(source),
        cases=cases, scope="Fixed neutral guide, no bending stiffness or friction. "
        "A 0.30 mm clearance is enforced at sampled cable stations. Full Boolean overlap is "
        "an independent check. Wire equilibrium is approximate; TPU deformation and dynamics are excluded.",
    ))


if __name__ == "__main__":
    main()
