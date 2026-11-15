"""Estimate static gravity balance along K7 using the audited local paths."""

import mujoco
import numpy as np
from scipy.optimize import linprog

from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from project_runtime import load_mujoco_model
from study_s1_local_force_geometry import load_reference


def pull_jacobian(model, data, route, origins):
    points = np.asarray(route["points_source_mm"])
    weight = np.asarray(route["distal_weights"])
    xyz = np.zeros_like(points)
    jacobian = np.zeros((len(points), 3, model.nv))
    for owner, fractions in zip(route["owners"], (1. - weight, weight), strict=True):
        body = data.body(owner)
        local = (points - origins[owner]) * .001
        world = local @ body.xmat.reshape(3, 3).T + body.xpos
        xyz += fractions[:, None] * world
        for index, fraction in enumerate(fractions):
            if fraction:
                jac = np.zeros((3, model.nv))
                mujoco.mj_jac(model, data, jac, None, world[index], body.id)
                jacobian[index] += fraction * jac
    segment = np.diff(xyz, axis=0)
    length = np.linalg.norm(segment, axis=1)
    direction = segment / length[:, None]
    return float(length.sum()), -np.einsum("ij,ijk->k", direction, np.diff(jacobian, axis=0))


def best_balance(matrix, demand, cap, extra_uncapped=None):
    """Minimize the worst absolute residual in rotational coordinates."""
    count = len(demand)
    cable_count = matrix.shape[1]
    bounds = [(0., cap)] * cable_count
    if extra_uncapped is not None:
        matrix = np.hstack((matrix, extra_uncapped))
        bounds.extend([(0., None)] * extra_uncapped.shape[1])
    slack = -np.ones((count, 1))
    result = linprog(np.r_[np.zeros(matrix.shape[1]), 1.],
                     A_ub=np.vstack((np.hstack((matrix, slack)), np.hstack((-matrix, slack)))),
                     b_ub=np.r_[demand, -demand], bounds=bounds + [(0., None)],
                     method="highs")
    if not result.success:
        raise ValueError(result.message)
    tensions = result.x[:-1]
    residual = demand - matrix @ tensions
    fit = dict(pull_N=tensions[:cable_count].tolist(), residual_Nm=residual.tolist(),
               max_abs_residual_Nm=float(max(abs(residual))))
    if extra_uncapped is not None:
        fit["pad_force_N"] = tensions[cable_count:].tolist()
    return fit


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    combined_ref = provenance(folder / "whole_arm_combined_cuff_k7.json")
    combined = load_reference(combined_ref)
    motion = load_reference(combined["motion"])
    scene = load_reference(motion["scene_check"])
    routed = load_reference(combined["routed"])
    routes = routed["routes"]
    route_samples = routed["paths"]["K7_SC1_pos"]["samples"]
    axes = load_reference(motion["anatomical_axes"])["axes_source"]
    axis_matrix = np.asarray([axes[name] for name in ("flexion", "abduction", "longitudinal")])
    model_path = ROOT / scene["model"]["file"]
    if provenance(model_path) != scene["model"]:
        raise ValueError("K7 geometry changed after the motion check")
    model = load_mujoco_model(model_path)
    data = mujoco.MjData(model)
    origins = {key: np.asarray(xyz) for key, xyz in scene["neutral_body_origins_mm"].items()}
    joints = ("SC", "AC", "GH", "elbow")
    indices = np.flatnonzero(np.isin(model.dof_jntid, [model.joint(name).id for name in joints]))
    rows = []
    samples = motion["paths"]["K7_SC1_pos"]["accepted"][:len(combined["combined"])]
    for index, (row, checked) in enumerate(zip(samples, combined["combined"], strict=True)):
        if row["input_deg"] != checked["input_deg"]:
            raise ValueError("K7 samples are out of phase")
        data.qpos[:] = row["qpos"]
        mujoco.mj_forward(model, data)
        entries = [pull_jacobian(model, data, route, origins) for route in routes.values()]
        np.testing.assert_allclose([entry[0] * 1000. for entry in entries],
                                   [route_samples[index]["routes"][key]["length_mm"] for key in routes],
                                   rtol=0., atol=1e-7)
        matrix = np.column_stack([entry[1] for entry in entries])
        demand = data.qfrc_bias.copy()
        left, right = max(0, index - 1), min(index + 1, len(samples) - 1)
        ds = np.deg2rad(samples[right]["input_deg"] - samples[left]["input_deg"])
        direction = np.empty(model.nv)
        mujoco.mj_differentiatePos(model, direction, ds,
                                   np.asarray(samples[left]["qpos"]), np.asarray(samples[right]["qpos"]))
        leverage = direction @ matrix
        # A one-dimensional estimate supposes that a separate guide or
        # constraint supplies all torques perpendicular to the path.
        gravity_path = float(direction @ demand)
        cuff_states = [next(s for s in combined["cases"][head]["K7_SC1_pos"][index]["states"]
                            if s["command_N"] == 10.) for head in ("FJ1500", "FJ1504")]
        # Apply the equal and opposite cuff couples to GH only. The local
        # solver reports projections on axes in the scapula reference frame.
        cuff_local = np.linalg.solve(axis_matrix, np.sum([s["humeral_torque_axis_Nm"] for s in cuff_states], axis=0))
        jac_h = np.zeros((3, model.nv))
        jac_s = np.zeros((3, model.nv))
        mujoco.mj_jacBody(model, data, None, jac_h, data.body("humerus").id)
        mujoco.mj_jacBody(model, data, None, jac_s, data.body("scapula").id)
        cuff_world = data.body("scapula").xmat.reshape(3, 3) @ cuff_local
        cuff = (jac_h - jac_s).T @ cuff_world
        rows.append(dict(
            input_deg=row["input_deg"], path_gravity_demand_Nm=gravity_path,
            ideal_path_cable_leverage_mm=(leverage * 1000.).tolist(),
            ideal_path_pull_at_best_channel_N=(gravity_path / max(leverage)).item()
            if gravity_path > 0 and max(leverage) > 0 else None,
            gravity_demand_Nm=demand[indices].tolist(),
            cable_leverage_mm=(matrix[indices] * 1000.).tolist(),
            cuff_assumed_couple_Nm=cuff[indices].tolist(),
            best_at_10N_per_line=best_balance(matrix[indices], demand[indices], 10.),
            best_with_10N_cuff_probes=best_balance(matrix[indices], (demand - cuff)[indices], 10.),
            best_at_30N_per_line=best_balance(matrix[indices], demand[indices], 30.),
        ))
    save(folder / "whole_arm_k7_gravity_balance.json", dict(
        combined=combined_ref, generator=provenance(ROOT / "tools/audit_whole_arm_k7_gravity_balance.py"),
        routes=list(routes), rotational_dofs=[model.joint(model.dof_jntid[i]).name for i in indices],
        samples=rows, scope="Quasi-static virtual-work estimate for prescribed samples, "
        "nominal mass model and frictionless tension-only lines. Joints held quasi-static with zero velocity. "
        "Eight routed cable derivatives use prescribed interpolation. The two cuff cases add their "
        "independently solved 10 N couples. Contact reactions, joint stop reactions, TPU external loads "
        "and carriage force balance are outside the calculation. 10 N and 30 N are probe caps."))
    print(f"{len(rows)} poses; {len(indices)} rotational coordinates; {len(routes)} variable-tension paths")


if __name__ == "__main__":
    main()
