"""Evaluate the existing unilateral scapulothoracic pad hypothesis on K7."""

import mujoco
import numpy as np
from scipy.optimize import linprog

from audit_whole_arm_k7_gravity_balance import best_balance
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from project_runtime import load_mujoco_model
from study_s1_local_force_geometry import load_reference


def pad_reaction(model, data, pads):
    scapula = data.body("scapula")
    rotation = scapula.xmat.reshape(3, 3)
    result = []
    total = np.zeros(model.nv)
    for pad in pads:
        point = scapula.xpos + rotation @ pad["local_point_m"]
        vector = point - pad["rib_point_m"]
        normal = np.asarray(pad["normal"])
        height = float(vector @ normal)
        lateral = float(np.linalg.norm(vector - height * normal))
        compression = pad["rest_separation_m"] - height
        engaged = lateral <= pad["radius_m"] and 0. < height < pad["rest_separation_m"]
        force = compression * pad["stiffness_N_per_m"] if engaged else 0.
        jac = np.zeros((3, model.nv))
        mujoco.mj_jac(model, data, jac, None, point, scapula.id)
        leverage = jac.T @ normal
        total += leverage * force
        result.append(dict(id=pad["candidate_index"], separation_mm=1000. * height,
                           lateral_offset_mm=1000. * lateral, compression_mm=1000. * max(compression, 0.),
                           engaged=engaged, bottomed=height <= 0., reaction_N=force,
                           leverage=leverage))
    return total, result


def fit_long_travel_pads(cable, pad, gravity, travel_mm, stiffness, pad_cap):
    """Fit one neutral preload per pad and independent cable pulls at each pose."""
    pose_count, axis_count, line_count = cable.shape
    pad_count = pad.shape[2]
    force_drop = travel_mm * stiffness
    preload_bounds = list(zip(force_drop.max(axis=0), pad_cap + force_drop.min(axis=0), strict=True))
    if any(low > high for low, high in preload_bounds):
        return dict(feasible=False, minimum_peak_force_N=np.ptp(force_drop, axis=0).tolist())
    torque = np.zeros((pose_count, axis_count, pad_count + pose_count * line_count))
    torque[:, :, :pad_count] = pad
    for i in range(pose_count):
        start = pad_count + i * line_count
        torque[i, :, start:start + line_count] = cable[i]
    compensated_gravity = gravity + np.einsum("ijk,ik->ij", pad, force_drop)
    matrix = torque.reshape(-1, torque.shape[-1])
    slack = -np.ones((len(matrix), 1))
    result = linprog(np.r_[np.zeros(matrix.shape[1]), 1.],
                     A_ub=np.vstack((np.hstack((matrix, slack)), np.hstack((-matrix, slack)))),
                     b_ub=np.r_[compensated_gravity.ravel(), -compensated_gravity.ravel()],
                     bounds=preload_bounds + [(0., 30.)] * (pose_count * line_count) + [(0., None)],
                     method="highs")
    if not result.success:
        raise ValueError(result.message)
    preload = result.x[:pad_count]
    pulls = result.x[pad_count:-1].reshape(pose_count, line_count)
    pad_force = preload - force_drop
    residual = gravity - np.einsum("ijk,ik->ij", pad, pad_force) - np.einsum("ijk,ik->ij", cable, pulls)
    np.testing.assert_allclose(np.max(abs(residual)), result.fun, rtol=0., atol=1e-8)
    return dict(feasible=True, stiffness_N_per_mm=stiffness, neutral_preload_N=preload.tolist(),
                neutral_compression_mm=(preload / stiffness).tolist(), pad_forces_N=pad_force.tolist(),
                pulls_N=pulls.tolist(), residual_Nm=residual.tolist(),
                max_abs_residual_Nm=np.max(abs(residual), axis=1).tolist())


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    gravity_ref = provenance(folder / "whole_arm_k7_gravity_balance.json")
    gravity = load_reference(gravity_ref)
    combined = load_reference(gravity["combined"])
    motion = load_reference(combined["motion"])
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    pad_ref = provenance(DEST / "unified_SC_AC/passive_pad_engagement_paths.json")
    pad_file = load_reference(pad_ref)
    old_model = load_reference(pad_file["endpoint_model"])
    original_chain = load_reference(old_model["chain_check"])
    for segment in ("clavicle", "scapula", "humerus"):
        np.testing.assert_allclose(original_chain["neutral_body_origins_m"][segment],
                                   np.asarray(scene["neutral_body_origins_mm"][segment]) * .001,
                                   atol=1e-8, rtol=0.)
    original_source = load_reference(original_chain["reference"])
    if any(original_source["sources"][key] != source for key, source in transfer["shared_sources"].items()
           if key.startswith("thorax_reference:")):
        raise ValueError("Pad data reference a different thorax")
    model_path = ROOT / scene["model"]["file"]
    if provenance(model_path) != scene["model"]:
        raise ValueError("K7 model has changed")
    model = load_mujoco_model(model_path)
    data = mujoco.MjData(model)
    indices = np.flatnonzero(np.isin(model.dof_jntid, [model.joint(name).id for name in ("SC", "AC", "GH", "elbow")]))
    poses = motion["paths"]["K7_SC1_pos"]["accepted"][:len(gravity["samples"])]
    rows = []
    for state, stored in zip(poses, gravity["samples"], strict=True):
        if state["input_deg"] != stored["input_deg"]:
            raise ValueError("Mismatched gravity sample")
        data.qpos[:] = state["qpos"]
        mujoco.mj_forward(model, data)
        reaction, pads = pad_reaction(model, data, pad_file["pads"])
        pad_matrix = np.column_stack([pad.pop("leverage")[indices] for pad in pads])
        demand = np.asarray(stored["gravity_demand_Nm"])
        matrix = np.asarray(stored["cable_leverage_mm"]) * .001
        np.testing.assert_allclose(demand, data.qfrc_bias[indices], atol=1e-9, rtol=0.)
        previous = stored["best_at_30N_per_line"]["max_abs_residual_Nm"]
        fit = best_balance(matrix, demand - reaction[indices], 30.)
        rows.append(dict(input_deg=state["input_deg"], humerus_deg=state["humerus_direction_change_deg"],
                         pads=pads, pad_coord_torque_Nm=reaction[indices].tolist(),
                         pad_leverage_mm=(1000. * pad_matrix).tolist(),
                         no_pad_max_residual_Nm=previous, with_pads=fit,
                         reduction_Nm=previous - fit["max_abs_residual_Nm"]))
    for expected, actual in zip(pad_file["pads"], rows[0]["pads"], strict=True):
        np.testing.assert_allclose(actual["reaction_N"], expected["neutral_force_N"], atol=1e-6, rtol=0.)
    travel = []
    for index, pad in enumerate(pad_file["pads"]):
        samples = [row["pads"][index] for row in rows]
        height_min = min(row["separation_mm"] for row in samples)
        height_max = max(row["separation_mm"] for row in samples)
        travel.append(dict(id=pad["candidate_index"], required_normal_travel_mm=height_max-height_min,
                           peak_lateral_offset_mm=max(row["lateral_offset_mm"] for row in samples),
                           additional_free_height_at_last_sample_mm=height_max - 1000. * pad["rest_separation_m"],
                           minimum_peak_force_to_just_touch_all_samples_N=(
                               (height_max - height_min) * pad["stiffness_N_per_m"] * .001)))
    # Keep the attachment points, normals and footprint from the audited model.
    pad = .001 * np.asarray([row["pad_leverage_mm"] for row in rows])
    cable = .001 * np.asarray([row["cable_leverage_mm"] for row in gravity["samples"]])
    demand = np.asarray([row["gravity_demand_Nm"] for row in gravity["samples"]])
    heights = np.asarray([[p["separation_mm"] for p in row["pads"]] for row in rows])
    if any(p["bottomed"] or p["lateral_offset_mm"] > 1000. * definition["radius_m"]
           for row in rows for p, definition in zip(row["pads"], pad_file["pads"], strict=True)):
        raise ValueError("Long-travel pad trial requires contact within the existing footprint")
    # Upper-bound the benefit of these normal directions by allowing each pad
    # reaction to vary independently, with no shared spring law or force cap.
    best_possible = []
    for cable_i, pad_i, demand_i in zip(cable, pad, demand, strict=True):
        best_possible.append(best_balance(cable_i, demand_i, 30., pad_i))
    spring_fits = {str(k): fit_long_travel_pads(cable, pad, demand, heights - heights[0], k, 30.)
                   for k in (1., 2., 3., 10.)}
    primary_index = next(i for i, definition in enumerate(pad_file["pads"])
                         if definition["candidate_index"] == 4)
    primary_only = fit_long_travel_pads(cable, pad[:, :, [primary_index]], demand,
                                       heights[:, [primary_index]] - heights[0, primary_index], 1., 30.)
    save(folder / "whole_arm_k7_scapular_support.json", dict(
        gravity=gravity_ref, pads=pad_ref, generator=provenance(ROOT / "tools/audit_whole_arm_k7_scapular_support.py"),
        samples=rows, travel_diagnostic=travel, ideal_independent_pad_balance=best_possible,
        long_travel_pad_trial=dict(pad_force_cap_N=30., cable_force_cap_N=30., common_stiffness_trials=spring_fits,
                                   pad_4_only=primary_only),
        scope="The existing three 12 mm candidate pads are modeled as unilateral normal springs. "
        "Contact requires positive compression inside the fixed rib patch perimeter, with no bottoming. "
        "Reuses the nominal gravity load and eight tension-only paths at 30 N per wire. "
        "This prescribed-pose calculation does not integrate motion, prove stable contact, validate pad material, "
        "or include TPU or cuff loading. The added long-travel fit tests one preload per existing pad "
        "across all poses, assuming the spring can physically provide the required travel. "
        "The ideal independent pad forces are uncapped to bound the benefit of the existing contact directions."))
    print("K7 pad engagement:", [sum(row["engaged"] for row in pose["pads"]) for pose in rows])


if __name__ == "__main__":
    main()
