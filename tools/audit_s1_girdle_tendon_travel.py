"""Measure shoulder endpoint travel on the rigid-checked SC/AC path."""

import argparse

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from simulate_elbow_antagonistic_load import write_csv
from study_s1_coupled_scapula import ShoulderAssembly
from study_s1_local_force_geometry import load_reference, span


def place(transform, point):
    return transform[:, :3] @ point + transform[:, 3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--combine-flexion", action="store_true",
                        help="Probe SC/AC with additional relative GH flexion")
    args = parser.parse_args()
    local = DEST.parent / "s1_local_muscle_load"
    path_ref = provenance(local / "FJ1500_cap_relief/FJ1504_AC_blind_coupled_motion.json")
    endpoints_ref = provenance(local / "gh_endpoint_geometry.json")
    rigid = load_reference(path_ref)
    endpoints = load_reference(endpoints_ref)
    origins = load_reference(rigid["pose_reference"])["neutral_origins_m"]
    model = ShoulderAssembly({}, {}, origins)
    np.testing.assert_allclose(model.gh, endpoints["gh_centre_source_mm"], atol=.001, rtol=0.)
    candidate = rigid["cases"]["short_right_screw"]
    forward = candidate["directions"]["1"]
    if candidate["neutral_hits"] or forward["reverse_replay_hits"] or len(forward["path"]) < 2:
        raise ValueError("No clear outward and return SC/AC path")
    path = forward["path"]
    output = local / "FJ1504_AC_blind_shoulder_span_travel.csv"
    if args.combine_flexion:
        combined_flexion(model, path, endpoints, path_ref, endpoints_ref,
                         output.with_name("FJ1468_girdle_flexion_span.csv"))
        return
    rows = []
    summary = {}
    for name, channel in endpoints["neutral"].items():
        origin, distal = (np.array(channel[key]) for key in ("proximal_port_source_mm", "distal_port_source_mm"))
        owner = channel["origin_body"]
        neutral = np.linalg.norm(origin - distal)
        samples = []
        for pose in path:
            angle = pose["scapular_coronal_deg"]
            sc = np.array(pose["SC_rotvec_deg"])
            clavicle, scapula, ac = model.kinematics(angle, sc)
            np.testing.assert_allclose(place(scapula, model.gh) - model.gh, pose["GH_displacement_mm"],
                                       atol=1e-6, rtol=0.)
            proximal_frame = {"clavicle": clavicle, "scapula": scapula}[owner]
            origin_world = place(proximal_frame, origin)
            distal_world = place(scapula, distal)
            length, direction, arm = span(origin_world, distal_world, place(scapula, model.gh))
            # Transform the reference GH axes with the scapula to report comparable components.
            components = np.asarray(endpoints["axes_source"]) @ (scapula[:, :3].T @ arm)
            ac_arm = (np.cross(distal_world - place(clavicle, model.ac), direction)
                      if owner == "clavicle" else np.zeros(3))
            samples.append([length, length - neutral, *components, *ac_arm])
            rows.append([name, owner, angle, *sc, np.rad2deg(ac.magnitude()), *samples[-1]])
        values = np.array(samples)
        summary[name] = dict(origin_body=owner, neutral_span_mm=float(neutral),
                             span_change_bounds_mm=[float(values[:, 1].min()), float(values[:, 1].max())],
                             final_span_change_mm=float(values[-1, 1]),
                             span_excursion_mm=float(np.ptp(values[:, 0])),
                             ac_source_y_moment_arm_mm=[float(values[:, -2].min()), float(values[:, -2].max())])
    if not rows or not np.isfinite(np.array([row[2:] for row in rows])).all():
        raise ValueError("Missing or invalid SC/AC span samples")
    write_csv(output, ["channel", "origin_body", "scapula_coronal_deg", "sc_x_deg", "sc_y_deg", "sc_z_deg",
                       "ac_relative_angle_deg", "span_mm", "delta_from_neutral_mm", "flexion_arm_mm",
                       "abduction_arm_mm", "longitudinal_arm_mm", "ac_source_x_arm_mm",
                       "ac_source_y_arm_mm", "ac_source_z_arm_mm"], rows)
    save(output.with_suffix(".json"), dict(
        generator=provenance(ROOT / "tools/audit_s1_girdle_tendon_travel.py"),
        rigid_path=path_ref, local_endpoints=endpoints_ref, table=provenance(output),
        sampled_poses=len(path), final_scapula_coronal_deg=path[-1]["scapular_coronal_deg"],
        final_sc_rotvec_deg=path[-1]["SC_rotvec_deg"], channels=summary,
        scope="Endpoint chords at fixed relative GH and elbow coordinates. Evaluated on the saved "
              "positive rigid SC/AC path. AC arms hold SC fixed and use global source axes. "
              "Return traverses the same poses in reverse. Uses baseline "
              "fixed-endpoint ports; Y free terminals and actuator-side routing are excluded. "
              "This is a geometric adjustment estimate, not flexible contact or driven motion.",
    ))
    print(f"Measured {len(path)} SC/AC poses for {len(summary)} shoulder channels")


def combined_flexion(model, path, endpoints, path_ref, endpoints_ref, output):
    channel = endpoints["neutral"]["FJ1468"]
    if channel["origin_body"] != "clavicle":
        raise ValueError("FJ1468 proximal port must follow the clavicle")
    origin, distal = (np.array(channel[key]) for key in ("proximal_port_source_mm", "distal_port_source_mm"))
    flexion_axis = np.array(endpoints["axes_source"][0])
    angles = np.arange(31.)
    distals = model.gh + Rotation.from_rotvec(np.deg2rad(angles[:, None]) * flexion_axis).apply(distal - model.gh)
    gh_only = np.linalg.norm(distals - origin, axis=1)
    rows = []
    by_flexion = []
    for pose in path:
        angle, sc = pose["scapular_coronal_deg"], np.array(pose["SC_rotvec_deg"])
        clavicle, scapula, _ = model.kinematics(angle, sc)
        np.testing.assert_allclose(place(scapula, model.gh) - model.gh, pose["GH_displacement_mm"],
                                   atol=1e-6, rtol=0.)
        origin_world = place(clavicle, origin)
        world = distals @ scapula[:, :3].T + scapula[:, 3]
        lengths = np.linalg.norm(world - origin_world, axis=1)
        additional = lengths - gh_only
        by_flexion.append(additional)
        for flexion, length, delta in zip(angles, lengths, additional, strict=True):
            rows.append([angle, *sc, flexion, length, length - gh_only[0], delta, delta - additional[0]])
    extra = np.asarray(by_flexion)
    np.testing.assert_allclose(extra[0], 0., atol=1e-9, rtol=0.)
    if not np.isfinite(extra).all():
        raise ValueError("Invalid combined endpoint geometry")
    peak = np.unravel_index(np.argmax(extra), extra.shape)
    write_csv(output, ["scapula_coronal_deg", "sc_x_deg", "sc_y_deg", "sc_z_deg", "gh_flexion_deg",
                       "span_mm", "delta_from_neutral_mm", "extra_vs_same_gh_fixed_girdle_mm",
                       "deviation_from_additive_mm"], rows)
    save(output.with_suffix(".json"), dict(
        generator=provenance(ROOT / "tools/audit_s1_girdle_tendon_travel.py"), rigid_path=path_ref,
        local_endpoints=endpoints_ref, table=provenance(output), sampled_girdle_poses=len(path),
        gh_flexion_grid_deg=angles.tolist(), sample_rows=len(rows),
        maximum_extra_vs_same_gh_fixed_girdle_mm=float(extra[peak]),
        maximum_extra_at_deg=[path[peak[0]]["scapular_coronal_deg"], float(angles[peak[1]])],
        final_girdle_extra_vs_same_gh_mm=extra[-1].tolist(),
        maximum_absolute_deviation_from_additive_mm=float(np.max(np.abs(extra - extra[:, :1]))),
        scope="Prescribed 0–30 degree GH flexion probes about the saved local flexion axis. "
              "Only the girdle path at fixed neutral GH was rigid-checked in the source. "
              "Combined poses, real cable arc, force and actuator payout are unverified.",
    ))
    print(f"Measured {len(rows)} combined endpoint probes")


if __name__ == "__main__":
    main()
