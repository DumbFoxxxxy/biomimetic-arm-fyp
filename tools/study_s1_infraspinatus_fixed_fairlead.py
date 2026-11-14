"""Screen a supported spherical fairlead for the local FJ1500 tendon."""

from itertools import product

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, rod, save
from check_shoulder_front_liner_keepers import is_compliant
from study_elbow_series_package import bounded, overlaps, union
from study_s1_local_force_geometry import load_reference


def wrap(a, b, centre, radius):
    """Shortest exterior string around a sphere; return lengths and tangent forces."""
    da, db = a-centre, b-centre
    la, lb = np.linalg.norm(da), np.linalg.norm(db, axis=1)
    if la <= radius or np.min(lb) <= radius:
        return None
    u, v = da/la, db/lb[:, None]
    cosine = np.clip(v@u, -1., 1.)
    angle = np.arccos(cosine)
    start_angle, end_angle = np.arccos(radius/la), np.arccos(radius/lb)
    wrapped = angle > start_angle+end_angle
    sine = np.linalg.norm(v-cosine[:, None]*u, axis=1)
    if np.any(wrapped & (sine < 1e-5)):
        return None
    direction = (v-cosine[:, None]*u)/np.maximum(sine[:, None], 1e-12)
    first = centre+radius*(np.cos(start_angle)*u+np.sin(start_angle)*direction)
    last_angle = angle-end_angle
    last = centre+radius*(np.cos(last_angle)[:, None]*u+np.sin(last_angle)[:, None]*direction)
    direct = np.linalg.norm(b-a, axis=1)
    lengths = np.where(wrapped, np.sqrt(la**2-radius**2)+np.sqrt(lb**2-radius**2)+radius*(last_angle-start_angle), direct)
    toward_b = np.where(wrapped[:, None], last-b, a-b)
    toward_a = np.where(wrapped[:, None], first-a, b-a)
    toward_b /= np.linalg.norm(toward_b, axis=1)[:, None]
    toward_a /= np.linalg.norm(toward_a, axis=1)[:, None]
    return lengths, wrapped, toward_b, toward_a, (u, direction, start_angle, last_angle)


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source_file = folder / "FJ1500_humeral_exit_screen.json"
    source = read(source_file)
    revision = source["selected"]
    path = load_reference(source["baseline"])
    smooth = load_reference(path["parent"])
    coupling = load_reference(smooth["coupling"])
    old = load_reference(coupling["prior_path"])
    rest = load_reference(old["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    reference = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(reference["local_geometry"])
    assembly = load_reference(path["shoulder_assembly"])
    model = next(row for row in load_reference(reference["kinematics"])["models"] if row["head"] == "FJ1500")
    gh = np.array(geometry["gh_centre_source_mm"])
    flexion, _, longitudinal = np.array(geometry["axes_source"])
    a = np.array(model["neutral_crimp_outer_exit_A_source_mm"])+rest["rest_offset_source_mm"]
    bone_exit = np.array(model["neutral_crimp_outer_exit_B_source_mm"])+revision["parameters"]["shift_source_mm"]
    angles = np.arange(97.)
    rotations = Rotation.from_rotvec(angles[:, None]*np.pi/180.*flexion)
    b = gh+rotations.apply(bone_exit-gh)
    direct = np.linalg.norm(a-b, axis=1)
    base_arm = np.cross(b-gh, (a-b)/direct[:, None]) @ longitudinal
    np.testing.assert_allclose(base_arm, revision["parameters"]["longitudinal_moment_arm_mm"], atol=1e-8)

    owners = assembly["intended_body_owners"]
    solids = {key: load_component(row) for key, row in {**assembly["sources"], **revision["meshes"]}.items()
              if not is_compliant(key) and not key.endswith("/tendon_wire") and key.startswith(("bone:", "GH_", "collar:"))}
    scapula = solids["bone:scapula"]
    verts, faces = mesh_geom.manifold_mesh_to_arrays(scapula.to_mesh64())
    triangles = verts[faces]
    probes = list(range(0, 97, 8))
    obstacles = [bounded({key: solid.transform(rigid_transform(rotations[i], gh, gh)) if owners[key] == "humerus_ideal_attachment" else solid
                          for key, solid in solids.items()}) for i in probes]
    mf = A.I._manifold()
    wire_radius = .30
    trials = []
    max_support_reach = 40.
    for xyz, surface_radius in product(product(range(-45, 16, 5), range(-30, 41, 5), range(-30, 26, 5)), (2.5, 4.)):
        centre = gh+xyz
        result = wrap(a, b, centre, surface_radius+wire_radius)
        if result is None:
            continue
        length, contacted, toward_b, toward_a, _ = result
        if not contacted.any():
            continue
        arm = np.cross(b-gh, toward_b) @ longitudinal
        trials.append(dict(centre=centre.tolist(), radius_mm=surface_radius, wrapped_count=int(contacted.sum()),
                           length_span_mm=float(np.ptp(length)), minimum_longitudinal_arm_mm=float(arm.min()),
                           peak_guide_reaction_at_10N=float((10.*np.linalg.norm(toward_a+toward_b, axis=1)).max())))
    direct_span = float(np.ptp(direct))
    improving = sorted((row for row in trials if row["length_span_mm"] < direct_span and row["minimum_longitudinal_arm_mm"] > 0.),
                       key=lambda row: row["length_span_mm"])
    verified = []
    for row in improving:
        centre = np.array(row["centre"])
        radius = row["radius_mm"]
        distances, nearest = mesh_geom._closest_all(centre, triangles)
        foot = nearest[distances.argmin()]
        distance = float(distances.min())
        row["scapula_distance_mm"] = distance
        if distance <= radius or distance > max_support_reach:
            row["geometry"] = "outside_screened_support_reach"
            continue
        sphere = mf.Manifold.sphere(radius, 64).translate(centre)
        if (hits := overlaps(sphere, obstacles[0])):
            row["geometry"] = "sphere_intersects_neutral_structure"
            row["neutral_hits_mm3"] = hits
            continue
        direction = (centre-foot)/distance
        stem = rod(foot-1.5*direction, centre, 1.3)
        if (stem^scapula).volume() < 1.:
            row["geometry"] = "insufficient_scapula_contact_stock"
            continue
        support = sphere+stem
        result = wrap(a, b, centre, radius+wire_radius+.02)
        if result is None:
            continue
        length, contacted, toward_b, _, arc = result
        u, transverse, begin, end = arc
        findings = []
        wires = []
        for j, i in enumerate(probes):
            pts = np.array([a, b[i]])
            if contacted[i]:
                theta = np.linspace(begin, end[i], 25)
                pts = np.vstack([a, centre+(radius+wire_radius+.02)*(np.cos(theta)[:, None]*u+
                                                                      np.sin(theta)[:, None]*transverse[i]), b[i]])
            wire = A.D.outer_envelope(pts, wire_radius)
            wires.append(wire)
            support_hits = overlaps(support, {key: solid for key, solid in obstacles[j].items() if key != "bone:scapula"})
            wire_hits = overlaps(wire, {**obstacles[j], **bounded({"support_stem": stem-sphere})})
            if support_hits or wire_hits:
                findings.append(dict(angle_deg=i, support_hits_mm3=support_hits, wire_hits_mm3=wire_hits))
        row["geometry"] = "sampled"
        row["findings"] = findings
        row["scapula_contact_source_mm"] = foot.tolist()
        row["geometric_stem_length_mm"] = distance+1.5
        row["path_length_mm"] = length.tolist()
        row["longitudinal_arm_mm"] = (np.cross(b-gh, toward_b)@longitudinal).tolist()
        # Bypass the direct stem obstruction with a supported two-leg arm.
        # The intermediate point specifies rigid stock, never a cable waypoint.
        mount_end = foot-1.5*direction
        side_a = A.D.unit(np.cross(direction, np.eye(3)[np.argmin(abs(direction))]))
        side_b = np.cross(direction, side_a)
        checked_stems = 0
        bent = None
        for fraction, clearance, bearing_angle in product((.35,.55,.75), (4.,8.,12.), np.arange(0., 360., 30.)):
            theta = np.deg2rad(bearing_angle)
            middle = foot+fraction*(centre-foot)+clearance*(np.cos(theta)*side_a+np.sin(theta)*side_b)
            stem_stock = union([rod(mount_end, middle, 1.3), rod(middle, centre, 1.3),
                                mf.Manifold.sphere(1.3, 24).translate(middle)])
            if (stem_stock^scapula).volume() < 1.:
                continue
            checked_stems += 1
            envelope = stem_stock+sphere
            if any(overlaps(envelope, {key: solid for key, solid in targets.items() if key != "bone:scapula"})
                   or (wire ^ (stem_stock-sphere)).volume() > .001
                   or overlaps(wire, targets) for targets, wire in zip(obstacles, wires, strict=True)):
                continue
            bent = dict(neutral_stem_centreline_source_mm=[mount_end.tolist(), middle.tolist(), centre.tolist()],
                        stem_length_mm=float(np.linalg.norm(middle-mount_end)+np.linalg.norm(centre-middle)),
                        geometric_support_root_contact_mm3=float((stem_stock^scapula).volume()))
            break
        row["bent_stem_checks"] = checked_stems
        row["bent_stem_clear_candidate"] = bent
        verified.append(row)
    dense_results = []
    for row in verified:
        stem = row["bent_stem_clear_candidate"]
        if stem is None:
            continue
        root, elbow, centre = np.array(stem["neutral_stem_centreline_source_mm"])
        radius = row["radius_mm"]
        sphere = mf.Manifold.sphere(radius, 64).translate(centre)
        stock = union([rod(root, elbow, 1.3), rod(elbow, centre, 1.3), mf.Manifold.sphere(1.3, 24).translate(elbow)])
        full = stock+sphere
        _, contact, _, _, (u, transverse, begin, end) = wrap(a, b, centre, radius+wire_radius+.02)
        hits = []
        for i, rotation in enumerate(rotations):
            arm = rigid_transform(rotation, gh, gh)
            targets = bounded({key: solid.transform(arm) if owners[key] == "humerus_ideal_attachment" else solid
                               for key, solid in solids.items()})
            theta = np.linspace(begin, end[i], 33)
            arc = centre+(radius+wire_radius+.02)*(np.cos(theta)[:, None]*u+np.sin(theta)[:, None]*transverse[i])
            wire = A.D.outer_envelope(np.vstack([a, arc, b[i]]) if contact[i] else np.array([a,b[i]]), wire_radius)
            support_hits = overlaps(full, {k:v for k,v in targets.items() if k != "bone:scapula"})
            wire_hits = overlaps(wire, {**targets, **bounded({"fairlead": full})})
            if support_hits or wire_hits:
                hits.append(dict(angle_deg=i, support_hits_mm3=support_hits, wire_hits_mm3=wire_hits))
        dense_results.append(dict(centre_source_mm=centre.tolist(), radius_mm=radius, tested_angles=len(angles), hits=hits,
                                  span_mm=row["length_span_mm"], stem_length_mm=stem["stem_length_mm"]))
    save(folder / "FJ1500_fixed_fairlead_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_fixed_fairlead.py"), source=provenance(source_file),
        direct_length_span_mm=direct_span, baseline_minimum_longitudinal_arm_mm=float(base_arm.min()),
        direct_chord_by_degree_mm=direct.tolist(), max_scapular_surface_reach_mm=max_support_reach,
        cap_exit_source_mm=a.tolist(), neutral_bone_exit_source_mm=bone_exit.tolist(), grid_trials=trials,
        improving_preserving_torque_count=len(improving), geometric_samples_deg=probes,
        sampled_support_count=len(verified), candidates_without_sampled_hits=sum(not row["findings"] for row in verified),
        bent_stem_candidates_without_sampled_hits=sum(row["bent_stem_clear_candidate"] is not None for row in verified),
        dense_bent_stem_screen=dense_results,
        scope="Scapula-fixed sphere with stem intersecting scapular stock as a support concept; no fastener design."
              " Frozen cap, straight exterior string and spherical tangency ignore crimp exit angles, TPU, friction and installed motion."))
    print("Direct", direct_span, "wrapped", len(trials), "improving", len(improving), "support checks", len(verified),
          "clean", sum(not row["findings"] for row in verified), flush=True)


if __name__ == "__main__":
    main()
