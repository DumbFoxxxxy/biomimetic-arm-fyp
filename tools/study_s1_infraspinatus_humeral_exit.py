"""Screen collar-supported FJ1500 exit relocation and shoulder leverage."""

from copy import deepcopy

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from check_shoulder_front_liner_keepers import is_compliant
from study_elbow_series_package import bounded, overlaps, union
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_smooth_tendon import centreline
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1500_outward_path_screen.json"
    path = read(path_file)
    smooth = load_reference(path["parent"])
    coupling = load_reference(smooth["coupling"])
    old = load_reference(coupling["prior_path"])
    rest = load_reference(old["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    reference = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(reference["local_geometry"])
    assembly = load_reference(path["shoulder_assembly"])
    model = next(row for row in load_reference(reference["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(path["input_route"])["routes"] if row["head"] == HEAD)
    sleeve_key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    sleeve = load_reference(old["sleeve"])
    if assembly["sources"][sleeve_key] != sleeve["original"]:
        raise ValueError("Mismatched assembly")
    sources = {**assembly["sources"], sleeve_key: sleeve["selected"]["source"]}
    parts = {key: load_component(source) for key, source in sources.items()
             if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
             and not key.endswith("/tendon_wire")}
    owners = assembly["intended_body_owners"]
    jaw_key = "collar:jaw_upper"
    jaw = parts[jaw_key]
    moving_keys = sorted(k for k in parts if k.startswith(f"link:{HEAD}/bone_"))
    if len(moving_keys) != 6:
        raise ValueError("Expected the complete bone-side eye and bolt package")
    cap_owned = {k for k in parts if owners[k] == f"free_terminal:{HEAD}" or k.startswith(f"link:{HEAD}/muscle_")}
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion, _, longitudinal = np.array(geometry["axes_source"])
    pin_a, pin_b = (np.array(model[f"neutral_{side}_pin_source_mm"]) for side in ("cap", "bone"))
    axis_a, axis_b = (np.array(model[f"{side}_pin_axis_source"]) for side in ("cap", "bone"))
    cap_bias = np.array(rest["rest_offset_source_mm"])
    adjusted, ports = rest_geometry(model, route["ports"], cap_bias)
    origin, stop = (np.array(route["ports"][k][face]) for k, face in
                    (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm")))
    coupling_axis = Rotation.from_rotvec(A.D.unit(origin-stop)*np.deg2rad(path["selected"]["coupling_roll_deg"])).apply(path["direction_source"])
    rest_exit_a = np.array(adjusted["neutral_crimp_outer_exit_A_source_mm"])
    rest_exit_b = np.array(model["neutral_crimp_outer_exit_B_source_mm"])
    angles = np.arange(97.)
    rot = Rotation.from_rotvec(angles[:, None]*flexion*np.pi/180.)

    def kinematics(shift):
        posed_b = centre + rot.apply(rest_exit_b+shift-centre)
        chord = rest_exit_a-posed_b
        lengths = np.linalg.norm(chord, axis=1)
        levers = np.cross(posed_b-centre, chord/lengths[:, None])
        return dict(neutral_chord_mm=float(lengths[0]), fixed_cap_chord_span_mm=float(np.ptp(lengths)),
                    longitudinal_moment_arm_mm=(levers@longitudinal).tolist(), flexion_moment_arm_mm=(levers@flexion).tolist())

    base = kinematics(np.zeros(3))
    stationary = {key: part.translate(cap_bias) if key in cap_owned else part for key, part in parts.items()
                  if key not in moving_keys and key != jaw_key}
    base_axis = A.D.unit(-longitudinal + (longitudinal @ axis_b)*axis_b)
    side = np.cross(axis_b, base_axis)
    cases = []
    geometries = {}
    for distance in (4., 8., 12., 16., 20.):
        for turn in (0., -20., 20., -40., 40.):
            angle = np.deg2rad(turn)
            shift = distance*(np.cos(angle)*base_axis + np.sin(angle)*side)
            p = pin_b + shift
            # Twin webs connect the new bearing ring to the retained jaw.
            bearing = rod(p-1.5*axis_b, p+1.5*axis_b, 3.2)
            webs = [rod(pin_b+s*.9*axis_b, p+s*.9*axis_b, 1.45) for s in (-1., 1.)]
            stock = union([bearing, *webs])
            if (stock ^ jaw).volume() < 1.:
                continue
            bore = rod(p-5.*axis_b, p+5.*axis_b, 1.10)
            moved = {key: parts[key].translate(shift) for key in moving_keys}
            # Remove interference from the added webs at washers, nut and
            # swivel, while retaining the entire new bearing annulus.
            guard = bearing-bore
            pad = A.I._manifold().Manifold.cube([.24]*3, True)
            clearance = union(part.minkowski_sum(pad) for part in moved.values())
            revised = (((jaw+stock)-bore)-((clearance-guard)-jaw)).as_original().simplify(1e-6)
            if (guard-revised).volume() > .001:
                raise ValueError("Relief cut entered the new bearing")
            external = bounded(stationary)
            hits = overlaps(revised-jaw, external)
            pack_hits = {key: overlap for key, part in moved.items()
                         if (overlap := overlaps(part, {**external, **bounded({jaw_key: revised})}))}
            stats = kinematics(shift)
            connected = sum(piece.volume() > 1e-5 for piece in revised.decompose()) == 1
            row = dict(offset_mm=distance, roll_deg=turn, shift_source_mm=shift.tolist(), new_pin_source_mm=p.tolist(),
                       connected=connected, added_support_volume_mm3=float((revised-jaw).volume()),
                       removed_old_jaw_volume_mm3=float((jaw-revised).volume()),
                       bracket_hits_mm3=hits, package_hits_mm3=pack_hits, **stats)
            cases.append(row)
            if connected and not hits and not pack_hits:
                geometries[distance, turn] = revised, moved
            print(distance, turn, round(stats["fixed_cap_chord_span_mm"], 2), len(hits), len(pack_hits), flush=True)
    viable = [row for row in cases if (row["offset_mm"], row["roll_deg"]) in geometries
              and row["longitudinal_moment_arm_mm"][0]*base["longitudinal_moment_arm_mm"][0] > 0.]
    verified = []
    for trial in sorted(viable, key=lambda row: row["fixed_cap_chord_span_mm"])[:3]:
        shift = np.array(trial["shift_source_mm"])
        candidate, moved = geometries[trial["offset_mm"], trial["roll_deg"]]
        moved_model = deepcopy(adjusted)
        for key in ("neutral_bone_pin_source_mm", "neutral_crimp_outer_exit_B_source_mm"):
            moved_model[key] = (np.array(moved_model[key])+shift).tolist()
        strand = np.array(moved_model["active_strand_centreline_source_mm"])
        strand[-2:] += shift
        moved_model["active_strand_centreline_source_mm"] = strand.tolist()
        length_offset = trial["neutral_chord_mm"]-base["neutral_chord_mm"]
        installed = {**parts, **moved, jaw_key: candidate}
        rows = []
        for angle in (0., 30., 60., 68., 75., 85., 96.):
            gh = Rotation.from_rotvec(flexion*np.deg2rad(angle))
            states = solve_branch(moved_model, ports, centre, gh, bias["stiffness"], [0., 5., 10.],
                                  retain_leads=True, installed_length_offset=length_offset,
                                  lateral_coupling=(coupling_axis, path["selected"]["gain"]), strand_bending_Nmm2=10.)
            arm_pose = rigid_transform(gh, centre, centre)
            static_obstacles = bounded({key: solid for key, solid in stationary.items()
                                        if owners[key] != "humerus_ideal_attachment" and key not in cap_owned})
            moving_support_hits = overlaps((candidate-jaw).transform(arm_pose), static_obstacles)
            for state in states:
                cap, hits = collision_state(installed, owners, state, cap_bias, gh, centre, pin_a, pin_b+shift, axis_a, axis_b)
                rows.append(dict(GH_deg=angle, cap_pin_source_mm=cap.tolist(), rigid_hits=hits,
                                 added_support_hits_mm3=moving_support_hits, **state))
        verified.append(dict(parameters=trial, length_adjustment_mm=length_offset, cases=rows,
                             maximum_cap_translation_mm=max(np.linalg.norm(row["cap_translation_source_mm"]) for row in rows),
                             failing_states=sum(bool(row["rigid_hits"] or row["added_support_hits_mm3"]) for row in rows)))
    selected = next((row for row in verified if not row["failing_states"]), None)
    if selected:
        p = selected["parameters"]
        revised, moved = geometries[p["offset_mm"], p["roll_deg"]]
        meshes = {}
        for key, part in {jaw_key: revised, **moved}.items():
            _, meshes[key] = export(part, folder / f"FJ1500_exit_{key.split('/')[-1].split(':')[-1]}.obj", pin_b)
        selected["meshes"] = meshes
        reloaded = {**parts, **{key: load_component(row) for key, row in meshes.items()}}
        shift = np.array(p["shift_source_mm"])
        revised_model = deepcopy(adjusted)
        for key in ("neutral_bone_pin_source_mm", "neutral_crimp_outer_exit_B_source_mm"):
            revised_model[key] = (np.array(revised_model[key])+shift).tolist()
        strand = np.array(revised_model["active_strand_centreline_source_mm"])
        strand[-2:] += shift
        revised_model["active_strand_centreline_source_mm"] = strand.tolist()
        records = []
        for angle in range(97):
            gh = Rotation.from_rotvec(flexion*np.deg2rad(angle))
            state = solve_branch(revised_model, ports, centre, gh, bias["stiffness"], [10.], retain_leads=True,
                                 installed_length_offset=selected["length_adjustment_mm"],
                                 lateral_coupling=(coupling_axis, path["selected"]["gain"]), strand_bending_Nmm2=10.)[0]
            pin, terminal_hits = collision_state(reloaded, owners, state, cap_bias, gh, centre, pin_a, pin_b+shift, axis_a, axis_b)
            arm = rigid_transform(gh, centre, centre)
            cap_rotation = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
            turn_a, turn_b = np.deg2rad(state["swivel_angles_deg"])
            cap_pose = rigid_transform(cap_rotation, pin_a, pin)
            eye_a = rigid_transform(cap_rotation*Rotation.from_rotvec(axis_a*turn_a), pin_a, pin)
            eye_b = rigid_transform(gh*Rotation.from_rotvec(axis_b*turn_b), pin_b+shift, centre+gh.apply(pin_b+shift-centre))
            posed = {}
            for key, part in reloaded.items():
                frame = None
                if key in (f"link:{HEAD}/muscle_groove_ring", sleeve_key):
                    frame = eye_a
                elif key in (f"link:{HEAD}/bone_groove_ring", f"link:{HEAD}/bone_uncrimped_sleeve"):
                    frame = eye_b
                elif owners[key] == f"free_terminal:{HEAD}":
                    frame = cap_pose
                elif owners[key] == "humerus_ideal_attachment":
                    frame = arm
                posed[key] = part.transform(frame) if frame is not None else part
            rigid_targets = bounded(posed)
            support_targets = {key: part for key, part in rigid_targets.items() if key != jaw_key}
            support_hits = overlaps((reloaded[jaw_key]-jaw).transform(arm), support_targets)
            fastener_hits = {key: hit for key in moving_keys
                             if (hit := overlaps(posed[key], {other: part for other, part in rigid_targets.items()
                                                             if other != key}))}
            exits, tips = (np.array(state[key]) for key in ("tendon_exits_source_mm", "lead_tips_source_mm"))
            h = np.linalg.norm(tips[1]-tips[0])/3.
            points = centreline([tips[0], tips[0]+h*A.D.unit(tips[0]-exits[0]),
                                 tips[1]+h*A.D.unit(tips[1]-exits[1]), tips[1]], 65)
            wire_hits = overlaps(A.D.outer_envelope(points, .30), rigid_targets)
            records.append(dict(GH_deg=angle, terminal_hits=terminal_hits, added_support_hits_mm3=support_hits,
                                fastener_hits_mm3=fastener_hits, central_wire_hits_mm3=wire_hits,
                                cap_translation_mm=float(np.linalg.norm(state["cap_translation_source_mm"])),
                                longitudinal_torque_Nm=float(longitudinal @ state["GH_torque_source_Nm"])))
        selected["readback_10N_flexion"] = records
    save(folder / "FJ1500_humeral_exit_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_humeral_exit.py"),
        solver=provenance(ROOT / "tools/study_s1_y_terminal_equilibrium.py"), baseline=provenance(path_file),
        original_pin_source_mm=pin_b.tolist(), original_geometry=sources[jaw_key], baseline_fixed_cap=base,
        trials=cases, checked_free_end_cases=verified, selected=selected,
        scope="Integral bearing and webs on original upper collar jaw, with complete moved bone-side eye package."
              " Fixed-cap chord is a geometry sensitivity, not actuator stroke."
              " Prescribed shoulder poses. Central strand checked at 10 N. TPU deformation, free drive and strength remain unverified."))
    print("Physically mounted neutral candidates", len(viable), "selected", selected is not None, flush=True)


if __name__ == "__main__":
    main()
