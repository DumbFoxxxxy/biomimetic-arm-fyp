"""Check load-dependent cuff terminal candidates on the shared whole-arm poses."""

import argparse
from copy import deepcopy
from itertools import combinations

import mujoco
import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_series_elbow_integration import transforms
from audit_s1_y_retained_lead_clearance import rigid_transform, screen
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, rod, save
from project_runtime import load_mujoco_model
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--offset-infraspinatus", action="store_true")
    modes.add_argument("--shared-infraspinatus-support", action="store_true")
    modes.add_argument("--combined-k7", action="store_true")
    parser.add_argument("--include-moving-hardware", action="store_true")
    args = parser.parse_args()
    if args.include_moving_hardware and not args.combined_k7:
        parser.error("--include-moving-hardware requires --combined-k7")
    folder = DEST.parent / "whole_arm_first_pass"
    motion_ref = provenance(folder / "whole_arm_shared_motion.json")
    motion = load_reference(motion_ref)
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    base_ref = provenance(folder / "tabletop_base_check.json")
    base = load_reference(base_ref)
    if base["scene"] != motion["scene_check"]:
        raise ValueError("Base and motion use different assemblies")
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
               "base:T12": base["source"]}
    rigid = {key: load_component(source) for key, source in sources.items()}
    owners = {**transfer["owner_by_source"], "base:T12": "thorax"}
    origins = {key: np.asarray(scene["neutral_body_origins_mm"][key]) for key in set(owners.values())}
    model_file = ROOT / scene["model"]["file"]
    if provenance(model_file) != scene["model"]:
        raise ValueError("Shared model changed")
    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    eq_ref = provenance(DEST.parent / "s1_local_muscle_load/y_terminal_retained_lead_equilibrium.json")
    equilibrium = load_reference(eq_ref)
    historical = load_reference(equilibrium["shoulder_assembly"])
    historical_sources = historical["sources"]
    geometry = load_reference(equilibrium["local_geometry"])
    centre = np.array(geometry["gh_centre_source_mm"])
    np.testing.assert_allclose(centre, origins["humerus"], atol=1e-6, rtol=0.)
    ports = {row["head"]: row["ports"] for row in load_reference(geometry["port_geometry"])["routes"]}
    axes = np.asarray(geometry["axes_source"])
    stiffness = (.2, .4, 50., 1.)
    load_sequence = (0., 5., 10., 5., 0.)
    specifications = load_reference(equilibrium["kinematics"])["models"]
    paths = ("K7_SC1_pos", "GH_axial_positive", "GH_axial_negative")
    revision_ref = None
    installed_length_adjustment = 0.
    fixed_revision_overlap = {}
    shared_support_ref = None
    if args.offset_infraspinatus or args.shared_infraspinatus_support or args.combined_k7:
        revision_ref = provenance(DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief/FJ1500_humeral_exit_screen.json")
        selected = load_reference(revision_ref)["selected"]
        specifications = [deepcopy(row) for row in specifications if args.combined_k7 or row["head"] == "FJ1500"]
        shifted = next(row for row in specifications if row["head"] == "FJ1500")
        delta = np.array(selected["parameters"]["shift_source_mm"])
        for key in ("neutral_bone_pin_source_mm", "neutral_crimp_outer_exit_B_source_mm"):
            shifted[key] = (np.asarray(shifted[key]) + delta).tolist()
        strand = np.asarray(shifted["active_strand_centreline_source_mm"])
        strand[-2:] += delta
        shifted["active_strand_centreline_source_mm"] = strand.tolist()
        historical_sources = {**historical_sources, **selected["meshes"]}
        for key, source in selected["meshes"].items():
            if key not in {"link:FJ1500/bone_groove_ring", "link:FJ1500/bone_uncrimped_sleeve"}:
                if key == "collar:jaw_upper":
                    # The upstream candidate locates the jaw directly on the
                    # humerus; it is absent from the shared fixture inventory.
                    owners[key] = "humerus"
                if owners.get(key) != "humerus":
                    raise ValueError(f"Missing humeral support in shared assembly: {key}")
                rigid[key] = load_component(source)
        changed_keys = set(selected["meshes"])
        if args.shared_infraspinatus_support or args.combined_k7:
            shared_support_ref = provenance(folder / "whole_arm_infraspinatus_support.json")
            support = load_reference(shared_support_ref)
            adapter = "elbow:FJ1477_bone_adapter"
            if (support["replaced_adapter"] != sources[adapter]
                    or support["neutral_hits_mm3"]
                    or any(row["incremental_hits_mm3"] for row in support["K7_added_material_samples"])):
                raise ValueError("Shared support is not cleared on the current chassis")
            rigid[adapter] = load_component(support["source"])
            del rigid["collar:jaw_upper"]
            changed_keys.add(adapter)
        for key in sorted(changed_keys):
            if key in rigid:
                hits = overlaps(rigid[key], bounded({name: solid for name, solid in rigid.items() if name != key}))
                if hits:
                    fixed_revision_overlap[key] = hits
        installed_length_adjustment = selected["length_adjustment_mm"]
    combined_cuffs = {}
    moving_sources = {}
    moving_frames = {}
    if args.combined_k7:
        routed_ref = provenance(folder / "whole_arm_posterior_detour_audit.json")
        routed = load_reference(routed_ref)
        if routed["motion"] != motion_ref or fixed_revision_overlap:
            raise ValueError("Combined K7 requires matched paths and clear fixed hardware")
        paths = ("K7_SC1_pos",)
        path_limit = min(row["initial_clear_samples"] for row in routed["paths"][paths[0]]["summary"].values())
    cases = {}
    for spec in specifications:
        head = spec["head"]
        pin_a = np.asarray(spec["neutral_cap_pin_source_mm"])
        pin_b = np.asarray(spec["neutral_bone_pin_source_mm"])
        axis_a = np.asarray(spec["cap_pin_axis_source"])
        axis_b = np.asarray(spec["bone_pin_axis_source"])
        cap_keys = [key for key, owner in historical["intended_body_owners"].items()
                    if owner == f"free_terminal:{head}"]
        eye_keys = [[f"link:{head}/{end}_{part}" for part in ("groove_ring", "uncrimped_sleeve")]
                    for end in ("muscle", "bone")]
        moving = {key: load_component(historical_sources[key])
                  for key in (*cap_keys, *eye_keys[0], *eye_keys[1])}
        if args.include_moving_hardware:
            moving_sources.update({key: historical_sources[key] for key in moving})
        records = {}
        for name in paths:
            samples = motion["paths"][name]["accepted"]
            stride = 3 if name == "K7_SC1_pos" else 10
            indices = sorted({*range(0, len(samples), stride), len(samples) - 1})
            if args.combined_k7:
                indices = range(path_limit)
            result = []
            for index in indices:
                pose = samples[index]
                data.qpos[:] = pose["qpos"]
                mujoco.mj_forward(model, data)
                frames = transforms(model, data, origins)
                scap = frames["scapula"]
                rotation = Rotation.from_matrix(scap[:, :3].T @ frames["humerus"][:, :3])
                placed = {key: solid.transform(frames[owners[key]]) for key, solid in rigid.items()}
                obstacles = bounded(placed)
                try:
                    states = solve_branch(spec, ports[head], centre, rotation, stiffness, load_sequence,
                                          retain_leads=True, installed_length_offset=(
                                              installed_length_adjustment if head == "FJ1500" else 0.))
                except RuntimeError as error:
                    result.append(dict(input_deg=pose["input_deg"], solver_error=str(error)))
                    continue
                rows = []
                for step, state in enumerate(states):
                    cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
                    angle_a, angle_b = np.deg2rad(state["swivel_angles_deg"])
                    moved_pin = pin_a + state["cap_translation_source_mm"]
                    cap_pose = rigid_transform(cap, pin_a, moved_pin)
                    eye_a = rigid_transform(cap * Rotation.from_rotvec(axis_a * angle_a), pin_a, moved_pin)
                    eye_b = rigid_transform(rotation * Rotation.from_rotvec(axis_b * angle_b), pin_b,
                                            centre + rotation.apply(pin_b - centre))
                    hardware = dict(placed)
                    for keys, frame in ((cap_keys, cap_pose), (eye_keys[0], eye_a), (eye_keys[1], eye_b)):
                        hardware.update({key: moving[key].transform(frame).transform(scap) for key in keys})
                        if args.include_moving_hardware and state["commanded_muscle_tension_N"] == 10.:
                            matrix = np.column_stack((scap[:, :3] @ frame[:, :3],
                                                      scap[:, :3] @ frame[:, 3] + scap[:, 3]))
                            moving_frames.setdefault(index, {}).update({key: matrix.tolist() for key in keys})
                    tips = np.asarray(state["lead_tips_source_mm"])
                    tips = tips @ scap[:, :3].T + scap[:, 3]
                    exit_a, exit_b = np.asarray(state["tendon_exits_source_mm"])
                    span = float(np.linalg.norm(tips[1] - tips[0]))
                    # Probe the tensioned central span only. The real cable also has lead bends.
                    wire_hits = (overlaps(A.D.outer_envelope(tips, .30), obstacles)
                                 if state["short_tendon_tension_N"] > 1e-5 else None)
                    if args.combined_k7:
                        combined_cuffs.setdefault((index, step), {})[head] = (
                            {key: hardware[key] for key in moving},
                            A.D.outer_envelope(tips, .30) if wire_hits is not None else None,
                            tips.tolist())
                    rows.append(dict(command_N=state["commanded_muscle_tension_N"],
                                     tendon_N=state["short_tendon_tension_N"], central_chord_mm=span,
                                     retained_lead_lengths_mm=[float(np.linalg.norm(x - y))
                                                               for x, y in zip((exit_a, exit_b), state["lead_tips_source_mm"])],
                                     humeral_torque_axis_Nm=(axes @ np.asarray(state["GH_torque_source_Nm"])).tolist(),
                                     terminal_shift_mm=float(np.linalg.norm(state["cap_translation_source_mm"])),
                                     cap_affine_source_mm=cap_pose.tolist(),
                                     cap_pin_neutral_source_mm=pin_a.tolist(),
                                     origin_port_source_mm=ports[head]["origin"]["inner_face_source_mm"],
                                     stationarity_residual_N=state["stationarity_residual_N"],
                                     central_wire_hits_mm3=wire_hits,
                                     **screen(hardware, set(moving))))
                result.append(dict(input_deg=pose["input_deg"], states=rows))
            records[name] = result
            print(head, name, len(result), "poses", flush=True)
        cases[head] = records
    if args.combined_k7:
        combined = []
        for index, pose in enumerate(motion["paths"][paths[0]]["accepted"][:path_limit]):
            data.qpos[:] = pose["qpos"]
            mujoco.mj_forward(model, data)
            frames = transforms(model, data, origins)
            obstacles = bounded({key: solid.transform(frames[owners[key]]) for key, solid in rigid.items()})
            wires = {}
            for key, wire in routed["routes"].items():
                points = np.asarray(wire["points_source_mm"])
                a, b = (frames[owner] for owner in wire["owners"])
                w = np.asarray(wire["distal_weights"])[:, None]
                points = (1. - w) * (points @ a[:, :3].T + a[:, 3]) + w * (points @ b[:, :3].T + b[:, 3])
                wires[key] = rod(*points, .30) if len(points) == 2 else A.D.outer_envelope(points, .30)
            baseline_hits = {key: hits for key, shape in wires.items() if (hits := overlaps(shape, obstacles))}
            cable_pair_hits = [dict(left=a, right=b, mm3=volume)
                               for (a, x), (b, y) in combinations(wires.items(), 2)
                               if (volume := A.I.intersection_volume(x, y)) > .001]
            steps = []
            for step, load in enumerate(load_sequence):
                cuffs = combined_cuffs.get((index, step), {})
                if len(cuffs) != 2:
                    steps.append(dict(load_N=load, missing_cuffs=sorted({"FJ1500", "FJ1504"} - cuffs.keys())))
                    continue
                first, second = (cuffs[key] for key in ("FJ1500", "FJ1504"))
                first_parts = {**first[0], **({"FJ1500": first[1]} if first[1] is not None else {})}
                second_parts = {**second[0], **({"FJ1504": second[1]} if second[1] is not None else {})}
                cross = {label: hits for label, shape in first_parts.items()
                         if (hits := overlaps(shape, bounded(second_parts)))}
                neighboring = {label: hits for label, shape in {**first_parts, **second_parts}.items()
                               if (hits := overlaps(shape, bounded(wires)))}
                own = [cases[head][paths[0]][index]["states"][step] for head in ("FJ1500", "FJ1504")]
                steps.append(dict(load_N=load, cuff_cross_hits_mm3=cross, cuff_to_neighbor_wires_mm3=neighboring,
                                  local_cuff_hits=any(row["overlaps"] or row["central_wire_hits_mm3"] for row in own),
                                  tensioned_span_count=sum(cuff[1] is not None for cuff in cuffs.values()),
                                  central_tendon_tips_world_mm={key: item[2] for key, item in cuffs.items()
                                                                if item[1] is not None}))
            combined.append(dict(input_deg=pose["input_deg"], humerus_deg=pose["humerus_direction_change_deg"],
                                 baseline_wire_hits_mm3=baseline_hits, baseline_wire_pair_hits=cable_pair_hits,
                                 steps=steps))
            print("Combined K7", pose["input_deg"], len(baseline_hits), flush=True)
        if args.include_moving_hardware:
            baseline_ref = provenance(folder / "whole_arm_combined_cuff_k7.json")
            baseline = load_reference(baseline_ref)
            if baseline["combined"] != combined or baseline["support"] != shared_support_ref:
                raise ValueError("Replayed cuff data does not match the checked K7 run")
            save(folder / "whole_arm_k7_moving_cuff_geometry.json", dict(
                combined=baseline_ref, moving_sources=moving_sources,
                frames=[dict(input_deg=row["input_deg"], matrices=moving_frames[index])
                        for index, row in enumerate(combined)],
                generator=provenance(ROOT / "tools/audit_whole_arm_floating_terminals.py")))
            return
        save(folder / "whole_arm_combined_cuff_k7.json", dict(
            motion=motion_ref, support=shared_support_ref, routed=routed_ref, cases=cases, combined=combined,
            generator=provenance(ROOT / "tools/audit_whole_arm_floating_terminals.py"),
            scope="Both independently solved cuff local equilibria overlaid on identical prescribed poses. "
                  "Checks moving cuff hardware and tensioned chord against neighboring cables and each other; "
                  "includes the shared support. Does not solve contact-coupled loading, slack cable, TPU or full-arm dynamics."))
        return
    filename = ("whole_arm_shared_infraspinatus_audit.json" if args.shared_infraspinatus_support else
                "whole_arm_offset_infraspinatus_audit.json" if args.offset_infraspinatus
                else "whole_arm_floating_terminal_audit.json")
    save(folder / filename, dict(
        motion=motion_ref, base=base_ref, original_kinematics=eq_ref, stiffness=list(stiffness),
        humeral_revision=revision_ref, installed_length_adjustment_mm=installed_length_adjustment,
        fixed_revision_overlap_mm3=fixed_revision_overlap,
        shared_support=shared_support_ref,
        loads_N=list(load_sequence), axes=geometry["axes"], cases=cases,
        generator=provenance(ROOT / "tools/audit_whole_arm_floating_terminals.py"),
        scope="One cuff terminal assessed at a time with its original free cap and swivel hardware "
              "reintroduced on the shared chassis. Conditional, uncalibrated local static equilibrium. "
              "Straight tensioned central tendon chord only; smooth lead bends, its TPU belly, mutual "
              "cuff interactions and the other cables are outside this screen. A slack chord has no "
              "validated shape. Repeated tension levels sample a load cycle at fixed joint pose."))


if __name__ == "__main__":
    main()
