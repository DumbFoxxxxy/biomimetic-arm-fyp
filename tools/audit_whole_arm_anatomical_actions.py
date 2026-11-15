"""Name the axial probes and check routed shoulder and scapular actions."""

import argparse
from itertools import combinations

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import find_hits, transforms
from audit_whole_arm_k7_gravity_balance import pull_jacobian
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import A, load_component, rod, save
from project_runtime import load_mujoco_model
from study_s1_local_force_geometry import load_reference


def unit(vector):
    return vector / np.linalg.norm(vector)


def centre(part):
    bounds = np.asarray(part.bounding_box())
    return (bounds[:3] + bounds[3:]) / 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shoulder-elevation", action="store_true", help="Check the current wires through K2 and K3")
    args = parser.parse_args()
    folder = DEST.parent / "whole_arm_first_pass"
    combined_ref = provenance(folder / "whole_arm_combined_cuff_k7.json")
    combined = load_reference(combined_ref)
    motion = load_reference(combined["motion"])
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    routes = load_reference(combined["routed"])["routes"]
    support = load_reference(combined["support"])
    hardware = load_reference(load_reference(support["audit"])["humeral_revision"])["selected"]["meshes"]
    base_ref = provenance(folder / "tabletop_base_check.json")
    base = load_reference(base_ref)
    if base["motion"] != combined["motion"]:
        raise ValueError("Base and joint poses do not match")
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
               "elbow:FJ1477_bone_adapter": support["source"], "base:T12": base["source"]}
    sources.update({key: source for key, source in hardware.items() if key.startswith("link:") and key in sources})
    bones = {key: load_component(source) for key, source in sources.items()}
    owners = {**transfer["owner_by_source"], "base:T12": "thorax"}
    bone_keys = list(bones)
    keys = bone_keys + list(routes)
    count = len(bones)
    pairs = np.asarray([(i, j) for i, j in combinations(range(len(keys)), 2)
                        if j >= count or owners[keys[i]] != owners[keys[j]]])
    girdle_lines = {"FJ1459_g3", "FJ1536_g2"}
    girdle_pair_indices = [(i, j) for i, j in pairs if keys[i] in girdle_lines or keys[j] in girdle_lines]
    origins = {key: np.asarray(value) for key, value in scene["neutral_body_origins_mm"].items()}
    model_path = ROOT / scene["model"]["file"]
    if provenance(model_path) != scene["model"]:
        raise ValueError("Shared arm model changed")
    model = load_mujoco_model(model_path)
    data = mujoco.MjData(model)
    axes = load_reference(motion["anatomical_axes"])["axes_source"]
    mujoco.mj_forward(model, data)
    sternum = centre(bones["thorax_reference:FJ3290"])
    spine = centre(bones["thorax_reference:FJ3166"])
    anterior_delta = sternum - spine
    anterior = np.array([0., np.sign(anterior_delta[1]), 0.])
    if abs(anterior_delta[1]) < 3. * np.linalg.norm(anterior_delta[[0, 2]]):
        raise ValueError("Source Y axis does not align with the sternum-to-spine direction")
    # The neutral path has a near-straight elbow. Use its flexion axis to
    # construct a transverse hand-direction marker for the axial sign test.
    humerus = np.asarray(axes["longitudinal"])
    elbow_axis = data.xaxis[model.joint("elbow").id].copy()
    hand = unit(np.cross(elbow_axis, humerus))
    medial = np.array([np.sign(sternum[0] - origins["humerus"][0]), 0., 0.])
    if hand @ anterior < .9:
        raise ValueError("Bent-elbow direction marker does not face anteriorly in neutral")
    axial = {}
    for name in ("GH_axial_positive", "GH_axial_negative"):
        path = motion["paths"][name]
        data.qpos[:] = path["accepted"][1]["qpos"]
        mujoco.mj_forward(model, data)
        turned_hand = data.body("humerus").xmat.reshape(3, 3) @ hand
        medial_change = float((turned_hand - hand) @ medial)
        action = "K4_internal" if medial_change > 0 else "K5_external"
        axial[action] = dict(source_path=name, marker_medial_change_at_first_step=medial_change,
                             rigid_range_deg=abs(path["accepted"][-1]["input_deg"]),
                             first_rigid_limit=path["first_limit"], checked_rigid_return=path["checked_return"])
    if len(axial) != 2:
        raise ValueError("Axial marker must identify opposing rotations")
    for action, head, filename, torque_sign in (
        ("K4_internal", "FJ1504", "whole_arm_floating_terminal_audit.json", -1.),
        ("K5_external", "FJ1500", "whole_arm_shared_infraspinatus_audit.json", 1.),
    ):
        entry = axial[action]
        ref = provenance(folder / filename)
        record = load_reference(ref)
        if record["motion"] != combined["motion"]:
            raise ValueError(f"Cuff poses for {action} belong to another shared arm")
        rows = record["cases"][head][entry["source_path"]]
        accepted, first = [], None
        for row in rows:
            high = next(state for state in row["states"] if state["command_N"] == 10.)
            if torque_sign * high["humeral_torque_axis_Nm"][2] <= 0. or any(
                state["overlaps"] or state["central_wire_hits_mm3"] for state in row["states"]
            ):
                first = row
                break
            accepted.append(dict(input_deg=row["input_deg"], torque_Nm=high["humeral_torque_axis_Nm"][2],
                                 minimum_gap_mm=min(state["min_gap_capped_1mm"] for state in row["states"])))
        entry["individual_cuff_test"] = dict(source=ref, id=head, accepted=accepted, first_limit=first)

    def probe(angle, delta, include_wires=True, gap_lines=girdle_lines):
        mujoco.mj_resetData(model, data)
        mujoco.mj_integratePos(model, data.qpos, delta * np.deg2rad(angle), 1.)
        mujoco.mj_forward(model, data)
        frames = transforms(model, data, origins)
        placed = [bones[key].transform(frames[owners[key]]) for key in bone_keys]
        records = {}
        if include_wires:
            for key, route in routes.items():
                a, b = (frames[owner] for owner in route["owners"])
                neutral = np.asarray(route["points_source_mm"])
                w = np.asarray(route["distal_weights"])[:, None]
                xyz = (1. - w) * (neutral @ a[:, :3].T + a[:, 3]) + w * (neutral @ b[:, :3].T + b[:, 3])
                placed.append(rod(*xyz, .30) if len(xyz) == 2 else A.D.outer_envelope(xyz, .30))
                length, leverage = pull_jacobian(model, data, route, origins)
                records[key] = dict(length_mm=1000. * length, signed_leverage_mm=float(1000. * (leverage @ delta)))
        hits = find_hits(placed, keys if include_wires else bone_keys,
                         pairs if include_wires else pairs[pairs[:, 1] < count])
        ac_point = data.body("scapula").xpos * 1000.
        gh_point = data.body("humerus").xpos * 1000.
        dx = gh_point - origins["humerus"]
        line = gh_point - ac_point
        old = origins["humerus"] - origins["scapula"]
        rotation = np.arctan2(line[2], -line[0]) - np.arctan2(old[2], -old[0])
        gap = None
        if include_wires and not hits:
            gap_pairs = girdle_pair_indices if gap_lines == girdle_lines else [
                (i, j) for i, j in pairs if keys[i] in gap_lines or keys[j] in gap_lines]
            gap = min(float(A.I.gap(placed[i], placed[j], 3.)) for i, j in gap_pairs)
        return dict(input_deg=angle, qpos=data.qpos.tolist(), anterior_displacement_mm=float(dx @ anterior),
                    AC_anterior_displacement_mm=float((ac_point - origins["scapula"]) @ anterior),
                    superior_displacement_mm=float(dx[2]), coronal_scapular_deg=float(np.rad2deg(rotation)),
                    routes=records, girdle_wire_gap_capped_3mm=gap, hits=hits)

    def sampled_path(name, delta, limit, step, gap_lines):
        accepted, first = [], None

        def forward(angle):
            return probe(angle, delta, gap_lines=gap_lines)

        for index in range(round(limit / step) + 1):
            row = forward(index * step)
            if row["hits"]:
                if index:
                    for tick in range(1, round(step * 10)):
                        fine = forward(round((index - 1) * step + tick * .1, 1))
                        if fine["hits"]:
                            row = fine
                            break
                        accepted.append(fine)
                first = row
                break
            accepted.append(row)
        for row in reversed(accepted):
            repeat = forward(row["input_deg"])
            if repeat["hits"] or repeat["routes"] != row["routes"]:
                raise ValueError(f"Could not reproduce the return samples for {name}")
        print(name, len(accepted), "last", accepted[-1]["input_deg"] if accepted else None,
              "first hit", first["hits"][0] if first else None, flush=True)
        return accepted, first

    def route_summary(rows, names):
        summaries = {}
        if rows:
            for key in names:
                lengths = np.asarray([row["routes"][key]["length_mm"] for row in rows])
                leverage = [row["routes"][key]["signed_leverage_mm"] for row in rows]
                summaries[key] = dict(excursion_mm=float(np.ptp(lengths)), net_mm=float(lengths[-1] - lengths[0]),
                                      signed_leverage_range_mm=[min(leverage), max(leverage)])
        return summaries

    if args.shoulder_elevation:
        shoulder = {}
        gh = model.joint("GH").dofadr.item()
        for name, axis, active in (("K2_GH_flexion", "flexion", "FJ1468"),
                                   ("K3_GH_abduction", "abduction", "FJ1467")):
            delta = np.zeros(model.nv)
            delta[gh:gh + 3] = axes[axis]
            rigid = motion["paths"][name]["accepted"]
            upper = rigid[-1]["input_deg"]
            accepted, first = sampled_path(name, delta, upper, 1., {active})
            for row in accepted:
                if row["input_deg"].is_integer():
                    np.testing.assert_allclose(row["qpos"], rigid[int(row["input_deg"])]["qpos"], atol=1e-10, rtol=0.)
                row["active_wire_gap_capped_3mm"] = row.pop("girdle_wire_gap_capped_3mm")
            direction_limit = next((i for i, row in enumerate(accepted)
                                    if row["routes"][active]["signed_leverage_mm"] <= 0.), len(accepted))
            aligned = accepted[:direction_limit]
            shoulder[name] = dict(active=active, upper_deg=upper, accepted=accepted, first_limit=first,
                                  checked_return=len(accepted), first_nonpositive_leverage_deg=(
                                      accepted[direction_limit]["input_deg"] if direction_limit < len(accepted) else None),
                                  active_pull_last_deg=aligned[-1]["input_deg"] if aligned else None,
                                  active_pull=route_summary(aligned, [active]).get(active),
                                  line_summary=route_summary(accepted, routes),
                                  minimum_active_wire_gap_mm=min((row["active_wire_gap_capped_3mm"] for row in accepted), default=None))
        save(folder / "whole_arm_shoulder_elevation.json", dict(
            combined=combined_ref, base=base_ref,
            generator=provenance(ROOT / "tools/audit_whole_arm_anatomical_actions.py"), paths=shoulder,
            scope="Go-and-return rigid and eight fixed-end wires on the current support and tabletop base, at 1 degree "
                  "with 0.1 degree first-contact refinement. Prescribed GH poses with SC, AC, elbow and slider held neutral. "
                  "Signed leverage is the exact derivative of the prescribed blended wire path. Gap is taken for "
                  "the active wire, capped at 3 mm. Floating cuff poses, TPU deformation and muscular dynamics are excluded."))
        return

    sc, ac = (model.joint(name).dofadr.item() for name in ("SC", "AC"))
    def direction(axis, sign, ratio):
        delta = np.zeros(model.nv)
        delta[sc + axis], delta[ac + axis] = sign, sign * ratio
        return delta

    def screen(name, axis, sign, ratio):
        delta = direction(axis, sign, ratio)
        accepted, first = sampled_path(name, delta, 10., .5, girdle_lines)
        if len(accepted) > 1 and name in ("protraction", "retraction"):
            if sign * accepted[-1]["anterior_displacement_mm"] <= 0.:
                raise ValueError(f"Motion does not match its anatomical direction: {name}")
        return dict(SC_axis=axis, SC_direction=sign, AC_per_SC_same_axis=ratio, accepted=accepted,
                    checked_return=len(accepted), first_limit=first, line_summary=route_summary(accepted, sorted(girdle_lines)),
                    minimum_girdle_gap_mm=min((row["girdle_wire_gap_capped_3mm"] for row in accepted), default=None))

    scapular = {name: screen(name, axis, sign, ratio)
                for name, axis, sign, ratio in (("protraction", 2, 1., .5), ("retraction", 2, -1., .5),
                                                ("upward_rotation", 1, 1., 0.))}
    retraction_trials = []
    if len(scapular["retraction"]["accepted"]) < 2:
        for ratio in (0., 1., -.5, 1.5, -1., -1.1, -1.2, -1.3, -1.4, 2., -2.):
            delta = direction(2, -1., ratio)
            trial = probe(.5, delta, include_wires=False)
            retraction_trials.append(dict(AC_per_SC_same_axis=ratio, input_deg=.5, hits=trial["hits"],
                                          anterior_displacement_mm=trial["anterior_displacement_mm"]))
            if trial["hits"] or trial["anterior_displacement_mm"] >= 0.:
                continue
            alternative = screen("retraction_alternative", 2, -1., ratio)
            if len(alternative["accepted"]) > 1:
                keeps_pull = all(row["routes"]["FJ1536_g2"]["signed_leverage_mm"] > 0.
                                 for row in alternative["accepted"][:2])
                retraction_trials[-1]["initial_rhomboid_pull_supported"] = keeps_pull
                retraction_trials[-1]["full_screen_last_clear_deg"] = alternative["accepted"][-1]["input_deg"]
                if "retraction_alternative" not in scapular or keeps_pull:
                    scapular["retraction_alternative"] = alternative
                if keeps_pull:
                    break
    save(folder / "whole_arm_anatomical_actions.json", dict(
        combined=combined_ref, base=base_ref, generator=provenance(ROOT / "tools/audit_whole_arm_anatomical_actions.py"),
        anterior_axis=anterior.tolist(), medial_axis=medial.tolist(), neutral_bent_elbow_marker=hand.tolist(),
        thorax_axis_landmarks_source_mm=dict(manubrium=sternum.tolist(), T4=spine.tolist()), axial=axial,
        scapular=scapular, retraction_AC_ratio_trials=retraction_trials,
        scope="SC/AC probes to 10 degrees, sampled at 0.5 degrees and at 0.1 degrees from the last "
                                "clear sample to the first blocked sample. Go/return exact rigid and eight fixed-end wire "
                                "envelopes with the shared support and tabletop base. No moving cuff or TPU is tested in "
                                "the new scapular probes. Anterior displacement is measured at the GH centre, with GH held "
                                "neutral relative to scapula. Axial labelling uses a synthetic bent elbow marker; rigid motion "
                                "limits are from the original straight-elbow probes, not the bent probe. Cuff tests "
                                "evaluate each channel independently, including the local load and unload sequence."))


if __name__ == "__main__":
    main()
