"""Compare local attachment geometry and passive compliance for the current elbow."""

from pathlib import Path

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from simulate_elbow_antagonistic_load import HEADS, build_local_model


def rotate(radial, axis, angle):
    return np.cos(angle)[..., None] * radial + np.sin(angle)[..., None] * np.cross(axis, radial)


def profile(origin, centre, axis, distal, angles):
    offset = distal - centre
    axial = axis * (axis @ offset)
    posed = centre + axial + rotate(offset - axial, axis, angles)
    delta = origin - posed
    length = np.linalg.norm(delta, axis=-1)
    arm = np.cross(posed - centre, delta / length[:, None]) @ axis
    return posed, length - length[0], arm


def candidates(head, origin, centre, axis, distal, angles):
    offset = distal - centre
    axial = axis * (axis @ offset)
    radial = offset - axial
    radius = np.linalg.norm(radial)
    shifts = (0, -5, -10, -15, -20, -25, -30, -35, -40) if head == HEADS[1] else (-5, 0, 5)
    options = {}
    for shift in shifts:
        for phase in range(-30, 31, 5):
            point = centre + axial + rotate(radial * (1. + shift / radius), axis, np.deg2rad(phase))
            poses, stroke, arm = profile(origin, centre, axis, point, angles)
            options[shift, phase] = dict(
                shift_mm=shift, phase_deg=phase, point_source_mm=point.tolist(),
                displacement_from_original_mm=float(np.linalg.norm(point - distal)),
                poses=poses, stroke=stroke, arm=arm,
            )
    return options


def summarize(option, active_sign):
    poses = (0, 30, 60, 90, 116, 145)
    stroke, arm = option["stroke"], option["arm"]
    return dict(
        **{key: value for key, value in option.items() if key not in ("poses", "stroke", "arm")},
        stroke_mm={str(angle): float(stroke[angle]) for angle in poses},
        signed_moment_arm_mm={str(angle): float(arm[angle]) for angle in poses},
        minimum_effective_arm_mm_to_116=float((active_sign * arm[:117]).min()),
        minimum_effective_arm_mm_to_145=float((active_sign * arm[:146]).min()),
        first_wrong_direction_deg=next((i for i, value in enumerate(active_sign * arm) if value <= 0.), None),
    )


def screen_geometry(option, origin, obstacles_by_pose):
    result = {}
    for angle in range(146):
        # A 0.5 mm solid envelope flags a candidate bore centre inside hardware.
        point = option["poses"][angle]
        probe = A.I._manifold().Manifold.sphere(.5, 24).translate(point)
        chord = A.D.outer_envelope(np.stack([origin, point]), .3)
        for key, target in obstacles_by_pose[angle].items():
            for label, envelope in (("point", probe), ("straight_span", chord)):
                name = f"{label}:{key}"
                if name not in result:
                    volume = float(A.I.intersection_volume(envelope, target))
                    if volume > .001:
                        result[name] = dict(elbow_deg=angle, overlap_mm3=volume)
    return result


def main():
    check_file = DEST / "elbow_passive_model/model_check.json"
    check = read(check_file)
    source = ROOT / check["model"]["file"]
    assembly_file = DEST / "elbow_extension_stop/assembly.json"
    assembly = read(assembly_file)
    if provenance(source) != check["model"]:
        raise RuntimeError("Passive elbow source changed")
    model = mujoco.MjModel.from_xml_string(build_local_model(source, check))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    axis = data.joint("elbow").xaxis.copy()
    centre = data.joint("elbow").xanchor.copy() * 1000.
    elbow = model.joint("elbow").dofadr.item()
    angles = np.deg2rad(np.arange(146.))
    gravity = []
    for angle in angles:
        data.joint("elbow").qpos[:] = angle
        mujoco.mj_forward(model, data)
        gravity.append(-float(data.qfrc_bias[elbow]))
    mujoco.mj_resetData(model, data)
    mujoco.mj_forward(model, data)
    owners = {row["source_key"]: row["owner"] for row in check["geom_owners"]}
    obstacles = {key: (owners[key], load_component(value)) for key, value in assembly["sources"].items()
                 if key in ("bone:humerus", "bone:forearm") or key.startswith("elbow:elbow_")}
    obstacles_by_pose = [
        {key: part.transform(A.I.rotation_about(centre, axis, angle)) if owner == "forearm" else part
         for key, (owner, part) in obstacles.items()}
        for angle in range(146)
    ]
    groups = {}
    tested = {}
    selected = {}
    for head, sign in zip(HEADS, (1., -1.), strict=True):
        origin, distal = [data.site(f"{head}_{end}").xpos.copy() * 1000. for end in ("origin", "distal")]
        options = candidates(head, origin, centre, axis, distal, angles)
        groups[head] = options
        if sign < 0.:
            representative = [(0, 0), (-10, 0), (-20, 0), (-30, 0), (-40, 0)]
            baseline_arm = min(sign * options[0, 0]["arm"])
            local_options = (key for key, option in options.items()
                             if option["displacement_from_original_mm"] <= 20.
                             and min(sign * option["arm"]) >= baseline_arm)
            representative.extend(sorted(local_options, key=lambda key: options[key]["stroke"][116])[:12])
        else:
            representative = [(0, 0), (0, -10), (0, 10), (5, 0)]
        # Search within a 15 mm endpoint move for a larger minimum arm over the full target.
        best = max((option for option in options.values()
                    if option["displacement_from_original_mm"] <= 15.),
                   key=lambda option: min(sign * option["arm"]))
        representative.append((best["shift_mm"], best["phase_deg"]))
        tested[head] = {}
        clear = []
        for key in dict.fromkeys(representative):
            option = options[key]
            row = summarize(option, sign)
            overlaps = screen_geometry(option, origin, obstacles_by_pose)
            row["straight_span_geometry_first_overlaps"] = overlaps
            tested[head][f"{key[0]},{key[1]}"] = row
            if not overlaps and option["displacement_from_original_mm"] <= 20. and min(sign * option["arm"]) > 0.:
                clear.append(key)
            print(head, key, "stroke at 116", round(row["stroke_mm"]["116"], 2),
                  "min arm to 145", round(row["minimum_effective_arm_mm_to_145"], 2),
                  "overlaps", list(row["straight_span_geometry_first_overlaps"]), flush=True)
        if sign < 0.:
            selected[head] = min(clear, key=lambda key: options[key]["stroke"][116])
        else:
            selected[head] = max(clear, key=lambda key: min(options[key]["arm"]))

    comparisons = []
    flexor_options, extensor_options = (groups[head] for head in HEADS)
    # Extension stiffness variants are hypotheses for the full extensor assembly.
    # Compression stiffness and flexor stiffness are kept at the current reference.
    for flex_key, ext_key in (((0, 0), (0, 0)), ((0, 0), (-20, 0)),
                              (selected[HEADS[0]], selected[HEADS[1]])):
        flex, ext = flexor_options[flex_key], extensor_options[ext_key]
        for extension_stiffness in (200., 100., 50., 0.):
            rows = {}
            for angle in (0, 30, 116, 145):
                f_arm, e_arm = flex["arm"][angle], ext["arm"][angle]
                f_stroke, e_stroke = flex["stroke"][angle], ext["stroke"][angle]
                torque = (gravity[angle] + 200. * f_stroke * f_arm * 1e-6
                          + extension_stiffness * e_stroke * e_arm * 1e-6)
                rows[str(angle)] = dict(
                    elbow_deg=angle, opposing_torque_Nm=float(-torque),
                    flexor_tension_to_balance_N=float(-torque / (f_arm * .001)) if f_arm > 0. else None)
            comparisons.append(dict(flexor_option=list(flex_key), extensor_option=list(ext_key),
                                    extensor_extension_stiffness_N_per_m=extension_stiffness,
                                    static_results=rows))
    output = DEST.parent / "elbow_local_muscle_load"
    save(output / "elbow_local_attachment_options.json", dict(
        generator=provenance(Path(__file__).resolve()), source=check["model"],
        rigid_geometry=provenance(assembly_file), step_deg=1., tested_options=tested,
        grid={head: [summarize(option, sign) for option in groups[head].values()]
              for head, sign in zip(HEADS, (1., -1.), strict=True)},
        local_probe_shortlist=selected, passive_comparisons=comparisons,
        scope="Hypothetical neutral-cap endpoint offsets, not installed anchors. Holds proximal origins and "
        "elbow axis fixed. Point and straight chord probes screen only both bones and central hinge hardware; "
        "they do not assess a curved muscle or installed connector. Quasi-static balance fixes the caps "
        "and neglects contact, preload and friction. Zero extensor stiffness is a lower-bound comparison."))


if __name__ == "__main__":
    main()
