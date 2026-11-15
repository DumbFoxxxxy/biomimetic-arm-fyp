"""Generate the first-pass cable handoff from the checked shared-arm records."""

import csv
from pathlib import Path

import numpy as np

from build_channel_first_inventory import DEST, provenance
from build_s1_teres_origin_mount import save
from study_s1_local_force_geometry import load_reference


CHANNELS = {
    "FJ1468": ("三角肌前束", "锁骨入口、肱骨终端、TPU 通道随动"),
    "FJ1467": ("三角肌中束", "肩胛入口、肱骨终端、2 mm 偏置与新版 TPU 内衬"),
    "FJ1513": ("三角肌后束", "肩胛入口、肱骨终端、TPU 通道随动"),
    "FJ1507": ("大圆肌", "肩胛与肱骨端盖、局部随钢丝通道运动"),
    "FJ1486": ("肘屈肌", "肱骨入口与前臂转座"),
    "FJ1477": ("肘伸肌", "肱骨侧滑块入口与前臂转座"),
    "FJ1459_g3": ("前锯肌选定束", "胸廓、肩胛导向眼，固定方式待核对"),
    "FJ1536_g2": ("大菱形肌选定束", "胸廓、肩胛导向眼，固定方式待核对"),
    "FJ1500": ("冈下肌", "肩胛入口、浮动端盖及肱骨侧上移转眼和共用支座"),
    "FJ1504": ("肩胛下肌", "肩胛入口、浮动端盖及肱骨侧转眼"),
}


def source(folder, name):
    ref = provenance(folder / name)
    return ref, load_reference(ref)


def xyz(point):
    return ", ".join(f"{value:.3f}" for value in point)


def save_csv(path: Path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    combined_ref, combined = source(folder, "whole_arm_combined_cuff_k7.json")
    routed = load_reference(combined["routed"])
    motion = load_reference(combined["motion"])
    scene = load_reference(motion["scene_check"])
    seats_ref, seats = source(folder, "girdle_seat_pair.json")
    _, tpu = source(folder, "whole_arm_k7_tpu_moving_cuff_audit.json")
    if tpu["combined"] != combined_ref:
        raise ValueError("TPU and wire data use different combined runs")
    branch = routed["paths"]["K7_SC1_pos"]
    verified = []
    for row in combined["combined"]:
        if (row["baseline_wire_hits_mm3"] or row["baseline_wire_pair_hits"]
                or any(step.get("missing_cuffs") or step["local_cuff_hits"]
                       or step["cuff_cross_hits_mm3"] or step["cuff_to_neighbor_wires_mm3"]
                       for step in row["steps"])):
            break
        verified.append(row)
    if len(verified) < 2 or any(a["input_deg"] != b["input_deg"]
                                for a, b in zip(verified, branch["samples"][:len(verified)], strict=True)):
        raise ValueError("Missing common K7 samples")
    labels = []
    csv_rows = []
    guide_points = {}
    for key, spec in CHANNELS.items():
        name, support = spec
        if key in routed["routes"]:
            item = routed["routes"][key]
            owners = item["owners"]
            curve = np.asarray(item["points_source_mm"])
            points = curve[[0, -1]]
            frames = [np.asarray(scene["neutral_body_origins_mm"][owner]) for owner in owners]
            weights = np.asarray(item["distal_weights"])
            if len(weights) != len(curve) or np.any((weights < -1e-10) | (weights > 1. + 1e-10)):
                raise ValueError(f"Invalid motion weights on {key}")
            guide_points[key] = dict(muscle=name, owners=owners, definition=item["type"], source=item["source"],
                                     centreline_source_mm=curve.tolist(),
                                     centreline_local_mm_by_owner=[(curve - origin).tolist() for origin in frames],
                                     distal_owner_blend_weights=weights.tolist())
            if key in ("FJ1459_g3", "FJ1536_g2"):
                defined_seats = [seats["seats"][f"{key}_{role}"] for role in ("origin", "insertion")]
                np.testing.assert_allclose(points, [seat["port_source_mm"] for seat in defined_seats], atol=1e-8, rtol=0.)
                guide_points[key]["physical_guide_eyes"] = defined_seats
            lengths = [row["routes"][key]["length_mm"] for row in branch["samples"][:len(verified)]]
            change = lengths[-1] - lengths[0]
            excursion = np.ptp(lengths)
            signed = -np.diff(lengths) / np.deg2rad(np.diff([row["input_deg"] for row in verified]))
            direction = "收线" if np.min(signed) > 0 else "放线" if np.max(signed) < 0 else "方向变号"
            metric = f"{change:+.2f} / {excursion:.2f} mm，{direction}，等效力臂 {min(signed):+.2f} 至 {max(signed):+.2f} mm"
            local = points - [scene["neutral_body_origins_mm"][owner] for owner in owners]
            local_text = [xyz(p) for p in local]
            point_type = "名义路径两端"
            limit = branch["summary"][key]["first_limit"]
            if limit:
                obstacle = "、".join(limit["rigid_hits_mm3"]) or "相邻钢丝"
                limit_text = f"SC {limit['input_deg']:.2f}°，{obstacle}"
            else:
                limit_text = "所查 K7 路径内未见超阈值相交"
        else:
            states = [next(state for state in case["states"] if state["command_N"] == 10.)
                      for case in combined["cases"][key]["K7_SC1_pos"][:len(verified)]]
            points = next(step for step in verified[0]["steps"] if step["load_N"] == 10.)[
                "central_tendon_tips_world_mm"][key]
            neutral_state = states[0]
            origin = np.asarray(neutral_state["origin_port_source_mm"])
            initial_cap_pin = np.asarray(neutral_state["cap_pin_neutral_source_mm"])
            cap_pose = np.asarray(neutral_state["cap_affine_source_mm"])
            posed_pin = cap_pose[:, :3] @ initial_cap_pin + cap_pose[:, 3]
            guide_points[key] = dict(muscle=name, origin_owner="scapula", distal_owner="humerus",
                                     origin_port_source_mm=origin.tolist(), origin_port_local_mm=(
                                         origin - scene["neutral_body_origins_mm"]["scapula"]).tolist(),
                                     unloaded_floating_cap_pin_source_mm=initial_cap_pin.tolist(),
                                     loaded_floating_cap_pin_source_mm=posed_pin.tolist(), load_probe_N=10.,
                                     loaded_neutral_tendon_tips_world_mm=points,
                                     distal_exit_local_mm=(
                                         np.asarray(points[1]) - scene["neutral_body_origins_mm"]["humerus"]).tolist(),
                                     retained_lead_lengths_mm=neutral_state["retained_lead_lengths_mm"])
            owners = ("浮动肌腹端", "肱骨侧转眼")
            local_text = ["", ""]
            moment = [state["humeral_torque_axis_Nm"][-1] for state in states]
            metric = f"10 N 时长轴力矩 {min(moment):+.3f} 至 {max(moment):+.3f} N·m"
            change = excursion = None
            point_type = "10 N 中立位受拉中央跨段的两出线端"
            limit_text = "中央跨段通过；浮动端盖与自身 TPU 的配合待核对"
        labels.append(f"| {key} {name} | {owners[0]} → {owners[1]} | {support} | {metric} | {limit_text} |")
        csv_rows.append(dict(id=key, muscle=name, origin_owner=owners[0], distal_owner=owners[1],
                             point_definition=point_type, origin_source_xyz_mm=xyz(points[0]),
                             distal_source_xyz_mm=xyz(points[1]), origin_body_local_xyz_mm=local_text[0],
                             distal_body_local_xyz_mm=local_text[1], K7_signed_delta_mm=change,
                             K7_excursion_mm=excursion, support=support, K7_limit=limit_text))
    save_csv(folder / "whole_arm_cable_handoff.csv", csv_rows)
    save(folder / "whole_arm_guide_points.json", dict(
        combined=combined_ref, seats=seats_ref, routes=combined["routed"], data=guide_points,
        scope="Intermediate samples describe the prescribed centreline and the blend between two bone frames. "
              "Only separately identified eyes, mounts and ports are physical guide interfaces. Floating cuff "
              "tips and pin positions are taken from the neutral posture under the individual 10 N input probe."))
    notes = [
        "# 整臂首版 · 钢丝与接口交接",
        "",
        f"K7 共 {len(verified)} 个检查样点，最后共同通过的十线样点为 SC {verified[-1]['input_deg']:.2f}°，"
        f"对应肱骨方向变化 {verified[-1]['humerus_deg']:.2f}°。下表的八条固定端路径只统计这个共同区间。",
        "",
        "线长列为终点减中立位的名义线长变化及区间峰峰值。沿协调 SC 输入的等效力臂为负线长梯度，"
        "正号表示该路径正向运动时需要收线。肩袖两行给出独立求解下的长轴力矩，"
        "中央短肌腱的固定长度约束与执行器实际收线量须分别核对。",
        "",
        "| 编号与功能肌群 | 有效区段两端 | 导向或支承 | K7 所查区间的局部响应 | 首限或待核对项 |",
        "|---|---|---|---|---|",
        *labels,
        "",
        "坐标表给出八条路径两端在中立位的解剖源坐标及骨段局部坐标。肩袖两行记录的是 10 N 条件下的中央跨段端点，"
        "浮动端随力变动，故这两行留空骨段局部坐标。各通道的完整折线和插值权重见路径报告。",
        "",
        "[坐标 CSV](whole_arm_cable_handoff.csv) · [完整折线、导向眼与肩袖端口坐标](whole_arm_guide_points.json)"
        " · [路径与首限](whole_arm_posterior_detour_audit.json)"
        " · [肩袖局部计算](whole_arm_combined_cuff_k7.json) · [TPU 和五金接触](whole_arm_k7_tpu_moving_cuff_audit.json)",
    ]
    (folder / "whole_arm_cable_handoff.md").write_text("\n".join(notes) + "\n", encoding="utf-8")
    print(f"Exported {len(csv_rows)} channel records for {len(verified)} checked samples")


if __name__ == "__main__":
    main()
