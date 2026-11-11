"""Consolidate the completed sensitivity cases without treating them as acceptance."""
import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elastic_conduit_bench'
    names = ['FJ1468_E5_T5_ds5_dt2e-05', 'FJ1468_E15_T5_ds5_dt2e-05',
             'FJ1468_E40_T5_ds5_dt2e-05', 'FJ1468_E15_T0_ds5_dt2e-05',
             'FJ1468_E15_T5_ds5_dt1e-05', 'FJ1468_E15_T5_ds2.5_dt1e-05']
    cases, rows = {}, []
    for name in names:
        path = folder/f'{name}.json'
        report = read(path)
        for field in ('source', 'generator', 'model'):
            assert provenance(ROOT/report[field]['file']) == report[field]
        result = report['result']
        assert result['finite_and_no_limit_or_warning']
        assert max(r['constraint_force_reconstruction_max_abs_mixed_units'] for r in result['trace']) < 1e-8
        hold = min(result['trace'], key=lambda r: abs(r['time_s']-.24))
        cases[name] = report
        rows.append(dict(name=name, source=provenance(path),
            specification=report['specification'],
            summary={k: v for k, v in result.items() if k != 'trace'},
            hold={k: v for k, v in hold.items() if k != 'nodes_local_m'}))
    nominal = cases[names[1]]['result']
    time_refined = cases[names[4]]['result']
    space_refined = cases[names[5]]['result']
    base_hold = min(nominal['trace'], key=lambda r: abs(r['time_s']-.24))
    gravity_hold = min(cases[names[3]]['result']['trace'], key=lambda r: abs(r['time_s']-.24))
    additional_hold_displacement = float(np.linalg.norm(
        np.array(base_hold['nodes_local_m'])-np.array(gravity_hold['nodes_local_m']), axis=1).max())*1000.
    differences = dict(
        time_refinement_peak_displacement_difference_mm=abs(time_refined['maximum_node_displacement_mm']-nominal['maximum_node_displacement_mm']),
        spatial_refinement_peak_displacement_difference_mm=abs(space_refined['maximum_node_displacement_mm']-time_refined['maximum_node_displacement_mm']),
        spatial_refinement_minimum_length_change_difference_mm=abs(space_refined['minimum_loaded_length_change_mm']-time_refined['minimum_loaded_length_change_mm']),
        loaded_vs_zero_tension_same_time_max_node_difference_mm=additional_hold_displacement,
        qualification='Two spatial resolutions and two timesteps only. Peak scalar agreement is not a convergence proof for every force, transient or shape.')
    contact_file = folder/'anatomy_snapshot_screen.json'
    contact = read(contact_file)
    assert provenance(ROOT/contact['loaded_report']['file']) == contact['loaded_report']
    assert provenance(ROOT/contact['gravity_report']['file']) == contact['gravity_report']
    save(folder/'sensitivity_summary.json', dict(cases=rows, differences=differences,
        anatomy_snapshot_screen=provenance(contact_file),
        decision='FJ1468 requires guide/support revision before moving-route acceptance: computed free-space wire crosses humerus. No current route, muscle or whole-arm dynamic approval.',
        next_required_work=['bone-attached end interfaces with revised exit directions',
                            'physically supported bone-side guide/contact or minimum lattice support',
                            'recompute beam reaction and shape with that support',
                            'add tension-only slack spool and friction',
                            'integrate supported route into moving thorax/shoulder/elbow assembly']))
    lines = [
        '# 弹性通道支承试验：FJ1468', '',
        '这是三角肌前束通道的隔离边界试验，不是整臂驱动验收。它用于检查：静态弯曲通道只固定两端时，拉线会不会大量消耗在通道自身变形上。', '',
        '模型以当前中空通道为起点，近端位置及首段方向固定，远端位置通过连接约束固定、方向自由。内部节点全部可动；轴向和转动弹簧分别依据 EA/l、EI/l 建立。截面外半径1.6 mm、内半径0.8 mm；假设密度1200 kg/m³，模量5/15/40 MPa均为未实测的探索参数，不是论文给出的性能保证。阻尼随段长反比缩放。', '',
        '拉力在0–0.05 s平滑升至指定值，保持至0.25 s，0.25–0.30 s卸载，观察至0.5 s。输入是非负钢丝张力，关节和内部节点位置均由动力学求解；本试验尚未计算卷线器、独立钢丝松弛、摩擦、骨面接触或完整肌腹。', '',
        '| 假设模量 MPa | 拉力 N | 段长 mm | 步长 μs | 最大节点偏移 mm | 最多缩短 mm | 最大端部约束误差 mm |',
        '|---:|---:|---:|---:|---:|---:|---:|']
    for row in rows:
        spec, result = row['specification'], row['summary']
        mean_length = np.mean(spec['segment_rest_lengths_m'])*1000.
        lines.append(f"| {spec['assumed_modulus_MPa']:g} | {result['peak_tension_N']:g} | {mean_length:.3f} | {spec['timestep_s']*1e6:g} | {result['maximum_node_displacement_mm']:.3f} | {-result['minimum_loaded_length_change_mm']:.3f} | {result['maximum_distal_support_error_mm']:.5f} |")
    lines += ['',
        f"15 MPa、5 N主例在0.24 s时相对零张力对照的最大节点位置差为{additional_hold_displacement:.3f} mm。零张力下也有重力下垂，不能把总位移全部说成钢丝造成。0.5 s的残余位移是有限时长下的动态状态，不能解释为永久变形。", '',
        f"步长减半后最大偏移差{differences['time_refinement_peak_displacement_difference_mm']:.6f} mm；段长再减半后差{differences['spatial_refinement_peak_displacement_difference_mm']:.3f} mm，最多缩短量差{differences['spatial_refinement_minimum_length_change_difference_mm']:.3f} mm。这支持偏移不是粗离散产生的假象，但不是所有物理量的全局收敛证明。六组均无MuJoCo数值警告或轴向保护限位触发。±50%轴向限位仅为数值保护，不是已设计的机械限位。", '',
        '## 回放到骨骼中的筛查', '',
        '将主例最大偏移、0.24 s加载、最后回弹姿态以及零张力0.24 s姿态放回64个骨件中检查。加载和最大偏移姿态的名义钢丝会进入肱骨，通道包络也与肩胛和肱骨相交；另两处检查姿态没有发现大于0.001 mm³的重叠。该检查仅覆盖四个姿态，不是连续碰撞验收。', '',
        '自由空间梁没有骨接触，因此上述穿骨不是实际可实现的运动结果，而是说明模型必须补入骨面导向、接触或有实体支承的最小晶格结构。TPU包络相交只作为接触待解标记，不能仅按体积判软结构失败；钢丝路径穿骨则不能直接批准。', '',
        '## 当前决定与复现', '',
        '保留四条中空通道为静态几何候选，FJ1468暂不标记为运动路径通过。下一步先补端部实际连接及骨面导向／最小晶格支承，重新求解通道反力与形状，再接入全臂模型及只受拉、允许松弛的卷线驱动。当前没有增加悬空固定点，也没有采用工程背板。', '',
        '生成单例：项目Python环境下运行 `tools/simulate_channel_first_elastic_conduits.py`，可用 `--modulus`、`--tension`、`--segment-length`、`--timestep` 指定表中对照。之后运行 `tools/audit_channel_first_elastic_conduit_snapshots.py` 及 `tools/summarize_channel_first_elastic_conduits.py`。原始JSON保存节点轨迹和反力，XML保存对应边界模型。', '',
        '早期双层梁模型的细化在加载阶段异常结束，未产生合格结果；已保留在historical_two_body。改为每节点一个刚体后，67段模型正常加载并通过轴向运动学检查。具体原生退出原因未定，不把猜测当根因。恒定关节阻尼的中间模型另存historical_constant_joint_damping，当前表只使用阻尼随段长缩放后重新计算的结果。', '']
    (folder/'README.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
    print('Summarized', len(rows), 'verified cases; route acceptance remains withheld')


if __name__ == '__main__':
    main()
