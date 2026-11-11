"""Summarize contact evidence while retaining geometry-basis and hardware gaps."""
from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elastic_conduit_bone_contact'
    stems = ['FJ1468_T5_k100000_mu0_u5', 'FJ1468_T5_k200000_mu0_u5',
             'FJ1468_T5_k400000_mu0_u5', 'FJ1468_T5_k200000_mu0.15_u5',
             'FJ1468_T5_k200000_mu0_u2']
    rows = []
    for stem in stems:
        source_file = folder/f'{stem}.json'
        audit_file = folder/f'{stem}_trajectory_audit.json'
        source, audit = read(source_file), read(audit_file)
        for key in ('source', 'generator', 'surface_query', 'beam_generator', 'model', 'anatomy_reference'):
            assert provenance(ROOT/source[key]['file']) == source[key]
        assert audit['source'] == provenance(source_file)
        assert provenance(ROOT/audit['anatomy_reference']['file']) == audit['anatomy_reference']
        assert source['result']['finite_no_warning_or_axial_limit']
        assert audit['states_with_wire_bone_overlap'] == 0
        rows.append(dict(name=stem, source=provenance(source_file), trajectory_audit=provenance(audit_file),
            specification=source['specification'],
            summary={key: value for key, value in source['result'].items() if key != 'trace'},
            minimum_whole_wire_gap_mm_capped1_5=audit['minimum_wire_gap_mm_capped1_5'],
            checked_states=audit['checked_states'], checked_parts=audit['checked_parts']))
    verification_file = folder/'surface_query_verification.json'
    verification = read(verification_file)
    assert provenance(ROOT/verification['surface_module']['file']) == verification['surface_module']
    comparison_file = folder.parent/'reference_basis_comparison.json'
    comparison = read(comparison_file)
    for key in ('reference', 'current_selection', 'selected_assembly'):
        assert provenance(ROOT/comparison[key]['file']) == comparison[key]
    nominal, refined = rows[1]['summary'], rows[4]['summary']
    differences = dict(
        force_update_refinement_peak_displacement_difference_mm=abs(nominal['peak_node_displacement_mm']-refined['peak_node_displacement_mm']),
        force_update_refinement_maximum_penetration_difference_mm=abs(nominal['maximum_guide_envelope_penetration_mm']-refined['maximum_guide_envelope_penetration_mm']),
        force_update_refinement_whole_wire_gap_difference_mm=abs(rows[1]['minimum_whole_wire_gap_mm_capped1_5']-rows[4]['minimum_whole_wire_gap_mm_capped1_5']))
    save(folder/'summary.json', dict(cases=rows, contact_query_verification=provenance(verification_file),
        geometry_basis_comparison=provenance(comparison_file), refinement_differences=differences,
        status='conditional_fixed_end_contact_prediction_not_final_route_approval',
        decision='Bone contact resolves the examined free-space wire penetration in this 5 N guide bench. Retain route candidate; engineer end interfaces and reconcile the unique shoulder/elbow bone reference before moving-arm validation.',
        outstanding=['one reconciled shoulder/elbow engineering bone assembly',
                     'manufactured end mounts and cable terminal reaction seats',
                     'moving-bone contact and unified arm dynamics',
                     'independent tension-only slack spool and internal cable friction',
                     'cross-section/lumen distortion and minimum lattice support',
                     'all other routes and target motions']))
    lines = ['# 骨面接触支承：FJ1468通道', '',
        '结论：在本轮固定端部、5 N张力的等效模型中，实际骨面反力能够阻止此前自由空间模型的钢丝穿骨。五组参数对照、每组250个保存姿态的整段钢丝筛查均未发现大于0.001 mm³的骨骼重叠。结果支持保留该路线候选，不代表完整肩肘模型、工程接口或实物性能通过。', '',
        '## 模型及支承来源', '',
        '沿用中空TPU通道对应的可压缩弹性梁。近端位置及首段方向固定，远端位置固定、方向自由；中间节点全部可动。肱骨和肩胛骨的真实三角网格作为接触表面，压缩反力经运动节点的点雅可比传到梁上。没有将中间节点固定到世界，也没有使用工程背板。端部固定仍是边界条件，不是已完成安装件。', '',
        '点—三角形距离在带内直接计算。空间缓存只保存有界候选面，不插值距离场；解析立方体500点及两骨各300点完成验证，真实模型带内距离与全局搜索一致，符号与实体射线奇偶一致。该验证检查算法抽样，不证明所有表面输入均适用。', '',
        '模量暂取15 MPa，5 N张力在0.05 s内升起、保持至0.25 s、0.30 s前卸载，总计算0.5 s。所有材料、接触刚度、阻尼及摩擦均未实测。接触力密度以每米通道的压入量积分，单位N/m²；不能把表中的数值直接当作单点N/mm刚度。', '',
        '| 接触密度 N/m² | TPU—骨摩擦系数 | 力更新 μs | 最大节点偏移 mm | 最大名义压入 mm | 整段钢丝最小骨面间隙 mm |',
        '|---:|---:|---:|---:|---:|---:|']
    for row in rows:
        spec, result = row['specification'], row['summary']
        lines.append(f"| {spec['contact_penalty_density_N_per_m2']:g} | {spec['bone_TPU_friction_coefficient']:g} | {spec['contact_update_interval_s']*1e6:g} | {result['peak_node_displacement_mm']:.3f} | {result['maximum_guide_envelope_penetration_mm']:.3f} | {row['minimum_whole_wire_gap_mm_capped1_5']:.3f} |")
    lines += ['',
        '骨面摩擦系数0.15仅代表TPU外壁与骨面的探索对照；钢丝在TPU通道内的摩擦尚未加入，不能据此宣称已满足驱动链的内部摩擦预测要求。每组均无MuJoCo数值警告或轴向保护限位触发。', '',
        f"名义组力更新由100 μs细化至40 μs后，最大偏移变化{differences['force_update_refinement_peak_displacement_difference_mm']:.6f} mm，整段最小间隙变化{differences['force_update_refinement_whole_wire_gap_difference_mm']:.6f} mm。动力学积分步长20 μs，钢丝实体检查的保存时间间隔约2 ms；不宣称全时间数学无干涉证明。", '',
        '压入量是固定截面包络和接触罚模型中的几何量，不是实际晶格或孔腔变形。相对名义0.8 mm壁厚，0.06–0.20 mm仍需局部变形／样件检查；不能仅因接触判失败，也不能据此证明不夹线。名义组仍有约8.95 mm通道缩短发生在端点位置固定时，后续必须区分通道变形消耗的收线与真正的骨段驱动行程。', '',
        '## 几何基准与后续工作', '',
        '当前接触计算依据thorax_arm_reference的共享肱骨装配。与S1_CURRENT_CANDIDATE中的实际骨骼求差，锁骨、肩胛骨及肱骨均存在实体差异；详见reference_basis_comparison.json。不得把本表移作S1-Y2合并装配通过证据，原文件和候选指针均未替换。', '',
        '下一步先收敛唯一肩肘骨骼基准及端部接口，明确驱动端、反力座和短肌腱终端，再将骨面接触随运动骨段带入统一模型。该结果尚不支持添加任意空间固定点；需要时才设计实际导向或最小晶格支承。完整肌腹、内部钢丝摩擦、松线与卷线器模型及其余肌群仍需继续。', '',
        '## 复现', '',
        '项目Python环境下运行tools/simulate_channel_first_conduit_bone_contact.py，参数为--contact-density、--friction和--update-steps。运行tools/verify_channel_first_bone_contact.py验证查询；对每个输出stem运行tools/audit_channel_first_contact_trajectory.py --stem STEM；最后运行tools/summarize_channel_first_bone_contact.py。XML只包含梁及驱动器，单独打开不会重现Python施加的骨接触力。', '']
    (folder/'README.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
    print('Verified five contact cases and 1250 complete-wire/anatomy states; conditional result only')


if __name__ == '__main__':
    main()
