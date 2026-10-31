"""M5 evidence consolidation with explicit local/whole-arm scope boundaries."""
from build_distal_swivels import *
def run():
    read=lambda name:json.loads((DEST/name).read_text(encoding='utf8'))
    assembly=read('assembly.json');neutral=read('neutral_checks.json');exports=read('export_checks.json');whole=read('whole_arm_STEP_check.json');routes=read('straight_route_screen.json');prints=read('print_transforms.json');animation=read('terminal_articulation.json')
    assert len(assembly['instances'])==neutral['instance_count']==131
    assert not neutral['rigid_interferences'] and not neutral['TPU_contact_review_flags']
    assert not routes['new_terminal_vs_fixed_bony_context']
    assert routes['assembly_sha256']==hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest()
    assert exports['pass'] and len(exports['parts'])==22 and len(exports['scene'])==131 and all(r['pass'] for r in exports['prescribed_terminal_pose_checks'])
    assert exports['neutral_scene_sha256']==hashlib.sha256((DEST/'neutral_view.xml').read_bytes()).hexdigest()
    assert whole['whole']['pass'] and whole['whole']['solid_count']==129 and whole['source_assembly_sha256']==hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest()
    for check in [whole['whole'],*whole['modules'].values()]:assert check['STEP_sha256']==hashlib.sha256((ROOT/check['STEP']).read_bytes()).hexdigest()
    for p in whole['parts'].values():
        assert p['STEP_vs_OBJ_vertex_error_mm']<.001
        assert p['OBJ_sha256']==hashlib.sha256((ROOT/p['source_OBJ']).read_bytes()).hexdigest()
        assert p['STEP_sha256']==hashlib.sha256((ROOT/p['source_STEP']).read_bytes()).hexdigest()
    summaries=[];evidence=['assembly.json','neutral_checks.json','export_checks.json','whole_arm_STEP_check.json','straight_route_screen.json','print_transforms.json','terminal_articulation.json']
    for fid in ['FJ1486','FJ1477']:
        design=read(fid+'/design.json');checks=read(fid+'/mechanism_checks.json');cad=read(fid+'/cad_roundtrip.json');inputs={p['name']:p for p in read(fid+'/cad_input.json')['parts']}
        assert set(design['parts'].values())=={1} and len(design['parts'])==11
        assert not design['selected']['hard_interferences']
        assert checks['minimum_sampled_articulation_prefix_deg']==25 and len(checks['articulation_paths'])==24
        assert checks['all_tested_pullout_directions_obstructed_by_3mm'] and len(checks['retention_directions'])==28
        assert checks['all_rigid_assembly_paths_clear'] and len(checks['assembly_paths'])==10 and not checks['assembly_soft_flags']
        if DISTAL_LOAD:assert design['integral_clamp_stops']['present'] and checks['all_integral_clamp_stops_pass']
        if DISTAL_LOAD:
            canonical=read(fid+'/canonical_stop_mesh.json');assert canonical['OBJ_sha256']==hashlib.sha256((DEST/fid/'distal_socket_base.obj').read_bytes()).hexdigest()
            assert canonical['positive_material_components']==1
            evidence+=[fid+'/canonical_stop_mesh.json']
        for path in checks['articulation_paths']:assert len(path['forward'])==len(path['return'])==26
        for name,digest in checks['source_mesh_sha256'].items():assert hashlib.sha256((DEST/fid/(name+'.obj')).read_bytes()).hexdigest()==digest
        assert len(cad['parts'])==11
        for p in cad['parts']:assert p['pass'] and abs(p['input_volume_mm3']-inputs[p['part']]['volume_mm3'])<1e-5
        obstructions=[r['first_obstruction']['offset_mm'] for r in checks['retention_directions']]
        summary={'muscle':fid,'assembly_paths':10,'sampled_articulation_prefix_deg':25,'radial_articulation_paths':24,'pullout_directions_at_neutral_terminal_orientation':28,'first_pullout_obstruction_range_mm':[min(obstructions),max(obstructions)],'ball_centre_display_mm':design['selected']['ball_centre_display_mm'],'retainer_roll_deg':design['selected']['retainer_roll_deg'],'housing_position_local_mm':design['selected']['position_local_mm'],'rear_nut_installation':[p['label'] for p in checks['assembly_paths'] if p['label'].startswith('M2_nut')]};summaries.append(summary)
        evidence += [fid+'/'+name for name in ['design.json','mechanism_checks.json','cad_roundtrip.json','wire_stub_check.json']]
    assert len(prints['parts'])==4
    for p in prints['parts']:
        folder=DEST/p['muscle'];assert p['source_OBJ_sha256']==hashlib.sha256((folder/(p['part']+'.obj')).read_bytes()).hexdigest();assert p['print_STL_sha256']==hashlib.sha256((folder/'print_candidates'/(p['part']+'.stl')).read_bytes()).hexdigest()
    assert animation['sha256']==hashlib.sha256((ROOT/animation['file']).read_bytes()).hexdigest() and animation['actual_frame_count']>100
    branch=None
    if DISTAL_LOAD:
        branch=read('first_branch_checks.json');assert branch['assembly_sha256']==hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest()
        assert all(p['opposition_to_neutral_active_pull']>.999 and not p['rigid_interferences'] and not p['static_TPU_contact_flags'] for p in branch['modules'])
        evidence+=['first_branch_checks.json']
        animation_inputs=read('terminal_animation_inputs.json');assert animation_inputs['GIF_sha256']==animation['sha256']
        for filename,digest in animation_inputs['input_sha256'].items():assert hashlib.sha256(__import__('pathlib').Path(filename).read_bytes()).hexdigest()==digest
        evidence+=['terminal_animation_inputs.json']
    status={'stage':'M5B distal-load spherical endpoints' if DISTAL_LOAD else 'M5A direct-reference spherical endpoints','neutral_manufactured_instances':131,'rigid_CAD_solids':129,'new_endpoint_parts':22,'local_assembly_paths_passed':20,'neutral_output_first_branch_opposes_active_pull':DISTAL_LOAD,'endpoint_summaries':summaries,'M5_full_installed_arm_motion_verified':False,'passive_steel_routing_verified':False,'crimp_or_strength_verified':False,'TPU_actuation_verified':False,'frozen':False,'evidence_sha256':{name:hashlib.sha256((DEST/name).read_bytes()).hexdigest() for name in evidence}}
    (DEST/'review_status.json').write_text(json.dumps(status,indent=2),encoding='utf8')
    text=['# M5 肌肉远端球头连接候选','', '新增实体端帽底座、分体球窝保持盖、球头/压接坯及M2紧固件。被动钢丝通过可转向端头连接至端帽，不再把端帽外的虚拟折点当作有支撑的导向。当前只确认端部组件的局部几何与安装，仍未完成整臂钢丝路线或TPU驱动。','']
    if DISTAL_LOAD:text+=['本版为M5B：球头中立出口已朝向原内部钢丝末段牵引方向的反向，零位参考下的反向阻抗投影为+1。两个输出支路各检查了距球心3.5–30mm的直线钢丝，未检出硬件或TPU相交。30mm端点只是空间检查范围，不是已经存在的导向点/支架。尚须安装真实转向，将这条远端支路接回前臂。原朝向旧力臂孔位的M5A保留在相邻目录，未覆盖。','']
    if DISTAL_LOAD:text+=['两处M2耳部新增一体环形止口，高0.3mm，控制球窝两半的间距；保持盖与止口接触后承接名义夹紧载荷，不靠一个悬空缝隙保持球头间隙。几何检查确认实体止口存在、与盖接触，并会阻挡继续闭合0.01mm；未验证螺钉预紧、材料压缩或蠕变。','']
    text+=['## 结构与尺寸','', '- 每块11件：端帽/球窝底座1、保持盖1、球头压接坯1、M2×10螺栓2、螺母2、垫圈4。原TPU端垫与两枚M3连接保留。','- 名义球头半径3.5mm、球窝内半径3.7mm/外半径5.2mm；颈半径1.2mm，压接坯外半径1.6mm、未压接孔半径0.45mm。未选定承载材料和实际压接工艺；螺钉仅按当前设计包络，实物头型/尺寸与工具空间仍须复核。','- 0.6mm钢丝端段在未压接孔内几何检查通过；孔有间隙，不能靠这个孔本身夹持钢丝。M1主动钢丝止挡仍在原端帽中心，和新的被动肌腱端头分开。','', '| 肌肉 | 球窝中心（端帽局部/mm） | 保持盖滚转/° | 刚性装入路径 | 转向采样 | 零位端头直拉阻挡/mm |','|---|---|---:|---:|---|---|']
    for s in summaries:text.append(f'| {s["muscle"]} | {s["housing_position_local_mm"]} | {s["retainer_roll_deg"]} | 10/10 | 24方向各0–25°往返 | {s["first_pullout_obstruction_range_mm"][0]:.1f}–{s["first_pullout_obstruction_range_mm"][1]:.1f} |')
    text+=['', '转向常规步长1°，干涉边界细化阈值0.1°；装入最后5mm步长0.1mm，其余不大于1mm。直拉防脱只检查零位端头朝向下28条直线、步长0.1mm，不证明所有姿态/曲线脱出、承载、防松或寿命。','', '伸肌侧原轴向螺母入口受端帽阻挡。已保留失败记录，并改为先从侧向放入后垫圈和螺母，再从正面装M2螺栓。需要手动保持未锁紧螺母，未暗示存在未建模的防转卡槽。20条局部安装路径均无刚性干涉，未检出TPU装入接触；这不是整臂完整装配通过。','', '## CAD、网格与演示','', '- `cad/swivel_arm_review.step`：129件中立位刚性实体；131件中立位场景另含两块TPU。所有实体按源OBJ与STEP顶点双向核对<0.001mm，整体回读有效。','- 两个模块各自`cad/swivel_terminal_assembly.step`为11件端部组件；分面BRep，不含原生特征树/装配配合。独立零件STL/STEP保留公共装配坐标。','- `terminal_articulation.gif`只用mocap规定球头的局部±25°转向，TPU保持零位；不是受力、收缩或整臂驱动仿真。12个代表姿态的编译网格旋转与源几何已核对。','- 近景图只为查看端部而隐藏其他骨骼；131件模型的零位相交另行检查，未以隐藏物体制造通过。','- 4件刚性外壳打印摆放候选位于各模块`print_candidates/`，带源哈希和坐标变换；未切片、未检查支撑/打印/实物公差。金属坯与紧固件不作为TPU承载打印件。','', '## 仍未通过：直连走线与输出方向','', '重新检查了真实球头中心至真实前臂止挡的直线，不在端帽外添加无支撑折点。这里只把端帽保持在中立位以筛查空间，未给TPU设置不存在的导轨或轴向运动规律。','', '| 肌肉 | 此筛查所需最大球头转向/° | 直线穿透硬件采样/151 | 球头金属与运动部件相交采样 | 反向阻抗投影（零位） |','|---|---:|---:|---:|---:|']
    for r in routes['muscles']:text.append(f'| {r["muscle"]} | {r["maximum_required_terminal_rotation_deg"]:.3f} | {r["straight_wire_penetration_samples"]} | {r["terminal_metal_interference_samples"]} | {r["output_work_sign_diagnostic"]["opposition_at_elbow_zero"]:.3f} |')
    text+=['', '直连旧前臂止挡/孔口的线路存在穿透，不能作为最终钢丝路线。表中端头相交是前臂转动、端帽仍保持零位的筛查结果，不能据此把局部转向通过转移成安装后的整臂通过。M5B的端口已面向远端输出，不把已拒绝的直连旧孔位当作球窝转角要求。','', '输出方向投影用原中立位内部钢丝末段的牵引方向作参考：正值表示外部肌腱反抗该方向，负值表示顺向。直连旧力臂线路的投影在151个姿态均为负，因此不能直接当作正确输出；M5B已把中立出口方向改为反向。此检查不涉及力的数值，也不证明TPU变形后的方向。下一步要建立有实体支撑的转向或重新布置端部/肌腹，同时检查输出方向与运动空间，之后再估算行程。','', 'M4的32.345/42.229mm来自已拒绝的轴向假设，不是本M5的新实测行程，也不是本端部的驱动能力。端头偏置还会给端帽/TPU端垫带来力矩，未验受力；不能把端帽强制锁在轴线上来隐藏这个问题。骨架本体沿用M4，但新增端部未做整个安装态的13路径/柔性运动验收，M4结论不自动转移。','', '## 复现','', '在项目根目录执行：','', '```powershell',"& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_distal_swivels.py --muscle FJ1486"+(' --distal-load' if DISTAL_LOAD else ''), "& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_distal_swivels.py --animate --cutaway"+(' --distal-load' if DISTAL_LOAD else ''), "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/review_distal_swivels.py"+(' --distal-load' if DISTAL_LOAD else ''),'```','', '完整重建顺序见rebuild_distal_swivels.py。新增零件、历史失败装配顺序、证据哈希均独立保存。肩部/肩胛肌群、整体机架、实际驱动和打印实测仍属于未完成目标。']
    (DEST/'SWIVEL_ENDPOINT_REVIEW.md').write_text('\n'.join(text)+'\n',encoding='utf8')
    if DISTAL_LOAD:
        out=OUT
        entry={'candidate':'M5B anatomical-arm geometry and distal-load endpoints','scope':'current review candidate, not frozen or full-arm actuation verified','candidate_directory':str(DEST.relative_to(ROOT)),'neutral_CAD':whole['whole']['STEP'],'neutral_MuJoCo':str((DEST/'neutral_view.xml').relative_to(ROOT)),'review':str((DEST/'SWIVEL_ENDPOINT_REVIEW.md').relative_to(ROOT)),'local_terminal_animation':animation['file'],'last_full_bone_path_evidence':str((BASE/'bone_motion_audit.json').relative_to(ROOT)),'M5B_installed_whole_arm_motion_verified':False,'TPU_actuation_verified':False,'physical_test_verified':False,'full_arm_muscle_layout_complete':False,'frame_and_motor_reaction_complete':False,'frozen':False}
        (out/'CURRENT_CANDIDATE.json').write_text(json.dumps(entry,ensure_ascii=False,indent=2),encoding='utf8')
        prefix=str(DEST.relative_to(out)).replace('\\','/')
        index=['# 当前人体骨骼仿生手臂候选入口','', '当前为M5B：真实解剖骨形＋工程关节候选＋两块异形TPU晶格肌肉＋远端球头端头。未定版，也未完成整臂驱动。','',f'- [本轮检查与限制]({prefix}/SWIVEL_ENDPOINT_REVIEW.md)',f'- [129件刚性中立位CAD]({prefix}/cad/swivel_arm_review.step)',f'- [整臂零位图]({prefix}/whole_neutral.png)',f'- [局部球头转向演示]({prefix}/terminal_articulation.gif)（规定姿态，TPU不动，不是驱动仿真）',f'- [装配与台架准备]({prefix}/ASSEMBLY_AND_BENCH_PROTOCOL.md)','', '131件中立位场景含两块TPU。新端头22件、20条局部安装路径与局部转向几何已核对；出口朝远端、反抗原中立位内部钢丝牵引。后续必须建立实体转向/接收端，再重新计算肌肉行程与变形。','', '最后完整骨架13路径检查为M4；新增端头不自动继承整臂运动通过。原直连旧力臂路线已拒绝，历史目录保留作对照。','', '在项目根目录查看当前整臂零位：','', '```powershell',"& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_distal_swivels.py --distal-load --whole",'```','', '待完成：实际钢丝走线与支撑、TPU非线性变形和实测、肩部/肩胛肌群总体驱动布局、整体机架与电机反力闭环、受力/制造验证。原始人体骨骼和肌肉来源文件均保留。']
        (out/'CURRENT_CANDIDATE.md').write_text('\n'.join(index)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in status.items() if k not in ['evidence_sha256','endpoint_summaries']},indent=2))
if __name__=='__main__':run()
