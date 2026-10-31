from build import *
import shutil

def run():
    a=json.loads((OUT/'source_audit.json').read_text(encoding='utf8'));r=json.loads((OUT/'range_audit.json').read_text(encoding='utf8'));muscles=json.loads((OUT/'muscle_parameters.json').read_text(encoding='utf8'));sc=json.loads((OUT/'scapula_task_audit.json').read_text(encoding='utf8'))
    checks=[];data=json.loads((OUT/'cad_input.json').read_text(encoding='utf8'))
    for p in data['parts']:
        raw=(OUT/'assets'/(p['name']+'.stl')).read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];expected=np.array(p['vertices'])[p['faces']];err=float(np.max(np.abs(tri-expected)));checks.append({'part':p['name'],'CAD_input_vs_export_STL_max_coordinate_error_mm':err,'pass':err<1e-4})
    for row in a['bones']:
        p=ROOT/row['source'];checks.append({'source':row['bone'],'source_stl_sha_unchanged':hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],'pass':hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']})
    dump(OUT/'geometry_consistency.json',{'checks':checks,'pass':all(x['pass'] for x in checks),'qualification':'CAD input versus exported STL, source preservation; separate STEP report validates BRep/volume/bounds. No FEM or printed validation.'})
    with (OUT/'anchors_REFERENCE_ONLY_mm.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['muscle','role','x_mm','y_mm','z_mm','owner','status'])
        for md in muscles:
            for role,owner in [('A_sheath_reaction_mm','humerus'),('C_wire_terminal_and_tendon_mm','moving_TPU_end'),('B_bone_insertion_mm','forearm')]:w.writerow([md['id'],role,*md[role],owner,'candidate interface; NOT drill-ready'])
    with (OUT/'wire_paths_REFERENCE_ONLY_mm.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['muscle','segment','index','x_mm','y_mm','z_mm'])
        for md in muscles:
            for i,p in enumerate(md['channel_centerline_mm']):w.writerow([md['id'],'active wire channel reference',i,*p])
            for i,p in enumerate([md['C_wire_terminal_and_tendon_mm'],md['B_bone_insertion_mm']]):w.writerow([md['id'],'passive tendon neutral reference',i,*p])
    paths=json.loads((OUT/'motion_paths.json').read_text(encoding='utf8'))
    with (OUT/'elbow_tendon_reference_NOT_TPU_STROKE.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['elbow_deg','flexor_C0_B_mm','extensor_C0_B_mm','qualification'])
        for p in paths['K1']:w.writerow([p['q_deg'][-1],*p['neutral_C_to_moving_B_distance_mm'],'C held at neutral on humerus; NOT actual TPU shortening'])
    txt=['# A1 人体骨骼与晶格肌肉首轮报告','','**状态：研究候选，未定版；不能按整臂打印装配或宣称肌肉已实现驱动。**','','## 已完成的数字原型','','- 以源文件尺寸建立五骨整臂，尺桡骨同一刚体，SC2/AC1/GH3/肘1。','- 源模型毫米制、OBJ/STL平移关系、哈希、闭合性和关节点拟合已登记。','- 原骨形未切削、未缩放，所有原始文件及V2-S1保留。',f"- GH球面拟合半径 {a['GH_fit']['radius_mm']:.3f} mm，径向RMS {a['GH_fit']['radial_rms_mm']:.3f} mm。",f"- 肘轴工程校准位移 {a['engineering_axis_calibration']['displacement_mm']:.3f} mm，轴向变化 {a['engineering_axis_calibration']['axis_change_deg']:.3f}°；骨面本身未修改。",'- 生成两块完整几何肌肉：变截面BCC晶格、外表面交叉网格、中心孔道和带孔端板。','- 全臂其他肌群提供候选体积和路径，明确标为 VOLUME_ONLY。','','## 骨骼路径结果','','以下结果仅针对骨面、理想关节及规定路径，不包含参考关节硬件或TPU驱动能力。','常规步长≤1°，首次干涉边界≤0.1°；布尔穿透报告阈值0.001 mm³。','最小间隙搜索上限5 mm，报告5代表至少5，不代表恰好5。','联合抬升测量实际肱骨方向，未把关节角简单相加。','','|动作|原目标|原底线|连续采样范围|首次限制|','|---|---:|---:|---:|---|']
    for row in r['actions'][:7]:
        limit='；'.join('/'.join(h['pair']) for h in row['first_limit']) or '所测路径未见骨面穿透，未证明最大范围'
        txt.append(f"|{row['action']}|{row['target']}°|{row['floor']}°|{row['actual_continuous_sampled_prefix']:.2f}°|{limit}|")
    txt+=['','## 肩胛额外动作的独立检查','','保留机构链后，SC转角不能代表纯平移。指定纯前后平移及绕肩胛法向旋转时，','位置与方向不能同时精确跟踪。允许误差为位置0.1 mm、方向0.1°，误差范围内的','微小前缀不视为已经覆盖验收动作；本轮四项均不登记整体通过。','','|动作|指定目标|误差容许下的路径前缀|验收|','|---|---:|---:|---|']
    for row in sc['actions']:txt.append(f"|{row['action']}|{row['target']} {row['units']}|{row['continuous_tracking_prefix']:.3f} {row['units']}|未通过：姿态跟踪／碰撞联合验证未完成|")
    txt+=['','## 肘部肌肉组件','','|组件|有效段端距|孔道|晶格杆径|连通性|','|---|---:|---:|---:|---|']
    for md in muscles:txt.append(f"|{md['id']}|{md['active_length_mm']:.2f} mm|2.4 mm|2 mm|单一材料连通体；零位未穿骨|")
    txt+=['','A 为导管反力与肌腹固定端，C 为钢丝终端和短肌腱连接端，B 为骨端候选插入点。','当前未设计完成骨端夹具和A反力支架，因此A/B不能作为钻孔定稿点。','C不被虚构滑块约束。零位C到运动B的距离仅用于路径观察，不直接换算TPU缩短率。','','## 尚未完成、必须继续的工作','','1. 参考销轴、轴套存在7处与骨面的零位重叠；需要真实孔腔、安装桥、SC两轴实体结构和GH可拆装保持器。当前不是完整的工程关节装配。','2. 必须在上述硬件落实后重新完成全部往返碰撞检查；当前骨面范围不能当成装配范围。','3. 肩胛纯滑移/内外旋与当前SC2/AC1链存在姿态跟踪矛盾，需要在保留链下重定义可实现的耦合动作，或另行明确机构调整；不私自增设自由度。','4. 肌肉附着夹具、通道大变形、肌肉间接触和全程骨面接触尚未验收。','5. 单肌肉打印、可重复行程、恢复、摩擦和端部损伤均待实测；没有25%行程能力保证。','','## CAD与演示','','`cad/` 为分面BRep STEP；`assets/*_print_mm.stl` 是摆正的单肌肉样件。','STEP回读结果见 cad_roundtrip.json；几何一致性见 geometry_consistency.json。','neutral_review.png 展示零位真实晶格；bone_motion.gif 只展示骨面已通过的路径。','查看脚本在移动时隐藏零位刚体肌肉，不能将画面当成柔性仿真。','生成、查看与手动测试说明见随附 README.md 和 PRINT_AND_TEST.md。']
    (OUT/'REVIEW.md').write_text('\n'.join(txt)+'\n',encoding='utf8')
    for n in ['README.md','PRINT_AND_TEST.md']:shutil.copyfile(Path(__file__).parent/n,OUT/n)
    # Track all delivered files, without claiming completion/freeze.
    manifest=json.loads((OUT/'manifest.json').read_text(encoding='utf8'));manifest['files']={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'};manifest['known_blockers']=['engineering joint hardware incomplete','scapular pose tasks fail tracking','muscle installation and actuation unvalidated'];dump(OUT/'manifest.json',manifest)
    print('consistency checks',len(checks),'pass',all(x['pass'] for x in checks),flush=True)

if __name__=='__main__':run()
