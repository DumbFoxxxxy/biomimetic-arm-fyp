"""Evidence-gated M7 review and candidate index; never an actuation pass."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_guide_ground_frames import *
from project_runtime import WorkspacePath as Path

def read(p):return json.loads(p.read_text(encoding='utf8'))
def magnitude(row):
    c=np.array(row['last_clear_controls_deg']);return float(max(abs(c[9]),*(np.linalg.norm(c[i:i+3]) for i in [0,3,6])))
def run():
    design=read(DEST7/'design.json');assembly=read(DEST7/'assembly_checks.json');exports=read(DEST7/'export_checks.json');cad=read(DEST7/'whole_arm_STEP_check.json');motion=read(DEST7/'bone_motion_audit.json');base=read(BASE/'bone_motion_audit.json');render=read(DEST7/'render_inputs.json');cables=read(DEST7/'neutral_supported_branch_checks.json')
    assert all(m['selected'] for m in design['modules'])
    assert len(assembly['paths'])==16 and assembly['all_rigid_paths_clear'] and assembly['TPU_contact_flag_count']==0 and assembly['all_seats_pass']
    for p,h in assembly['source_sha256'].items():assert digest(Path(p))==h
    assert exports['all_pass'] and len(exports['registration'])==109 and exports['assembly_sha256']==digest(DEST7/'assembly.json')
    assert cad['whole']['pass'] and cad['whole']['solid_count']==151 and cad['changed_subassemblies']['solid_count']==24 and cad['fresh_changed_vertex_checks']==24 and cad['unchanged_parent_proofs_reused']==127 and cad['assembly_sha256']==digest(DEST7/'assembly.json')
    for r in cad['parts'].values():assert digest(ROOT/r['source_OBJ'])==r['OBJ_sha256'] and digest(ROOT/r['source_STEP'])==r['STEP_sha256']
    assert cad['parent_manifest_sha256']==digest(CURRENT/'whole_arm_STEP_check.json')
    for part in ['whole','changed_subassemblies']:assert digest(ROOT/cad[part]['path'])==cad[part]['sha256']
    assert motion['complete_prescribed_suite'] and len(motion['actions'])==13 and not motion['neutral']['hits']
    for p,h in motion['mesh_evidence_sha256'].items():assert digest(Path(p))==h
    for p,h in motion['input_evidence_sha256'].items():assert digest(DEST7/p)==h
    for p,h in render['input_sha256'].items():assert digest(Path(p))==h
    assert render['meshes']==155 and render['assembly_sha256']==digest(DEST7/'assembly.json')
    assert cables['all_clear'] and len(cables['checks'])==4 and cables['source_assembly_sha256']==digest(DEST7/'assembly.json')
    for p,h in cables['source_sha256'].items():assert digest(Path(p))==h
    reference={r['probe']:r for r in base['actions']};rows=[];regressions=[]
    for r in motion['actions']:
        old=magnitude(reference[r['probe']]);actual=magnitude(r);regressed=actual<old-.1000001
        if regressed:regressions.append(r['probe'])
        rows.append({'action':r['probe'],'M4_reference_amplitude_deg':old,'M7_continuous_amplitude_deg':actual,'actual_humerus_elevation_deg':r['actual_humerus_elevation_deg'],'first_limit':r['first_limit'],'range_regression':regressed})
    result={'scope':'M7 current geometry candidate audit:only determinate rigid bodies and local changed installation, no installed compliant actuation proof','rigid_suite_complete':True,'range_comparison':rows,'regressions_vs_M4':regressions,'rigid_instances':109,'neutral_instances':153,'free_TPU_or_distal_instances_omitted_from_motion':44,'CAD_neutral_rigid_solids':151,'new_or_revised_parts':24,'complete_arm_assembly_verified':False,'TPU_actuation_verified':False,'physical_test_verified':False,'full_arm_muscle_layout_complete':False,'downstream_transmission_complete':False,'motor_reaction_complete':False,'frozen':False}
    (DEST7/'review_status.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    table=['| 动作 | M4 对照幅值 ° | 本轮连续幅值 ° | 实际肱骨方向变化 ° | 首个限制 |','|---|---:|---:|---:|---|']
    for r in rows:
        limits='；'.join(' / '.join(v['pair']) for v in r['first_limit']) or '该规定路径内未检出穿透'
        table.append(f"| {r['action']} | {r['M4_reference_amplitude_deg']:.4f} | {r['M7_continuous_amplitude_deg']:.4f} | {r['actual_humerus_elevation_deg']:.4f} | {limits} |")
    text='''# M7 后方绕行的实体导向轮支架候选

本轮把 M6 首段槽轮连接到既有两枚肱骨适配件：适配件、支架和键式接收座为连通实体，
沿用原来的两组肱骨 M3 安装孔。骨长、骨骼关节中心和骨形均未修改。
原模型、M5B 及两轮失败支架保留。本候选未冻结，也不是完整肌肉驱动手臂。

## 实际结构和装配

每块肌肉：原适配件的低位支撑枝延伸到肘后方，再连接槽轮叉架的键式接收座。
接收座承靠面位于槽轮局部 Z=-5 mm，与键后面名义接触；实测几何间隙为 0，
再向里移动 0.01 mm 会发生实体重叠，避免用悬空间隙假装夹紧止口。
键的侧隙每边 0.15 mm，尚未考虑打印和安装公差。

固定钢丝力的名义路径：槽轮 → 轴套/轴螺栓 → 叉架 → 键式座/M3 紧固件 →
连通支架 → 原肱骨双螺栓座。该路径存在实际实体连接，不是仿真中凭空固定的点。
轴承摩擦、实际预紧、杆件刚度、屈曲、疲劳和可制造性仍需材料与打印试验。

原正面 M3×20 方案的安装入口被前臂挡住；先装导向轮又挡住原前臂插入路径。
现改为名义 18 mm 轴长的 M3 背面螺栓，先从侧面放入前侧松螺母和垫圈，
再装背面螺栓与后垫圈。具体商品、螺纹和工具空间未选定，不冒充完成紧固验证。
屈肌适配件先在安装面前 5 mm 的位置侧向移入，再沿法向就位；不是直线穿过前臂。

16 条改变后的局部刚性安装路径通过，末 5 mm 采样 0.1 mm，其余不超过 1 mm。
两套导向轮的插入检查保留全部中立位肌肉周边件；这些路径的 TPU 相交标记为 0。
裸适配件在肌肉装入前安装；本轮并未解决所有 TPU 肌腹的整机安装和变形问题。

## 本轮证据

- 10 件新制造/修订零件：2 枚含支架的适配件及 8 件键式安装紧固件。
- 再包含 M6 两套各 7 件槽轮，共 24 件新增或修订的刚性几何。
- 153 件中立位制造几何；其中 2 件是 TPU，151 件刚性中立位 CAD。
- 动作检查含 109 件确定随骨骼运动的刚性几何，另有胸廓等骨性参考参与碰撞。
  44 件 TPU/自由远端及旧的未压接止挡没有强加运动约束，未转移通过结论。
- 对本轮 24 件逐项复核 CAD 与源 OBJ 顶点；127 件不变部分的既有配准证据经源文件
  哈希和变换一致性核验后复用。整臂 STEP 重新读回，实体数与体积均检查。
- 109 件 MuJoCo 中立位网格配准误差最大约 0.0000204 mm；10 件新导出 STL/OBJ 复核。

下表幅值是对应控制路径中最大的关节转角幅值；联合路径的实际抬升看单列肱骨方向，
不以关节角相加。动作按不大于 1° 采样，首个碰撞区间细化到不大于 0.1°，
从零位连续计入，回程重新计算。碰撞阈值 0.001 mm³；这是有限路径检查，不是全姿态证明。
名义接触面会使最小间隙为 0；完整逐姿态间隙和首个限制见 motion_paths.json，
间隙大于 3 mm 的数值被截到 3 mm，不表示所有非接触区域有制造余量。

'''+ '\n'.join(table)+'''

## 仍需继续

槽轮后的 40 mm 钢丝检查段仍为自由末端，没有接上前臂接收机构，也没有产生驱动。
CAD 中立位组合保留 M5B 球头的名义姿态；静态走线图将这两个自由球头分别转向
10.7873°/15.4759°，使实际出口对准第一段切线。两者是同一零件几何的不同允许姿态，
不是两套不同孔道；该图示姿态及两条钢丝检查管已对全部153件中立位制造几何复核。
TPU 肌腹、自由端盖在整臂图中保持零位，未证明收缩量、拮抗伸长、碰撞后的变形或恢复。
24 件接口及 109 件刚性路径的通过不能推广到完整的柔性肌肉运动。
肩胛/肩部肌肉总体驱动布局、电机与导管反力闭环、材料标定及实物试验仍在范围内。

## 文件和复现

- `cad/grounded_guides_arm_review.step`：151 件刚性中立位组合，含自由端头/未压接坯。
- `cad/guide_frame_subassemblies.step`：24 件本轮接口组合，方便单独审查。
- `whole_arm.xml` / `whole_arm.png`：155 个实际网格的静态 MuJoCo 包络，含 2 段钢丝检查管。
- `arm.xml`：109 件确定刚体的关节模型；不含未验证的柔性运动。
- `assembly_checks.json`、`bone_motion_audit.json`、`export_checks.json`、`whole_arm_STEP_check.json`。

从项目根目录运行生成/验证脚本时均使用 `--swept-root`：build_guide_ground_frames.py、
prepare_grounded_guide_model.py、check_grounded_guide_assembly.py、audit_grounded_guides.py、
check_grounded_guide_exports.py、export_guide_ground_frames.py、export_grounded_guide_arm.py。
其中 CAD 脚本使用 D:/Python3.10.6/python.exe，其余几何脚本使用 D:/MujocoEnv/Scripts/python.exe；
渲染使用 C:/Python314/python.exe render_first_guides.py --grounded --swept-root。
脚本均位于 mujoco_arm/anatomical_concept/，生成在独立目录，不覆盖原始骨肌模型。
'''
    (DEST7/'GROUNDED_GUIDE_REVIEW.md').write_text(text,encoding='utf8')
    if regressions:print('M7 not adopted; actual range regressions',regressions,flush=True);return
    prefix=str(DEST7.relative_to(ROOT));index={'candidate':'M7 anatomical arm, M5B endpoints and posterior first-guide frames','scope':'current geometry candidate, not complete drive or frozen','candidate_directory':prefix,'neutral_CAD':str((DEST7/'cad/grounded_guides_arm_review.step').relative_to(ROOT)),'neutral_MuJoCo':str((DEST7/'whole_arm.xml').relative_to(ROOT)),'review':str((DEST7/'GROUNDED_GUIDE_REVIEW.md').relative_to(ROOT)),'last_full_bone_path_evidence':str((DEST7/'bone_motion_audit.json').relative_to(ROOT)),**{k:result[k] for k in ['complete_arm_assembly_verified','TPU_actuation_verified','physical_test_verified','full_arm_muscle_layout_complete','downstream_transmission_complete','motor_reaction_complete','frozen']}}
    (OUT/'CURRENT_CANDIDATE.json').write_text(json.dumps(index,indent=2),encoding='utf8');relative=DEST7.relative_to(OUT).as_posix();(OUT/'CURRENT_CANDIDATE.md').write_text(f'''# 当前人体骨骼仿生手臂候选入口

当前为 M7 后方绕行的首段导向轮支架候选，沿用人体骨骼与两块异形 TPU 肌肉、M5B 球头输出端。
支架具有实体骨侧连接。未定版，完整钢丝传动、TPU 变形与整臂驱动未完成。

- [本轮检查及限制]({relative}/GROUNDED_GUIDE_REVIEW.md)
- [151件刚性中立位CAD]({relative}/cad/grounded_guides_arm_review.step)
- [24件接口CAD]({relative}/cad/guide_frame_subassemblies.step)
- [MuJoCo整臂静态图]({relative}/whole_arm.png)
- [109件确定刚体动作检查]({relative}/bone_motion_audit.json)

中立位153件，动作模型109件；44件TPU/自由远端没有施加虚构运动约束。
槽轮后的钢丝仍是自由检查段，尚无前臂接收端或实际驱动行程。
原始模型、M5B及失败支架候选均保留。肩部/肩胛肌群、电机反力闭环和实测仍待完成。
''',encoding='utf8');print('M7 evidence checks passed; candidate index updated, goal remains incomplete',flush=True)

if __name__=='__main__':run()
