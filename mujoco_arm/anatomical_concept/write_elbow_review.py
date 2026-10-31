"""Evidence-bound E1 and selected whole-skeleton review; not a drive acceptance."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
def read(p):return json.loads(p.read_text(encoding='utf8'))

def run():
    a=read(DEST/'assembly.json');motion=read(DEST/'bone_motion_audit.json');geometry=read(DEST/'geometry_checks.json');scene=read(DEST/'audit_checks.json');retention=read(DEST/'assembly_and_retention_checks.json');cad=read(DEST/'cad_roundtrip.json');compound=read(DEST/'assembly_STEP_check.json');whole=read(DEST/'whole_arm_STEP_check.json')
    checks={'nine_connected_parts':len(a['components'])==9 and set(a['components'].values())=={1},'neutral_all_pairs':not a['neutral_overlaps'] and not motion['neutral']['hits'],'OBJ_STL':geometry['all_geometry_checks_pass'],'nine_individual_STEP':len(cad['parts'])==9 and all(r['pass'] for r in cad['parts']),'nine_body_STEP':compound['pass'],'whole_42_body_STEP':whole['pass'],'42_scene_geometry':len(scene['zero_pose_geometry'])==42 and all(r['pass'] for r in scene['zero_pose_geometry']),'sample_step_deg':scene['maximum_sampled_joint_increment_deg']<=1+1e-6,'13_motion_paths':len(motion['actions'])==13,'forearm_straight_retention':retention['all_forearm_directions_blocked_within_2mm'],'pin_axial_retention':retention['both_pin_directions_blocked_within_2mm'],'sequential_assembly':retention['sequential_installation_geometry_pass']}
    if not all(checks.values()):raise RuntimeError(checks)
    evidence=['assembly.json','cad_input.json','arm.xml','bone_motion_audit.json','geometry_checks.json','audit_checks.json','cad_roundtrip.json','assembly_STEP_check.json','whole_arm_STEP_check.json','assembly_and_retention_checks.json']
    status={'scope':'finite rigid geometry, nominal assembly and CAD review only; full musculoskeletal goal incomplete','checks':checks,'frozen':False,'sha256':{n:hashlib.sha256((DEST/n).read_bytes()).hexdigest() for n in evidence},'not_verified':['steel/TPU muscle geometry and actuation','whole-arm support fixture','printed tolerances and preload','general retention escape paths, loads and strength','arbitrary combined motion']}
    (DEST/'verification_status.json').write_text(json.dumps(status,indent=2),encoding='utf8')
    lines=['# E1 肘实体连接及整臂骨架候选','',
      '当前文件夹为新解剖路线的工程候选，不是定版或驱动通过。SC、AC、GH和肘现在均有真实连接几何；MuJoCo仍沿用理想关节约束驱动预设姿态，不能据此证明实物无松动、能够承载或受TPU肌肉驱动。', '',
      '## 机构与实际连接','',
      '保留独立拟合的肘心和转轴。肱骨侧为一体式双耳叉架，前臂侧为带轴套的中心座；两侧限位套伸入叉架孔，提供名义径向定位，销轴、两片大垫圈及M6螺母形成名义轴向保持。轴套随前臂，限位套/销轴/防脱附件随肱骨。', '',
      '尺骨与桡骨保留源位置，用两条直径6mm的实际连接桥组成一个连通打印件。FJ3349几何并未删除：桡骨实体包含在forearm_elbow_hub/FJ3391中，避免重复计数和重复显示。首轮不具备旋前旋后；以后需要该自由度时必须重新分体和设计关节。', '',
      '| 参数 | 名义尺寸 |', '|---|---:|', f"| 叉架两耳内侧间距 | {2*a['fork_half_gap_mm']:.1f}mm |", f"| 单耳厚度 | {a['ear_width_mm']:.1f}mm |", '| 前臂中心座宽度 | 16mm |','| 销轴直径 | 6mm |','| 轴套外径/内径 | 8 / 6.3mm |','| 限位套外径/内径 | 11 / 6.6mm |','| 叉架限位套孔直径 | 11.7mm |','| 大垫圈外径/厚度 | 18 / 1mm |',f"| 自定义销轴杆长 | {a['custom_pin_shank_length_mm']:.1f}mm |", '| 前臂中心座与两侧限位套名义单侧轴向间隙 | 0.2mm |', '',
      '销轴/螺母/垫圈是名义硬件包络，螺纹未建模；没有选定可直接购买的标准销轴。轴套、骨侧打印件材料和公差须经试装确定，不把全部零件都当成同一打印材料。', '',
      '## 可装配性修正','',
      '初版存在骨侧与限位套/附件交叠，已保存在iterations/initial_overlap。孔腔和端面修正后，虽然能屈伸，前臂仍无法装入；该中间版本保存在iterations/before_assembly_opening。', '',
      f"修订版增加局部侧方装配开口，工具只作用在肘心39mm半径内，并明确保留叉架耳部、侧梁和上端连接。该开口另移除肱骨侧{a['assembly_opening']['humeral_material_removed_mm3']:.1f}mm³材料；不是删除连接或关闭碰撞对。限位套孔改为贯通，允许分步装入。", '',
      '检查通过的安装顺序：','', '1. 将轴套放入前臂中心座。','2. 销轴和限位套尚未安装时，前臂连同轴套沿记录的侧方方向就位。','3. 两侧限位套分别从轴向穿过叉架孔装入。','4. 带头侧大垫圈插入销轴。','5. 放入另一片大垫圈并安装螺母。', '',
      '检查是理想刚体几何：前臂接近路径1mm采样，最后5mm细化到0.1mm；销轴和套筒路径1mm采样。螺母仅检查包络就位，不证明实际旋拧、工具可达性或预紧合适。整组前臂+套筒同时塞入的旧路径仍失败，当前采用分步顺序，不把失败项隐藏。', '',
      '## 运动与保持检查','',
      f"完整模型包含{motion['moving_component_count']}个活动实体和{motion['fixed_reference_count']}个固定骨参考/附件，保留所有活动实体对的碰撞检查。源模型5块上肢骨在本候选中为4个活动骨侧组件，尺桡骨已实体合并。13条规定路径往返采样≤1°，边界细化≤0.1°，实体交叠体积阈值0.001mm³。", '',
      '| 探针 | 请求最大输入° | 连续前缀最大输入° | 实际肱骨方向变化° | 首个限制 |', '|---|---:|---:|---:|---|']
    for row in motion['actions']:
        def angle(v):return max(abs(v[9]),*(np.linalg.norm(v[i:i+3]) for i in [0,3,6]))
        limit='；'.join(' / '.join(h['pair']) for h in row['first_limit']) or '规定路径内未检出'
        lines.append(f"| {row['probe']} | {angle(row['requested_controls_deg']):.3f} | {angle(row['last_clear_controls_deg']):.3f} | {row['actual_humerus_elevation_deg']:.3f} | {limit} |")
    lines+=['', '前臂六方向直移以及销轴组件正反轴向直移均在2mm内遇到保持结构。仅证明这些采样方向，不证明任意脱出路径、打印后的承载或轴向游隙。源SC/AC/GH处仍有零间隙界面；全模型最小间隙不构成正的制造余量证明，骨架未冻结。7个肋软骨没有按刚体骨验收，柔性接触未知。', '',
      '## 肌肉接口影响','',
      '全局候选附着点偏差见geometry_checks.json，其中肩袖8–12mm偏差主要来自前阶段GH修订。肘部新增修订使FJ1486肱肌起始区的33个候选点距骨面超过1mm、最大约6.45mm；FJ1480三头肌内侧头有1点约1.60mm。肱肌工程连接端需重新选面，不能直接冻结原源模型区域。当前三头肌/肱肌尺骨侧及肱二头肌桡骨侧候选点未出现超过1mm的骨面偏差。候选接近面不等于已经确认的解剖标志或螺钉位置。', '',
      '## CAD与查看','',
      '9个独立肘部STEP和9实体复合STEP均通过回读；另按层级替换旧骨件，导出整臂42个独立实体的anatomical_arm_review.step，并重新导入检查实体数、有效性和总体积。每个选定实体的来源及哈希见whole_arm_STEP_check.json，未重复叠加被替换的骨模型。STEP是三角面BRep，不含SolidWorks原生特征树/装配配合。生成脚本及JSON参数仍可编辑。', '',
      '整臂STEP只包含选定上肢/SC骨侧基座及硬件，不包含全部胸廓参考、独立机架或TPU肌肉。实际固定支架、电机与导管反力安装仍需设计。', '',
      '```powershell', "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/build_elbow_connection.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/audit_joint_model.py --elbow", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_SC_assembly.py --elbow", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_elbow_retention.py", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_cad.py --elbow", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_assembly.py --elbow", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_arm_review_cad.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_joint_audit.py --elbow", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/write_elbow_review.py", "& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_reference.py --elbow --supports", '```','',
      '查看器加--detail看肘部；按1–5播放已检查前缀，运动时隐藏静态肌肉。render_elbow_assembly.py生成MuJoCo装配动画；render_joint_motion.py --elbow生成整臂动作演示。二者均不模拟实际肌肉驱动。', '',
      '## 下一阶段','',
      '在本骨架候选上先完成肘部前后拮抗TPU肌肉：依据真实肌腹外形和工程端部定义中心线/截面，生成受控曲面网格、内部晶格、实际钢丝导向通道及端部过渡。明确钢丝终端、导管反力和骨侧连接，不把普通拉索驱动或肌肉渲染占位当成完成。先计算实际所需行程和力臂，再设计打印样件验证；25%压缩假设不能当成通过能力。整臂肌群扩展、支架/电机安装及物理验证仍属于完整目标。']
    (DEST/'ELBOW_AND_ARM_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf8');print(checks)

if __name__=='__main__':run()
