"""Evidence-bound review of G1; ideal elbow and muscle drive remain open."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact'
DEST=BASE/'GH_connection'
def read(p):return json.loads(p.read_text(encoding='utf8'))

def run():
    assembly=read(DEST/'assembly.json');audit=read(DEST/'bone_motion_audit.json');scene=read(DEST/'audit_checks.json')
    geometry=read(DEST/'geometry_checks.json');cad=read(DEST/'cad_roundtrip.json');compound=read(DEST/'assembly_STEP_check.json');capture=read(DEST/'capture_probes.json');paths=read(DEST/'motion_paths.json')
    checks={'12_connected_parts':len(assembly['components'])==12 and set(assembly['components'].values())=={1},'neutral_assembly':not assembly['neutral_overlaps'] and not audit['neutral']['hits'],'OBJ_STL':geometry['all_geometry_checks_pass'],'individual_STEP':len(cad['parts'])==12 and all(p['pass'] for p in cad['parts']),'compound_STEP':compound['pass'],'36_scene_geometry':len(scene['zero_pose_geometry'])==36 and all(r['pass'] for r in scene['zero_pose_geometry']),'sample_step':scene['maximum_sampled_joint_increment_deg']<1+1e-6,'13_paths':len(audit['actions'])==13,'26_pullout_directions':capture['all_26_directions_blocked_within_3mm']}
    if not all(checks.values()):raise RuntimeError(checks)
    evidence=[DEST/'assembly.json',DEST/'cad_input.json',DEST/'arm.xml',DEST/'bone_motion_audit.json',DEST/'geometry_checks.json',DEST/'cad_roundtrip.json',DEST/'assembly_STEP_check.json',DEST/'audit_checks.json',DEST/'capture_probes.json']
    status={'scope':'G1 finite-path geometry and roundtrip checks only; NOT whole-arm completion','checks':checks,'frozen':False,'evidence_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in evidence},'not_verified':['elbow physical connection','muscle interfaces/channels','steel+TPU deformation and drive','arbitrary escape path','manufacturing fit, load, strength and life']}
    (DEST/'verification_status.json').write_text(json.dumps(status,indent=2),encoding='utf8')
    lines=['# G1 盂肱实体连接候选（未定版）','',
      '以SC/AC紧凑候选为基础，GH中心保持源拟合位置。将肱骨头局部改为半径16mm的工程球头，颈部半径4.5mm；球窝内半径16.4mm、外半径19.5mm。两条颈部运动让位槽使用半径5.3mm的扫掠圆柱构造，保留两个独立防脱盖，每盖采用名义M3×12螺栓、螺母和两片垫圈。螺纹未建模。', '',
      '骨侧过渡、球窝、防脱盖及附件共12个连通实体。盖板/附件随肩胛（AC刚体）运动，球头及肱骨随GH刚体运动。全模型25个活动实体对42个固定骨参考/SC附件检查；没有排除相邻骨对。肘仍为理想铰链，未提供实体销轴和防脱。', '',
      '## 修改与装配修正','',
      '初版新增连接与既有锁骨及AC保持盖发生零位相交，已保存于iterations/initial_neutral_conflict/。修订版在真实相交区域开局部间隙，没有关闭碰撞对；之后重新检查零位与全路径。球窝曾有约529mm³的悬空块，加入明确的外侧连接支撑后成为单一实体。局部间隙修正产生0.00406mm³细碎残片，通过记录中心/半径的局部轮廓修整去除；没有删掉实际连接件。具体工具区域、连接端点和几何增减量见assembly.json。', '',
      f"肱骨相对原骨形移除{assembly['geometry_delta']['humerus_removed_mm3']:.1f}mm³，肩胛相对AC基线移除{assembly['geometry_delta']['scapula_removed_mm3']:.1f}mm³并增加{assembly['geometry_delta']['scapula_added_mm3']:.1f}mm³。此方案明显重做了关节近端，不能把它称为未改动的人体关节。原输入均保留。", '',
      '## 指定路径的连续前缀','',
      '每个探针从零位开始，正反向采样≤1°，首次干涉边界细化≤0.1°。未检出是指交叠体积不超过0.001mm³；间隙查询封顶3mm，不是数学上的全姿态证明或打印性能验收。组合运动测量实际肱骨方向，不把关节角相加。表中间隙为已接受前缀内完整1°采样点的最小值，不包含未存入路径数组的细化端点。', '',
      '| 探针 | 请求输入° | 连续前缀输入° | 肱骨方向变化° | 采样最小间隙mm（封顶3） | 首个限制 |', '|---|---:|---:|---:|---:|---|']
    for row in audit['actions']:
        c=np.array(row['requested_controls_deg']);a=np.array(row['last_clear_controls_deg'])
        angle=lambda v:max(abs(v[9]),*(np.linalg.norm(v[i:i+3]) for i in [0,3,6]))
        name=row['probe'];limit='；'.join(' / '.join(h['pair']) for h in row['first_limit']) or '请求路径内未检出'
        safe=[p for p in paths[name] if angle(np.array(p['controls_deg']))<=angle(a)+1e-8]
        gap=min(p['min_bone_gap_mm_capped_3'] for p in safe)
        lines.append(f"| {name} | {angle(c):.3f} | {angle(a):.3f} | {row['actual_humerus_elevation_deg']:.3f} | {gap:.4f} | {limit} |")
    lines+=['', '前屈与外展首次完整角度碰撞样本已分区：均发生在新增工程几何，不在保留的原肩胛骨实体中；详见limit_localization.json。后续仍可通过连接外形/支撑路由优化，当前未把范围降低改写为“满足原目标”。', '',
      '注意：全体间隙最小值为0，零位最近对为工程肩胛/FJ3384与AC_cap_right（同属AC刚体）。虽然未检出超过体积阈值的穿透，这仍不是正的制造间隙证明，不能据此直接冻结公差或打印装配。完整模型仍保留该对的碰撞检查。', '',
      '## 肌肉附着影响','',
      '候选解剖骨面不是冻结的工程锚点。以下区域的原候选点距新骨面超过1mm；其中肩袖近端偏差较大，必须设计新的骨侧连接端、过渡及短肌腱，不能把源肌肉原位叠加当成已经连上球头。', '',
      '| 肌肉编号 | 骨面 | 超过1mm点数 | 最大偏差mm |','|---|---|---:|---:|']
    for row in geometry['old_anatomical_patch_impact']:
        if row['points_more_than_1mm_from_new_surface']:lines.append(f"| {row['muscle']} | {row['structure']} | {row['points_more_than_1mm_from_new_surface']} | {row['maximum_deviation_mm']:.3f} |")
    lines+=['', '## 验证范围','',
      '12个独立STEP和12实体复合STEP均通过回读，网格与场景的36项位置检查通过。STEP为三角面BRep，保留独立实体，不含SolidWorks原生特征树/装配配合。球头26方向直拉探针在0.9mm内均遇到球窝/盖；这是有限直线检查，不是所有曲线路径的防脱证明。真实螺栓预紧、打印间隙、工具可达性、强度及寿命尚未验证。', '',
      '保留更多骨形的GH_preserve_head比较候选有独立文件夹。它保留了一个原肱骨局部，并以弧形连接连续接回骨干，但在很小的前屈/外展范围即相交；肱骨件STEP导出还报Invalid BRep。该候选未被选用，不继承G1的CAD或运动通过结论。', '',
      'TPU肌腹、钢丝通道、端部载荷路径和驱动尚未建立，本阶段演示隐藏肌肉。7个肋软骨不作为刚性骨验收，尚未验证其柔性接触。下一步是肘的实体连接，以及骨侧肌肉接口和成形肌腹/通道；不是继续把理想约束称为已制造关节。', '',
      '## 复现','', '```powershell', "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/build_GH_connection.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/audit_joint_model.py --GH", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_SC_assembly.py --GH", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_GH_capture.py", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_cad.py --GH", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_assembly.py --GH", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_joint_audit.py --GH", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/diagnose_GH_limits.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/write_GH_review.py", "& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_reference.py --GH --supports", '```','', '加--detail显示关节细节；按1–5播放已采样前缀。渲染用--render，动作GIF由render_joint_motion.py --GH生成。']
    (DEST/'GH_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf8');print(checks)
    print('GH review',DEST/'GH_REVIEW.md')

if __name__=='__main__':run()
