"""Record the integral-rotor experiment, including its failed assembly gate."""
from check_rotor_pocket import ROOT,BASE,DEST
import json,hashlib

def run():
    read=lambda name:json.loads((DEST/name).read_text(encoding='utf8'))
    motion=read('bone_motion_audit.json');cad=read('rotor_cad_roundtrip.json');neutral=read('neutral_checks.json');geometry=read('export_checks.json');joint=read('joint_assembly_checks.json');mount=read('mount_assembly_checks.json');revision=read('canonical_mesh_revision.json');rejected=read('rejected_insertion_relief.json')
    assert cad['pass'] and len(cad['parts'])==2
    assert cad['cad_input_sha256']==hashlib.sha256((DEST/'cad_input.json').read_bytes()).hexdigest()
    for path,digest in motion['mesh_evidence_sha256'].items():assert hashlib.sha256(__import__('pathlib').Path(path).read_bytes()).hexdigest()==digest
    for name,digest in motion['input_evidence_sha256'].items():assert hashlib.sha256((DEST/name).read_bytes()).hexdigest()==digest
    assert len(motion['actions'])==1 and motion['actions'][0]['probe']=='elbow_flexion'
    assert not motion['complete_prescribed_suite']
    assert neutral['neutral_rigid_geometry_clear'] and geometry['all_generated_roundtrips_pass'] and geometry['all_scene_registration_pass']
    assert mount['all_rigid_paths_clear']
    for name,digest in joint['geometry_sha256'].items():assert hashlib.sha256((DEST/name).read_bytes()).hexdigest()==digest
    names=['assembly.json','cad_input.json','selection.json','canonical_mesh_revision.json','rotor_cad_roundtrip.json','bone_motion_audit.json','motion_paths.json','joint_assembly_checks.json','mount_assembly_checks.json','neutral_checks.json','export_checks.json','rejected_insertion_relief.json']
    status={'stage':'M3 integral rotor experiment, not adopted','elbow_continuous_rigid_prefix_deg':150*motion['actions'][0]['continuous_prefix_fraction'],'whole_arm_13_path_verified':False,'two_rigid_CAD_parts_valid':True,'neutral_rigid_geometry_clear':True,'bare_elbow_installation_pass':joint['sequence_clear'],'muscle_mount_rigid_paths_clear':True,'TPU_installation_contact_flags':mount['soft_flags'],'further_straight_slot_rejected_for_disconnection':not rejected['accepted'],'tendon_route_verified':False,'muscle_drive_verified':False,'frozen':False,'evidence_sha256':{name:hashlib.sha256((DEST/name).read_bytes()).hexdigest() for name in names}}
    (DEST/'candidate_status.json').write_text(json.dumps(status,indent=2),encoding='utf8')
    failed=next((p for p in joint['paths'] if not p['clear']),None)
    native=next(p for p in cad['parts'] if p['part']=='forearm_with_rotor')
    text=['# M3 同轴转子实验：未采用','', '此实验独立于M2安装版本。实际构造10mm名义钢丝中心半径、1mm槽宽、4mm总宽的环槽转子，轴向站位+12mm，空心连接颈融合进前臂；肱骨局部让位约756.904mm³，关节中心和叉架耳部不移动。','', '| 检查 | 当前结果 |','|---|---|','| 两件实体CAD | 有效；肱骨为分面实体，前臂由已验证原STEP与同参数24边形转子原生布尔合并 |',f'| CAD/制造网格一致性 | 前臂相对体积差 {native["relative_volume_error"]:.3g}；512个网格顶点至CAD表面最大 {native["sampled_mesh_vertex_to_STEP_surface_max_mm"]:.6f}mm，非全表面证明 |','| 全部中立位100组件 | 显著硬件/骨性背景相交0；TPU接触另记，不是变形验证 |',f'| 肘往返路径 | 0–{status["elbow_continuous_rigid_prefix_deg"]:.3f}°刚性几何连续可达；78个已归属组件，1°常规采样/0.1°边界阈值 |','| 整臂肩带13路径 | 本实验尚未重新执行，M2结果不可直接转移 |',f'| 肌肉近端安装 | 16条刚性路径通过；仍有{mount["soft_flags"]}个TPU装入采样接触标记 |',f'| 一体转子前臂装入 | **未通过**：原直线装入路径首个阻挡发生在偏移{failed["first_limit"]["offset_mm"] if failed else "无"}mm |','| 扩大该直线入口 | **拒绝**：会留下两个正体积骨骼分量，不删除碎片凑通过 |','| 肌腱端部/缠绕/导向/驱动 | 尚未接入，不计算实际肌肉驱动通过 |','', '## CAD修复依据','',f'直接缝合MF布尔网格得到无效BRep。改以原生24边形棱柱构造转子并与原前臂STEP融合，随后把有效CAD重新导出为制造/仿真网格。相对原试件网格顶点至另一表面最大修正约{revision["mesh_vertex_to_other_mesh_surface_max_mm"]:.6f}mm，相对体积变化{revision["relative_volume_change_from_rejected_mesh"]:.3g}。修订后重新执行了零位、导出配准和肘运动检查；未放宽一致性阈值来宣布通过。旧试件保留于before_native_canonical_mesh/。','', '24边形环槽离散后并非数学意义的任意角度旋转不变，因此肘路径使用实际网格旋转复核，不能仅以名义圆环对称性省略检查。','', '## 设计结论','', '当前一体转子不并入主方案。下一步比较真实可拆转子及扭矩连接，或检查其他装入方向；选定后再设计端帽至转子的肌腱、支架和钢丝固定。即便运动范围通过，可装配性也是独立门槛。','', '完整证据哈希见candidate_status.json。M2可查看CAD仍位于上两级cad/mounted_arm_rigid_review.step。']
    (DEST/'ROTOR_EXPERIMENT_REVIEW.md').write_text('\n'.join(text)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in status.items() if k!='evidence_sha256'},indent=2))
if __name__=='__main__':run()
