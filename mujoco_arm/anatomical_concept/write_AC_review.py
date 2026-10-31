"""Compare actual AC candidates; repair obsolete audit labels, never pass results."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection'

def read(path):return json.loads(path.read_text(encoding='utf8'))

def run():
    rows={};checks={}
    definition=read(BASE.parent/'definition.json')
    for label,folder in [('large','AC_connection'),('compact','AC_compact')]:
        dest=BASE/folder;report=read(dest/'bone_motion_audit.json');assembly=read(dest/'assembly.json')
        owners=dict(definition['bone_owner'])
        for name,owner in assembly['part_owners'].items():
            owners[{'clavicle_SC_AC_balls':'FJ3362','scapula_AC_socket':'FJ3384'}.get(name,name)]=owner
        report.update(scope='Independent prescribed-path solid audit; five moving bones plus ten AC accessories against fixed torso bones and SC accessories; no moving pairs excluded.',moving_bone_count=5,moving_component_count=len(owners),moving_owners=owners,included_physical_connections=['SC','AC'],not_included=['GH/elbow retention hardware','muscle deformation','actuation','manufacturing tolerance and loads'])
        (dest/'bone_motion_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
        rows[label]={r['probe']:r for r in report['actions']}
        cad=read(dest/'cad_roundtrip.json')['parts'];scene=read(dest/'audit_checks.json')
        checks[label]={'neutral_clear':not report['neutral']['hits'],'mesh':read(dest/'geometry_checks.json')['all_geometry_checks_pass'],'CAD':len(cad)==len(assembly['components']) and all(p['pass'] for p in cad),'compound_STEP':read(dest/'assembly_STEP_check.json')['pass'],'scene':len(scene['zero_pose_geometry'])==26 and all(r['pass'] for r in scene['zero_pose_geometry']),'sampling':scene['maximum_sampled_joint_increment_deg']<1+1e-6,'capture':read(dest/'capture_probes.json')['all_six_straight_directions_blocked_within_2mm']}
        if not all(checks[label].values()):raise RuntimeError('Candidate checks incomplete '+label)
    lines=['# 肩锁实体连接候选对照','',
      '两方案均保留原拟合AC中心，不移动GH中心。大方案球半径10mm、球窝外半径13mm；紧凑方案分别6mm、9mm。两者保留完整球头、骨侧连续过渡、球窝、分体保持盖、M3螺栓/螺母/垫圈，均为12个连通实体。紧凑方案螺栓名义M3×12，大方案M3×16；螺纹未建模。',
      '', '当前后续设计优先使用AC_compact，尚未冻结。它恢复了大方案损失的盂肱抬升路径，代价是孤立AC上回旋探针的范围更小。不是在全部动作上优于大方案，也尚未验证缩小球头的强度。', '',
      '## 运动采样（数字几何）','',
      '下表为从零位连续可达的指定路径前缀；角度为输入旋转向量的模，联合抬升列为实际肱骨方向相对零位的夹角。不是人体临床ROM或已承诺的新验收标准。', '',
      '| 探针 | 大方案° | 紧凑方案° | 紧凑方案首个限制 |', '|---|---:|---:|---|']
    for name,a in rows['large'].items():
        b=rows['compact'][name]
        def angle(r):
            if name=='coordinated_raise':return r['actual_humerus_elevation_deg']
            c=np.array(r['last_clear_controls_deg']);return max(abs(c[9]),*(np.linalg.norm(c[i:i+3]) for i in [0,3,6]))
        limit='；'.join(' / '.join(h['pair']) for h in b['first_limit']) or '请求路径内未检出'
        lines.append(f'| {name} | {angle(a):.3f} | {angle(b):.3f} | {limit} |')
    compact=BASE/'AC_compact';impact=read(compact/'geometry_checks.json')['old_anatomical_patch_impact'];affected=[r for r in impact if r['points_more_than_1mm_from_new_surface']]
    lines+=['','## 肌肉接口影响','', '候选骨面区域不是已冻结锚点。以下区域有采样点距工程修改表面超过1mm，必须重映射；不能直接复用原解剖附着面：', '', '| 肌肉编号 | 骨面 | 超过1mm点数 | 最大偏差mm |', '|---|---|---:|---:|']
    for r in affected:lines.append(f"| {r['muscle']} | {r['structure']} | {r['points_more_than_1mm_from_new_surface']} | {r['maximum_deviation_mm']:.3f} |")
    lines+=['','## 验证边界','',
      '全部活动零件依其所属刚体运动：肩锁球窝、保持盖和附件随肩胛运动，没有误固定到世界坐标。15个活动实体与42个固定实体检查，邻接骨没有整对排除。角度采样≤1°、边界细化≤0.1°，实体交叠阈值0.001mm³，最小间隙查询封顶3mm。7个肋软骨没有作为刚性实体验收；柔性接触仍未知。回程重新计算，演示仅播放首次碰撞之前的连续前缀。', '',
      'OBJ/STL回读、各12个独立STEP实体和复合STEP分别检查。STEP为三角面BRep，可导入CAD编辑，但不含原生特征树和装配配合。六方向球体直拉均被固定球窝/保持盖阻挡：仅是有限方向几何检查，不能证明所有逃脱路径或抗拉承载。', '',
      'SC/AC实体候选已加入模型；GH与肘仍是理想约束，不能据此宣称完整工程骨架已完成。肌肉在运动中隐藏，未模拟TPU、钢丝收缩或力传递。实际机架、装配工具可达性、打印公差、紧固预紧、防脱强度及寿命仍待验证。', '',
      '## 复现入口','', '从项目根目录运行（Python环境见已有SC说明）：', '', '```powershell', "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/build_AC_connection.py --compact", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/audit_joint_model.py --AC --compact", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_SC_assembly.py --AC --compact", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_AC_capture.py --compact", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_cad.py --AC --compact", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_assembly.py --AC --compact", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_joint_audit.py --AC --compact", "& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_reference.py --AC --compact --supports", '```','', '去掉--compact可复现大方案。视窗按1–5播放已检查路径，M切换肌肉、T切换胸廓参考、R回到零位；运动肌肉隐藏避免误示驱动。']
    (BASE/'AC_CANDIDATE_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf8');print(checks);print('review',BASE/'AC_CANDIDATE_REVIEW.md')

if __name__=='__main__':run()
