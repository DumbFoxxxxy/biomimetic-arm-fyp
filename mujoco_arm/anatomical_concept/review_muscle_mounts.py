"""Consolidate M2 evidence without converting rigid checks into muscle claims."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate'

def run():
    read=lambda name:json.loads((DEST/name).read_text(encoding='utf8'))
    audit=read('bone_motion_audit.json');paths=read('motion_paths.json')
    assembly=read('assembly.json');neutral=read('neutral_checks.json')
    fit=read('mount_assembly_checks.json');exports=read('export_checks.json')
    cad=read('cad_roundtrip.json');cad_input=read('cad_input.json')
    expected={'elbow_flexion','GH_flexion','GH_abduction','GH_axial_positive','GH_axial_negative','SC_elevation','SC_protraction','SC_retraction','SC_rotation_X','AC_upward','AC_rotation_X','AC_rotation_Z','coordinated_raise'}
    assert audit['complete_prescribed_suite'] and set(paths)==expected
    assert {a['probe'] for a in audit['actions']}==expected
    for name,digest in audit['input_evidence_sha256'].items():
        assert hashlib.sha256((DEST/name).read_bytes()).hexdigest()==digest,(name,'stale motion evidence')
    max_step=0.;norm_error=0.
    for name,path in paths.items():
        q=np.array([p['qpos'] for p in path]);assert q.shape[1]==13
        assert np.allclose(q[0],[1,0,0,0]*3+[0])
        for offset in (0,4,8):
            part=q[:,offset:offset+4]
            norm_error=max(norm_error,float(abs(np.linalg.norm(part,axis=1)-1).max()))
            dot=np.sum(part[:-1]*part[1:],axis=1)
            max_step=max(max_step,float(np.degrees(2*np.arccos(np.clip(abs(dot),0,1))).max()))
        max_step=max(max_step,float(abs(np.degrees(np.diff(q[:,-1]))).max()))
        action=next(a for a in audit['actions'] if a['probe']==name)
        assert action['forward_samples']==action['return_samples']==len(path)
        safe=action['continuous_prefix_fraction']
        end=np.array(action['requested_controls_deg']);nonzero=np.flatnonzero(abs(end)>1e-8)
        for p in path:
            fraction=np.array(p['controls_deg'])[nonzero[0]]/end[nonzero[0]]
            if fraction<=safe+1e-8:assert not p['hits'],(name,'collision inside accepted prefix')
    assert max_step<=1.000001 and norm_error<1e-8
    assert not audit['neutral']['hits'] and neutral['neutral_rigid_geometry_clear']
    assert exports['all_generated_roundtrips_pass'] and exports['all_scene_registration_pass']
    assert len(exports['parts'])==len(assembly['generated_components'])==27
    assert len(exports['scenes']['neutral_view.xml'])==100 and len(exports['scenes']['arm.xml'])==78
    assert fit['all_rigid_paths_clear'] and len(fit['paths'])==16
    inputs={p['name']:p for p in cad_input['parts']}
    assert len(cad['parts'])==len(inputs)==27
    for p in cad['parts']:
        assert p['pass'] and p['valid']
        source=inputs[p['part']]
        # STEP exporter records the source used for readback; reject stale reports.
        assert abs(p['input_volume_mm3']-source['volume_mm3'])<1e-5
        assert p['input_vertex_count']==len(source['vertices'])
    combined=read('mounted_arm_STEP_check.json')
    assert combined['pass'] and combined['solid_count']==98 and combined['BRep_valid']
    assert combined['assembly_sha256']==hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest()
    assert combined['STEP_sha256']==hashlib.sha256((ROOT/combined['STEP']).read_bytes()).hexdigest()
    for r in combined['parts'].values():
        assert r['STEP_vs_OBJ_vertex_error_mm']<.001
        assert r['OBJ_sha256']==hashlib.sha256((ROOT/r['source_OBJ']).read_bytes()).hexdigest()
        assert r['STEP_sha256']==hashlib.sha256((ROOT/r['source_STEP']).read_bytes()).hexdigest()
    evidence=['assembly.json','cad_input.json','cad_roundtrip.json','mounted_arm_STEP_check.json','arm.xml','neutral_view.xml','bone_motion_audit.json','motion_paths.json','export_checks.json','neutral_checks.json','mount_assembly_checks.json']
    result={'stage':'M2 proximal mounting candidate','manufactured_instances':100,'motion_included_rigid_instances':78,'motion_excluded_TPU_and_free_distal_instances':22,'maximum_actual_joint_sample_step_deg':max_step,'maximum_quaternion_norm_error':norm_error,'rigid_insertion_paths':16,'soft_insertion_contact_flags':fit['soft_flags'],'complete_assembly_verified':False,'whole_module_motion_verified':False,'drive_verified':False,'physical_test_verified':False,'frozen':False,'evidence_sha256':{name:hashlib.sha256((DEST/name).read_bytes()).hexdigest() for name in evidence}}
    (DEST/'review_status.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    lines=['# M2 肘部肌肉近端安装候选','', '保留真实骨形与原关节中心，新增肱骨安装座、两件可拆转接支架及紧固件。此阶段只解决近端固定与反力接口，远端传力尚未连接。','', '## 已核对的几何证据','',f'- 100 个中立位组件；27 个新增实体的 OBJ/STL 与 STEP 回读通过。',f'- 已归属骨段的 78 个刚性组件参加 13 条往返路径检查；22 个 TPU 或自由远端组件不计入刚体运动通过。',f'- 实际四元数/肘轴采样最大步长 {max_step:.8f}°；边界细化阈值 0.1°，交叠阈值 0.001 mm³，间隙查询上限 3 mm。',f'- 16 条刚性安装路径无干涉；仍有 {fit["soft_flags"]} 条采样接触标记涉及 TPU，完整装配未验证。','', '## 路径结果','', '| 路径 | 请求控制量模长/° | 从零位连续可达控制量/° | 实际肱骨抬升/° | 首个限制对 |','|---|---:|---:|---:|---|']
    for a in audit['actions']:
        end=np.array(a['requested_controls_deg']);requested=max(*(np.linalg.norm(end[i:i+3]) for i in (0,3,6)),abs(end[-1]))
        limits='; '.join(' / '.join(h['pair']) for h in a['first_limit']) or '采样路径内未检出'
        lines.append(f'| {a["probe"]} | {requested:.3f} | {requested*a["continuous_prefix_fraction"]:.3f} | {a["actual_humerus_elevation_deg"]:.3f} | {limits} |')
    lines+=['','联合路径的控制量模长不是实际抬升；实际抬升单独根据肱骨方向计算。完整请求路径含碰撞后的诊断姿态，只有从零位连续无碰撞前缀可计入可达范围。','', '## 当前 CAD 与复现','', '- `cad/mounted_arm_rigid_review.step`：98个独立刚性实体，含尚未接入传动的自由远端附件；两块TPU未转成STEP。所有实体分别按装配OBJ顶点与回读STEP顶点双向匹配，误差小于0.001mm；整体实体数和总体积复核通过。','- 两块TPU在各自模块目录，装配平移见assembly.json；打印摆放候选另有坐标，不能直接混入装配。','- 在项目根目录执行 `C:/Python314/python.exe mujoco_arm/anatomical_concept/view_elbow_muscle_mounts.py` 查看中立位100组件。加 `--render` 重建正反视图。','- 执行 `D:/MujocoEnv/Scripts/python.exe mujoco_arm/anatomical_concept/review_muscle_mounts.py` 检查当前证据一致性并重建本报告。','', '## 设计修订与限制','', '- 屈肌/伸肌整体向前/后平移 16/14 mm，安装站位写在 assembly.json；未移动关节中心。','- U 形支架开 6.6 mm 侧装槽，保留两枚共用 M3 螺钉；骨侧穿栓在台面预装肌肉后安装。','- 伸肌法兰局部削边并采用低头螺钉包络，使正向 GH 轴向范围由约 63.19° 提升至当前记录值；未通过更改原请求角度制造通过。','- 法兰最窄名义垫圈边余量仅 0.1 mm，未做强度、公差或打印验收。','- 最小间隙可能为零；上述无显著交叠并不保证制造间隙。','- STEP 为可导入的分面实体，不含原生特征树、螺纹或装配配合。','- 未模拟 TPU 变形、主动收缩、钢丝摩擦、压接、疲劳；禁止把中立位整体展示当成肌肉驱动演示。','', '## 下一步','', '设计真实的远端肌腱输出、肘侧导向/转子与反力支架，计算连续行程和力臂；随后打印单肌肉测量可重复行程。肩部肌群及整体固定支架仍属于总目标。','', '证据摘要及文件哈希见 review_status.json。']
    (DEST/'MOUNT_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='evidence_sha256'},indent=2))

if __name__=='__main__':run()
