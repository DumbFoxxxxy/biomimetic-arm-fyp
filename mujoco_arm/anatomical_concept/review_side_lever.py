"""M4 evidence audit: distinguish installable rigid lever from failed tendons."""
from build_elbow_side_lever import *

def run():
    folder=DEST/'axis_negative_r15';read=lambda name:json.loads((folder/name).read_text(encoding='utf8'))
    design=read('design.json');assembly=read('assembly.json');fit=read('assembly_checks.json');motion=read('bone_motion_audit.json');paths=read('motion_paths.json');cad=read('cad_roundtrip.json');combined=read('side_lever_arm_STEP_check.json');exports=read('export_checks.json');neutral=read('neutral_checks.json');diagnostic=read('tendon_axial_diagnostic.json');contact=read('tendon_family_contacts.json')
    assert len(assembly['instances'])==111 and len(design['parts'])==12
    assert fit['all_rigid_paths_clear'] and len(fit['paths'])==8
    for name,digest in fit['source_mesh_sha256'].items():assert hashlib.sha256((folder/(name+'.obj')).read_bytes()).hexdigest()==digest
    assert motion['complete_prescribed_suite'] and len(motion['actions'])==13
    for name,digest in motion['input_evidence_sha256'].items():assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest
    for path,digest in motion['mesh_evidence_sha256'].items():assert hashlib.sha256(__import__('pathlib').Path(path).read_bytes()).hexdigest()==digest
    max_step=0.;norm_error=0.
    for name,path in paths.items():
        q=np.array([p['qpos'] for p in path]);assert q.shape[1]==13
        for start in [0,4,8]:
            values=q[:,start:start+4];norm_error=max(norm_error,float(abs(np.linalg.norm(values,axis=1)-1).max()));dot=np.sum(values[:-1]*values[1:],axis=1);max_step=max(max_step,float(np.degrees(2*np.arccos(np.clip(abs(dot),0,1))).max()))
        max_step=max(max_step,float(abs(np.degrees(np.diff(q[:,-1]))).max()))
        action=next(r for r in motion['actions'] if r['probe']==name);end=np.array(action['requested_controls_deg']);index=np.flatnonzero(abs(end)>1e-8)[0]
        assert action['forward_samples']==action['return_samples']==len(path)
        for p in path:
            if p['controls_deg'][index]/end[index]<=action['continuous_prefix_fraction']+1e-8:assert not p['hits']
    assert max_step<=1.000001 and norm_error<1e-8 and not motion['neutral']['hits']
    inputs={p['name']:p for p in read('cad_input.json')['parts']};assert len(cad['parts'])==12
    for p in cad['parts']:assert p['pass'] and abs(p['input_volume_mm3']-inputs[p['part']]['volume_mm3'])<1e-5
    assert combined['pass'] and combined['solid_count']==109
    assert combined['assembly_sha256']==hashlib.sha256((folder/'assembly.json').read_bytes()).hexdigest()
    assert combined['STEP_sha256']==hashlib.sha256((ROOT/combined['STEP']).read_bytes()).hexdigest()
    assert all(r['STEP_vs_OBJ_vertex_error_mm']<.001 for r in combined['parts'].values())
    assert exports['all_generated_roundtrips_pass'] and exports['all_scene_registration_pass'] and len(exports['scenes']['neutral_view.xml'])==111 and len(exports['scenes']['arm.xml'])==89
    assert neutral['neutral_rigid_geometry_clear'] and neutral['manufactured_instance_count']==111
    for name,digest in read('geometry_checks.json')['source_mesh_sha256'].items():assert hashlib.sha256((folder/(name+'.obj')).read_bytes()).hexdigest()==digest
    assert diagnostic['source_lever_design_sha256']==hashlib.sha256((folder/'design.json').read_bytes()).hexdigest()
    # Clarify inherited generic report labels, without changing any geometric
    # result: uncrimped stoppers were prescribed packing envelopes, not fixed.
    scope=read('rigid_motion_scope.json');scope.update(scope='M4 prescribed rigid packing:87 nominal bone/bolted pieces plus2 uncrimped passive stop positions. No grounding/drive proof for stoppers.',prescribed_uncrimped_stops=['flexor_tendon_stop','extensor_tendon_stop'],nominally_bolted_or_bone_instances=87)
    (folder/'rigid_motion_scope.json').write_text(json.dumps(scope,indent=2),encoding='utf8')
    motion.update(scope='M4 prescribed all-pair rigid geometry:87 nominal bone/bolted pieces plus2 uncrimped stop envelopes; not whole-module motion or drive proof',rigid_motion_scope=scope)
    (folder/'bone_motion_audit.json').write_text(json.dumps(motion,ensure_ascii=False,indent=2),encoding='utf8')
    neutral['scope']='all111 manufactured neutral instances; hard/soft classification only, no drive validation'
    (folder/'neutral_checks.json').write_text(json.dumps(neutral,indent=2),encoding='utf8')
    desired={r['muscle']:r for r in diagnostic['modules'] if (r['muscle'],r['lever_eyelet']) in [('FJ1486','flexor'),('FJ1477','extensor')]};assert all(r['complete_sign_consistent_axial_family'] for r in desired.values())
    assert contact['source_design_sha256']==hashlib.sha256((folder/'design.json').read_bytes()).hexdigest()
    first_cap=next((r for r in contact['samples'] if r['cap_rigid_interferences']),None)
    evidence=['design.json','assembly.json','geometry_checks.json','assembly_checks.json','cad_roundtrip.json','side_lever_arm_STEP_check.json','bone_motion_audit.json','motion_paths.json','neutral_checks.json','export_checks.json','tendon_axial_diagnostic.json','tendon_family_contacts.json']
    status={'stage':'M4 detachable rigid lever candidate; direct tendon family rejected','manufactured_neutral_instances':111,'neutral_rigid_CAD_solids':109,'prescribed_rigid_scene_instances':89,'nominally_bolted_or_bone_instances':87,'prescribed_uncrimped_stops':['flexor_tendon_stop','extensor_tendon_stop'],'unvalidated_TPU_and_free_distal_instances':22,'new_connected_components':12,'finite_rigid_installation_paths_passed':8,'installation_TPU_contact_flags':fit['soft_flags'],'maximum_actual_joint_sample_step_deg':max_step,'quaternion_norm_error':norm_error,'full13_rigid_path_suite_completed':True,'direct_tendon_family_verified':False,'muscle_actuation_verified':False,'frozen':False,'evidence_sha256':{name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in evidence}}
    (folder/'review_status.json').write_text(json.dumps(status,indent=2),encoding='utf8')
    lines=['# M4 可拆肘部钢丝力臂候选','', '以解剖骨形和当前独立工程关节为基础，在尺骨42/58mm名义截面设置两处可拆M3安装位。关节中心不移动，尺桡骨仍作为同一前臂刚体。零位按公共坐标对齐。此交付解决刚性力臂及其安装，并没有完成钢丝/TPU驱动。','', '## 结构与装配','',f'- 首轮参考候选：15mm名义力臂、轴向站位−56mm；屈肌/伸肌两孔相位分别为{design["eyelet_phase_diagnostic"]["FJ1486"]["phase_deg"]}°/{design["eyelet_phase_diagnostic"]["FJ1477"]["phase_deg"]}°。这是未验力学的中间半径候选，不是完成优化的唯一参数。',f'- 12个新实体：修订前臂、可拆力臂、两套M3螺栓/螺母/双垫圈、两件尚未压接的肌腱止挡。每件正材料分量为1。',f'- 前臂实际移除{design["forearm_removed_mm3"]:.3f}mm³，含穿孔、局部端面座及力臂轮廓让位；原始骨形、M2和失败候选均保留。没有整体删除骨段或关闭相邻碰撞对。','- 先安装前臂/轴套，再安装可拆力臂、螺栓/垫圈/螺母，最后放入止挡；8条名义刚性装入路径通过，末5mm细化0.1mm，其余≤1mm。此检查不证明完整整臂装配、工具空间或打印/螺纹性能。','- 111件中立位展示，109件刚性STEP；两块TPU另为网格。89个刚性网格参加规定路径，其中87件有名义骨段/螺栓归属，2件未压接止挡仅按力臂携带的位置做包装检查，不代表已固定。22件TPU/自由远端组件未参加骨段刚体通过结论。','', '## CAD 与运动复核','', '- `cad/side_lever_arm_review.step`：109个独立刚性实体；逐件STEP与装配OBJ顶点双向核对<0.001mm，整体回读有效。分面BRep，不含原生特征树或装配配合。',f'- 实际关节采样最大{max_step:.8f}°，边界细化阈值0.1°；显著刚性交叠阈值0.001mm³。最小间隙包含合法名义接触，不能据此声称正制造间隙。','', '| 路径 | 请求控制模长/° | 连续可达控制量/° | 实际肱骨方向变化/° | 首个限制 |','|---|---:|---:|---:|---|']
    for r in motion['actions']:
        end=np.array(r['requested_controls_deg']);requested=max(*(np.linalg.norm(end[i:i+3]) for i in [0,3,6]),abs(end[-1]));limits='; '.join(' / '.join(h['pair']) for h in r['first_limit']) or '采样内未检出'
        lines.append(f'| {r["probe"]} | {requested:.3f} | {requested*r["continuous_prefix_fraction"]:.3f} | {r["actual_humerus_elevation_deg"]:.3f} | {limits} |')
    lines+=['', '联合抬升按肱骨实际方向测量，不把肩带角度简单相加。规定的完整诊断路径可以包含碰撞后的姿态，只有从零位连续无碰撞的前缀计入可达范围。','', '## 钢丝与TPU：尚未通过的部分','', '孔位相位搜索只采用一个明确的假设：远端端帽沿原肌腹轴向移动，且所画肌腱自由段长度不变。没有在MuJoCo中加入不存在的导轨/导向约束。端帽出口参考点的弯折缺少实际导向/终端解，不能作为真实受拉绳索的已验证路径。','', '| 肌肉 | 假设轴向行程/mm | 所画钢丝与硬件相交采样数/151 | 首个钢丝限制 |','|---|---:|---:|---|']
    for c in contact['muscles']:
        d=desired[c['muscle']];first=c['first_rigid_penetration'];limit='无' if first is None else str(first['angle_deg'])+'°: '+', '.join(h['part'] for h in first['hits'])
        lines.append(f'| {c["muscle"]} | {d["axial_excursion_mm"]:.3f} | {c["rigid_penetration_samples"]} | {limit} |')
    lines += ['',f'预测端帽姿态在{contact["cap_interference_samples"]}/151个角度与刚性结构相交；首个出现于{first_cap["angle_deg"] if first_cap else "无"}°。这些交叠来自假设轨迹，说明此轨迹不能直接采用，并不宣称实际柔性系统最多只能到这个角度。', '', '实际TPU应变、横向运动、接触/摩擦、恢复、拉力、钢丝弯曲半径、压接和防松均未验证。零位无相交、相位方向一致和规定骨架路径通过，都不代表TPU可以驱动到150°。','', '## 下一步','', '保留此可拆力臂作为接口候选，设计真实的端帽肌腱连接和有实体支撑的导向；重新规划端帽运动与肌肉布置以避免刚性穿透。随后再估算所需肌肉行程，打印单肌肉并测量可重复行程。肩部/肩胛肌群、机架和整臂驱动仍未完成。','', '## 复现','', '在项目根目录运行：','', '```powershell',"& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_side_lever.py --detail", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/review_side_lever.py",'```','', '完整生成顺序见rebuild_side_lever.py。当前未冻结；证据哈希见review_status.json。']
    (folder/'SIDE_LEVER_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in status.items() if k!='evidence_sha256'},indent=2))
if __name__=='__main__':run()
