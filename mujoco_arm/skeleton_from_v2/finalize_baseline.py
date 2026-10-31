"""Freeze a bounded geometry baseline only after the explicit evidence gates."""
import json,hashlib,shutil,csv,xml.etree.ElementTree as ET
import numpy as np
from verify_joint_candidate import ROOT,OUT,mujoco

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run():
    report=json.loads((OUT/'retained_shape_rom.json').read_text());paths=json.loads((OUT/'retained_shape_paths.json').read_text())
    cad=json.loads((OUT/'cad_roundtrip.json').read_text());vertices=json.loads((OUT/'cad_vertex_correspondence.json').read_text());baseline=json.loads((OUT/'v2_bone_geometry_baseline.json').read_text())
    fingerprint=report['mesh_cache_sha256'];assert digest(OUT/report['mesh_cache'])==fingerprint
    assert cad['input_cache_sha256']==vertices['input_cache_sha256']==fingerprint and cad['pass'] and cad['assembly_valid'] and vertices['pass']
    assert len(report['actions'])==len(paths)==12 and all(a['reverse_collision_count']==0 for a in report['actions'])
    source=ROOT/'mujoco_arm/cable_design_v2/arm_cable_design_v2.xml';assert digest(source)==baseline['source_xml_sha256']
    model=mujoco.MjModel.from_xml_path(str(OUT/'skeleton_review.xml'))
    violations=[]
    for name,seq in paths.items():
        for q in seq:
            angles=np.radians(q)
            for j in range(model.njnt):
                value=angles[model.jnt_qposadr[j]]
                if model.jnt_limited[j] and not model.jnt_range[j,0]-1e-8<=value<=model.jnt_range[j,1]+1e-8:violations.append([name,j,float(value)])
    assert not violations
    dest=OUT/'V2-S1';(dest/'assets').mkdir(exist_ok=True,parents=True);(dest/'cad').mkdir(exist_ok=True)
    tree=ET.parse(OUT/'skeleton_review.xml');root=tree.getroot();root.set('model','V2_S1_bone_geometry_baseline');root.find('compiler').set('meshdir','.')
    for geom in root.findall('.//geom'):
        name=geom.get('name');mesh=geom.get('mesh')
        if mesh:
            sourceasset=OUT/(name+'_review_local.stl');target=dest/'assets'/(name+'.stl');shutil.copy2(sourceasset,target)
            root.find(f'.//asset/mesh[@name="{mesh}"]').set('file','assets/'+target.name)
    ET.indent(tree,space='  ');tree.write(dest/'skeleton.xml',encoding='utf-8')
    portable=mujoco.MjModel.from_xml_path(str(dest/'skeleton.xml'))
    for attr in ['body_pos','body_quat','jnt_pos','jnt_axis','jnt_range','mesh_vert','mesh_face']:
        assert np.array_equal(getattr(model,attr),getattr(portable,attr)),attr
    for part in (OUT/'cad').glob('*.step'):shutil.copy2(part,dest/'cad'/part.name)
    for src,target in [('retained_shape_rom.json','rom_audit.json'),('retained_shape_paths.json','motion_paths.json'),('cad_roundtrip.json','cad_roundtrip.json'),('cad_vertex_correspondence.json','cad_geometry_comparison.json'),('skeleton_motion_review.gif','motion.gif')]:shutil.copy2(OUT/src,dest/target)
    if (OUT/'before_after.png').exists():shutil.copy2(OUT/'before_after.png',dest/'before_after.png')
    labels={'K1':'肘屈曲','K2':'GH前屈','K3':'GH外展','K4':'GH内旋','K5':'GH外旋','K6':'肩胛上回旋','K7':'联合抬升','SLIDE_FWD':'肩胛前滑','SLIDE_BACK':'肩胛后滑','SCAP_EXT':'肩胛外旋','SCAP_INT':'肩胛内旋','PECTORAL':'抱合组合姿态'}
    reasons={'K2':'锁骨—肱骨首次碰撞约113.75°','K6':'已验证路径至19°；20°的局部路径搜索未满足间隙，不是已证明的结构极限','K7':'沿已验证肩带路径达到146.02°，未覆盖170°；不是全姿态最大值','SLIDE_BACK':'锁骨—胸廓干涉','SCAP_INT':'肩胛—胸廓干涉；继续轮廓切除超过保留形态的调整预算','PECTORAL':'仅验证V2既定组合姿态，不等同于解剖抱合角度验收'}
    table=[]
    for a in report['actions']:
        value=a['actual_slide_mm'] if a.get('actual_slide_mm') is not None else a['actual_range_deg'];unit='mm' if a.get('actual_slide_mm') is not None else '°'
        actual=f'{value:.2f}{unit}' if value is not None else '25°前屈 / 8°外展 / −20°轴转'
        if a['action'] in ['K6','K7']:actual='已验证至 '+actual
        row=[a['action'],labels[a['action']],a['target'],a['floor'],actual,reasons.get(a['action'],'完整目标路径通过')];table.append(row)
        for gap in a['end_pair_gaps']:
            key='/'.join(gap['bones']);a['minimum_sampled_pair_gaps_mm_capped_at_2'][key]=min(a['minimum_sampled_pair_gaps_mm_capped_at_2'].get(key,2),gap['gap_mm_capped_at_2'])
    report['frozen']=True;report['editable_STEP_roundtrip_verified']=True;report['native_solidworks_feature_tree_updated']=False
    (dest/'rom_audit.json').write_text(json.dumps(report,indent=2))
    with (dest/'acceptance.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['id','动作','原目标','原底线','验证范围','限制或说明']);w.writerows(table)
    normal=np.cross(np.array([44.55,112.87,106.69])-np.array([72.76,268.45,82.21]),np.array([159.06,187.54,7.3])-np.array([72.76,268.45,82.21]));normal/=np.linalg.norm(normal)
    if normal[2]<0:normal=-normal
    delta={'fixed_joint_centers_design_mm':{'SC':[40,235,-106],'AC':[156.86,233.24,-11.97],'GH':[159.06,187.54,7.3],'E':[184.34,-101.35,7.3]},'joint_center_displacements_mm':0,'SC_AC_length_mm':150,'GH_E_length_mm':290,'main_plate_standoff_vector_design_mm':(normal*8).tolist(),'SUP_INF_P':'nominal outline/reference controls displaced by stand-off away from protected joint regions; P remains an old reference marker, not an asserted surface contact','AC_ball_nominal_radius_mm':10,'final_SC_slot':json.loads((OUT/'sc_mount_repair.json').read_text()),'final_reverse_relief':json.loads((OUT/'reverse_clearance_repair.json').read_text()),'anchors_frozen':False}
    (dest/'geometry_delta.json').write_text(json.dumps(delta,indent=2))
    lines=['# V2-S1 骨架几何基线（2026-10-29）','','按用户选择，以保留形态、局部修正和记录实际范围为定版标准。本版锁定五个骨骼的几何和已验证路径，保留 V2 的 SC2 / AC1 / GH3 / 肘1 机构链；没有采用 GH 外移80mm或同心三轴。','','## 动作结果','','|编号|动作|原目标|原底线|已验证范围|限制/说明|','|---|---|---|---|---|---|']
    lines.extend('|'+ '|'.join('—' if v is None else str(v) for v in row)+'|' for row in table)
    lines.extend(['','K6、K7 的数值是已找到并验证的连续路径范围，不宣称为全局最大活动度。K6 仍低于原验收底线；K2、K7 未覆盖完整目标。原指标保留。','','## 几何变化','','- SC、AC、GH、肘心不移动，锁骨和肱骨的关节中心距离保留。肩胛主体沿法向渐进让位8mm，AC/GH附近保留连接过渡。SUP、INF及原P所在主体区域随轮廓调整；P不能继续当成已满足的贴胸接触点。','- AC球头增至名义R10。肩关节颈部、肘部及肩胛局部轮廓按实际运动包络避让；保留骨骼空腔。没有删除整对关节碰撞检查。','- 胸锁支座开局部让位槽，保护轴承座区域。五个骨骼均保持单一材料实体；内腔不算断开的骨骼。','- 肩胛材料净变化约−11.2%，肱骨约−7.0%，其余变化见 geometry_delta.json。主体让位不是整块肩胛带着关节中心平移。','','## 验证和交付','','- 12条动作路径，常规相邻姿态各关节变化不大于1°，碰撞边界细化至0.1°；前后滑边界采用0.02°细化，对应标记点变化小于0.1mm。验证了已接收路径的去程和回程，关节限位无越界。','- 五骨之间逐对实体布尔检查，重叠体积数值阈值0.001mm³；零位及接收路径不保留实体穿透。各对最小采样间隙见 rom_audit.json，超过2mm的结果按2mm封顶。没有作连续全姿态证明。','- cad/ 包含五个零位零件STEP及整体STEP。它们是可直接编辑的分面BRep实体，不包含原SolidWorks参数特征历史。STEP回读合法性、实体数量、体积、包围盒、顶点对应及重新网格化后的差集体积检查均通过。','- skeleton.xml 与经检查的网格一致；motion.gif 是规定姿态的骨骼运动演示，不是肌肉驱动仿真。原V2与原始SLDPRT/SLDASM保留。','','## 本次定版的边界','','仅定骨架几何基线。关节采用V2的理想运动学约束；未验证SC交叉销等辅助硬件的完整实体装配、防脱、受力、强度、制造公差或实物性能。胸廓沿用原实心几何包络，任务书的3mm壳体尚未落实。GH凹面经包络避让，不宣称已按原R10.2/深3mm解析球窝完成制造设计。','钢丝锚点与TPU肌肉未定版，旧V2 sites只保留为后续设计种子，部分点位可能已不落在修正骨面上。原生SolidWorks特征树未回写，本版可编辑CAD交付为STEP。'])
    text='\n'.join(lines)+'\n';(dest/'README.md').write_text(text,encoding='utf-8');(ROOT/'design/V2_SKELETON_FREEZE.md').write_text(text+'\n交付目录：`mujoco_arm/exports/skeleton_from_v2/V2-S1/`。\n',encoding='utf-8')
    brief=ROOT/'design/DESIGN_BRIEF.md';content=brief.read_text(encoding='utf-8');marker='<!-- V2-S1 BASELINE START -->';end='<!-- V2-S1 BASELINE END -->'
    if marker in content:content=content.split(marker)[0]+content.split(end,1)[1].lstrip('\n')
    heading,rest=content.split('\n',1)
    update=f'''\n{marker}\n## 0A 最新骨架几何基线 V2-S1（2026-10-29）\n\n本节优先于下文历史版本速查及相互冲突的旧表述。已按“保留形态、记录实际范围”的标准锁定骨架几何，详见 `design/V2_SKELETON_FREEZE.md`；交付目录为 `mujoco_arm/exports/skeleton_from_v2/V2-S1/`。\n\n- 保留V2的SC2 / AC1 / GH3 / 肘1链；SC、AC、GH和肘心坐标不变。肩胛主体渐进让位8mm，SUP/INF/P主体控制区域相应调整；原P仅作参考标记，不再宣称零位贴胸。\n- 当前验证范围：肘145°、GH前屈113.69°、外展120°、内旋70°、外旋90°；肩胛上回旋验证至19°、联合抬升验证至146.02°。K6低于原底线；原K1–K7目标和底线不降低。\n- AC球头名义R10；肩、肘与SC支座有局部实体避让。精确轮廓以本版STEP/网格为准，旧GH解析球窝尺寸未作为本版制造验收通过项。\n- 本版为五骨几何与规定动作路径的基线。胸廓3mm壳体、辅助关节硬件、防脱、实物公差及承载未完成；理想关节假设沿用。钢丝+TPU肌肉及其锚点属于下一阶段。\n- STEP为已回读核验的可编辑分面BRep；原生SolidWorks参数特征树未回写，原V2和原始CAD保留。\n{end}\n\n'''
    brief.write_text(heading+'\n'+update+rest.lstrip('\n'),encoding='utf-8')
    status={'version':'V2-S1','geometry_baseline_frozen':True,'freeze_standard':'preserve form and document verified paths; user permits below-floor ranges','all_original_targets_met':False,'all_original_floors_met':False,'editable_STEP_roundtrip_verified':True,'native_solidworks_feature_tree_updated':False,'hardware_ready':False,'TPU_drive_verified':False,'wire_anchors_frozen':False,'source_V2_preserved':True,'joint_limit_violations':violations,'accepted_forward_pose_count':sum(len(q) for q in paths.values()),'accepted_reverse_pose_count':sum(len(q) for q in paths.values()),'input_mesh_cache_sha256':fingerprint,'files_sha256':{str(p.relative_to(dest)):digest(p) for p in dest.rglob('*') if p.is_file() and p.name!='freeze_manifest.json'}}
    (dest/'freeze_manifest.json').write_text(json.dumps(status,indent=2));print(json.dumps({k:v for k,v in status.items() if k!='files_sha256'}),flush=True)
if __name__=='__main__':run()
