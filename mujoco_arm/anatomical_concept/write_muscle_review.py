"""Evidence-bound elbow module review, explicitly separate from installed drive."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles'
def read(p):return json.loads(p.read_text(encoding='utf8'))

def run():
    excursion=read(DEST/'excursion_screen.json');placement=read(DEST/'neutral_placement_screen.json');rows=[];checks={};hashes={}
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;design=read(folder/'design.json');geometry=read(folder/'geometry_checks.json');assembly=read(folder/'end_assembly_checks.json');cad=read(folder/'hard_caps/cad_roundtrip.json')
        source=read(folder/'cad_input.json')['parts'];cached=read(folder/'hard_caps/cad_input.json')['parts'];caps=[p for p in source if p['name'].endswith('rigid_cap')]
        prints=read(folder/'print_candidates/print_transforms.json')['parts']
        result={'20_parts_single_material_each':len(design['components'])==20 and set(design['components'].values())=={1} and geometry['single_connected_material_each'],'OBJ_STL':geometry['all_roundtrips_pass'],'module_neutral_clear':not design['neutral_module_overlaps'],'nominal_wire_threading_geometry':geometry['nominal_steel_threading_geometry_pass'],'nominal_channel_clear':design['wire_channel_TPU_overlap_mm3']<1e-6,'11_end_fitting_paths':len(assembly['paths'])==11 and assembly['all_paths_collision_free'],'rigid_caps_STEP':len(cad['parts'])==2 and all(p['pass'] for p in cad['parts']),'rigid_cap_CAD_input_current':caps==cached,'nominal_guide_wall':min(p['nominal_tube_wall_lower_bound_mm'] for p in design['actual_short_tube_guides'])>=.6-1e-6}
        if not all(result.values()):raise RuntimeError({fid:result})
        result['three_print_files_current']=len(prints)==3 and all(p['print_STL_roundtrip_pass'] and p['source_OBJ_sha256']==hashlib.sha256((folder/(p['part']+'.obj')).read_bytes()).hexdigest() and p['print_STL_sha256']==hashlib.sha256((folder/'print_candidates'/(p['part']+'.stl')).read_bytes()).hexdigest() for p in prints)
        if not result['three_print_files_current']:raise RuntimeError('Stale print candidate '+fid)
        checks[fid]=result
        for name in ['design.json','cad_input.json','geometry_checks.json','end_assembly_checks.json','hard_caps/cad_roundtrip.json','print_candidates/print_transforms.json']:
            hashes[fid+'/'+name]=hashlib.sha256((folder/name).read_bytes()).hexdigest()
        rows.append((fid,design,geometry,assembly))
    status={'scope':'two unmounted manufacturing geometry candidates; not whole-arm completion','checks':checks,'evidence_sha256':hashes,'frozen':False,'installed_drive_pass':False,'not_verified':['bone-side adapters and grounded reaction path','passive tendon guides/drum and full-range moment arms','nonlinear deformation and soft contact','printing, crimp grip, preload, fatigue and loads','shoulder/scapular muscle manufacture and whole-arm actuation']}
    (DEST/'verification_status.json').write_text(json.dumps(status,indent=2),encoding='utf8')
    lines=['# 肘部异形TPU肌肉制造候选','',
      '已生成肱肌参考FJ1486和肱三头肌外侧头参考FJ1477两块候选。它们包含实际三维晶格、参数化外表面菱形网、多个短管状导向单元、穿线孔、TPU端垫、刚性端板和金属附件包络。每块20个组件，TPU肌腹自身为单连通实体；不是长方体占位或纯渲染外壳。当前模型尚未安装到骨架，未验证收缩、驱动或制造性能。', '',
      '## 与论文方法的对应','',
      '依据本地论文《Accelerating Musculoskeletal Robotics Through Parametric Design and 3D‐Printed Flexible Structures》第2.3、4.2节及图4，并对照其PSM前作图4：保留肌腹体积→晶格/外网格→管状导向→端部的组成思路。参考论文p=2mm及15/10/20mm尺寸建立候选；本脚本局部轴依次为横向15、非主要受载方向20、工程收缩方向10mm。', '',
      '截面来自真实OBJ沿肱骨工程方向的三角面截交线，17个截面拟合周期三次B样条，轴向也用三次B样条；相当于权重为1的曲面参数化。内网格采用规则八分支体单元后裁切，外网格沿参数面生成。尚未逐式复现论文式(3)的体积参数坐标映射，不能直接继承论文的应变、刚度或耐久数值。', '',
      '钢丝路径通过固定端点、截面内候选点的短步长路径确定，段间在真实导向管处折转。该路径是工程设计，未经解剖肌纤维标定；不采用PCA代替肌纤维。管轴沿相邻直线跨度的角平分线，保留名义管壁。压接后的实际摩擦和受压变形仍需试验。', '',
      '| 项目 | FJ1486 | FJ1477 |','|---|---:|---:|']
    metrics=[('设计肌腹轴向跨度mm',lambda d,g:d['body_span_mm']),('TPU材料体积mm³',lambda d,g:d['TPU_material_volume_mm3']),('肌腹包络内材料比例',lambda d,g:d['material_inside_loft_fraction']),('拟合包络位于原肌肉外的体积比例',lambda d,g:d['loft_outside_original_anatomy_mm3']/d['loft_volume_mm3']),('最大导向折转角°',lambda d,g:max(p['incoming_outgoing_turn_deg'] for p in d['actual_short_tube_guides'])),('导向管名义壁厚下界mm',lambda d,g:min(p['nominal_tube_wall_lower_bound_mm'] for p in d['actual_short_tube_guides'])),('0.6mm钢丝名义最小间隙mm',lambda d,g:g['nominal_steel_min_clearance_mm_capped_2'])]
    for label,fn in metrics:lines.append('| '+label+' | '+' | '.join(f'{fn(d,g):.4f}' for fid,d,g,a in rows)+' |')
    lines+=['', '拟合前裁去源肌肉轴向两端各10%，用工程端垫和过渡取代。上表“位于原肌肉外的比例”只衡量拟合包络的外溢，不是总体积保留率或完整表面误差证明；局部端垫加宽、表面杆件和连接加强材料另有设计变化。', '',
      '## 端部与传力','',
      '主动功来自电机与钢丝。近端端板固定后提供导管反力接口；主动钢丝经过短导向管，在远端金属压接止挡处终止，拉近两端并压缩晶格。远端端板的横向孔连接被动短钢丝肌腱，再向关节输出。TPU提供可变形体积、弹性和缓冲，当前没有假定其能自行主动收缩。', '',
      '每块模块：1件TPU肌腹、2件刚性端板、4个名义M3×12螺栓、4个M3螺母、8片垫圈、1个压接止挡，共20件。端垫直径22mm、厚4mm，M3孔距14mm；TPU孔道切削半径0.85mm，检查的保留目标半径0.8mm；0.6mm钢丝仅为首轮候选规格。刚性端板含导管孔/台阶、止挡座及被动肌腱横孔。金属止挡必须实际压接后才可能夹住钢丝，几何包络的通孔不能证明夹持力。', '',
      '新检查曾发现螺母终态可放置但不能装入；已保留iterations/before_nut_access，并在TPU端部开侧向装入槽。现已检查两端板、四组螺母/内垫圈、四螺栓及止挡共11条1mm采样路径。螺纹啮合、工具可达性、预紧及TPU端部蠕变未验证。', '',
      '局部分离的生成枝条以记录端点的真实短连接修复，导向孔和螺母装入槽保留为空腔。只对生成过程产生的微小封闭孔腔和松散碎片做记录的局部修整，未删除源解剖模型。各制造STL/OBJ回读均通过；四个刚性端板STEP回读通过。TPU仍以制造网格和可复现生成参数交付。', '',
      '## 接入骨架仍需完成','',
      '解剖候选端点直线牵引仅作为对照筛查，不代表上述弯曲肌腹的实际驱动路线，也未给不存在的滑轮施加约束：', '',
      '| 参考肌肉 | 直线端点总行程mm | 首次力臂符号不符° | 与骨架相交的姿态数/151 |', '|---|---:|---:|---:|']
    for r in excursion['muscles']:lines.append(f"| {r['id']} | {r['direct_route_geometric_excursion_mm']:.3f} | {r['first_wrong_or_zero_torque_sign_deg'] if r['first_wrong_or_zero_torque_sign_deg'] is not None else '此探针未出现'} | {r['colliding_samples']} |")
    lines+=['', '直拉路径不能冻结为最终走线。理想半径8/10/12/14mm卷绕的150°行程分别约20.94/26.18/31.42/36.65mm，只是后续尺寸对照；没有据此生成“已验证的滑轮驱动”结论，也不把行程比例当成TPU可用压缩率。', '',
      '| 模块 | 原位置刚性干涉对数 | 沿候选外移方向首次采样刚性清空mm |', '|---|---:|---:|']
    for r in placement['muscles']:lines.append(f"| {r['id']} | {len(r['neutral_hard_interferences'])} | {r['first_sampled_rigid_clear_offset_mm']} |")
    lines+=['', '外移结果仅用于设计骨侧支架起点，尚未采用为装配位置；检查只覆盖零位。TPU与骨接触列为变形/摩擦待查，既不按刚性穿透一律判失败，也未当作已通过。近端支架、导管/电机反力固定、远端实际肌腱导向和关节范围内的力臂/行程仍待实现。肩部与肩胛肌群的制造、驱动和真实机架仍属于整臂目标。', '',
      '## 文件和复现','',
      'FJ1486/、FJ1477/内的TPU_muscle.stl为制造几何；proximal_rigid_cap.stl和distal_rigid_cap.stl为刚性端板。其余金属STL是硬件包络，不作为TPU打印件。hard_caps/cad/含端板STEP。print_candidates/为另外旋转/落床的候选文件，变换有记录；这些文件只用于打印准备，装配仍使用父目录公共坐标。切片支撑、材料/孔配合和打印质量未验证。', '',
      'MuJoCo模块查看器读取与制造STL同源的OBJ，以支持较多三角面；静态钢丝仅示意真实导向之间的名义路径。图片不包含已验证的收缩变形或驱动动画。', '',
      '```powershell', "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/parameterize_elbow_muscles.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/build_shaped_elbow_muscles.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_shaped_muscles.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_muscle_end_assembly.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/screen_elbow_excursion.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_muscle_placement.py", "& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_muscle_caps.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/prepare_muscle_prints.py", "& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/write_muscle_review.py", "& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_shaped_muscles.py --render", "& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_shaped_muscles.py --id FJ1486", '```', '',
      '先完成骨侧适配、实际关节走线及打印样件检查，再决定肌腹尺寸和可用范围。测试记录模板中的实测项保持null，未从论文或数字形状填入“实测”结果。']
    (DEST/'MUSCLE_MODULE_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    template={'scope':'real bench measurements only; all entries below are unmeasured','modules':[]}
    for fid,d,g,a in rows:template['modules'].append({'id':fid,'printer':None,'nozzle_mm':None,'TPU_grade_and_batch':None,'print_orientation':None,'active_wire_diameter_mm_measured':None,'reference_cap_spacing_mm_measured':None,'tests':[],'sample_test_schema':{'cycle':None,'commanded_wire_draw_mm':None,'measured_cap_spacing_mm':None,'wire_input_force_N':None,'distal_output_force_N':None,'terminal_slip_mm':None,'recovery_after_release_mm':None,'guide_or_end_damage_notes':None}})
    target=DEST/'bench_measurements_template.json'
    if not target.exists():target.write_text(json.dumps(template,indent=2),encoding='utf8')
    print(checks)

if __name__=='__main__':run()
