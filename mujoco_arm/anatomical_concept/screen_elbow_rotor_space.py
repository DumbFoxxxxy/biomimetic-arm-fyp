"""Actual coaxial grooved rotor solids, not ideal-drum constraints.

Screen rotationally symmetric rotor + hub-collar additions against the current
neutral manufactured arm. No tendon route, load or complete drive is accepted.
"""
from build_elbow_muscle_mounts import ROOT,OUT,DEST as MOUNT,load,rod,unit,mf,np
from build_SC_connection import write
import json,hashlib

DEST=MOUNT/'distal_rotor_space'
def run():
    DEST.mkdir(exist_ok=True)
    a=json.loads((MOUNT/'assembly.json').read_text(encoding='utf8'))
    elbow=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
    e=json.loads((elbow/'assembly.json').read_text(encoding='utf8'))
    from build_reference import SHIFT
    E=np.array(e['centre_source_mm'])-SHIFT;axis=unit(np.array(e['axis_source']))
    solids={r['key']:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for r in a['instances']}
    records={r['key']:r for r in a['instances']};boxes={k:np.array(s.bounding_box()) for k,s in solids.items()}
    forearm=solids['FJ3391'];humerus=solids['FJ3368'];rows=[];relieved=[]
    def cyl(low,high,radius):return rod(E+axis*low,E+axis*high,radius)
    for sign in [-1,1]:
        for station in [12,16,20,24,28]:
            centre=sign*station
            for radius in [8,10,12,14]:
                # 4mm wide ring: two1.5mm flanges and1mm wide groove. Steel
                # centre radius is r for nominal0.6mm wire, groove root r-.3.
                lo,hi=centre-2,centre+2
                ring=cyl(lo,lo+1.5,radius+1)+cyl(lo+1.5,hi-1.5,radius-.3)+cyl(hi-1.5,hi,radius+1)
                start,end=sorted([sign*7.6,centre]);collar=cyl(start,end,8)
                added=(ring+collar)-cyl(min(lo,start)-1,max(hi,end)+1,6)
                bb=np.array(added.bounding_box());hard=[];soft=[]
                for key,s in solids.items():
                    if key=='FJ3391':continue # this is the part receiving actual fused material
                    other=boxes[key]
                    if np.any(np.maximum(bb[:3]-other[3:],other[:3]-bb[3:])>=0):continue
                    volume=max(0,float((added^s).volume()))
                    if volume>.001:
                        (soft if records[key]['kind']=='TPU' else hard).append({'part':key,'volume_mm3':volume})
                join_volume=max(0,float((added^forearm).volume()))
                connected=join_volume>.001 and len([c for c in added.decompose() if c.volume()>.001])==1
                name=f'rotor_{"positive" if sign>0 else "negative"}_{station}_r{radius}'
                row={'candidate':name,'axis_station_mm':centre,'nominal_steel_centre_radius_mm':radius,'groove_width_mm':1.,'flange_outer_radius_mm':radius+1,'centre_bore_radius_mm':6,'collar_outer_radius_mm':8,'forearm_fusion_overlap_mm3':join_volume,'addition_connected':connected,'hard_interferences':hard,'soft_contact_review_flags':soft,'neutral_rigid_clear':connected and not hard,'physical_tendon_route_exists':False,'moment_arm_verified':False,'frozen':False}
                # A real bearing/rotor pocket, not a collision exclusion: retain
                # the fork ears at +/-34mm and keep pin/sleeves unchanged.
                cut=cyl(lo-.35,hi+.35,radius+1.35)+cyl(start-.35,end+.35,8.35)
                revised=(humerus-cut).as_original().simplify(.001)
                positive=[c for c in revised.decompose() if c.volume()>.001]
                remainder=max(0,float((revised^added).volume()))
                removed=max(0,float((humerus-revised).volume()))
                relief={'humerus_positive_components':len(positive),'removed_humerus_mm3':removed,'remaining_rotor_humerus_overlap_mm3':remainder,'radial_and_axial_nominal_clearance_mm':.35,'fork_ears_preserved_by_station_extent':max(abs(lo-.35),abs(hi+.35),abs(start-.35),abs(end+.35))<e['fork_half_gap_mm'],'all_other_hard_clear':not any(h['part']!='FJ3368' for h in hard),'full_motion_verified':False}
                relief['candidate_neutral_eligible']=connected and len(positive)==1 and remainder<=.001 and relief['all_other_hard_clear'] and relief['fork_ears_preserved_by_station_extent']
                row['explicit_humeral_rotor_pocket']=relief
                if relief['candidate_neutral_eligible'] and radius>=10:
                    relieved.append((removed,name,revised,(forearm+added).as_original().simplify(.001),row,added))
                rows.append(row)
                if row['neutral_rigid_clear']:
                    write(DEST/(name+'.stl'),added.as_original().simplify(.001))
                print(name,'hard',len(hard),'soft',len(soft),'connected',connected,flush=True)
    selected=None
    if relieved:
        removed,name,hum,fore,row,added=min(relieved,key=lambda r:(r[0],r[1]))
        choice=DEST/'pocket_candidate';choice.mkdir(exist_ok=True)
        # 0.001mm boolean ancestry leaves microfacets rejected by OpenCascade.
        # Use the same0.01mm simplification as the verified parent skeleton.
        hum=hum.as_original().simplify(.01);fore=fore.as_original().simplify(.01)
        parts=[write(choice/'humerus_rotor_pocket.stl',hum),write(choice/'forearm_with_rotor.stl',fore)]
        rotor_folder=choice/'rotor_addition';rotor_folder.mkdir(exist_ok=True)
        rotor=write(rotor_folder/'rotor_addition.stl',added.as_original().simplify(.001))
        (rotor_folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[rotor]}),encoding='utf8')
        assert all(p['components']==1 for p in parts)
        (choice/'cad_input.json').write_text(json.dumps({'units':'mm','parts':parts}),encoding='utf8')
        selected={'candidate':name,'selection_rule':'among eligible radius>=10mm candidates, least removed humeral material; no strength optimization','source_humerus_sha256':hashlib.sha256((ROOT/records['FJ3368']['file']).read_bytes()).hexdigest(),'removed_humerus_mm3':removed,'geometry':row,'export_simplification_tolerance_mm':.01,'parts':[p['name'] for p in parts],'motion_verified':False,'tendon_and_drive_verified':False}
        (choice/'selection.json').write_text(json.dumps(selected,indent=2),encoding='utf8')
    result={'scope':'actual axisymmetric rotor/collar space screen; no physical tendon, anchors, output constraint or actuation proof','source_assembly_sha256':hashlib.sha256((MOUNT/'assembly.json').read_bytes()).hexdigest(),'source_forearm_OBJ_sha256':hashlib.sha256((ROOT/records['FJ3391']['file']).read_bytes()).hexdigest(),'centre_display_mm':E.tolist(),'axis':axis.tolist(),'rigid_threshold_mm3':.001,'relative_transform_note':'axisymmetric addition is invariant under ideal elbow-axis rotation, so its geometry relative to humerus is constant; coupled arm paths and tendon movement still need checking','candidates':rows,'neutral_rigid_clear_count':sum(r['neutral_rigid_clear'] for r in rows),'selected_explicit_pocket_candidate':selected,'drive_pass':False}
    result['relative_transform_note']='Nominal circular rotor is axisymmetric; actual24-sided facets are not invariant under arbitrary rotation. This is a neutral screen only. The selected pocket experiment separately samples the actual rotating geometry.'
    (DEST/'space_screen.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    lines=['# 肘轴同轴肌腱转子空间筛查','','构造了带实际环槽、侧挡边及连接至前臂轮毂的空心轴颈实体。所有相邻硬件均参加检查；仅接收新增材料的前臂本体作真实融合处理。','',f'40 个候选中 {result["neutral_rigid_clear_count"]} 个在当前零位没有显著刚性相交。TPU 接触单独记录。', '', '环槽内钢丝半径为8/10/12/14mm，轴向站位为±12/16/20/24/28mm；槽宽1mm，名义钢丝0.6mm。只有真实穿线、缠绕固定、端帽至转子的导向和力臂确认后，才可用半径计算所需肌肉行程。','', '| 候选 | 首个刚性限制 | TPU 接触数 | 与前臂融合体积/mm³ |','|---|---|---:|---:|']
    for r in rows:
        limit=r['hard_interferences'][0]['part'] if r['hard_interferences'] else '零位未检出'
        lines.append(f'| {r["candidate"]} | {limit} | {len(r["soft_contact_review_flags"])} | {r["forearm_fusion_overlap_mm3"]:.3f} |')
    if selected:lines+=['',f'显式局部间隙槽试件：{selected["candidate"]}，移除肱骨材料 {selected["removed_humerus_mm3"]:.3f} mm³，保留单连通肱骨、叉架两侧耳及原销轴/限位套。新增转子实际融合到前臂；两件输出在 pocket_candidate/。这是待全路径检查的候选，尚未替换 M2 结构。']
    lines+=['','尚未选定最终转子。槽底棱边、弯曲半径、端部固定和扭矩传递强度未验证；不能将此筛查计作肌肉运动通过。']
    (DEST/'ROTOR_SPACE_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    print('Neutral rigid clear candidates',result['neutral_rigid_clear_count'],flush=True)

if __name__=='__main__':run()
