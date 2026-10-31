from build import *
import shutil

def run():
    for path in [OUT/'elbow_connection',OUT/'elbow_connection/muscle_mounts']:
        report=path/'full_arm_audit.json'
        if not report.exists():continue
        data=json.loads(report.read_text(encoding='utf8'));paths=json.loads((path/'full_arm_motion_paths.json').read_text(encoding='utf8'));txt=['# E1整臂几何检查及肌肉安装候选','','本报告检查真实三角网格和E1全部实体附件。SC/AC/GH尚未安装的参考件、','机架占位没有冒充验收对象，整臂工程验收仍为未通过。','','|路径|去程点数|回程点数|整条规定路径是否清晰|首个限制|','|---|---:|---:|---|---|']
        for row in data['actions']:
            seq=np.array([p['q_deg'] for p in paths[row['action']]]);row['tested_joint_min_deg']=seq.min(0).tolist();row['tested_joint_max_deg']=seq.max(0).tolist()
            limit='无检测到的实体穿透' if row['first_limit'] is None else '; '.join('/'.join(h['pair']) for h in row['first_limit']['hits'])
            txt.append(f"|{row['action']}|{row['forward_samples']}|{row['return_samples']}|{'是' if row['whole_path_component_clear'] else '否'}|{limit}|")
        txt+=['','闭合路径的末态是零位，不表示只测试了0度；各关节实际测试最小/最大值','已记录在JSON中。RAISE_THEN_FLEX_ELBOW先前屈96度，再屈肘145度，','随后保持屈肘放下，再伸肘；外展路径同理。抱合路径同时包括前屈25度、','外展8度、旋转−20度及屈肘116度。去程/回程中的碰撞区没有被算作连续可达。','','检查精度：角度步长不大于1度，首次边界细化至0.1度；体积报告阈值0.001mm³。','所有28对被包含零件均检查，包括同一刚体中的零件。间隙搜索上限5mm。','不包含动态受力、软肌肉变形和任意联合姿态的完整证明。']
        if (path/'interfaces.json').exists():
            interface=json.loads((path/'interfaces.json').read_text(encoding='utf8'));txt+=['','## A/B安装接口','','在肱骨上增加两块反力安装板和支撑，在前臂增加两只带3.4mm通孔的肌腱孔耳。','安装板的两只M3孔与TPU端板孔同轴；名义孔径3.4mm，导管入口4.2mm。','安装板端面与TPU端面有2.5mm间隔，后续装配需要相应的隔套/夹持压板。','4.2mm通孔本身不能固定PTFE导管；法兰接头或夹持件仍必须设计与验证，','不能宣称已经具备导管反力保持能力。C保持自由，不添加虚构滑块。',f"零位TPU与新骨端几何检查：{sum(x['pass'] for x in interface['neutral_TPU_checks'])}/{len(interface['neutral_TPU_checks'])}通过。",'本轮只是零位软硬件几何检查，运动中的肌肉变形、接触和端部固定仍待验证。']
        dump(report,data);(path/'CONNECTED_REVIEW.md').write_text('\n'.join(txt)+'\n',encoding='utf8')
    print('connected reports updated')

if __name__=='__main__':run()
