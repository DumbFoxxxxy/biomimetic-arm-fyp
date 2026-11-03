"""Synchronize technical entries with the selected S1-Y1 evidence, preserving history."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / '02_人体骨骼晶格肌肉方案/输出/shoulder_design'
Y = D / 's1_y_muscle_pair'

def edit(path, replacements):
    text = path.read_text(encoding='utf-8')
    for old, new in replacements:
        if old not in text:
            raise ValueError(f'Missing expected text in {path}: {old[:80]}')
        text = text.replace(old, new, 1)
    path.write_text(text, encoding='utf-8')

section = '''### 2.15 当前S1-Y1：两块共形Y晶格肌腹

当前唯一选择见[候选指针](../../../design/S1_CURRENT_CANDIDATE.json)，结构与复现见[S1-Y1说明](s1_y_muscle_pair/README.md)。以2.14节的紧凑191对象为固定几何父级，只替换冈下肌与肩胛下肌两件TPU；其余189个源文件及装配变换不变，未移动骨骼、GH中心、锚点或钢丝路径。三角肌与大圆肌仍保留旧晶格。

沿内部通道中心线建立局部Y切向，横截面薄向为X，Z由右手系得到。把共享节点的单元映射到实际肌腹截面，再形成中心至八角的内部杆、边界菱形表面网格、连续通道与端部过渡。内部杆名义直径2 mm、表面杆1.5 mm；单元尺寸是生成密度输入，不是变形后处处相同的晶胞尺寸。

| 参数／证据 | 冈下肌FJ1500 | 肩胛下肌FJ1504 |
|---|---:|---:|
| 名义X/Y/Z尺寸 | 15/10/30 mm | 6/8/10 mm |
| 单元／轴向分层／表面面片 | 10／10／40 | 28／7／70 |
| 实际轴向站距 | 9.20577 mm | 7.73130 mm |
| 实际横向X投影间距范围 | 6.56500–18.23850 mm | 2.64230–11.67730 mm |
| 实际横向Z投影间距范围 | 13.02568–51.15164 mm | 1.13960–14.23512 mm |
| TPU体积 | 8412.55559 mm³ | 7726.39466 mm³ |

[独立回读报告](s1_y_muscle_pair/assembly_check.json)核对191对象、18145对零位配对和六主动线：无超0.001 mm³判据的刚性／主动线相交，各源均为单正材料连通、零边界及非流形边。每个轴向分层均有实际材料支撑的内部单元和外表面杆；这验证几何存在，不验证收缩或保护性能。四个近端M3孔路清，两端盖有限法向探针遇材料；端垫／钩／通道保护区的导出差值分别约0.03212／0.00701 mm³，不称接口严格零变化。

连续通道名义孔径1.72 mm、钢丝0.6 mm。除显式端接口外，0.75 mm壁厚包络的原始缺失量分别0.01958／0.00928 mm³；同宽包络径向外移0.002 mm以避开共面数值残差后，分别约0.00000165／0.00000024 mm³。两种结果均保留；0.002 mm不是打印公差或壁厚放宽。未验证弯折、摩擦、变形后的通畅与壁厚。

骨骼未改，但针对Y1清单重新执行了[裸骨路径](s1_y_muscle_pair/bare_motion_review/bare_joint_motion_checks.json)与[保持／11步放置](s1_y_muscle_pair/readback_keeper_checks.json)，范围与2.14节一致。裸骨仍排除89对象，不能作为带肌肉运动通过。完整刚体扫掠因自由末端缺少耦合约束而被[明确拒绝](s1_y_muscle_pair/rigid_constraint_guard_check.json)。

固定间距Y格即使网格闭合，也曾仅覆盖4/10与5/7轴向分层；该方案已拒绝为完整肌腹，记录保留在审查目录。当前共形方案解决了所查分层的材料覆盖，尚未解决打印工艺、实际变形或驱动力验证。局部截面压缩及端部融合使论文均匀晶胞性能不能直接转移。`stage_pass=false`，预测与实测驱动范围均未知。

'''

edit(D / 'SHOULDER_DESIGN.md', [
 ('版本：v0.9.11', '版本：v0.9.12'),
 ('当前研究起点为191对象紧凑埋入式GH保持件版', '当前研究起点为S1-Y1：紧凑GH骨架＋两块共形Y晶格肌腹，共191对象'),
 ('详见2.14节。', '骨架见2.14节，两块肌腹更新见2.15节。'),
 ('- [当前191对象研究起点](s1_GH_assembly_review/compact_recessed_keepers/README.md)、[源清单](s1_GH_assembly_review/compact_recessed_keepers/assembly_check.json)、[裸骨范围](s1_GH_assembly_review/compact_recessed_keepers/bare_motion_review/bare_joint_motion_checks.json)。', '- [当前S1-Y1研究起点](s1_y_muscle_pair/README.md)、[源清单](s1_y_muscle_pair/assembly_check.json)、[裸骨范围](s1_y_muscle_pair/bare_motion_review/bare_joint_motion_checks.json)。'),
 ('当前查看默认已指187对象', '历史187查看已另存，默认整肩查看以候选指针为准'),
 ('### 2.14 当前紧凑191对象研究起点及裸骨范围', '### 2.14 S1-Y1沿用的紧凑191对象骨架父级及裸骨范围'),
 ('[当前紧凑版](s1_GH_assembly_review/compact_recessed_keepers/README.md)', '[紧凑骨架父级](s1_GH_assembly_review/compact_recessed_keepers/README.md)'),
 ('## 3. 肌肉／近端支承父级已验证内容与继承边界', section + '## 3. 肌肉／近端支承父级已验证内容与继承边界'),
 ('当前选择见2.14节', '当前选择见2.15节'),
 ('最新装配仍使用原肌腹轴向体系：一般单元10/8/10 mm、三角肌12/10/12 mm，杆半径1 mm，局部Z沿中心线。生成器已增加独立Y沿中心线模式，默认Z在受检案例中与原脚本几何对称差为0；**增加Y模式没有自动替换最新装配的肌腹**，不能宣称已复制论文收缩行为。', '当前S1-Y1已替换冈下肌、肩胛下肌为局部Y沿中心线的共形晶格，参数与实际变尺度见2.15节。三角肌、大圆肌仍采用原Z轴体系；三角肌名义12/10/12 mm、其余原一般单元10/8/10 mm、杆半径1 mm。生成器保留旧Z默认，既有方向回归两项通过；Y方向一致及完整几何覆盖不代表已复制论文收缩行为。'),
 ('191对象研究起点实际驱动范围', 'S1-Y1实际驱动范围'),
 ('| 1a | 修正GH轴向捕获，协调保持盖与肩胛下肌空间 | 实际肱骨与诊断球闭盖受阻、开盖可装；开槽／五金路径及变化配对复验；新盖柔性冲突有处理依据 |', '| 1a | 核对实际GH紧固件、工具、预紧与座壁／横销承载 | 已有有限捕获及名义放置；下一步提交真实选件／工具进入、预紧与承载依据，不能重复列为尚无几何方案 |'),
 ('| 2 | 完成通道壁检查并确定可试验的工艺与终端 | 壁厚／孔路记录、实际切片与试样、压接／卡扣／抱箍试验；未测项保持未知 |', '| 2 | 对共形Y肌腹确定工艺、切片、通道与终端 | 已有零位壁厚包络／孔路记录；需实际切片、样件、穿线及压接／卡扣／抱箍试验，局部细格不能只凭STL闭合判可打印 |'),
])
with (D / 'SHOULDER_DESIGN.md').open('a', encoding='utf-8') as f:
    f.write('\n2026-11-02 v0.9.12：采用S1-Y1两块共形Y肌腹并复核源清单；其余189源不变。补实际晶格尺寸、覆盖、通道与接口证据，保留稀疏格失败及全部原指标。驱动／打印未验收。[更新前快照](../../../reviews/S1_2026-11-02/Y_muscle_pair/before_SHA256.json)。\n')

edit(ROOT / 'PROJECT_STATUS.md', [
 ('技术文件v0.9.11', '技术文件v0.9.12'),
 ('[191对象紧凑埋入式GH保持件](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_GH_assembly_review/compact_recessed_keepers/README.md)', '[S1-Y1：紧凑GH骨架＋两块共形Y肌腹](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y_muscle_pair/README.md)'),
 ('[肩部技术文件v0.9.11]', '[肩部技术文件v0.9.12]'),
 ('## 最新进展与限制\n', '## 最新进展与限制\n\n当前S1-Y1只替换冈下肌、肩胛下肌两件TPU，其余189源与变换不变。两肌腹为Y沿中心线的变截面内部晶格＋菱形表面网格＋连续通道，分别覆盖10/10、7/7轴向分层；三角肌、大圆肌保留旧格。名义晶胞输入不是处处一致的实际尺寸，不能转移论文材料性能。新清单已重做零位、裸骨和保持检查；固定间距稀疏格失败保留。\n'),
 ('| 肌肉与钢丝 | 四模块六主动线布局；两肩袖肌自由末端、肌腱及近端支承 |', '| 肌肉与钢丝 | 四模块六主动线；两块共形Y肌腹完整分层几何、连续通道及自由末端接口 |'),
 ('先核对191对象紧凑方案', '先核对S1-Y1骨架紧固与两块Y肌腹的切片／穿线可行性；核对191对象紧凑方案'),
])

edit(ROOT / 'design/CURRENT_BASELINE.md', [
 ('191对象紧凑埋入式GH保持件；后续几何与肌肉集成使用此版，不是制造冻结', 'S1-Y1：紧凑GH骨架＋冈下肌／肩胛下肌共形Y肌腹，共191对象；不是制造冻结'),
 ('技术文件2.14节', '技术文件2.14–2.15节'),
 ('现有装配肌腹仍为局部Z轴沿中心线；生成器Y模式及方向样件未自动替换整肌腹。', 'S1-Y1中冈下肌／肩胛下肌已采用Y沿中心线的共形晶格，含连续通道和菱形表面网格；三角肌／大圆肌仍为旧Z轴体系。名义尺寸仅为密度输入，实际晶胞变尺度，未标定的刚度与收缩率不能从论文或方向样件直接转移。'),
])

brief = ROOT / 'design/DESIGN_BRIEF.md'
old = brief.read_bytes()
marker = '## 历史记录阅读说明'.encode('utf-8')
head, history = old.split(marker, 1)
s = head.decode('utf-8').replace('191对象紧凑埋入式GH保持件版', 'S1-Y1：紧凑GH骨架＋两块共形Y肌腹，共191对象').replace('技术文件v0.9.11', '技术文件v0.9.12')
s += '冈下肌、肩胛下肌已替换为Y方向共形晶格并补菱形表面网格，其他189源不变。三角肌、大圆肌仍保留旧格；局部晶胞尺寸变化，打印及收缩／张力待实测。当前几何与复现见[S1-Y1说明](../02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y_muscle_pair/README.md)。\n\n'
brief.write_bytes(s.encode('utf-8') + marker + history)

for path in [ROOT / 'README.md', ROOT / '02_人体骨骼晶格肌肉方案/README.md']:
    edit(path, [
        ('191对象紧凑埋入式GH保持件版', 'S1-Y1：紧凑GH骨架＋两块共形Y肌腹，共191对象'),
        ('详见技术文件2.14节。', '两块肩袖肌已采用共形Y格与菱形表面网格，其余189源不变；实际变形与打印待验。详见技术文件2.14–2.15节。'),
    ])
scheme = ROOT / '02_人体骨骼晶格肌肉方案/README.md'
s = scheme.read_text(encoding='utf-8').replace('v0.9.11', 'v0.9.12').replace('187对象集成候选及GH试改', 'S1-Y1共形肌腹、紧凑GH骨架及历史比较')
scheme.write_text(s, encoding='utf-8')

edit(ROOT / 'RUNNING.md', [
 ('当前选择为191对象紧凑埋入式GH保持件', '当前选择为S1-Y1：紧凑GH骨架＋两块共形Y肌腹，共191对象'),
 ('s1_GH_assembly_review/compact_recessed_keepers/README.md', 's1_y_muscle_pair/README.md'),
])
edit(D / 's1_GH_assembly_review/compact_recessed_keepers/README.md', [
 ('# 当前S1研究起点：191对象紧凑埋入式保持件', '# S1-Y1的骨架父级：191对象紧凑埋入式保持件'),
 ('更新：2026-11-02。已选为后续几何与肌肉集成的研究起点，', '更新：2026-11-02。本目录保留骨架父级及旧肌腹；当前默认为[S1-Y1两块共形Y肌腹版](../../s1_y_muscle_pair/README.md)，骨骼沿用本目录，'),
])

print('Updated technical document v0.9.12 and seven entry documents; historical brief bytes preserved.')
