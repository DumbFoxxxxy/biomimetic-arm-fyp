"""Update active entries after adopting checked S1-Y2; preserve task history."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'

def update(p,pairs):
    s=p.read_text(encoding='utf8')
    for old,new in pairs:
        if old not in s:raise ValueError(f'Missing expected text: {p}: {old[:60]}')
        s=s.replace(old,new,1)
    p.write_text(s,encoding='utf8')

addition='''### 2.16 当前S1-Y2：封闭过渡空隙的明确填实

[S1-Y2结构与复现](s1_y2_transition_stock/README.md)为当前唯一几何研究起点。Y1保留为父级；只换冈下肌、肩胛下肌两件TPU，其余189源及变换、骨骼／钢丝／锚点不变。对已确认位于R3.3保留过渡区且不属于完整通道或端接口的四处闭腔显式填实，分别增加约13.11811与2.79255 mm³材料，约原肌腹0.156%与0.036%。无外部材料设计切除，但局部刚度／折叠影响未知，不能据体积比例推断性能不变。

这是分类后的设计变更，逐空腔位置／体积、通道与接口排除记录在[比较生成清单](s1_y_muscle_pair/transition_stock_comparison/transition_trials.json)。没有用导出清理自动填全部负壳，也不允许把名义通道填成实体。原开口搜索后伸过长的问题已查明；限制后伸能开敞肩胛下肌两处，但保留最小闭腔。与开口＋剩余填实方案比较后，本阶段选择保留外部轮廓的全闭腔填实，尚未确认力学优劣。

[独立单肌肉复核](s1_y_muscle_pair/transition_stock_comparison/independent_checks.json)重新检查覆盖、有限接触、原保护材料、177刚体及六主动线；[新191对象报告](s1_y2_transition_stock/assembly_check.json)重新计算18145对零位与六主动线，无超原0.001 mm³刚性／主动线相交。这两块修订肌腹无封闭负壳，不代表其他对象无闭腔。

[新裸骨路径](s1_y2_transition_stock/bare_motion_review/bare_joint_motion_checks.json)及[保持／放置](s1_y2_transition_stock/readback_keeper_checks.json)针对新源清单重新执行，裸骨范围、11步与有限捕获同2.14节；仍不作带肌肉或TPU驱动范围。当前选择、源哈希、钢丝引用和空白试验记录已更新，旧记录均保留。`stage_pass=false`，实际驱动范围未知。

'''
update(D/'SHOULDER_DESIGN.md',[
 ('版本：v0.9.12','版本：v0.9.13'),
 ('当前研究起点为S1-Y1：','当前研究起点为S1-Y2：'),
 ('两块肌腹更新见2.15节。','两块肌腹形态见2.15节，当前过渡区修订见2.16节。'),
 ('[当前S1-Y1研究起点](s1_y_muscle_pair/README.md)、[源清单](s1_y_muscle_pair/assembly_check.json)、[裸骨范围](s1_y_muscle_pair/bare_motion_review/bare_joint_motion_checks.json)', '[当前S1-Y2研究起点](s1_y2_transition_stock/README.md)、[源清单](s1_y2_transition_stock/assembly_check.json)、[裸骨范围](s1_y2_transition_stock/bare_motion_review/bare_joint_motion_checks.json)'),
 ('### 2.15 当前S1-Y1：两块共形Y晶格肌腹','### 2.15 S1-Y1父级：两块共形Y晶格肌腹'),
 ('当前唯一选择见[候选指针]','本节保留Y1父级的几何记录；当前Y2见2.16节，唯一选择见[候选指针]'),
 ('## 3. 肌肉／近端支承父级已验证内容与继承边界',addition+'## 3. 肌肉／近端支承父级已验证内容与继承边界'),
 ('当前选择见2.15节','当前选择见2.16节'),
 ('当前S1-Y1已替换','S1-Y1父级已替换'),
 ('| S1-Y1实际驱动范围 |','| S1-Y2实际驱动范围 |'),
])
with (D/'SHOULDER_DESIGN.md').open('a',encoding='utf8') as f:f.write('\n2026-11-02 v0.9.13：当前研究起点切换S1-Y2，记录4处非通道闭腔的明确填实、有限开口搜索修正与完整零位／裸骨／保持复核。外形／锚点不移，刚度影响及驱动／制造未验收。[更新前快照](../../../reviews/S1_2026-11-02/Y2_transition_stock/before_SHA256.json)。\n')

update(ROOT/'PROJECT_STATUS.md',[
 ('技术文件v0.9.12','技术文件v0.9.13'),
 ('[S1-Y1：紧凑GH骨架＋两块共形Y肌腹](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y_muscle_pair/README.md)', '[S1-Y2：紧凑GH骨架＋过渡区修订Y肌腹](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/README.md)'),
 ('[肩部技术文件v0.9.12]','[肩部技术文件v0.9.13]'),
 ('当前S1-Y1只替换','Y1父级只替换'),
 ('## 最新进展与限制\n','## 最新进展与限制\n\n当前S1-Y2在保留两肌腹外形的前提下，显式填实4处已排除通道／接口功能的封闭过渡空隙；分别增料约13.11811／2.79255 mm³。其余189源／变换不变。18145对零位、六主动线、裸骨及保持路径重新检查；局部刚度变化未测，打印／驱动未验收。\n'),
 ('新增试制检查：两块真实肌腹','历史Y1试制检查：两块真实肌腹'),
 ('当前S1-Y1指针不变，尚待实际切片选择及接口／工艺验证。','上述为历史比较；随后限制刀具后伸已开敞两处。当前选S1-Y2分类填实版，避免依赖细开口，实际工艺与力学仍待验证。'),
 ('先核对S1-Y1骨架','先核对S1-Y2骨架'),
])
update(ROOT/'design/CURRENT_BASELINE.md',[
 ('S1-Y1：紧凑GH骨架＋冈下肌／肩胛下肌共形Y肌腹','S1-Y2：紧凑GH骨架＋两块共形Y肌腹／封闭过渡空隙填实'),
 ('技术文件2.14–2.15节','技术文件2.14–2.16节'),
 ('S1-Y1中冈下肌／肩胛下肌已采用','S1-Y2沿用Y1的冈下肌／肩胛下肌Y格，两件肌腹已采用'),
 ('详细参数与证据保留在技术文件','Y2仅对已排除通道／端接口功能的4处闭腔显式填实，外形不变，局部刚度影响未测。详细参数与证据保留在技术文件'),
])
brief=ROOT/'design/DESIGN_BRIEF.md';raw=brief.read_bytes();marker='## 历史记录阅读说明'.encode('utf8');head,tail=raw.split(marker,1)
s=head.decode('utf8').replace('S1-Y1：','S1-Y2：').replace('技术文件v0.9.12','技术文件v0.9.13').replace('[S1-Y1说明]','[S1-Y2说明]').replace('s1_y_muscle_pair/README.md','s1_y2_transition_stock/README.md')
s+='本轮只对两块肌腹4处已分类的非通道／非接口闭腔显式增料填实，保留外部轮廓及其他189源；原Y1及失败搜索保留。当前191清单已重做零位、裸骨与保持检查，局部刚度／制造与驱动待验。\n\n'
brief.write_bytes(s.encode('utf8')+marker+tail)

for path in [ROOT/'README.md',ROOT/'02_人体骨骼晶格肌肉方案/README.md',ROOT/'RUNNING.md']:
    s=path.read_text(encoding='utf8').replace('S1-Y1：','S1-Y2：').replace('技术文件2.14–2.15节','技术文件2.14–2.16节')
    if path.name=='RUNNING.md':s=s.replace('s1_y_muscle_pair/README.md','s1_y2_transition_stock/README.md')
    if path.parent.name=='02_人体骨骼晶格肌肉方案':s=s.replace('v0.9.12','v0.9.13').replace('S1-Y1共形肌腹','S1-Y2过渡区修订肌腹').replace('及历史比较、历史比较及验证记录','及历史比较／验证记录')
    path.write_text(s,encoding='utf8')
update(D/'s1_y_muscle_pair/README.md',[
 ('当前研究候选，**未制造定版，未验证实际驱动范围**。','Y1历史研究父级，**未制造定版，未验证实际驱动范围**。当前默认已切换[S1-Y2过渡区修订](../s1_y2_transition_stock/README.md)。'),
 ('SHOULDER_DESIGN v0.9.12','SHOULDER_DESIGN v0.9.13'),
 ('第一条显示当前191对象实际零位网格','第一条显示候选指针选择的当前191对象实际零位网格；本目录Y1需显式指定`--report`'),
])
print('Technical v0.9.13 active entries updated, historical brief bytes preserved.')
