"""Record complete prescribed-eye strand witnesses without changing S1-Y2."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'
def edit(p,pairs):
    s=p.read_text(encoding='utf8')
    for a,b in pairs:
        if a not in s:raise ValueError(f'Missing expected text in {p}: {a[:65]}')
        s=s.replace(a,b,1)
    p.write_text(s,encoding='utf8')

section='''### 2.18 规定端环角的完整肌腱路由候选

[中央段候选与完整网格](s1_y2_transition_stock/free_strand_routes/README.md)补上2.17节局部转动时移除的中央线段。端盖、肱骨、TPU与夹具仍为零位，只规定两端眼环角；中央材料长度明确取未测的名义39.92927／36.64097 mm。保持原约1.6 mm出口直出段，调整匹配两端切向的曲线柄，使整根折返肌腱长度保持在0.0001 mm数值离散条件内；不添加物理导向或假定端盖随肱骨刚性运动。

每块肌肉的单端、同号双端和反号双端各做正／负搜索，共16条1°／0.1°候选路径，合格采样独立重建回程。新整根肌腱与两个同时放置的端环／套筒一起检查，实际钢丝中央段在场。首限区分假设长度的必要条件、当前曲线族限制、自接触和刚性配合冲突；当前曲线族失败不能当作所有路由不可行。表格及完整记录见分支说明，不计作GH或实际驱动范围。

[独立验证](s1_y2_transition_stock/free_strand_routes/verification_checks.json)用32768段密采样核对18例曲线弧长／曲率和出口切向，4个完整191对象上下文诊断场景重新回读配对。候选不是全局最短路由或张紧平衡；没有弯曲刚度／外载／接触等作用时，弯曲自由线不能直接作为正张力理想绳的平衡解。安装长度、压接、预紧及自由端盖—TPU平衡仍未知。

新增源几何发现：肩胛下肌两出口沿销轴错位约10.55503 mm，出口连线偏出端环平面约16.97720°；冈下肌对应约−0.25531 mm／0.36638°。需测试斜向出线、侧向接触、弯折及疲劳；单靠固定端点的销轴转动不能消除轴向错位。该细节不直接判为失效，亦不以名义2.55 mm端眼半径冒充厂商允许半径。

[MuJoCo端部查看](../../../mujoco_arm/tools/view_shoulder_free_strand_route.py)只显示实际整根钢丝及4个眼环／套筒，隐藏骨与TPU，静态DOF=0。规定角示例不是变形或驱动动画。当前S1-Y2几何及选择指针未改，后续需让肱骨、自由端盖与TPU共同求解，实际驱动范围未知。

'''
edit(D/'SHOULDER_DESIGN.md',[
 ('版本：v0.9.14','版本：v0.9.15'),
 ('## 3. 肌肉／近端支承父级已验证内容与继承边界',section+'## 3. 肌肉／近端支承父级已验证内容与继承边界'),
 ('| 4 | 完善已建立的自由末端／肌腱端点与张力单边约束，接入TPU平衡与变姿态路由 | 当前仅有源端点、有效段长和局部转动检查；需安装长度／压接／预紧、力矩与行程、柔性平衡及实测对照，不借用裸骨或旧刚体范围 |',
  '| 4 | 完善已建立的自由末端／肌腱端点与张力单边约束，接入TPU平衡与变姿态路由 | 已有规定端环角的完整肌腱候选；仍需安装长度／压接／预紧、最短路由或线弯曲模型、力矩与行程、柔性平衡及实测，不借用裸骨或旧刚体范围 |'),
])
with (D/'SHOULDER_DESIGN.md').open('a',encoding='utf8') as f:f.write('\n2026-11-02 v0.9.15：补规定端环角的中央段等长候选，整根钢丝在场；18例密曲线／4导出场景复核及实际MuJoCo五件细节图。记录肩胛下肌约16.98°斜向出线。几何与指针不变，仍非安装长度、张紧平衡或驱动范围。[更新前快照](../../../reviews/S1_2026-11-02/free_strand_routes/before_SHA256.json)。\n')

edit(ROOT/'PROJECT_STATUS.md',[
 ('技术文件v0.9.14','技术文件v0.9.15'),
 ('[肩部技术文件v0.9.14]','[肩部技术文件v0.9.15]'),
 ('## 最新进展与限制\n','## 最新进展与限制\n\n新增[中央段完整肌腱候选](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/free_strand_routes/README.md)：规定端环角、保留零位端盖／肱骨，补上中央线并检查整根等名义长度的钢丝与两端硬件；16条采样路径、18例独立密曲线和4个实际网格诊断场景复核。肩胛下肌零位斜向出线约16.98°，需弯折／侧向接触试验。候选不是张紧平衡或全局最短路由，完整TPU驱动仍未求解。\n'),
 ('完整TPU平衡和中央段变姿态路由仍未解。','现已补规定端环角的中央段候选；自由端盖／肱骨变姿态下的物理路由和完整TPU平衡仍未解。'),
])
edit(ROOT/'design/CURRENT_BASELINE.md',[
 ('技术文件2.14–2.17节','技术文件2.14–2.18节'),
 ('## 文档分工','新增完整肌腱候选用于规定端环角的局部诊断：保留两出口切向和名义长度，中央线在场；它不证明最短路由、张紧平衡或实际长度。肩胛下肌约16.98°斜向出线需测试，源骨架／锚点不改。\n\n## 文档分工'),
])
brief=ROOT/'design/DESIGN_BRIEF.md';raw=brief.read_bytes();marker='## 历史记录阅读说明'.encode('utf8');head,tail=raw.split(marker,1)
s=head.decode('utf8').replace('技术文件v0.9.14','技术文件v0.9.15')
s+='[中央段完整候选](../02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/free_strand_routes/README.md)补上规定端环角下的整根等名义长度肌腱及端部硬件；端盖／肱骨保持零位，未替代柔性平衡或实际驱动。当前S1-Y2源与指针不变。\n\n'
brief.write_bytes(s.encode('utf8')+marker+tail)

for p in [ROOT/'README.md',ROOT/'02_人体骨骼晶格肌肉方案/README.md']:
    s=p.read_text(encoding='utf8').replace('技术文件2.14–2.17节','技术文件2.14–2.18节')
    if p.parent.name=='02_人体骨骼晶格肌肉方案':s=s.replace('v0.9.14','v0.9.15')
    p.write_text(s,encoding='utf8')
edit(D/'s1_y2_transition_stock/README.md',[
 ('## 当前数字证据','另补[规定端环角的整根肌腱候选](free_strand_routes/README.md)，中央段与两端硬件均在场，含MuJoCo细节查看。此项保留端盖／肱骨零位和未测的名义长度，不作张紧平衡或GH范围。本装配源及选择指针保持不变。\n\n## 当前数字证据'),
])
print('Technical v0.9.15 records complete diagnostic strands, while preserving current S1-Y2 sources.')
