"""Record source-based free-terminal constraints without changing geometry."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'

def edit(p,pairs):
    s=p.read_text(encoding='utf8')
    for a,b in pairs:
        if a not in s:raise ValueError(f'Missing expected text in {p}: {a[:60]}')
        s=s.replace(a,b,1)
    p.write_text(s,encoding='utf8')

section='''### 2.17 自由末端／肌腱约束与局部端环转动

[约束与局部检查](s1_y2_transition_stock/free_terminal_kinematics/README.md)建立真实销轴和两压接套外出口的端点模型，未改变S1-Y2几何。整根含端环／尾段的95.41817／92.12987 mm路径，不能当作中央肌腱自由长度；自由段名义中心线实际39.92927／36.64097 mm，零位出口直线39.92595／36.14845 mm。曲线差值不是实测松弛，安装材料长度与预紧仍未知。零位直线钢丝探针对其余190对象无超判据相交，不证明变姿态路由已求解。

对规定肱骨姿态，每块肌肉包含自由端盖6变量和两端眼环角2变量。理想无滑移、不可伸长、无阻直线张紧条件下，长度约束秩1，仍余7未知量；松弛时不激活，不能用负张力推动端盖。必须接入TPU变形、端盖接触／力矩平衡、眼环摩擦／止挡、主动线分布反力与真实绕障路由。代码拒绝未知安装长度或未核查直线，并保留假设／实测来源，不输出虚构张力或已解平衡。

| 端环局部组 | 负向最后未超判据角 | 正向最后未超判据角 | 首个刚性对象 |
|---|---:|---:|---|
| 冈下肌肌肉侧 | −180°检查上限 | +10.9° | 正向自由端盖 |
| 冈下肌骨侧 | −5.1° | +152.9° | 抱箍上半夹体 |
| 肩胛下肌肌肉侧 | −14.5° | +180°检查上限 | 负向自由端盖 |
| 肩胛下肌骨侧 | −165.1° | +19.3° | 抱箍下半夹体 |

每次仅转一端实际槽环／套筒／局部端线，其他对象零位；自己的中央段路由未求解、原整根肌腱源从步骤中移除。故表格既不是完整肌腱路径，也不是GH或TPU驱动角度。1°常规、0.1°首界、独立回程，14个边界姿态另行重放一致；若干最后采样间隙为0，0.001 mm³体积判据不是正制造间隙。

柔性相交另行保留：冈下肌肌肉侧负向约−123.8°、骨侧正向约+75.1°首次细化采样接触未变形三角肌，不以刚性范围取代柔性接触验证。雅可比的16组合成代数状态中央差分一致，但这些不是可达平衡姿态。[测量空白记录](s1_y2_transition_stock/free_terminal_kinematics/measurement_template.json)要求安装自由长度／4压接套保持、预紧、端盖姿态、眼环摩擦及拉线过程；全部未测。完整驱动范围继续未知。

'''
edit(D/'SHOULDER_DESIGN.md',[
 ('版本：v0.9.13','版本：v0.9.14'),
 ('## 3. 肌肉／近端支承父级已验证内容与继承边界',section+'## 3. 肌肉／近端支承父级已验证内容与继承边界'),
 ('双孔套及主动线止挡均为未压接候选，不能视为可靠终端固定。','双孔套及主动线止挡均为未压接候选，不能视为可靠终端固定。两套筒外出口间的中央自由段与整根长度不同，见2.17节。'),
 ('| 4 | 建立自由末端、固定长度肌腱及柔性接触模型 | 真实约束归属、力矩与行程、姿态平衡、与实测对照；不得借用旧刚体范围 |','| 4 | 完善已建立的自由末端／肌腱端点与张力单边约束，接入TPU平衡与变姿态路由 | 当前仅有源端点、有效段长和局部转动检查；需安装长度／压接／预紧、力矩与行程、柔性平衡及实测对照，不借用裸骨或旧刚体范围 |'),
])
with (D/'SHOULDER_DESIGN.md').open('a',encoding='utf8') as f:f.write('\n2026-11-02 v0.9.14：补真实压接出口／中央段长、张力单边约束及4端环局部转动检查；16组合成状态数值导数、14个边界独立重放与零位直线探针。几何与指针不变，柔性平衡／实际驱动未解。[更新前快照](../../../reviews/S1_2026-11-02/free_terminal_kinematics/before_SHA256.json)。\n')

edit(ROOT/'PROJECT_STATUS.md',[
 ('技术文件v0.9.13','技术文件v0.9.14'),
 ('[肩部技术文件v0.9.13]','[肩部技术文件v0.9.14]'),
 ('## 最新进展与限制\n','## 最新进展与限制\n\n新增[自由末端／肌腱约束与局部端环检查](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/free_terminal_kinematics/README.md)：两中央自由段名义路径39.92927／36.64097 mm，不能用整根含端环95／92 mm代替。已建立端点／张力单边约束，4端环1°／0.1°局部往返检查；给定肱骨时一条张紧长度仅约束8变量中的1个，完整TPU平衡和中央段变姿态路由仍未解。当前几何与指针不变。\n'),
 ('3. **真实约束与驱动未知。**','3. **真实约束与驱动未知。** 已有端点及条件长度约束、局部转角证据；安装长度、预紧、压接保持和柔性平衡仍缺数据。'),
])
edit(ROOT/'design/CURRENT_BASELINE.md',[
 ('技术文件2.14–2.16节','技术文件2.14–2.17节'),
 ('## 文档分工','当前自由端／肌腱模型按真实销轴、两压接外出口建模；中央自由段名义39.92927／36.64097 mm与整根钢丝长度分开。安装长度、预紧、压接保持与TPU平衡未测；局部眼环范围不是GH范围，端盖仍不刚性绑定空间或肱骨。\n\n## 文档分工'),
])
brief=ROOT/'design/DESIGN_BRIEF.md';raw=brief.read_bytes();marker='## 历史记录阅读说明'.encode('utf8');head,tail=raw.split(marker,1)
s=head.decode('utf8').replace('技术文件v0.9.13','技术文件v0.9.14')
s+='已新增[自由末端／肌腱约束](../02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/free_terminal_kinematics/README.md)，区分整根钢丝与两压接出口间自由段。实际安装长度／预紧／压接、TPU平衡及变姿态路由待验证，局部端环角不当作肩部驱动范围。当前几何与指针不变。\n\n'
brief.write_bytes(s.encode('utf8')+marker+tail)

for p in [ROOT/'README.md',ROOT/'02_人体骨骼晶格肌肉方案/README.md']:
    s=p.read_text(encoding='utf8').replace('技术文件2.14–2.16节','技术文件2.14–2.17节')
    if p.parent.name=='02_人体骨骼晶格肌肉方案':s=s.replace('v0.9.13','v0.9.14')
    p.write_text(s,encoding='utf8')
edit(D/'s1_y2_transition_stock/README.md',[
 ('## 当前数字证据','新增[自由末端／肌腱约束及4端环局部运动](free_terminal_kinematics/README.md)，明确两压接出口间的中央段长度、张力单边约束和未解的7个端部变量。局部端环扫描、零位直线探针均不作肩部驱动范围；本轮不改本装配源或指针。\n\n## 当前数字证据'),
])
print('Technical v0.9.14 updated; S1-Y2 source geometry and original historical brief unchanged.')
