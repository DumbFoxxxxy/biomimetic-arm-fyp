# 整臂首版 · 共享装配与局部钢丝

## 查看最新装配

[首版整臂 MJCF](whole_arm_first_pass_assembled.xml) · [零件来源与回读检查](whole_arm_first_pass_assembly.json) · [十条钢丝与接口交接表](whole_arm_cable_handoff.md) · [锚点坐标 CSV](whole_arm_cable_handoff.csv) · [钢丝折线及导向接口坐标](whole_arm_guide_points.json)

现有整臂的十条局部钢丝与 OpenSim 肌束完成了逐束对照。肘部沿用已有导向与转座，肩胛路径保留附着区和作用方向差距记录。[对照报告、逐束数据与 MuJoCo 回放](OPENSIM_COMPARISON.md)

首版装配包含桌面胸廓、肩胛、肩部 S1、肘部、八块 TPU、修订后的两处弹簧端座和中立位回位弹簧。用 `--first-pass` 可旋转查看 239 个几何体及八条固定端钢丝。肩袖活动件和两条中央短肌腱按不同载荷改变位姿，在下方 K7 回放中显示。

### K7 整臂演示

[播放首版双视角回放](whole_arm_k7_first_pass.mp4) · [查看终点截图](whole_arm_k7_first_pass.png) · [渲染来源和帧数](whole_arm_k7_first_pass_media.json)

MuJoCo 回放覆盖 14 个去程样点及同一路径的回程，共 27 帧。桌面底座、胸廓、肩胛、肩关节、肘部及回位弹簧同屏，八块 TPU 使用几何审核中的插值，20 件肩袖活动件使用局部 10 N 工况的位姿。十条线显示已布置的局部传力段或肩袖中央短肌腱。终点肱骨方向相对中立位改变约 17.33°。动画使用预设关节姿态，整臂独立驱动和 TPU 受力变形仍待验证。

### K7 重力平衡估算

在十线共同通过几何检查的 14 个姿态中，沿预设协调方向的重力需求为 3.94–6.76 N·m。此力矩以 SC 路径输入为广义坐标，包含随该输入联动的肩胛、盂肱和肘部贡献。若另有约束负责其余方向，仅由最有利的一条钢丝沿该坐标抗重力，名义拉力约 59.7–105.8 N。两端运动和分段线长与原几何结果逐姿态交叉核对。

再让八条固定端路径分担十个转动坐标的重力。在每条线 0–10 N 的假设下，最小化最大关节力矩余差后仍剩 2.30–2.95 N·m。加入两条肩袖各 10 N 局部加载工况的等效 GH 力偶后，余差区间相同。八线每路上限提高到 30 N，仍有 2.25–2.89 N·m 的余差，主缺口集中在胸锁和肩锁的对应方向。10 N、30 N 为比较工况。肩胛—胸廓接触反力、底座以外的支撑或关节约束、TPU 载荷和滑块受力均未纳入，实际独立驱动能力仍待验证。[逐姿态受力估算](whole_arm_k7_gravity_balance.json)

将既有三块单边支承垫的压缩反力加入每路 30 N 的估算。中立位三块均接触，SC 2.67° 样点起全部脱开，原因是法向间距超过垫的自由高度，横向偏移仍在原 12 mm 半径内。SC 0.67° 样点的最大关节力矩余差由 2.35 降到 1.28 N·m，终点仍为 2.89 N·m。三块垫的法向行程需求约为 6.58、2.35、7.70 mm。按原假设刚度 10 N/mm，单靠加高来保持全程接触，三块垫的峰值压力至少约 65.8、23.5、77.0 N。该下限用于比较支承结构，承压、材料刚度和受力动态仍待核对。[三垫接触与力矩对照](whole_arm_k7_scapular_support.json)

沿相同的 14 个姿态比较长行程支承。每路钢丝仍以 30 N 为计算上限。若只保留第 4 号接触点并假定 1 N/mm 的线性弹簧，拟合需要中立位预压约 26.1 N、法向工作行程 6.58 mm，终点反力约 19.5 N。全路径最大关节力矩余差为 1.33 N·m，低于三垫同为 1 N/mm 的 1.61 N·m。即使允许原三个接触点的反力在各姿态独立变化、取消垫的力上限，终点仍至少剩 1.33 N·m，主要涉及 SC 和 AC。由此将第 4 号点保留为被动支承候选，余差记录为当前八线与原接触法向的受力缺口。拟合出的预压和行程尚未落实到具体可装配件。[长行程支承与理想反力对照](whole_arm_k7_scapular_support.json)

回位弹簧让位前的对照装配可用 `--latest-tpu` 查看。[先前中立位 MJCF](whole_arm_latest_neutral.xml)包含 238 个几何体，作为修订端座和晶格前的参考。

共享参照采用[带实体胸锁与肩锁连接候选的 98 对象装配](../channel_first/unified_SC_AC/assembly.json)。在此基础上整合肩肘固定件、肘部两转座与伸肌滑块、肩胛带四件导向眼。整臂模型含 229 件刚性结构、13 个速度自由度和 6 件既有中立位 TPU 参考几何。查看时可载入另行检查的桌面底座，胸廓在仿真中固定。[组合模型](whole_arm_shared.xml) · [重新加载检查](whole_arm_shared_scene_check.json)

## 胸廓桌面底座

底座采用 310 × 262 × 10 mm 台脚和短立柱，立柱顶端为按第十二胸椎外形生成的承座，预留约 0.2 mm 名义胶层空间。承座、立柱和台脚导出为单个连通实体。几何检查覆盖固定胸廓及 710 个运动采样姿态，移动刚性件距底座至少 15 mm（检查距离上限）。台脚到胸廓的支承方案采用 T12 粘接，胸廓骨件之间的实物连接、胶层承载和抗倾覆尚待核对。现有回放姿态及重力边界沿用共同模型。[底座网格](tabletop_T12_socket.obj) · [底座尺寸及间隙](tabletop_base_check.json)

## 肩胛带两块 TPU 晶格

沿 FJ1459_g3 前锯肌和 FJ1536_g2 大菱形肌的已选钢丝生成了两块变截面晶格。横截面方向和宽度取自人体肌肉网格，沿肌腹收窄至两端加厚过渡。每块均有直径 1.6 mm 的通道和两端各 2 mm 的裸露短线。中立位按相邻骨面局部让位，导出回读后两块均为单个连通实体，且无超阈值刚体相交。与可读入的五件既有 TPU 网格在中立位无超阈值相交；肘屈肌参考存在非流形边，无法参与本次实体求交。[前锯肌晶格](FJ1459_g3_TPU_lattice.obj) · [大菱形肌晶格](FJ1536_g2_TPU_lattice.obj)

在 11 个肩胛和 31 个联合抬升姿态中，用固定截面、沿钢丝轴向伸缩的几何近似检查新肌腹。两块之间、各自通道与钢丝之间以及与对侧钢丝均未测得相交。前锯肌需要约 32.1% 的轴向伸长，最大肋骨重叠约 1.07 mm³；大菱形肌需要约 2.36% 轴向缩短，最大肩胛骨重叠约 0.070 mm³。以上是几何输入与局部接触信号，实际晶格行程、接触变形及端部固定待验证。[参数及逐姿态结果](girdle_tpu_geometry.json)

## 新增导向眼与肩肘接口整合

- 前锯肌 FJ1459_g3 与大菱形肌 FJ1536_g2 各增加两个随骨端部导向眼，名义孔径 1.3 mm，基座按骨面修形。四件 OBJ 回读后，与两条直径 0.60 mm 钢丝一起在原 98 件装配中通过 21 个去程姿态和 21 个回程姿态。四件的连接方式、钢丝端部止挡和受载能力待确认。[导向眼装配与去回程记录](girdle_seat_pair.json)
- 将既有肩肘接口转到共享骨架时，中立位出现 12 个刚性重叠。保留共享胸锁、肩锁和盂肱结构；在来源坐标 Z=1180 mm 处接入已建立止挡的肱骨远段，并恢复 4 个肘肌支架孔和 2 个三角肌 FJ1467 盲孔。修订件重新导出后，新增或修订对象相关配对在 0.001 mm³ 阈值下无重叠。全部刚体跨骨段配对亦在 MuJoCo 回读零位复核。[移植前后对照](interface_transfer_audit.json)
- 六件 TPU 网格用于中立位查看，回放时隐藏。整合模型中的局部钢丝过孔、止挡和转向支撑，以及 TPU 变形和组合空间仍需验证。

## 整臂动作范围与首限

以下范围来自共享装配的规定姿态。刚体按所属骨段运动；肘部两转座和滑块保持中立，采样 1°，单轴首碰边界加密至 0.1°。数值是模型几何范围，实测范围待补。原目标和底线保留。[逐姿态记录](whole_arm_shared_motion.json)

首版件的弹簧与 TPU 组合复核覆盖 K1 和 K7。K2、K3、K6 及肩胛前后运动的八线范围对应相同骨架与端口的刚体和钢丝检查，尚待加上装配件随动与弹簧、柔性肌腹的空间复核。K4、K5 的肩袖结果为单通道加载，整臂同时受力仍按缺口登记。

| 指标 | 原目标 | 原底线 | 首轮刚体结果 | 限制或待验证 |
|---|---:|---:|---|---|
| K1 肘屈曲 | 145° | 116° | 原端盖八线去回程 0–8.16°；端盖让位候选 0–31.31°；仅两条肘线的加长导轴方案 0–126.01°；固定转座及滑块的刚体探针 0–145° | 让位候选在 31.41° 首次超过 FJ1467 钢丝相交阈值。加入弹簧几何后中立位即有接触，屈肘后有肌腹挤压，含回位弹簧的整臂可用范围待确认 |
| K2 GH 屈曲 | 120° | 96° | 0–120° 刚体探针；0–41.7° 八线几何去回程 | 三角肌前束 FJ1468 持续为屈曲方向，收线 21.85 mm。41.4° 样点出现零间隙，41.8° 钢丝与肘屈肌转接座相交超阈值 |
| K3 GH 外展 | 120° | 96° | 0–120° 刚体探针；0–13.6° 八线几何去回程 | 三角肌中束 FJ1467 持续为外展方向，收线 5.16 mm。11° 样点出现零间隙，13.7° 钢丝与肩胛骨相交超阈值 |
| K4 GH 内旋 | 70° | 56° | 源长轴负向 0–90°，达到刚体探针上限；同范围内 FJ1504 每 10° 的独立局部加载通过 | 肩胛下肌的 10 N 输入产生 −0.095 至 −0.007 N·m 长轴力矩，卸载序列亦复核；整臂同时加载、其他线及 TPU 仍待复查 |
| K5 GH 外旋 | 90° | 72° | 源长轴正向 0–82.2° 刚体探针；FJ1500 独立加载通过 0° 和 10° 两个样点 | 刚体 82.3° 首碰肩胛骨与伸肌近端螺栓。独立局部加载在 20° 样点首次出现冈下肌短肌腱碰共用支座 |
| K6 肩胛上回旋 | 50° | 36° | 正向候选的肩胛参考线冠状面转动 0–10°，含八线去回程 | 到达本轮上限。以 AC 到 GH 连线的冠状面转角近似上回旋。大菱形肌线需收 2.45 mm，前锯肌线放 9.77 mm。活动肩袖及 TPU 仍待复查 |
| K7 联合抬升 | 170° | 132° | 肱骨方向相对中立位为 0–40.00° 刚体几何探针；带十线、TPU 插值和肩袖活动件的组合演示为 0–17.33° | 刚体达到探针上限。组合回放止于末个共同检查样点，余下样点未同时复核。TPU 仍有局部几何重叠，独立驱动待验证 |

K1 核对了 230 个刚体和桌面支承的源文件、所属骨段及安装位置，使用与 K7 终版预览相同的 MJCF。复用历史动作中的肘轴、两转座和伸肌滑块位移，肘屈肌钢丝收线 2.26 mm、沿路径等效力臂 +9.22 至 +17.85 mm，伸肌放线 3.79 mm、等效力臂 −35.42 至 −17.08 mm。FJ1467 与伸肌滑动端盖在最后有正间隙的 8.08° 样点相距约 0.0024 mm。边界是相邻原始姿态插值，按约 0.086° 步距细化。K1 回程沿相同姿态逆序复算。加长导轴额外材料在 K7 的 14 个样点与固定刚体、八条固定端钢丝及肩袖 10 N 活动件、短肌腱未出现超阈值相交。当前共用回放采用原长导轴；加长方案与整臂 TPU、弹簧和装配工序的兼容性继续保留为待验证事项。[统一装配的 K1 八线与 K7 导轴对照](whole_arm_k1_shared.json)

端盖让位候选以 30° 的短程演示为设计目标。在滑动端盖坐标中扫掠 FJ1467 的路径，从原端盖局部去除约 25.90 mm³，回读后保留单一连通实体。设计检查未切入两根导轴、紧固件、弹簧移动销及中立位伸肌肌腹的保护区。八线去回程中屈肌收线 11.53 mm，伸肌放线 12.04 mm，31.31° 样点的最小间隙约 0.019 mm。K7 的 14 个样点中，修订端盖与所查刚体、钢丝、10 N 肩袖活动件及短肌腱无超阈值相交，端盖与共用支座最近约 0.15 mm。该候选保留原长导轴与钢丝路径，可在同一个查看器回放 K1 和 K7。端盖承载、加工公差、回位弹簧和肌腹行程需后续校核。[让位端盖 OBJ](FJ1477_shared_cap_wire_relief.obj) · [让位与原方案对照](whole_arm_k1_cap_relief.json)

在让位后的 K1 轨迹上加回回位弹簧和八块 TPU 肌腹。弹簧按两个原有销心构造，销距由 26.00 增至 39.90 mm。59 个去程样点中，弹簧未碰八条固定端钢丝，但从中立位起与固定支座、滑动端盖分别重叠约 0.47、0.48 mm³，与三角肌晶格的名义重叠为 10.34–13.65 mm³。历史弹簧销止退螺钉方案尚未合入这套支座和端盖。按约 5° 间隔及终点复核八个 TPU 姿态，两块随动肘肌与所查刚体、钢丝和弹簧无超阈值相交；三角肌与肘伸肌从 4.98° 样点开始出现约 14.61 mm³ 的名义重叠，31.31° 时约 113.43 mm³。肩部保持中立，肘部肌腹沿端点插值；实际挤压量及阻力待验证。[弹簧与 TPU 组合空间数据](whole_arm_k1_combined_space.json)

在原两销心周围做固定、滑动端座让位，分别去除约 2.56 和 2.17 mm³。以弹簧扫掠位置为依据局部修剪三角肌晶格，第一轮与第二轮合计去除约 144.47 mm³，保留一个正体积连通件和三条钢丝的 0.72 mm 通道保护半径。第二轮复核 K1 的 59 个姿态，弹簧与刚件、固定端钢丝和静止 TPU 均无超阈值相交，8 个随动肘肌样点中弹簧也未发生超阈值相交。此处弹簧端部用原销模拟，轴向止退和装配顺序待验证。三角肌与伸肌依然在屈肘中接触，所查最大名义重叠为 19.66° 样点的 147.23 mm³，需单独校核局部柔性行程。[固定端座](K1_spring_relief_fixed_seat.obj) · [滑动端座](K1_spring_relief_moving_seat.obj) · [三角肌晶格](K1_spring_relief_deltoid_refined.obj) · [K1 复核](whole_arm_k1_spring_relief_refined.json)

把同样的端座、晶格和弹簧放入 K7 的 14 个姿态复算。弹簧与刚件、钢丝及所查肩袖活动件无超阈值相交，三角肌随肩部运动时仍有 12 个样点出现弹簧局部重叠，最高约 1.54 mm³。三角肌晶格与所查钢丝保持无超阈值相交，三角肌与伸肌的峰值名义重叠仍约 14.90 mm³。弹簧接触量来自规定姿态及 TPU 插值，接触力和磨擦还需验证。[同一让位方案的 K7 复算](whole_arm_k7_spring_relief.json)

K2 与 K3 按 1° 采样、首碰加密至 0.1°，使用八条当前固定端钢丝、肩袖共用支座和桌面底座。K2 的 FJ1468 等效力臂为 +22.29 至 +33.93 mm，K3 的 FJ1467 为 +21.22 至 +22.03 mm。两条路线的末段名义间隙已为零。最后仍有正间隙的已测样点分别是 K2 的 41.3° 和 K3 的 10.0°，对应约 0.0009 mm 和 0.0014 mm 的几何间隙，装配公差待验证。当前回放规定关节姿态，活动肩袖及 TPU 变形的组合空间仍待复核。[屈曲外展逐姿态、钢丝行程与首限](whole_arm_shoulder_elevation.json)

### 肩胛前后运动与肩部旋转

以胸骨柄和第四胸椎确定胸廓前方。肩胛前后运动用 GH 中心相对胸廓的前后位移表示，GH 相对肩胛保持中立。用肘屈曲方向构造一个垂直于肱骨长轴的参考向量，观察它在长轴探针的第一步朝向或远离胸廓中线，分别标定 K4 内旋与 K5 外旋。内外旋表中的刚体范围仍来自原伸肘附近的探针姿态。

| 动作 | 八线去回程的几何范围 | 局部钢丝方向与首限 |
|---|---|---|
| 肩胛前伸 | SC 正向 0–5.5°，GH 中心向前 0–15.49 mm，AC 中心向前 0–13.15 mm | 前锯肌线正向力臂保持到 4.0° 样点，对应前移 11.23 mm、收线 1.69 mm，大菱形肌放线 6.43 mm。4.5° 的前锯肌力臂已变号；5.6° 首次超过钢丝与前锯肌肩胛侧导向眼的相交阈值，后段名义间隙已降至零 |
| 肩胛后缩 | SC 负向 0–5.7°，AC 同轴反向补偿比例 1.1，GH 中心向后 0–12.88 mm | 大菱形肌线正向力臂保持到 5.0° 样点，对应后移 11.34 mm、收线 0.61 mm，前锯肌放线 4.75 mm。5.5° 的大菱形肌力臂已变号；5.8° 首次超过前锯肌钢丝与胸廓侧导向眼的相交阈值，后段名义间隙已降至零 |
| 肩胛上回旋 | SC 轴 1 正向 0–10°，对应冠状面参考线转动 0–10° | 大菱形肌线等效力臂 +12.04 至 +16.04 mm，前锯肌线 −57.85 至 −53.19 mm，肩胛线参与的最小名义间隙约 0.033 mm |

镜像后缩在 SC 0.1° 样点即出现肩胛下肌近端盖与 FJ3340 肋骨相交。表中后缩改用 AC 与 SC 反向协调，从小范围比例试探中选出第一个保留大菱形肌初段正向力臂的分支。前伸与后缩分别从中立位回放，实际协调与近零间隙的装配公差仍待验证。

肩胛局部路径用 0.5° 步长检查，首碰用 0.1° 细化，包含共享支座、桌面底座及八条固定端钢丝。肩袖单通道加载记录随内外旋方向重标。肩胛动作的活动肩袖和 TPU 随动、全臂受力和实物回程仍待验证。[逐帧方向、行程、间隙及回程](whole_arm_anatomical_actions.json)

[端点等效力臂和行程](whole_arm_anchor_motion.json)包含四条肩部、两条肘部及两条肩胛带钢丝的骨段归属、坐标和逐姿态结果。以肘部为例，FJ1486 的直线代理全程缩短，FJ1477 全程伸长，峰峰线长分别为 53.94 mm 和 111.28 mm。它们是两转座和滑块固定时的端点距离变化。冈下肌 FJ1500 与肩胛下肌 FJ1504 的浮动端单独登记。整臂可用钢丝行程及驱动能力需在实际导向和 TPU 结构上重算。

### 联合抬升回放

比较了四条 SC/AC 路径，每条分别测肩胛带单独运动和同时抬肩屈肘。首选回放 `K7_SC1_pos` 有 31 个去程采样，已按同一姿态集复查回程。在该回放中，两条肩胛钢丝同时与所有刚性件及彼此求交；线材与刚性件最小间隙约 **0.033 mm**，仍需设计公差和受载验证。SC 轴 2 正向分支首先出现前锯肌钢丝碰自身导向眼，SC 轴 1 负向分支首先是肩胛骨碰 FJ3342 肋骨。

沿 `K7_SC1_pos`，前锯肌 FJ1459_g3 名义线长增加 9.766 mm，大菱形肌 FJ1536_g2 名义线长减少 2.453 mm。这两组方向只适用于该采样分支。该分支由规定姿态生成，模型尚未求解肩胛胸廓支承接触与肌力平衡。

### 既有肩肘钢丝路径的装配复核

把四条肩部、两条肘部既有折线和两条肩胛钢丝放入同一模型。肩部和肘部按各线的端部所属骨段插值。联合抬升检查全部 31 个去程采样，回程使用相同姿态；单轴为每隔约 5° 的初筛。下表是联合抬升从中立位到首个受阻采样之前的区间，角度为 SC 轴 1 输入。肱骨方向变化在该路径约为 SC 输入的四倍。钢丝直径按 0.60 mm 检查，仍需结合实际导向和张紧状态验证。[折线、逐姿态行程及干涉数据](whole_arm_local_wire_audit.json) · [首碰线段](whole_arm_local_wire_limits.json)

| 通道 | 已清 SC 采样上限 | 区间峰峰线长 | 首个超阈值相交 |
|---|---:|---:|---|
| FJ1467 中束 | 0.67° | 0.823 mm | SC 1.00°，肘伸肌 A 侧导轴，距近端约 121 mm |
| FJ1513 后束 | 4.67° | 0.924 mm | SC 5.00°，同一导轴，距近端约 109–123 mm |
| FJ1468 前束 | 7.67° | 0.266 mm | SC 8.00°，屈肌骨侧转接座及螺栓，距近端约 147–153 mm |
| FJ1507 大圆肌 | 10.00° | 19.138 mm | 31 个样点内未发现 |
| FJ1486 肘屈肌 | 10.00° | 1.915 mm | 31 个样点内未发现 |
| FJ1477 肘伸肌 | 10.00° | 17.643 mm | 31 个样点内未发现 |

八条纳入的钢丝之间在该联合路径未出现超过阈值的相交。冈下肌 FJ1500、肩胛下肌 FJ1504 因浮动端缺可复现的受拉姿态而单独登记。肘部单轴初筛在转座和滑块均保持中立时，屈肌钢丝到 35° 样点清楚，40° 时在距近端约 107–117 mm 处遇到肱骨。145° 刚体范围须与此钢丝结果分别使用。后续调整优先针对靠远端的三角肌导向，以及肘屈肌的转座随动。

### 肩袖浮动端局部受力

在同一共享骨架上分别检查冈下肌 FJ1500 和肩胛下肌 FJ1504。复用原有浮动端、转眼和压接套，给定关节姿态后，以未标定的弹性参数求解 0 → 5 → 10 → 5 → 0 N 局部加载。先检查 K7 的 11 个姿态和 GH 长轴正、负向各 10 个姿态。求解包含浮动端平移与转动及两端转眼，检查硬件和受拉的中央直跨段，保留松线形态、压接出线弯折和 TPU 整体变形待验证。

| 通道 | K7 局部筛查 | 10 N 输入时 GH 长轴力矩分量 | 主要限制 |
|---|---|---:|---|
| FJ1500 原位 | 中立位在加载中即有相交 | +0.073 至 +0.158 N·m | 骨侧压接套与肘伸肌转接座，以及中央线段的刚体相交 |
| FJ1504 | 11 个姿态的加载和卸载均无所查相交 | −0.135 至 −0.094 N·m | K7 的最大预测端盖平移约 26.3 mm，肌腹是否能支撑该位移尚待核对 |

FJ1504 的源长轴负向共 10 个姿态均通过局部筛查，达到 90° 检查上限。正向到 30° 的采样通过，40° 样点开始出现局部干涉。原位 FJ1500 仍无从中立出发的完整无干涉加载区间。此前设计的 16 mm 上移方案令冈下肌浮动件和中央短肌腱在 K7 的 10 N 样点避开原相交，但上移支座本身与伸肌转接座及螺栓分别重叠约 122.57 与 17.65 mm³，因此保留为候选。所查结果用于确定下一处工程接口，实际驱动及整臂同时受拉仍待验证。[原位浮动端](whole_arm_floating_terminal_audit.json) · [16 mm 上移候选](whole_arm_offset_infraspinatus_audit.json)

整合版将冈下肌的上移支座与肘伸肌转接座合为一件，在旧支座局部保留螺栓让位。保留原转接座及冈下肌转轴的承压区域，两个部分合并前的共用材料体积为 122.32 mm³。连通网格回读后，中立位和 K7 的 31 个刚性姿态未发现新增相交。配套冈下肌浮动端在 K7 的 11 个姿态均完成 0 → 5 → 10 → 5 → 0 N 的局部加载与卸载检查。长轴正向在 20° 样点中央线碰新支座，负向在 40° 样点活动套碰肘伸肌 B 侧固定螺钉，前一采样分别为 10° 和 30°。支座装入和紧固顺序、两束肩袖同时受力、与其余钢丝和 TPU 的联合检查待完成。[整合支座网格](FJ1500_FJ1477_shared_support.obj) · [支座配合](whole_arm_infraspinatus_support.json) · [浮动端复核](whole_arm_shared_infraspinatus_audit.json)

两条肩袖的浮动端与八条已过初筛的局部钢丝已在 K7 同一装配中复核。14 个姿态共 70 组加载状态中，未发现所查肩袖活动件、受拉中央跨段、相邻钢丝或新增固定支座的超阈值相交。最后共同通过的样点为 SC 4.33°，实际肱骨方向变化约 17.33°。两条肩袖分别求静力平衡后放入同一姿态比较。MuJoCo 的 27 帧回放显示 10 N 工况中受拉的中央跨段，并使用新版固定支座及骨侧紧固件。浮动硬件仅参与求交，回放只显示跨段；完整 TPU 接触、压接出线弯折与整臂动力学待验证。[十线路径组合检查](whole_arm_combined_cuff_k7.json)

### K7 八块 TPU 肌腹组合检查

在十线回放的 14 个姿态上，按钢丝和端部运动插值摆放八块 TPU 肌腹，分别检查肌腹与肌腹、钢丝及固定五金的重叠。该计算沿用肩袖 10 N 的局部预测位姿，交叠量用于定位空间冲突。肩袖活动五金已经纳入十线检查，尚未与这批变形肌腹逐对求交；弹簧与肘部晶格的检查仍对应单独的屈肘回放。

| 主要配对 | 出现样点 | 最大名义重叠 |
|---|---:|---:|
| 三角肌肌腹与肘伸肌肌腹 | 14 / 14 | 464.83 mm³ |
| 三角肌肌腹与肘伸肌近端盖 | 14 / 14 | 249.15 mm³ |
| 三角肌肌腹与肘伸肌 A 导轴 | 14 / 14 | 210.68 mm³ |
| 三角肌肌腹与冈下肌肌腹 | 9 / 14 | 3.50 mm³ |
| 大圆肌肌腹与自身钢丝 | 10 / 14 | 14.28 mm³ |
| 三角肌肌腹与中束 FJ1467 钢丝 | 14 / 14 | 1.78 mm³ |

以三角肌中立位相交体积为切削输入，保留三条钢丝中心线周围半径 0.95 mm 的保护区。试算去除约 1627.61 mm³ 材料，晶格分成七个正体积实体，近端盖仍残留约 0.0051 mm³ 重叠，故未采用该切削。后续三角肌与肘伸肌相邻区需要调整晶格连接或肌腹局部布置，并让完整 TPU 通道跟随对应钢丝。K7 约 17.33° 仍用于十线几何回放，TPU 布置的关节运动范围待复核。[八块肌腹逐姿态结果](whole_arm_k7_tpu_audit.json) · [保留通道的切削对照](whole_arm_k7_deltoid_relief_trial.json)

补做一根半径 0.65 mm 的连接肋，并去除五个没有进入通道保护区的残片，总体积约 7.25 mm³。对近端盖旁残留相交处做了局部修边，保留距通道中心线 0.72 mm 的保护范围。新网格回读为一个连通实体，与所查肘部构件和伸肌肌腹在中立位无超阈值重叠。沿原 K7 的 14 个姿态复算，三角肌与伸肌的最大名义重叠由 464.83 降至 257.20 mm³，伸肌近端盖由 249.15 降至 90.22 mm³，A 侧导轴由 210.68 降至 12.33 mm³。FJ1467 钢丝和变形肌腹的名义相交峰值仍为 1.78 mm³，需复查通道随动。当前保留此网格作为中立装配候选，K7 受载可用范围待确认。[重连后的晶格](K7_deltoid_reconnected.obj) · [逐姿态对照](whole_arm_k7_deltoid_reconnected_audit.json)

中束 FJ1467 的钢丝与晶格在中立位已有约 1.78 mm³ 交叠，因此将三角肌局部变形绑定到钢丝近邻的随动补偿后，交叠仍存在。按已选中心线，仅在相交处重开半径 0.52 mm 的孔道，共去除 5.42 mm³ 材料。导出回读保持连通，前、后束的通道保护区未损失材料。重新检查 14 个 K7 姿态，三条三角肌钢丝均未与插值后的三角肌实体发生超阈值相交。三角肌与肘伸肌的肌腹重叠仍待处理。[随动对照](whole_arm_k7_deltoid_channel_tracking.json) · [修孔后网格](K7_deltoid_reconnected_bore.obj) · [修孔后逐姿态检查](whole_arm_k7_deltoid_bore_audit.json)

运动让位版从 4 个 K7 姿态的相交体积估算中立位切削区域，并保留中心线半径 0.72 mm 的通道保护体积。初次全量切削形成 26 个实体。选用单块相交体积超过 2 mm³ 的区域生成切削后，另外清除约 87.90 mm³ 不与通道保护区相交的断片。总移除约 870.36 mm³，网格重新导入后为单个连通实体。

沿完整的 14 个既有姿态复算，三角肌与肘伸肌肌腹的最大名义重叠从 257.01 降到 14.90 mm³，与近端盖由 90.07 降到 0.43 mm³，与 A 导轴由 11.22 降到 1.65 mm³，与伸肌骨侧转接座由 41.97 降到 10.19 mm³。所查姿态中三角肌仍无钢丝实体侵入。几何插值、截面承载和运动中的柔性接触仍需复核。[运动让位网格](K7_deltoid_swept_relief.obj) · [生成参数](whole_arm_k7_deltoid_swept_relief.json) · [复算结果](whole_arm_k7_deltoid_swept_audit.json)

大圆肌 FJ1507 的通道在中立位未测得干涉。按两端直线轴向插值时，钢丝自 SC 1.33° 样点起穿入 TPU，最大名义重叠 14.28 mm³。改为让通道邻域跟随同一钢丝路径的局部随动分布后，重新计算全部 14 个姿态，未出现超阈值的 FJ1507 钢丝与 TPU 相交。大圆肌与自身近端盖的最大名义重叠从 59.77 降至 8.33 mm³，仍可能影响端部装配和实际活动。本轮沿用原有大圆肌网格，验证的是几何随动假设，尚需实测通道如何随 TPU 变形。[大圆肌通道随动对照](whole_arm_k7_teres_channel_tracking.json)

以大圆肌通道随动模型为基准，将四个姿态的端盖和肩胛接触区映射回中立位，局部去除约 36.50 mm³ 材料，保留沿原钢丝中心线半径 0.85 mm 的保护体积。新网格导出回读后为一个连通实体。重新检查全部 14 个姿态，近端盖名义重叠峰值由 8.33 降至约 0.0046 mm³，远端盖无超阈值相交，肩胛骨处仍记录约 0.0023 mm³。钢丝在所有样点均未与大圆肌实体相交。更新后的中立位 MJCF 引用新的大圆肌网格。实际加载及端盖装入仍需核对。[大圆肌让位网格](K7_teres_end_relief.obj) · [逐姿态对照](whole_arm_k7_teres_end_relief.json)

最新八块 TPU 网格也与两条肩袖的 20 件活动五金进行了组合空间检查。复用已核对的 10 N 局部加载状态，逐姿态放置端盖、转眼及压接件，覆盖 K7 的 14 个姿态。超过 0.001 mm³ 的相交集中于各自的肌腹—浮动端盖接口，肩胛下肌最大约 20.37 mm³，冈下肌最大约 13.85 mm³。其余肌腹与这些活动五金未发现超阈值相交。该结果来自肌腹端点插值，接口中的贴合、预压与装入净空仍需按连接结构区分，尚不能用这些体积判断实际运动阻力。[活动五金位姿](whole_arm_k7_moving_cuff_geometry.json) · [八块 TPU 逐姿态对照](whole_arm_k7_tpu_moving_cuff_audit.json)

### 三角肌导向复核

三角肌名称按输入模型标注。FJ1467 为肩峰部（中束），FJ1513 为肩胛冈部（后束）。历史脚本或文件名中的 posterior 指当时的走线标签，复现时以编号定位。

沿同一 31 姿态 K7 轨迹，复用前束导向的变形分布，同时在后束和中束近端保留 3 mm 固定区及 12 mm 平滑过渡。此前前束与屈肌支架、后束与伸肌 A 侧导轴的相交在本轮采样中消失。仅调整变形分布时，FJ1467 仍在 SC 2.67° 碰 A 侧导轴。[导向分布对照](whole_arm_documented_deltoid_audit.json)

将 FJ1467 远端节点 66–75 在最大处偏置 2 mm 后，八条被检查的钢丝共同达到 SC 4.33° 的样点，对应该分支约 17.33° 的肱骨方向变化。下一采样 SC 4.67° 时，FJ1467 首先碰肩胛骨，重叠约 0.00126 mm³。该方案保留原两端锚点，偏置段仍需配套导向、与 TPU 通道复核，并检查小曲率半径处的折线精度和实际可弯性。[FJ1467 偏置与八线复核](whole_arm_posterior_detour_audit.json)

已在既有三角肌晶格上按这条偏置线补做 FJ1467 内衬，外半径 0.95 mm，通道半径 0.52 mm。导出回读为单个连通实体，前束和 FJ1513 的通道保护区无材料损失。FJ1467 受检区段在 14 个姿态中均无钢丝穿入 TPU；中立位与肘伸肌近端盖有 0.0051 mm³ 重叠，其余 13 个样点未发现超阈值的局部刚体重叠。保留该中立位接触记录，完整肌腹变形和通道受载能力待验证。[FJ1467 TPU 模型](deltoid_posterior_detour_TPU.obj) · [局部通道检查](whole_arm_posterior_liner.json)

可用 `--routed-k7` 播放带底座的八线回放，范围自动截到所有纳入检查的钢丝共同通过的最后样点。该回放共 27 帧，显示局部走线。冈下肌与肩胛下肌的浮动端、偏置 TPU 的变形和受力需另行验证。

### 肘部转座与滑块随动回放

将已有肘部回放的屈曲、两端转座和滑块坐标映射至共享整臂，胸锁、肩锁及盂肱保持中立，并沿用与该轨迹匹配的加长导轴。203 个去回程样点中，移动刚性件、两条肘部钢丝和桌面底座无超阈值相交。肘屈曲峰值 126.01°，滑块行程峰值 59.25 mm，两条名义折线的峰峰线长变化分别为屈肌 68.74 mm、伸肌 40.67 mm。旧短导轴在 103.98° 样点已与滑块侧端盖相交。此结果来自对既有运动坐标的移植和几何复核，整体肩肘受力、回位弹簧净空和柔性接触仍待检查。[肘部随动记录](whole_arm_passive_elbow_audit.json)

在此回放中进一步装入回位弹簧的端钩、锁止件及两块既有肘部晶格的局部让位版。TPU 按两端转座与滑块的实际姿态插值，按源编号读入以保留通道拓扑。43 个样点包含屈曲峰值、滑块峰值、去回程及接触区加密。回位弹簧、锁止件、两条内线和两块 TPU 之间的所查配对均无超过 0.001 mm³ 的相交。另有 5 个回程样点出现屈肌外侧与肱骨相交，最大为屈肘 15.68° 时的 0.518 mm³，所查内线保有通道空间。该接触按插值外形登记，需在受载模型或样件中核验。MuJoCo 回放显示锁止五金，弹簧及变形晶格的结果在几何报告中查看。[弹簧、肘肌和内线检查](whole_arm_elbow_flexible_audit.json)

## 锚点与路线

按人体肌肉附着区各取中央代表点及沿走向的两侧点，共 72 个直拉组合。使用直径 0.60 mm 的钢丝包络对 98 个刚体求交，阈值为 0.001 mm³。表中的原始点来自解剖邻近区域，属于附着候选。

72 个未偏置的组合均在中立位发生实体重叠。对每组最接近通线的三个组合，筛查了长度 2、4、6 mm 的局部偏置。偏置后以直径 1.50 mm 的探针检查从骨面候选到出线点的接头轴线，保留所属骨段的重叠信息供实体设计。7 组得到可继续检查的候选。其坐标、所属骨段、偏置方向和全部冲突配对见[初筛数据](girdle_route_screen.json)。

下表对每条路线列出本轮绝对 SC 角最大的连续采样分支。行程为该分支的最大与最小名义线长之差；间隙测量去掉钢丝两端各 2 mm，以避免将端点座间隙混入线段间隙。各分支按原协调路径采样，并以相反顺序复查。

| 路线 | 局部偏置 | SC 局部轴及连续区间 | 名义行程 | 线段最小间隙 | 首个限制 |
|---|---:|---|---:|---:|---|
| FJ1456_g1 胸小肌 | 2 mm | 轴 0，0 至 −8.703° | 5.081 mm | 1.110 mm | 锁骨—肋骨 |
| FJ1459_g2 前锯肌 | 6 mm | 轴 1，0 至 −5.898° | 0.768 mm | 0.044 mm | 肩胛支承点离开垫面 |
| FJ1459_g3 前锯肌 | 2 mm | 轴 0，0 至 −8.703° | 3.597 mm | 0.744 mm | 锁骨—肋骨 |
| FJ1520_g1 斜方肌下部 | 6 mm | 轴 0，0 至 −8.000° | 4.997 mm | 0.021 mm | 钢丝与刚体相交 |
| FJ1536_g1 大菱形肌 | 2 mm | 轴 0，0 至 −8.703° | 1.025 mm | 0.869 mm | 锁骨—肋骨 |
| FJ1536_g2 大菱形肌 | 2 mm | 轴 0，0 至 −8.703° | 2.668 mm | 1.411 mm | 锁骨—肋骨 |
| FJ1460_g1 锁骨下肌 | 4 mm | 轴 0，0 至 −8.703° | 0.614 mm | 0.829 mm | 锁骨—肋骨 |

这些是规定 SC 和 AC 姿态下的几何结果。例：沿轴 0 负向的协调分支，前锯肌 FJ1459_g3 的钢丝缩短、FJ1536_g2 大菱形肌的钢丝伸长，可作为同一路径上的对拉候选。FJ1536_g1 的等效力臂沿该分支变号，暂不作为该分支的主动线。最终作用方向以及肩胛上回旋 K6 仍需按实际肩胛姿态确认。

FJ1459_g1 仍需要局部绕行或重新取点。上斜方肌、肩胛提肌和小菱形肌三条路径尚无颈段实体反力支撑。当前偏置端点是供设计安装座的位置，安装座实体、通道、转向摩擦、夹线、TPU 肌腹和多线组合间隙均待验证。

## 复现与查看

在项目根目录，使用已安装 MuJoCo、SciPy 和 Manifold3D 的 Python 环境。当前机器可用以下命令。

```powershell
$env:PYTHONPATH = '.deps'
& 'D:/MujocoEnv/Scripts/python.exe' tools/screen_whole_arm_girdle_routes.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_girdle_routes.py --route FJ1459_g3
& 'D:/MujocoEnv/Scripts/python.exe' tools/build_whole_arm_girdle_seats.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_interface_transfer.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/build_whole_arm_shared_scene.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/screen_whole_arm_shared_motion.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/measure_whole_arm_anchor_motion.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_local_wires.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_local_wires.py --localize
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_local_wires.py --documented-deltoid
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_local_wires.py --posterior-detour
& 'D:/MujocoEnv/Scripts/python.exe' tools/build_whole_arm_posterior_liner.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_floating_terminals.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_floating_terminals.py --offset-infraspinatus
& 'D:/MujocoEnv/Scripts/python.exe' tools/build_whole_arm_infraspinatus_support.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_floating_terminals.py --shared-infraspinatus-support
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_floating_terminals.py --combined-k7
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --deltoid-relief
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --reconnect-deltoid
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --track-deltoid-channels
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --repair-deltoid-bore
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --sweep-deltoid-relief --sweep-elbow-min-volume 2
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --verify-swept-deltoid
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --track-teres-channel
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --relieve-teres-ends
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_floating_terminals.py --combined-k7 --include-moving-hardware
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --check-moving-cuffs
& 'D:/MujocoEnv/Scripts/python.exe' tools/summarize_whole_arm_first_pass.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_gravity_balance.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_scapular_support.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_anatomical_actions.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_anatomical_actions.py --shoulder-elevation
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k1_shared.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k1_shared.py --cap-relief
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k1_combined_space.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k1_combined_space.py --spring-relief
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k1_combined_space.py --refine-spring-relief
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --check-moving-cuffs --k1-spring-relief
& 'D:/MujocoEnv/Scripts/python.exe' tools/build_whole_arm_first_pass_assembly.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --check-moving-cuffs --k1-spring-relief --render-replay
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --first-pass
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --base --cuff-support
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --base --motion K7_SC1_pos --routed-k7 --combined-cuffs
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_local_wires.py --passive-elbow
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_elbow_flexibles.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/build_whole_arm_table_base.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/build_whole_arm_girdle_tpu.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --latest-tpu
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --latest-tpu --write-assembly
& 'D:/MujocoEnv/Scripts/python.exe' tools/audit_whole_arm_k7_tpu.py --check-moving-cuffs --render-replay
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --base --girdle-tpu
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --motion K2_GH_flexion
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --shoulder-action K2_GH_flexion
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --shoulder-action K3_GH_abduction
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --scapular-action protraction
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --scapular-action retraction_alternative
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --scapular-action upward_rotation
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --base --motion K7_SC1_pos
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --base --motion K7_SC1_pos --routed-k7
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --base --passive-elbow
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --shared-elbow
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --shared-elbow --cap-relief
& 'D:/MujocoEnv/Scripts/python.exe' tools/view_whole_arm_shared_scene.py --motion K7_SC1_pos --routed-k7 --combined-cuffs --cap-relief
```

查看器可以围绕整臂旋转，橙线表示当前选中路线。加 `--check` 可逐帧复算钢丝长度。

录制 K7 双视角视频还需要在用于运行脚本的 Python 环境中提供 Pillow 和 imageio-ffmpeg，以及可用于 MuJoCo 离屏渲染的 OpenGL 环境。
