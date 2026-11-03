# 骨干侧肌肉接口比较：两片式抱箍

2026-11-02。依据当前 S1 技术方案继续比较“骨干侧接口＋肌腱”的可能性。此处交付接口候选及配合几何，不是已连通的新肌肉驱动装配，也未替代 S1 基线。实际驱动范围未知。

## 接口结构

以独立 GH 候选的真实工程肱骨网格为基准，沿已标定肱骨轴比较距 GH 40、50、60 mm 三个站位。每个站位包含 12 个实际独立实体：两半刚性夹体、两半 TPU 衬垫、两枚名义 M3×16 螺栓、两枚螺母及四片垫圈。

夹体轴向长 10 mm；按 0.5 mm 站距截取该段骨截面、求联合凸包并外偏 1 mm 形成衬垫外包络，再外偏 2.5 mm 形成刚性壁。两半分缝为 1 mm，夹紧孔半径 1.7 mm，肌腱眼孔半径 1.1 mm。衬垫内面按真实骨体布尔扣除；这是名义贴骨接触，不是正的制造间隙。最终实际网格仍逐件对骨检查，截面采样本身不充当全高度证明。

夹紧桥最初伸过垫圈座，造成约 1–1.5 mm³ 相交，已在实际夹体上加工座面让位，没有排除五金碰撞。初版保存在 `iterations/initial_station50/`。40 mm 站位下半夹体首次回读有两条非流形边，原件保留；0.001 mm 有记录简化后完成回读。所有站位最终 12 件均为单正实体、边界与非流形边为 0。

夹紧力、接触摩擦、衬垫压缩、防转容量及打印强度尚未验证。名义螺栓／螺母包络不含实际螺纹啮合。

## 三站位比较

| 站位 | 冈下肌端面至眼孔直线距离 | 肩胛下肌端面至眼孔直线距离 | 与现有刚性件的零位冲突 | 结论 |
|---|---:|---:|---|---|
| 40 mm | 42.39 mm | 38.84 mm | 夹体／五金与大圆肌远端盖；最大约 10.52 mm³ | 不纳入当前配合候选 |
| 50 mm | 52.26 mm | 48.77 mm | 同上；最大约 15.49 mm³ | 不纳入当前配合候选 |
| 60 mm | 62.15 mm | 58.69 mm | 未检出超过 0.001 mm³ 的刚性相交 | 继续做局部配合比较，未采纳为新 S1 |

40 mm 站位衬垫／骨还存在约 0.001384 mm³ 贴骨残差，未转为制造通过。上述长度只是旧肌腹端面到新眼孔的直线参考距离，绕行及终端会增加长度；不能当作已经设计好的“短肌腱”。

把两条力线端点直接替换为眼孔的零位直线代理均保持正张力方向覆盖，但归一化正余量从原约 0.005683 降至约 0.002208／0.001748／0.001601。此指标不是力容量，也不是实物必需张力比；真实的肌腹端面、肌腱与接触传力尚未建立。

## 接口自身的条件运动检查

1° 常规步长、边界 0.1°，相交阈值 0.001 mm³；间隙查询封顶 2 mm。以下仅检查新增 12 件接口与固定肩胛／锁骨／GH 附件等的相对运动；假设接口及未变形衬垫理想随肱骨运动。未验证夹紧是否足以维持该约束。

| 站位 | 屈曲＋／－ | 外展＋／－ | 长轴＋／－ |
|---|---|---|---|
| 40 mm | +120°／−84.2° | +120°／−42.3° | +120°／−120° |
| 50 mm | +120°／−100.6° | +120°／−90.1° | +120°／−120° |
| 60 mm | +120°／−120° | +120°／−109.8° | +120°／−120° |

±120° 为检查上限，不是求得的极限。该表不包含原远端盖等其他活动件的首限，不能用作整肩运动通过结论；完整配合组结果见下节。

## 60 mm 站位的真实 TPU 配合

新增接口占用了一部分原三角肌空间。按实际 12 件接口及源坐标三轴正负 0.2 mm 平移包络局部让位；这一构造不保证任意方向都有 0.2 mm 间隙。

直接让位产生约 17.876 mm³ 的断连晶格，未删除该功能片段。最初仅搜索最近切面端点找不到有效补肋；改为先筛选能容纳完整 R1 肋的自由端点，再用约 3.27 mm 的真实晶格连接肋接回主体，增加约 12.213 mm³ 材料。该肋不穿过骨、接口或钢丝保护区。

重新导出曾在远处出现 0.026 mm 长的四面共边，定位为改变 STL 局部原点导致旧顶点重新量化；沿用父模型局部原点后，实际回读通过：单正材料实体、14 个封闭负壳保留，边界／非流形均为 0。负壳表示数字空腔，不表示打印可行性已经确认。

| 联合检查 | 结果 | 仍未解决 |
|---|---|---|
| 46 件零位配合 | 原 34 件参考结构＋12 件接口 | 原远端盖仍保留，未连接新肌腱 |
| 三角肌修改 | 移除约 212.911 mm³，原体外增加约 12.214 mm³ | 刚度、恢复、打印及接触变形未验 |
| 六路旧钢丝及 1.61 mm 通道保护区 | 与接口让位区无超阈值冲突 | 不是新肌腱或新驱动路径 |
| 六路钢丝对 46 件实体 | 未检出超过 0.001 mm³ 的相交 | 弯曲半径、摩擦与运动未验 |
| 大于 0.05 mm³ 的剩余相交 | 三角肌／工程肩胛约 1.84044 mm³ | 原有柔性接触，未判通过 |

46 件重新进行条件刚性路径检查，把 `collar:` 组件明确分配给肱骨活动组；三角肌等柔性肌腹不作为刚体扫掠。结果仍为屈曲 +114.2°／−64.2°、外展 +61.1°／−42.6°、长轴 +68.1°／−47.0°。原远端盖未被删除，因此主要首限不变；不存在由本轮接口获得的整肩范围提升。

## 力链尚缺的部分

眼孔是骨侧固定接口候选；还没有肌腹远端的独立肌腱终端、有效肌腱路径及受力平衡。端面到眼孔的参考直线会穿过现有端盖或眼孔边缘，不是可直接穿线的设计。

后续必须联合设计肌腹形状与终端，使内部钢丝止挡载荷通过真实肌腱传给接口。自由肌腹端盖不能凭空固定在世界坐标，也不能把它继续当作刚性骨附件来证明新布局的驱动范围。需要重新核对论文晶格工作方向、有效行程、接触与力臂。

## 复现与证据

在项目根目录、已安装运行依赖的环境中执行：

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) '.deps')
python mujoco_arm/tools/build_shoulder_shaft_collar.py --station 40
python mujoco_arm/tools/build_shoulder_shaft_collar.py --station 50
python mujoco_arm/tools/build_shoulder_shaft_collar.py --station 60
python mujoco_arm/tools/check_shoulder_shaft_collar.py
python mujoco_arm/tools/fit_shoulder_deltoid_to_collar.py --station 60
python mujoco_arm/tools/check_shoulder_rigid_motion.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/deltoid_relief/assembly_check.json' --moving-prefix 'collar:'
python mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/deltoid_relief/assembly_check.json'
```

- 站位比较：[40 mm](station_40/interface_check.json)、[50 mm](station_50/interface_check.json)、[60 mm](station_60/interface_check.json)。
- 接口运动／接触：[40 mm](station_40/motion_integration_check.json)、[50 mm](station_50/motion_integration_check.json)、[60 mm](station_60/motion_integration_check.json)。
- [46 件联合配合](station_60/deltoid_relief/assembly_check.json)、[重新计算的刚性路径](station_60/deltoid_relief/rigid_motion_check.json)。
- [完整参考结构静态图](station_60/deltoid_relief/available_parts.png)；[只显示接口及骨／GH 上下文的图](station_60/available_parts.png)。均为 DOF=0，不是变形或驱动演示。
- [原 GH 独立候选](../GH_only_candidate/README.md)、[肩部技术审查稿](../../SHOULDER_DESIGN.md)。

历史拒绝结果和早期几何位于 `iterations/` 及配合分支的 `iterations/`。当前以 `assembly_check.json` 和其输入哈希为证据。所有报告仍为 `stage_pass=false`，整臂目标继续推进。
