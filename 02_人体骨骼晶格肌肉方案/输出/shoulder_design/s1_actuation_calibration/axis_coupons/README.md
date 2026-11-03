# 晶格方向与连续通道对照样件

2026-11-02。几何样件已生成并回读；材料、打印及力／位移数据均未实测。这是推进单肌肉标定的准备，不代表肩部或整臂驱动通过。

## 要解决的问题

AIS2026 §2.3 将 x 定义为外部压缩方向、y 为收缩方向、z 为非承载单元尺寸。现有肩部生成器将 z 沿肌腹中心线布置，不能直接认为实现了论文的 y 向行为。论文还明确采用实际样件的经验映射，因此没有把其收缩比例、材料模量或旧 25% 假设写成当前肌肉能力。

本组有意比较 Y 和当前 Z 中心线候选；它不是完整的论文 Ex／Ey 重现实验，也不是作者原始 CAD。核心采用项目既有的体心到八角分支单元，名义杆径 2 mm。端部、连续通道、外表面及整块肌腹形状都会影响响应，不能把裸格结果直接转移到完整肌肉。

## 独立两端夹持

每一对使用相同 2×2×2 核心、相同 XYZ 摆放和打印轴，仅在受测方向的两个相对端面加反力板。没有另外两面的连续侧板跨接加载端。

旧 z 标定件的四面夹持板可能提供旁路载荷并限制横向变形，不能将其双向测试结果直接解释为未约束晶格的各向异性。这组样件单独存放，旧几何与记录保留；旧 D7 刚度比的预测窗口不作为当前实测通过条件。

P0 裸格使用平行压板压缩；P1 在同一核心上增加当前候选连续通道（外半径 1.8 mm、孔半径 0.86 mm），并给出实际压力片、未压接止挡和 0.6 mm 钢丝路由。P1 先验证手拉压缩、通道和回弹，不能直接把点加载压力片下的刚度当作材料弹性模量。

| 单元 X/Y/Z | 加载轴 | 核心名义标距 | 核心名义截面积 | 几何版本 |
|---|---|---:|---:|---|
| 15/10/30 mm | Y | 20 mm | 1800 mm² | P0 裸格、P1 通道 |
| 15/10/30 mm | Z | 60 mm | 600 mm² | P0 裸格、P1 通道 |
| 10/8/10 mm | Y | 16 mm | 400 mm² | P0 裸格、P1 通道 |
| 10/8/10 mm | Z | 20 mm | 320 mm² | P0 裸格、P1 通道 |

由于标距和截面积不同，不能直接比较同一位移下的力。P0 可登记名义应变 δ/L 与名义应力 F/A，注明边界及板摩擦后报告样件等效响应；不是 TPU 本体模量。P1 同时测输入钢丝位移与实际板间缩短量，不预设两者严格 1:1。

## 数字检查与制造边界

8 个几何版本共 20 个实际对象，全部 STL 回读满足单正材料实体、零边界及非流形边。P1 各件对之间未检出超过 0.001 mm³ 的材料相交，并记录止挡的有限几何保持检查。它不证明压接握持或承载。

[MuJoCo 静态对照图](comparison.png)显示实际网格，DOF=0，没有绘制虚构力曲线或变形。红色为参考尺寸组，蓝色为当前尺寸组；每组依次为 Y裸格／Y通道／Z裸格／Z通道。

TPU95A 仅为知识库给出的材料起始候选。实际材料、打印机、喷嘴、切片、支撑和工艺方向尚未确定；本文不指定论文未提供的温度／速度。预定同一 XYZ 方向打印，若因悬垂改变摆放，须记录实际变换并使成对样件保持一致。几何回读不等于已经切片或免支撑。

钢丝路由 `nominal_wire_NOT_print.stl` 不是 TPU 打印件；压力片与止挡是独立刚性／金属候选，毛坯止挡未经压接不能当作固定终端。P1 压力片依靠钢丝预张力压住远端板，仅用于压缩／恢复试验，不是完整拉伸样件。

## 测试记录

建议每类 3 份独立试样，共 24 份，使用同一材料批次和工艺。`measurements_BLANK.json` 提供 24 条全空记录，重复生成不覆盖已有文件。

1. 记录真实打印／干燥／支撑配置、杆径、标距、孔径与自然长度；先检查孔是否通线。
2. P0 使用只有两相对端面接触的平行压板，记录实际位移和可获得的力；P1 固定近端反力板、安装止挡及压力片，分别记录拉线和板间位移。
3. 从小位移开始，每步登记力、缩短、卡线及端部状态；发生卡线、端部损伤或不可恢复变形时停止该件并保留失败记录。
4. 做重复往返，登记卸载、恢复和 60 秒残余；没有力传感器时力保持 null，不用手感或电机命令代填。
5. 先区分方向、通道壁和夹持端部的影响，再决定真实肌腹参数；仍需在异形肌腹、壳网和实际端部上验证。

## 生成器方向修正

`build_shoulder_tpu.py` 的结构网格新增 `--path-axis y`：以单元 Y 为目标间距沿中心线均匀排站，横向用 X/Z；总路径长度会重新分配实际站距；论文局部基为 `[transport_e1, tangent, -transport_e2]`，保持右手系。默认 Z 仍是历史行为，没有自动旋转或替换现有 72 对象。

受检 60 mm 直线路径上，15/10/30 单元的 Y 模式形成 6 个轴向单元，Z 模式为 2 个；实际横向尺寸亦随映射变化。当前默认输出与整理前脚本的受检几何对称差为 0，证据见 [方向比较](../../../../../reviews/S1_2026-11-02/lattice_axis_mapping/comparison.json)。两个相关回归测试通过。

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) '.deps')
python mujoco_arm/tools/build_shoulder_axis_coupons.py
python mujoco_arm/tools/render_shoulder_axis_coupons.py
python mujoco_arm/tools/render_shoulder_axis_coupons.py --view
python -m unittest discover -s tests -p test_shoulder_lattice_axis.py -v
```

基础异形晶格的 Y 映射实验入口如下，必须选择独立输出目录：

```powershell
python mujoco_arm/tools/build_shoulder_tpu.py --module infraspinatus --mode structured --path-axis y --cell 15 10 30 --out '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_actuation_calibration/y_base_candidate' --json '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_actuation_calibration/y_base_candidate/report.json'
```

这条命令只生成基础异形晶格候选，不携带 72 对象候选的新端盖、卡扣和肌腱；未在本轮执行整肌腹替换，不能直接拿其 STL 覆盖已检查的装配。

## 证据与下一步

- [样件清单、源参数及网格哈希](manifest.json)
- [空白测量记录](measurements_BLANK.json)
- [静态查看检查](view_check.json)
- [近端反力轴线审查](../../s1_retained_joint_comparison/shaft_collar/station_60/origin_fixation_review/README.md)
- [肩部技术方案](../../SHOULDER_DESIGN.md)
- [打印参考](../../../../../design/TPU_PRINTING_REFERENCE.md)

近端前支承及背侧分压板已形成[92对象名义配合候选](../../s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/README.md)，未验证实物夹持或承载。下一步以实际样件确定晶格方向、通道影响与工艺参数，结合台架及自由末端约束推进；本组样件没有替换最新整肌腹。全部实测项保持未知，`stage_pass=false`。
