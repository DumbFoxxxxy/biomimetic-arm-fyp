# 仿生右臂 FYP 项目交接说明

更新时间：2026-10-26

## 1. 项目当前目标

本项目是本科 FYP 的右臂仿生桌面台架样机。当前阶段的目标不是完成最终实物，
而是：

- 用 SolidWorks 维持一套可继续修改的零位骨骼装配；
- 用 MuJoCo 展示 7 个受控自由度、16 根钢丝的拉动效果；
- 向导师证明拉动不同钢丝可以得到 K1-K7、肩胛滑动/旋转和胸大肌辅助动作；
- 用仿真结果辅助下一版 SolidWorks 接触结构设计。

当前不建模手、腕、桡骨和旋前旋后；暂不考虑 TPU 力、重力、钢丝柔性、摩擦细节
和实物球铰防脱。丝杠与电机位于手臂外，仿真中用理想拉线执行器代替。

## 2. 已确定的设计基准

坐标系：胸廓底面中轴为原点，`x` 向右、`y` 向上、`z` 向后，单位 mm。

```text
SC  = (40.00, 235.00, -106.00)
AC  = (156.86, 233.24, -11.97)
GH  = (159.06, 187.54, 7.30)
E   = (184.34, -101.35, 7.30)
P   = (62.18, 210.10, 91.39)
SUP = (72.76, 268.45, 82.21)
INF = (44.55, 112.87, 106.69)
```

主要尺寸：

| 项目 | 当前值 |
| --- | --- |
| 人体尺度 | 男 18-60 岁第 50 百分位，身高 1678 mm，1:1 |
| 胸廓 | 280 宽 x 212 深 x 280 高；R_pad 约 400 mm |
| 锁骨中心距 SC-AC | 150 mm |
| 肱骨 GH-E | 290 mm |
| 尺骨 E 到远端 | 237 mm |
| 肩胛板厚 | 5 mm |
| 肩胛关节盂平垫 | 40 x 40 mm |
| 肩峰端块 | 24 x 24 x 16 mm，叶片厚 5 mm |

肩胛骨的五个建模锚点为 SUP、INF、GH、AC、P。

## 3. 当前机构链

SolidWorks 和 MuJoCo 均采用简化的 7 自由度模型：

```text
胸廓
 -> SC 抬降
 -> SC 前后移
 -> 锁骨
 -> AC 内外旋
 -> 肩胛骨
 -> GH 屈伸
 -> GH 外展内收
 -> GH 内外旋
 -> 肱骨
 -> 肘屈伸
 -> 尺骨
```

受控自由度：

```text
肩胛滑动 + 肩胛上下回旋 + 肩胛内外旋
+ GH 屈伸 + GH 外展内收 + GH 内外旋
+ 肘屈伸
= 7 DOF
```

需要注意：

- SC 在 CAD 中是两轴十字轴万向节，在 MuJoCo 中由两个明确铰链表示；
- AC 在当前 MuJoCo 中不是完整球铰，而是一个 `ac_ext` 铰链；
- GH 用三个串联铰链模拟球铰；
- 当前 DAE/STL 视觉网格来自 SolidWorks 零件；
- 锁骨和肩胛骨的完整网格仍作为碰撞体，这是当前主要遗留问题。

当前 MuJoCo 关节限位：

```text
sc_elevation   -5..30 deg
sc_protraction -15..15 deg
ac_ext         -15..15 deg
gh_flex        -150..150 deg
gh_abd         -150..150 deg
gh_rot         -120..120 deg
elbow           0..145 deg
```

## 4. 16 根钢丝

原方案为 14 根，后来增加两根胸大肌钢丝，当前共 16 根：

```text
E01/E02   肘屈伸
G01/G02   GH 屈伸
G03/G04   GH 外展内收
G05/G06   GH 内外旋
S01/S02   肩胛上下回旋
S03/S04   肩胛前后滑动
S05/S06   肩胛内外旋
P01/P02   胸大肌锁骨部 / 胸肋部
```

权威锚点表：

```text
design/CABLE_ANCHORS_FINAL.md
tmp/cable_anchors_final.tsv
tmp/cable_anchors_final.csv
tmp/cable_anchors_final.py
```

当前理解必须保留：肌肉模拟不限定为“一块晶格 + 一根钢丝”。同一块扇形晶格可以由
两根钢丝控制，也可以使用 Y 形晶格和两根钢丝。因此钢丝走线数不等于肌肉块数，
也不等于自由度数。

P01/P02 当前起点在锁骨或胸骨前侧，胸前外侧导点固定在胸廓，止点在肱骨近端；
两根同时拉动主要产生肩前屈、内收和内旋辅助。

## 5. SolidWorks 当前状态

当前主文件：

```text
final_model/Thorax.SLDPRT
final_model/clavicle.SLDPRT
final_model/scapula.SLDPRT
final_model/humerus.SLDPRT
final_model/ulna.SLDPRT
final_model/universal_joint_pin.SLDPRT
final_model/assembly_zero.SLDASM
final_model/assembly_motion_gh.SLDASM
```

`assembly_zero.SLDASM` 是当前干净零位装配，**不要覆盖**。下一版另存为：

```text
final_model/assembly_zero_v2.SLDASM
```

历史装配、旧 14 路径、早期 Motion 实验和损坏副本在：

```text
final_model/archive/
```

不要恢复 `assembly_zero_fixed.SLDASM` 或其他历史损坏装配。

注意：`design/DESIGN_BRIEF.md` 的“当前版本速查”已经把 AC/GH 球铰设计值更新为
R10 级，并给出 GH 浅球窝 R10.2；但当前 CAD 中仍有较早的 R8 球头/球窝。下一版
CAD 应统一球径，不能直接把文档数值当作已建成的实物尺寸。

## 6. MuJoCo 当前状态

主文件：

```text
mujoco_arm/arm_tendon_demo.xml
mujoco_arm/config/tendon_demo_16.json
mujoco_arm/build_tendon_demo_model.py
mujoco_arm/check_tendon_demo.py
mujoco_arm/run_tendon_demo_gui.py
mujoco_arm/run_tendon_demo.py
mujoco_arm/render_tendon_demo.py
mujoco_arm/prepare_tendon_demo_meshes.py
mujoco_arm/export_tendon_demo_meshes.ps1
```

启动中文交互界面：

```text
mujoco_arm/open_tendon_demo_gui.bat
```

当前模型有：

```text
nq = 7
nu = 16
ntendon = 16
```

16 个执行器都是 `position` 型理想拉线执行器。减小目标长度等于拉入钢丝；
手臂不是由关节角马达直接摆姿势。GUI 支持中文钢丝名称、拉动毫米滑块、
姿势预设和当前预设状态显示。

## 7. 最近一次验证

以下命令在当前环境通过：

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\build_tendon_demo_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_tendon_demo.py
```

最近结果：

```text
PASS: all demo groups produce finite, visible motion.
```

K1-K7、滑动、肩胛旋转、P01/P02 都能产生可见运动。当前检查报告中的最大接触
穿透约 0.32 mm，但该数字来自完整网格碰撞代理，只能说明仿真没有明显数值失稳，
不能证明肩胛胸壁接触结构正确。

报告位置：

```text
mujoco_arm/exports/tendon_demo/check_report.json
```

## 8. 当前首要问题

K6/K7 时肩胛骨仍可能出现视觉上的胸廓穿透，或肩胛冈、背面锚定凸台先于预期的
肩胛胸壁接触面发生碰撞。根因不是简单的铰链限位，而是当前 MuJoCo 使用了
肩胛骨和锁骨的完整网格作为碰撞体：

```text
thorax_collision   椭球包络
clavicle_collision 完整网格
scapula_collision  完整网格
```

紧限位只是暂时掩盖问题。正确方向是让胸廓上的 R400 接触带与肩胛骨前面的独立
接触垫参与碰撞，完整骨骼网格只负责视觉。

## 9. 下一版 SolidWorks 必须完成的结构

### 9.1 肩胛骨接触垫

新增独立、可导出的接触几何：

```text
SCAP_CONTACT_PAD
约 50 x 130 mm，厚度 3-5 mm
边缘 R3-R5
与胸廓 R400 面在零位留 0.3-0.6 mm 间隙
```

要求：

- 接触垫位于肩胛骨前面；
- 不让肩胛冈、锚定平垫、螺钉孔、TPU 晶格安装区进入接触垫；
- 接触垫可以局部加厚，但不得改变已有五个锚点坐标；
- 接触垫应能单独导出，不依赖整块肩胛骨网格。

### 9.2 胸廓接触带

在胸廓背外侧建立明确接触几何：

```text
THORAX_SCAPULA_CONTACT_BAND
R400 单半径弧面
x 约 55-140 mm
y 覆盖肩胛骨实际滑动范围
```

接触带与外观壳体分开导出。不要继续用椭球包络代替真实接触面。

### 9.3 网格导出

下一版至少需要：

```text
thorax_visual.stl
scapula_visual.stl
thorax_contact_band.stl
scapula_contact_pad.stl
```

可选保留锁骨视觉网格，但锁骨碰撞不应再与胸廓完整网格默认互相碰撞。

建议新建 `assembly_zero_v2.SLDASM`，保持零位坐标不变，再导出网格。

## 10. 下一版 MuJoCo 修改顺序

1. 去掉完整 `scapula_collision` 与 `clavicle_collision` 的默认碰撞；
2. 加入 `thorax_contact_band` 和 `scapula_contact_pad`；
3. 先检查零位间隙，目标 0.3-0.6 mm；
4. 再检查 K6/K7 最大穿透，目标不超过 0.5 mm；
5. 确认肩胛骨背面特征不再参与胸廓碰撞；
6. 接触通过后，再重新核调 S01-S06 的锚点和走线；
7. 最后调整 GUI 预设和渲染视角，避免把透明胸廓外观误判为实际穿透。

验收标准：

```text
零位间隙                 0.3-0.6 mm
K6/K7 最大穿透           <= 0.5 mm
完整肩胛骨网格            不参与碰撞
所有 16 根钢丝             仍能产生有限、可见运动
```

## 11. 重要文档

当前应优先阅读：

```text
PROJECT_STRUCTURE.md
design/DESIGN_BRIEF.md
design/CABLE_ANCHORS_FINAL.md
design/MOTION_IMPLEMENTATION_ORDER.md
design/CABLE_ROUTING.md
design/IDEAS_backlog.md
mujoco_arm/TENDON_DEMO.md
```

`DESIGN_BRIEF.md` 很长，内部仍保留历史讨论和已作废表述。发生冲突时按以下优先级：

```text
当前模型文件
> DESIGN_BRIEF.md 开头的“当前版本速查”
> CABLE_ANCHORS_FINAL.md 的当前锚点表
> 正文中的历史段落
```

历史文档：

```text
design/archive/
archive/
tmp/archive/
final_model/archive/
```

## 12. 工作约束

- 用户希望以中文交流；
- 当前重点是 SolidWorks + MuJoCo，不扩展到 Unity；
- 不修改原始参考 PDF、知识库和参考骨骼；
- 项目目录不是 Git 仓库，没有版本控制兜底；
- 新设计结论先写入 `design/`，再改模型；
- 新 SolidWorks 版本不要覆盖已验收的 `assembly_zero.SLDASM`；
- 临时结果放 `tmp/`，MuJoCo 实验脚本放 `mujoco_arm/tools/`，结果放
  `mujoco_arm/exports/`。
