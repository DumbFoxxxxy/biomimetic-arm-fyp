---
name: muscle-to-tpu
description: Converts an anatomical muscle mesh (OBJ) into a printable TPU lattice muscle — moment-sign precheck, x/y/z/p lattice parameter selection, lattice plus continuous wire-channel generation, end-cap interfaces, and mandatory validation (connected components, wire-material clearance, volume roundtrip, print placement). Use when converting a human muscle model into a TPU actuator, adding a muscle to the musculoskeletal arm, or when the user mentions TPU 肌肉、晶格肌肉、PSM、钢丝通道、muscle module、muscle-to-TPU.
disable-model-invocation: true
---

# 人体肌肉 → TPU 晶格肌肉

把一块解剖肌肉网格（OBJ）变成可打印、可穿线、可受力的 TPU 晶格肌肉。
本项目的实现位于 `mujoco_arm/anatomical_concept/`，本 skill 固化的是**流程、参数语义与判据**，
不是重写代码。

## 不可协商的四条

先读这四条，它们来自肘部 11 轮迭代的真实返工。

### 1. 力矩符号必须在建模之前算

**不要先做几何、后查力矩。** 肘部 M1 漏了这一步，到 M5/M9 才发现：
伸肌力矩在 85° 反号、直线弦在 132° 碰肱骨，整轮重做。

对每一条力线，在目标角度区间上算：

```
moment(θ) = (B(θ) − C) × unit(A − B(θ)) · a
```

`C` 旋转中心、`a` 转轴、`A` 近端固定点、`B` 远端随动点。
**判据：区间内符号恒定，且与期望方向一致。** 中途过零即停止建模，改起止点或力线。

同时检查**直线弦是否穿骨**——若穿，说明实际走线必须绕行，弦长估算不能当行程。

### 2. 外形照抄，机构设计

外轮廓应当贴合解剖（肘部实测包络面 96%+ 落在原始肌肉体积内），但下面四项**必须工程化**：

| 保留 | 让位 |
|---|---|
| 肌腹外轮廓 | 端部（换刚性端盖 + 压接座） |
| 起止**位置**（真实骨面） | 通道走向（工程路径，见 Step 3，**不是肌纤维**） |
| 整体尺度与体积量级 | 与骨/相邻件的局部让位 |
| | 晶格密度分布（按需要的软硬） |

**薄处不要硬套**：单元尺寸 15×20×10、通道外径 3.2 mm，肌腹局部厚度不足 8 mm 就放不下，
应简化而不是硬塞。肩胛下肌、冈下肌这类贴骨薄肌要让位。

### 3. 参数语义 x / y / z / p（禁止按位置对表）

| 参数 | 含义 | 本项目取值 |
|---|---|---|
| **x** | **受压**方向（外力/接触）单元尺寸 | 15 mm |
| **y** | **收缩**方向（肌腹长轴）单元尺寸 | 10 mm |
| **z** | **非承载**方向单元尺寸 | 按标定曲线查（不要写死 20） |
| **p** | 杆直径 | 2 mm |

**代码里的变量名必须直接写语义**（`x_compression_mm` / `y_contraction_mm` / `z_free_mm`），
不要用 `cell_spacing_transverse/axial` 这类方位命名——肘部就是因此差点把 x 和 z 对调。

装到骨上之前必须确认**哪一轴真正受压**，并把它设为 x。放错等于把 x/z 互换。

### 4. 结论边界

以下关系**不成立**，不要写进任何结论：

- 单一连通实体 + 钢丝不相交 ≠ 能驱动
- 网格回读通过 ≠ 能打印
- 刚性路径无干涉 ≠ 肌肉能达到该角度

几何检查只能写成"几何检查通过"。

---

## 工作流

### Step 0 — 前置判断（不产出几何）

**已实现为脚本，不要手写**：

```bash
python mujoco_arm/tools/frames.py                                  # 坐标系自检（dE=0, dF=0）
python mujoco_arm/tools/force_line.py --self-test                  # 解析解自检
python mujoco_arm/tools/force_line.py --muscle FJ1486 --to 150 --step 1 --bone
```

| 输入 | 输出 | 判据 |
|---|---|---|
| 关节中心 C、轴 a、起止点 A/B | 力矩符号表（每角度）、符号翻转点、弦长变化 | 区间内符号恒定 |
| 骨骼网格 | 直线弦首次接触角度、最小间隙 | 目标区间内不接触，或已知接触角 |

**失败即停**，不要带着反号的力线进入建模。

`frames.py` 是坐标换算的唯一入口（`SHIFT / E / F`），不要再在脚本里手写那三行——
少减一个 `SHIFT` 会让整块肌肉落到 1000 mm 之外。

### Step 1 — 参数选择

1. `p = 2 mm`（可打印性与柔性的平衡点，论文实测 p 1.5→3 刚度差 28 倍）。
2. `x`、`y` 取 15 / 10，分别对齐受压轴与收缩轴。
3. `z` 由 z 标定曲线按目标刚度反查；**无曲线时不得外推**，只能取曲线已覆盖区间内的值。

### Step 2 — 几何生成

断面提取 → 周期 B 样条包络 → 八杆中心—顶点晶格 → 表面菱形杆网。

**八杆单元**：从单元中心向 8 个长方体的顶点各伸一根杆；在任一正交中截面上读到 4 根、
呈 X 形。这就是论文的 unit element，不是另一种拓扑。

**判据**：包络面 ≥95% 落在原始解剖体积内（允许局部超出，但不能脱形）。

### Step 3 — 通道与端部

- 通道用**连续管壁**（R0.85 内径 / R1.6 外径 / 0.75 壁厚），明显硬于周围稀疏晶格，
  否则变形时会被压扁。
- 端部沿用肘部 M10-R1 形式（刚性端盖 + 压接止挡座），只按尺寸重算。
- 多通道肌肉（如合并的三角肌）用 Y 形分叉：一端单点汇聚、一端分叉，避开"多块解剖体能否布尔"的问题。

**判据**：通道壁与晶格在同一个连通实体里。

#### 通道路径是怎么来的（**当前不是按肌纤维做的**）

| 环节 | 现在的确定方式 | 算"按发力"吗 |
|---|---|---|
| 两端附着点 | 扫相位 / 轴向 / 径向，取力矩符号正确且力臂足够大的解 | ✅ 是 |
| 中间路径 | 各断面极点上的**固定端最短步进路径**（`route_engineering_poles`） | ❌ 纯几何 |
| 导引段 | 沿入/出切线角平分线，半径按转角放大 `0.85 / cos α` | 几何 |
| 肌纤维方向 | **无数据**（源模型只有表面；`fiber_direction_verified: False`） | — |

代码原文注释：`Shortest-step interior path with fixed end poles; NOT anatomical fibers.`

**所以：要真正做到"按发力设计通道"，必须把"要求的作用线"作为约束加进路径求解**——

- 目标：偏离要求作用线最小 + 转弯次数/角度最小；
- 约束：全程在肌腹内、通道间隙 ≥ 钢丝半径、曲率满足钢丝可弯半径。

目前只对**端点**做了这件事，中段仍是几何解。**不要**用 PCA 主轴或几何中心线冒充纤维方向，
项目原则明令禁止；也不要把"端点选对了"说成"通道按发力设计"。

### Step 4 — 强制校验

```bash
python .cursor/skills/muscle-to-tpu/scripts/check_tpu_muscle.py <肌肉目录> --json report.json
```

四项全过才算通过，缺一项即停：

| 检查 | 判据 |
|---|---|
| 连通分量 | 主分量唯一；额外分量若 ≤32 面且 \|体积\| < 1e-3 mm³ 记为**碎屑**（报告但不判失败） |
| 主分量闭合性 | 每条无向边恰好 2 面（边界边 0、非流形边 0） |
| 钢丝—材料最小间隙 | > 钢丝半径（0.6 mm 钢丝 → 0.3 mm） |
| OBJ/STL 相对体积误差 | < 1e-6 |

**默认不做顶点焊接**（`--weld-tol 0`）：manifold3d 导出的 OBJ 已共享顶点索引，
再焊接会把相邻面错误合并，凭空造出非流形边（实测 FJ1477 会从 0 条变成 29 条）。
只有当外部来源的 OBJ 存在重复接缝顶点时才开启焊接。

### Step 5 — 交付与记录

```
<muscle_dir>/
├── TPU_muscle.obj / .stl        # 制造网格
├── continuous_channel_wall.obj  # 通道壁单独导出，供检查
├── design.json                  # 参数、坐标、线、判据
└── print_candidate/
    ├── TPU_muscle.obj / .stl    # 打印摆放
    └── placement.json           # 摆放变换 + 体积回读
```

`design.json` 必填字段见 [reference.md](reference.md)。

**打印方向必须与 z 标定样件一致**（同一台机、同一喷嘴/层高/批次；仅方向差异就造成约 1.5 倍模量差）。
不一致则该肌肉的实测数据不能与标定曲线对照。

---

## 常见失败与处置

| 症状 | 原因 | 处置 |
|---|---|---|
| 连通分量 > 1 | 肌腹被骨/让位切断，或分叉处没连上 | 检查让位是否切透；加桥接杆，不要靠"看起来连着" |
| 钢丝间隙 < 0.3 mm | 通道路径贴边，或通道被晶格杆侵占 | 加大通道外径或局部加密走线；不要缩小钢丝直径来凑 |
| 体积回读误差大 | OBJ/STL 朝向或单位不一致（mm vs m） | 统一 mm、检查面朝向 |
| 通道打不出来 | 壁厚 0.75 mm 低于喷嘴能力，或孔内生成支撑 | 先打短片段验证孔径，再打整件 |
| 装到骨上干涉 | 用了原始解剖外轮廓，没做让位 | 见"外形照抄、机构设计"的让位表 |

---

## 参考

- 参数依据、几何细节、`design.json` 字段、失败目录：[reference.md](reference.md)
- 现有实现：`mujoco_arm/anatomical_concept/build_shaped_elbow_muscles.py`
  （断面/B样条/晶格/表面网）、`build_continuous_muscle_channels.py`（连续通道）
- 方法论文：Yoshimura, Kawaharazuka, Okada, *Advanced Intelligent Systems* 2026
- 上游边界：`02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md`
