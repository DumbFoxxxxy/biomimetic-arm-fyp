# 仿生右臂 · 人体骨骼 + 工程关节 + 异形 TPU 晶格肌肉

用真实人体骨骼与解剖肌肉几何做一条可用钢丝驱动的右臂，再和人体对比——**差在哪、最大能到多少运动范围**。

**当前阶段**：肘部设计候选已交付（未实测）；肩部完成设计分解，准备进入 **L1（GH 单关节）**。
**最新周报**：[2026-11-01](reports/2026-11-01.md) · [全部周报](reports/README.md)

![肘部零位：真实骨骼 + 两块连续通道 TPU 肌肉](02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers/continuous_channel_candidate/elbow_neutral.png)

## 进度

| | 状态 |
|---|---|
| **肘部** | ✅ 设计候选已交付：127 件零位装配、125 实体 STEP、刚性 0–150° 266 条件姿态。**未实测** |
| **肩部** | 🔶 L0（GH 中心与三轴）已定；L1 选型中，尚未建模；L2 起未开始 |
| **工具链** | ✅ 肌肉→TPU 的 skill、坐标 / 力线 / 干涉三个模块，均带自检与回归 |
| **实物** | ⬜ 未开始：打印、压接、台架 |

## 已完成 / 未验证

这一栏是本项目最重要的纪律：**几何通过 ≠ 能用**。

| 已完成（有几何或脚本证据） | 未验证（不要当成结论） |
|---|---|
| 127 件零位装配，125 实体 STEP 回读有效 | TPU 实际可驱动范围（**未知**） |
| 两块肌肉各为单一连通实体，钢丝名义相交 0 | 通道能否打印、弯曲后是否闭塞 |
| 刚性肘 0–150°，266 条件姿态 + 独立回程无碰撞 | 连续管壁带来的收缩阻力与卡线 |
| 销轴捕获：外移 0.09 mm 不碰、0.11 mm 挡住 | 防脱力、防松、振动寿命 |
| 四类打印对比样件 + 12 行空白测量记录 | 压接夹持力（压接坯仍是名义几何） |
| 坐标链自检 dE = 0、dF = 0 | 台架夹具柔度（唯一必须实物标定的项） |
| 力线回归：−71.309 mm / +95.804 mm 与旧记录一致 | 完整装配（TPU + 钢丝 + 电机 + 导管） |
| 干涉判据统一后可复现存档（0.0250 / 0.0873 mm³） | 肩部：**只有文档，没有任何几何** |

## 两方案思路

两条路线**不互相继承尺寸、锚点或通过结果**。方案一并未放弃，只是目前没有达到设想，现处于并行评估中。

### [方案一 · 工程骨架 + 钢丝驱动](01_工程骨架方案/README.md)（并行评估中）

人工设计骨架（R400 胸廓、150 mm 锁骨、5 mm 肩胛板），预设 7 自由度链
（SC 十字轴 → AC → GH 三铰链 → 肘），先冻结锚点表，再用 16 根钢丝在 MuJoCo 里复现 K1–K7。

做到的部分不差：K1 肘 145°、K3 外展 120°、K4 内旋 70°、K5 外旋 90° 都达到目标，16 钢丝演示通过。

发现两个问题：

1. **肩胛骨与各肌肉的锚点很难确定。** 真实肩胛与胸廓之间没有骨性关节，靠肌肉张力和接触面撑着；
   人工骨架上只能用一条解析曲面近似，结果 K6 上回旋只到 19°（目标 50°）、K7 联合抬升只到 146°（目标 170°），
   再增加轮廓切除就超出保留形态的预算。
2. **即使强行确定锚点，这条手臂与"仿生"的关系也不大了。** 骨架形状是人为设定的、
   肌肉附着位置是设计者指定的——它验证的是机构逻辑，不是解剖结构。

### [方案二 · 人体骨骼 + 工程关节 + TPU 晶格肌肉](02_人体骨骼晶格肌肉方案/README.md)（当前推进）

起因是 Yoshimura 等 2026《Accelerating Musculoskeletal Robotics Through Parametric Design
and 3D-Printed Flexible Structures》：把肌肉体积参数化成 TPU 晶格 + 连续钢丝通道，
用同一套工艺做肌肉、肌腱、韧带与缓冲件。

目标不是"证明能造出来"，而是：**用真实解剖几何强行做一遍，再和人体手臂比——
差别在哪、最大能到多少运动范围。**

做法：真实 BodyParts3D 骨骼与肌肉外形 → 工程关节实体化 → 异形 TPU 肌腹内含连续钢丝通道 →
先做到"手拉钢丝可重复往返"，再逐层加自由度（L0–L7）。

### 两条路线的差别

| | 方案一 | 方案二 |
|---|---|---|
| 骨骼几何 | 人工设计 | 真实 BodyParts3D |
| 肌肉锚点 | 预先设定并冻结 | 由解剖结构决定 |
| 驱动方式 | 16 根钢丝直接拉关节 | 肌肉本体传力（钢丝在肌腹内） |
| 要回答的问题 | 机构能否复现指定动作 | 与人体差多少、极限在哪 |

## 肩部扩充路线

![肩部 L0–L7 扩充路线](docs/roadmap.svg)

原则：**每层放开一个关节、加一组肌肉、用台架验证一次**，上一层没验完不进入下一层。
详见 [SHOULDER_DESIGN.md](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)。

## 仓库导航

```
02_人体骨骼晶格肌肉方案/    当前主线
  ├─ 输入/                  人体骨骼 + 已筛选肌肉（BodyParts3D, CC BY-SA 2.1 JP）
  ├─ 输出/anatomical_concept/ 肘部候选、接口、试验样件、诊断
  └─ 输出/shoulder_design/    肩部技术文档
01_工程骨架方案/            方案一（并行评估中）
00_共用资料/                文献知识库、设计历史总账
mujoco_arm/                 源码（场景、建模、审计脚本）
.cursor/skills/             可复用流程（肌肉 → TPU 的完整做法）
mujoco_arm/tools/           frames / force_line / interference（带自检与回归）
docs/                       首页用的图
reports/                    每周汇报
```

## 快速复现

```powershell
git lfs install && git lfs pull          # 模型与图片走 LFS

python mujoco_arm/tools/frames.py                              # 坐标自检
python mujoco_arm/tools/force_line.py --self-test              # 力线解析自检
python mujoco_arm/tools/force_line.py --muscle FJ1486 --to 150 --step 1 --bone
python mujoco_arm/tools/regression_interference.py             # 干涉模块回归

python mujoco_arm/anatomical_concept/view_continuous_elbow.py  # 肘部零位
python mujoco_arm/anatomical_concept/view_paired_swivel.py --keepers --motion  # 刚性动画
```

环境与依赖说明见 [RUNNING.md](RUNNING.md)。MuJoCo 资源由 Python 读入内存后交给原生接口，
因此项目可以放在含中文的路径下。

## 工作原则

1. **结论边界**：几何与数值检查只能写成"几何通过"，不等于能驱动、能打印、能装配。
2. **分步验证**：每层先设可复算的门，上一层没验完不进入下一层；失败记录保留不覆盖。
3. **单一事实源**：坐标只在 `frames.py`、力线只在 `force_line.py`、干涉判据只在 `interference.py`，
   避免同一件事在十几个脚本里各写一遍。

---

<sub>文献依据与方法来源见 [knowledge_base](knowledge_base/README.md)；人体模型来自 BodyParts3D（CC BY-SA 2.1 Japan）。</sub>
