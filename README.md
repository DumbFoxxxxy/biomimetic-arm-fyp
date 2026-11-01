# 仿生右臂 · 人体骨骼 + 工程关节 + 异形 TPU 晶格肌肉

用真实人体骨骼几何与解剖肌肉外形，做一条**能用钢丝牵引、能在台架上被实测**的仿生右臂。

**当前阶段**：肘部设计候选已交付（未实测）；肩部完成设计分解，准备进入 **L1（GH 单关节）**。
**最新周报**：[2026-11-01](reports/2026-11-01.md) · [全部周报](reports/README.md)

![肘部零位：真实骨骼 + 两块连续通道 TPU 肌肉](02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers/continuous_channel_candidate/elbow_neutral.png)

## 一周一屏

| | |
|---|---|
| **本周完成** | 肘部封板；肩部拆成 L0–L7；肌肉→TPU 流程固化为可复用 skill 与三个工具模块 |
| **本周最关键的数字** | 连续通道增量 **+652.655 mm³**，与设计记录 652.657 差 **0.002 mm³** |
| **最大的未验证项** | TPU 实际驱动范围 —— 仍是未知，不能把 0–150° 当成它 |
| **需要导师定的** | ①中段通道要不要按力线重构 ②肱骨近端 6 个止点怎么排 |

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

## 技术路线

两条独立路线，**不互相继承尺寸、锚点或通过结果**。

| 路线 | 内容 | 状态 |
|---|---|---|
| [方案二 · 当前主线](02_人体骨骼晶格肌肉方案/README.md) | 人体骨形 + 工程关节 + 异形 TPU 晶格肌肉 | 肘部候选 → 肩部 L1 |
| [方案一 · 历史](01_工程骨架方案/README.md) | 工程骨架、V2–V8 探索、旧钢丝/PSM 路线 | 归档，结论不转移 |

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
01_工程骨架方案/            历史路线与存档
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
