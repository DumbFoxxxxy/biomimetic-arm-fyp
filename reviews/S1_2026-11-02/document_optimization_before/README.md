# 仿生右臂：人体骨骼＋工程关节＋异形 TPU 晶格肌肉

以真实解剖骨骼为基础，通过工程关节保证运动，以钢丝和可变形 TPU 晶格肌肉传力。允许短钢丝肌腱连接骨侧接口。

**当前状态（2026-11-02）：目标已恢复，正在按技术方案推进肩部 S1 关节与肌肉接口比较。** 肘部已有设计候选，肩部保留紧凑端盖候选，并新增独立 GH 保持关节对照；实际 TPU 驱动范围尚未实测。已允许重设计 S1，但尚未选定新结构。

| 入口 | 内容 |
|---|---|
| [PROJECT_STATUS](PROJECT_STATUS.md) | 唯一当前状态入口 |
| [肩部 S1 技术方案](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md) | 传力链、证据、限制及待讨论方向 |
| [当前有效规格](design/CURRENT_BASELINE.md) | 用户要求、比较指标与验收边界 |
| [设计任务书](design/DESIGN_BRIEF.md) | 当前入口与历史设计记录 |
| [运行说明](RUNNING.md) | 环境、依赖及查看方式 |

## 方案与目录

| 目录 | 用途 |
|---|---|
| [02_人体骨骼晶格肌肉方案](02_人体骨骼晶格肌肉方案/README.md) | 当前主线；输入骨骼／肌肉、肘部候选、肩部审查 |
| [01_工程骨架方案](01_工程骨架方案/README.md) | 历史路线对照；通过结果不转移到新骨骼 |
| `mujoco_arm/` | 场景、建模、查看及检查脚本 |
| [knowledge_base](knowledge_base/README.md) | 文献与方法参考 |
| `reviews/` | 版本审查、证据和文档快照 |
| `design/` | 有效规格、任务书和打印参考 |

肩部紧凑候选的条件刚性探针到达屈曲＋120°、外展＋61.1°、长轴＋79.3°／−47.0°；柔性肌肉及真实关节保持尚未验证。零位六路直线代理方向覆盖通过，也不能解释为可驱动。新增 GH 独立对照的条件刚性范围与该组结果不同，尚未采用为新 S1。完整输入、采样及未包含项见肩部技术方案。

本页旧进度与两路线比较保存在 [整理前快照](reviews/S1_2026-11-02/document_cleanup_before/README.md)。以下命令用于已有候选复现，不代表运行后即完成实体验收。

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
