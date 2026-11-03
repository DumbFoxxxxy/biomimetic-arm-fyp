# 仿生右臂：人体骨骼＋工程关节＋异形 TPU 晶格肌肉

以真实解剖骨骼为基础，通过工程关节保证运动，以钢丝和可变形 TPU 晶格肌肉传力。允许短钢丝肌腱连接骨侧接口。

**当前状态（2026-11-02）：正在推进肩部S1，最新集成几何候选为122对象台架接口版（父92对象保留），尚未定版。** 肘部候选独立保留；最新肩部结构的预测及实测驱动范围均未知。

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

最新肩部候选包含GH机械保持关节、60 mm骨干抱箍、两块肩袖肌自由末端／肌腱及近端支承。已有网格、孔路、零位钢丝和有限接触检查，台架固定、通道壁、压接、柔性运动及实物性能仍待验证。24／34／46对象的旧刚性范围不适用于72／92对象自由末端结构。最新几何入口、证据边界和原K1–K7指标见肩部技术方案。

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
