# 项目结构

更新时间：2026-10-25

当前工作交接与下一步优先级见 `HANDOFF.md`。新 agent 接手时先读该文件，
再读本结构与 `design/DESIGN_BRIEF.md` 的“当前版本速查”。

## 当前主链

| 目录 | 用途 |
| --- | --- |
| `final_model/` | 当前 SolidWorks 零件和装配 |
| `design/` | 设计任务书、锚点表和 Motion 顺序 |
| `mujoco_arm/` | 当前 16 钢丝 MuJoCo 演示 |
| `knowledge_base/` | 文献索引、笔记和全文 |
| `骨骼/`、`model/`、`muscle/` | 原始骨骼、肌肉和 CAD 参考资料 |
| `tmp/` | 当前设计计算和锚点生成源文件 |
| `archive/` | 已退出主流程的历史规划和临时数据 |

## SolidWorks

当前保留在 `final_model/` 的文件：

```text
Thorax.SLDPRT
clavicle.SLDPRT
scapula.SLDPRT
humerus.SLDPRT
ulna.SLDPRT
universal_joint_pin.SLDPRT
assembly_zero.SLDASM          # 干净零位装配
assembly_motion_gh.SLDASM     # 最新 Motion 实验装配
```

已移到 `final_model/archive/`：

- 旧的 14 路径装配和路径参考零件；
- 旧的 `装配体1.SLDASM`；
- 早期 `assembly_motion*.SLDASM`；
- 旧的 `cable_motion_work` 工作副本；
- 旧审计 STL 和报告；
- 一次性审计宏。

## MuJoCo 当前演示

当前主文件：

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
mujoco_arm/stl_io.py
```

启动：

```text
mujoco_arm/open_tendon_demo_gui.bat
```

SolidWorks 骨骼改变后重建：

```text
mujoco_arm/refresh_tendon_demo.bat
```

旧 BP3D、旧 cable-drive、旧模型和旧结果已移到：

```text
archive/mujoco_legacy/
```

当前肩胛分析工具在：

```text
mujoco_arm/tools/
```

K6/S01 独立工程基线在：

```text
mujoco_arm/engineering/
```

它使用解析 R400 接触带、肩胛接触垫和真实圆柱钢丝包角，不与
`arm_tendon_demo.xml` 共用碰撞或执行器配置。

面向实物硬点导出的 V2 钢丝设计链在：

```text
mujoco_arm/cable_design_v2/
```

它包含 16 根钢丝的动作分组、准静态张力评估、导轮生成和 SolidWorks 坐标导出。

V3 路线修复版本在：

```text
mujoco_arm/cable_design_v3/
```

V3 修复 G02、S03、S05/S06 的胸廓净空，并将原版本完整保留用于对比。

## 设计文档

当前权威文档：

```text
design/DESIGN_BRIEF.md
design/CABLE_ANCHORS_FINAL.md
design/MOTION_IMPLEMENTATION_ORDER.md
design/CABLE_ROUTING.md
design/IDEAS_backlog.md
```

旧 14 钢丝分析、旧 Motion 步骤和旧进度记录已移到：

```text
design/archive/
```

## 当前锚点计算源

保留：

```text
tmp/cable_anchors_final.py
tmp/cable_anchors_final.tsv
tmp/cable_anchors_final.csv
```

其他当前设计计算工具仍在 `tmp/`。旧的 14 钢丝长度分析、一次性夹具和测试文件已移到：

```text
tmp/archive/
archive/tmp_legacy/
```

## 后续规则

1. 新的正式设计结论先写进 `design/`，再改模型。
2. MuJoCo 新实验不要继续堆在 `mujoco_arm/` 根目录；实验脚本放 `mujoco_arm/tools/`，结果放 `mujoco_arm/exports/`。
3. SolidWorks 新版本使用新的文件名，不覆盖已验收版本；旧版本移到 `final_model/archive/`。
4. 临时文件放 `tmp/`，确认无用后再归档，不直接覆盖原始参考文件。
