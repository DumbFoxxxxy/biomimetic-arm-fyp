# 方案一：工程骨架与钢丝/PSM探索

以工程骨架、关节和预设锚点为基础，包含V2–V8探索、V2-S1和旧钢丝/TPU接口。

| 实际目录 | 内容 |
|---|---|
| SolidWorks | 原生零件和装配，未自动与后续STEP同步 |
| 早期模型 | 旧模型资料 |
| 输出 | 各阶段CAD、场景、报告及失败诊断 |
| 历史归档 | 退出主流程的旧资料 |
| 历史规划与交接 | 整理前规划快照 |

- [V2-S1骨架](输出/skeleton_from_v2/V2-S1/README.md)
- [PSM-W1](输出/wire_psm_from_v2s1/PSM-W1/README.md)
- [旧设计文档索引](设计文档索引.md)

源码实际位于根目录mujoco_arm，涉及cable_design_v2–v8、skeleton_from_v2、wire_psm_from_v2s1等子目录。
原“脚本”和“旧演示代码入口”联接已移除。运行说明见 [RUNNING.md](../RUNNING.md)。
这里的尺寸、锚点和通过结果不转移到人体骨骼方案。SolidWorks GUI外部引用没有在本次逐件核对。
