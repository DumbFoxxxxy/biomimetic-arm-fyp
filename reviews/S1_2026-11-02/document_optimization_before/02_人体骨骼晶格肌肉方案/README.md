# 方案二：人体骨骼＋工程关节＋异形TPU晶格肌肉

当前状态（2026-11-02）：目标已恢复，继续按技术方案推进 S1。肘部保留候选与试验准备；肩部已有紧凑端盖与 GH 独立保持关节几何对照，重设计尚未定案。均未制造定版，实际 TPU 驱动范围未知。当前入口见 [PROJECT_STATUS](../PROJECT_STATUS.md)。

| 实际目录 | 内容 |
|---|---|
| 输入/骨骼 | 人体骨骼与支撑参考 |
| 输入/muscle | 已筛选的肌肉参考模型 |
| 输出/anatomical_concept | 当前肘部候选、接口、试验样件及诊断 |
| 输出/shoulder_design | 肩部技术审查稿 v0.9.0、暂停前 S1 候选及历史记录 |
| 输出/anatomical_arm | 早期A1历史探索 |

- [肘部交付说明](输出/anatomical_concept/ELBOW_HANDOFF.md)
- [设计思路](输出/anatomical_concept/肘部设计思路.md)
- [当前候选](输出/anatomical_concept/CURRENT_CANDIDATE.md)
- [试验样件](输出/anatomical_concept/elbow_test_coupons/README.md)
- [肩部 S1 技术审查稿 v0.9.0](输出/shoulder_design/SHOULDER_DESIGN.md)

源码实际位于根目录mujoco_arm/anatomical_concept和mujoco_arm/anatomical_arm，原脚本目录联接已移除。
请按 [RUNNING.md](../RUNNING.md) 查看和复现，不需要恢复目录。
历史模型清单地址在读取时转换到当前实际目录；MuJoCo场景通过内存资源加载，支持中文项目路径。
完整外部BP3D源库仍需另行提供，当前已筛选输入保存在本方案中。
