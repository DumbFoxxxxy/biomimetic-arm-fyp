# 方案二：人体骨骼＋工程关节＋异形TPU晶格肌肉

当前状态（2026-11-02）：当前推进肩部S1。最新集成几何候选为187对象前衬垫定位保持版（原175对象及父级保留），尚未定版；肘部候选与试验准备独立保留。预测及实测TPU驱动范围未知，当前入口见[PROJECT_STATUS](../PROJECT_STATUS.md)。

**当前研究进展：195对象独立GH键槽／横向保持件版已完成有限防转、止退和13步裸骨放置检查；真实工具、螺纹预紧、强度及完整肌肉安装尚未验证，未替代187对象集成起点。** 详见技术文件2.13节。


| 实际目录 | 内容 |
|---|---|
| 输入/骨骼 | 人体骨骼与支撑参考 |
| 输入/muscle | 已筛选的肌肉参考模型 |
| 输出/anatomical_concept | 当前肘部候选、接口、试验样件及诊断 |
| 输出/shoulder_design | 肩部技术审查稿 v0.9.10、187对象集成候选及GH试改、历史比较及验证记录 |
| 输出/anatomical_arm | 早期A1历史探索 |

- [肘部交付说明](输出/anatomical_concept/ELBOW_HANDOFF.md)
- [设计思路](输出/anatomical_concept/肘部设计思路.md)
- [当前候选](输出/anatomical_concept/CURRENT_CANDIDATE.md)
- [试验样件](输出/anatomical_concept/elbow_test_coupons/README.md)
- [肩部 S1 技术审查稿 v0.9.10](输出/shoulder_design/SHOULDER_DESIGN.md)

源码实际位于根目录mujoco_arm/anatomical_concept和mujoco_arm/anatomical_arm，原脚本目录联接已移除。
请按 [RUNNING.md](../RUNNING.md) 查看和复现，不需要恢复目录。
历史模型清单地址在读取时转换到当前实际目录；MuJoCo场景通过内存资源加载，支持中文项目路径。
完整外部BP3D源库仍需另行提供，当前已筛选输入保存在本方案中。
