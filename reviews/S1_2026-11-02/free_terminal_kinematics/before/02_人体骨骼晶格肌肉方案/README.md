# 方案二：人体骨骼＋工程关节＋异形TPU晶格肌肉

当前状态（2026-11-02）：当前推进肩部S1。当前研究起点为S1-Y2：紧凑GH骨架＋两块共形Y肌腹，共191对象（历史方案保留），尚未定版；肘部候选与试验准备独立保留。预测及实测TPU驱动范围未知，当前入口见[PROJECT_STATUS](../PROJECT_STATUS.md)。

**当前研究进展：191对象版裸骨正屈曲／外展达到120°检查上限，11步名义放置和有限保持有证据；实际紧固、强度和钢丝／TPU驱动仍未知。** 两块肩袖肌已采用共形Y格与菱形表面网格，其余189源不变；实际变形与打印待验。详见技术文件2.14–2.16节。


| 实际目录 | 内容 |
|---|---|
| 输入/骨骼 | 人体骨骼与支撑参考 |
| 输入/muscle | 已筛选的肌肉参考模型 |
| 输出/anatomical_concept | 当前肘部候选、接口、试验样件及诊断 |
| 输出/shoulder_design | 肩部技术审查稿 v0.9.13、S1-Y2过渡区修订肌腹、紧凑GH骨架及历史比较／验证记录 |
| 输出/anatomical_arm | 早期A1历史探索 |

- [肘部交付说明](输出/anatomical_concept/ELBOW_HANDOFF.md)
- [设计思路](输出/anatomical_concept/肘部设计思路.md)
- [当前候选](输出/anatomical_concept/CURRENT_CANDIDATE.md)
- [试验样件](输出/anatomical_concept/elbow_test_coupons/README.md)
- [肩部 S1 技术审查稿 v0.9.13](输出/shoulder_design/SHOULDER_DESIGN.md)

源码实际位于根目录mujoco_arm/anatomical_concept和mujoco_arm/anatomical_arm，原脚本目录联接已移除。
请按 [RUNNING.md](../RUNNING.md) 查看和复现，不需要恢复目录。
历史模型清单地址在读取时转换到当前实际目录；MuJoCo场景通过内存资源加载，支持中文项目路径。
完整外部BP3D源库仍需另行提供，当前已筛选输入保存在本方案中。
