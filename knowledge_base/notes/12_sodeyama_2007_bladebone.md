# 12 · sodeyama_2007_bladebone

- **中文标题**：采用肩胛骨-锁骨结构的宽活动范围肩部设计与运动
- **原题**：The Designs and Motions of a Shoulder Structure with a Wide Range of Movement Using Bladebone-Collarbone Structures
- **作者**：Yoshinao Sodeyama, Tomoaki Yoshikai, Tamaki Nishino, Ikuo Mizuuchi, Masayuki Inaba
- **出处**：Proceedings of the 2007 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), pp. 3629-（2007）
- **类型**：会议论文　**等级**：重要　**分组**：肩带机构与机器人肩部实现
- **标签**：`肩胛骨` `锁骨` `几何约束` `运动生成` `仿真`
- **源文件**：`D:\XJTLU\PSP\The_designs_and_motions_of_a_shoulder_structure_with_a_wide_range_of_movement_using_bladebone-collarbone_structures.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\The_designs_and_motions_of_a_shoulder_structure_with_a_wide_range_of_movement_using_bladeb__11d6b9.txt`
- **规模**：6 页 / 提取 23500 字符

## 摘要

给出肩胛骨-锁骨结构的三条优点：活动范围大、结构外置腾出胸腔空间、柔性接近人；但代价是自由度多且受奇异几何约束、难以控制。因此作者设计了肩胛骨-锁骨结构的搭建流程，并开发了一个考虑几何约束的仿真系统，最后基于该仿真在真机上做运动实验。

## 要点

- 肩胛骨-锁骨结构的核心难点是几何约束下的运动生成，而不是机械强度。
- 必须先有能表达几何约束的仿真系统，才能规划真机动作。
- 这是先仿真、再真机这一工作流的早期范例。
- 与你目前的路线完全同构：先在 MuJoCo 里把约束跑通再做实物。

## 与本项目的关系

给你现在的工作流提供正当性：几何约束下动作难规划，所以先建仿真再上真机是被验证过的顺序。
