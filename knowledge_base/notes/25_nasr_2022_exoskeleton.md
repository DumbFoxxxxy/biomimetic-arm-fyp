# 25 · nasr_2022_exoskeleton

- **中文标题**：主动-被动肩部外骨骼的优化设计：人机交互的计算建模
- **原题**：Optimal design of active-passive shoulder exoskeletons: a computational modeling of human-robot interaction
- **作者**：Ali Nasr, Sydney Bell, John McPhee
- **出处**：Multibody System Dynamics, Vol. 57, pp. 73-106（2022）
- **类型**：期刊论文　**等级**：重要　**分组**：仿真、优化与联合建模
- **标签**：`多体动力学` `肌肉骨骼模型` `关节主动力矩` `代谢能量` `优化设计`
- **DOI**：10.1007/s11044-022-09855-8
- **源文件**：`D:\XJTLU\PSP\Optimal design of active-passive shoulder exoskeletons a computational modeling of human-robot interaction.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Optimal_design_of_active_passive_shoulder_exoskeletons_a_computational_modeling_of_human_r__7e4432.txt`
- **规模**：34 页 / 提取 98435 字符

## 摘要

外骨骼设计因为存在闭链、复杂交互和传动特征之间的耦合选择而困难。本文提出一套主动与被动肩部外骨骼的动态综合流程：先用六部分组成多体模型（上肢肌肉骨骼模型、最优控制器、外骨骼刚体、被动机构、动力驱动器和助力模型），再以六个实测任务（前伸、左右横伸、过顶伸、矢状面搬物、额状面搬物、过顶搬物）为期望运动，优化被动机构与电机参数，使人体关节主动力矩、功率、肌肉代谢能量和驱动耗电同时最小。

## 要点

- 多体模型里人体用肌肉骨骼模型表示，可与外骨骼联合求解，这是人机耦合仿真的标准做法。
- 优化目标不是单一性能，而是人体力矩、功率、代谢能量与能耗的组合。
- 期望运动来自实测的六个日常任务，而不是凭空指定的关节角曲线。
- 被动机构与主动电机需要联合优化，分开设计会次优。

## 与本项目的关系

对你后续用 Optuna 做多目标优化有直接参考：目标函数应该是组合型的（工作空间、力矩、质量、能量），约束来自真实运动轨迹。
