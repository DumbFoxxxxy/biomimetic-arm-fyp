# 43 · asaoka_2012_muscle_arrangement

- **中文标题**：用动力学仿真确定特定运动的最优多关节肌肉布置
- **原题**：Determining the optimal multiarticular muscle arrangement of a musculoskeletal robot for a specific motion using dynamics simulation
- **作者**：T. Asaoka, M. Kawamura, S. Kumakura, I. Mizuuchi
- **出处**：Proceedings of the 12th IEEE-RAS International Conference on Humanoid Robots, pp. 216-221（2012）
- **类型**：会议论文　**等级**：重要　**分组**：钢丝张力分配与肩胸接触优化
- **标签**：`多关节肌肉` `肌肉布置` `动力学优化` `钢丝张力分配` `K7`
- **DOI**：10.1109/humanoids.2012.6651523
- **源文件**：`D:\XJTLU\PSP\literature_inbox\43_determining_the_optimal_multiarticular_muscle_arrangement_of_a_musculoskeletal_robot_for_a_specific_.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\43_determining_the_optimal_multiarticular_muscle_arrangement_of_a_musculoskeletal_robot_fo__d81a8b.txt`
- **规模**：6 页 / 提取 26243 字符

## 摘要

论文不照搬人体肌肉布置，而是针对一个具体动态动作枚举三关节机器人上的多关节肌肉组合。关节力矩由肌肉张力和肌肉长度 Jacobian 耦合得到，作者在动力学仿真中评价产生的冲量特性，从 736281 种六肌肉布置中搜索最优方案，并用实体人工肌肉机器人复现该动作。

## 要点

- 多关节肌肉通过长度 Jacobian 同时影响多个关节，不能按单关节肌肉逐一设计。
- 评价指标来自完整动作的动力学结果，而不是零位静态力臂。
- 研究先枚举全部可行肌肉交叉模式，再对候选做动力学仿真。
- 共搜索 736281 种六肌肉、三关节布置。
- 搜索得到的最优布置不必与人体解剖布置相同。

## 与本项目的关系

K7 是典型的多关节协同动作。可借鉴其方法，对 S01-S04、G02、G07 的关节交叉模式做离散枚举，再对少量拓扑执行完整五点动力学/准静态校核，而不是连续微调一个错误拓扑。
