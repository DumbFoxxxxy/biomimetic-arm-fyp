# 42 · benhamida_2021_nsga2

- **中文标题**：面向康复任务的索驱并联机器人多目标优化设计
- **原题**：Multi-Objective optimal design of a cable driven parallel robot for rehabilitation tasks
- **作者**：Ines Ben Hamida, Med Amine Laribi, Abdelfattah Mlika, Lotfi Romdhane, Said Zeghloul, Giuseppe Carbone
- **出处**：Mechanism and Machine Theory, Vol. 156, 104141（2021）
- **类型**：期刊论文　**等级**：重要　**分组**：钢丝张力分配与肩胸接触优化
- **标签**：`钢丝张力分配` `正张成` `NSGA-II` `多目标优化` `绳索并联机器人` `Pareto` `可变绕行点`
- **DOI**：10.1016/j.mechmachtheory.2020.104141
- **源文件**：`D:\XJTLU\PSP\literature_inbox\42_multi_objective_optimal_design_of_a_cable_driven_parallel_robot_for_rehabilitation_tasks.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\42_multi_objective_optimal_design_of_a_cable_driven_parallel_robot_for_rehabilitation_task__74605f.txt`
- **规模**：24 页 / 提取 60316 字符

## 摘要

论文把索驱并联机器人的拓扑选择与尺寸综合放进同一个 NSGA-II 多目标优化。目标是最小化设备占用体积和轨迹上的最大钢丝张力，约束是目标轨迹必须位于工作空间内且所有钢丝张力保持正。张力用 Moore-Penrose 特解加零空间齐次解构造，最终从 Pareto 前沿按安全性等指标选择单一方案。

## 要点

- 决策向量同时包含离散拓扑类型和连续尺寸参数。
- 优化目标采用占用体积和最大钢丝张力，而不是只优化一个标量。
- 工作空间约束和正张力约束必须作为可行性条件处理。
- 冗余张力解写成 tau = W+ Wp + (I - W+ W)s，再用零空间项保证正张力。
- Pareto 前沿保留体积、张力和安全性之间的不同折中，不强行给出唯一最优。

## 与本项目的关系

V4 应把 G07 的路线拓扑、锚点和导轮数量一起优化，并同时报告收缩率、最大张力和净空，而不是用一个加权分数掩盖约束失败。该文的离散加连续编码也适合“0/1/2 个导轮”这种变量。
