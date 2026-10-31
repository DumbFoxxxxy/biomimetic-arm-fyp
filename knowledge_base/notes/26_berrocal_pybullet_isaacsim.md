# 26 · berrocal_pybullet_isaacsim

- **中文标题**：面向机器人与强化学习的 PyBullet 与 Isaac Sim 评测
- **原题**：Evaluating PyBullet and Isaac Sim in the Scope of Robotics and Reinforcement Learning
- **作者**：Enaitz Berrocal, Basilio Sierra, Hector Herrero
- **出处**：IEEE 会议论文（TECNALIA / UPV-EHU）（2024）
- **类型**：会议论文　**等级**：参考　**分组**：仿真、优化与联合建模
- **标签**：`仿真器选型` `保真度` `计算效率` `PPO` `Stable-Baselines3` `Benchmark`
- **源文件**：`D:\XJTLU\PSP\Evaluating_PyBullet_and_Isaac_Sim_in_the_Scope_of_Robotics_and_Reinforcement_Learning.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Evaluating_PyBullet_and_Isaac_Sim_in_the_Scope_of_Robotics_and_Reinforcement_Learning__38289a.txt`
- **规模**：6 页 / 提取 34003 字符

## 摘要

物理仿真器的选择显著影响强化学习训练的效率与效果。本文在保真度与计算效率两个维度上对比 PyBullet 与 Isaac Sim，用机械臂到达任务作为基准，评估两者对 RL 算法性能的影响，为研究者选择仿真器提供依据。

## 要点

- 仿真器选型要在保真度与计算效率之间权衡，不存在全面更优的一方。
- 评测使用统一的 RL 算法与任务，保证可比性。
- 结论是任务导向的。
- 仅讨论两仿真器，未涉及 MuJoCo。

## 与本项目的关系

你选 MuJoCo 主要因为柔性体 flex 与肌腱 tendon 支持。这篇可作为为何不用 PyBullet 或 Isaac Sim 的对照引用：二者关注 RL 训练吞吐而非柔性接触精度，且都不在本文评测范围内。
