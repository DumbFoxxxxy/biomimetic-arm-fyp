# 27 · du_2025_coadaptation

- **中文标题**：基于代理模型引导优化的仿人机器人设计与运动控制协同适应
- **原题**：Efficient co-adaptation of humanoid robot design and locomotion control using surrogate-guided optimization
- **作者**：Yidong Du, Xuechao Chen, Zhangguo Yu, Fei Meng, Zishun Zhou, Yuanxi Zhang, Qingqing Li, Qiang Huang
- **出处**：Biomimetic Intelligence and Robotics, Vol. 5, 100255（2025）
- **类型**：期刊论文　**等级**：参考　**分组**：仿真、优化与联合建模
- **标签**：`协同设计` `强化学习` `代理模型` `高维设计空间` `优化效率`
- **源文件**：`D:\XJTLU\PSP\Efficient co-adaptation of humanoid robot design and locomotion.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Efficient_co_adaptation_of_humanoid_robot_design_and_locomotion__63f268.txt`
- **规模**：11 页 / 提取 60576 字符

## 摘要

强化学习与算力的进步使数据驱动的运动控制与机械设计优化变得可行，可替代人工设计。但现有协同设计方法面临一个关键困难：高维设计空间加上昂贵的物理仿真导致搜索在算力上不可行。本文用代理模型引导优化来提升协同适应的效率。

## 要点

- 设计变量与控制策略联合优化，能避免先定结构再调控制带来的次优。
- 高维设计空间加昂贵仿真是协同设计的核心瓶颈。
- 代理模型（surrogate）用来降低对高保真仿真的调用次数。
- 与你用 Optuna NSGA-II 做大范围晶格参数搜索面临的问题完全同类。

## 与本项目的关系

若后面要把肌肉数量、走线、晶格参数一起优化，搜索空间会到几十维。这篇的代理模型引导思路比纯随机采样更省算力。
