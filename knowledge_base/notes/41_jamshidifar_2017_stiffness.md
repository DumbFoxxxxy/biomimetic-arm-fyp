# 41 · jamshidifar_2017_stiffness

- **中文标题**：运动学约束冗余索驱并联机器人：建模、冗余分析与刚度优化
- **原题**：Kinematically-Constrained Redundant Cable-Driven Parallel Robots: Modeling, Redundancy Analysis, and Stiffness Optimization
- **作者**：Hamed Jamshidifar, Amir Khajepour, Baris Fidan, Mitchell Rushton
- **出处**：IEEE/ASME Transactions on Mechatronics, Vol. 22, No. 2, pp. 921-930（2017）
- **类型**：期刊论文　**等级**：重要　**分组**：钢丝张力分配与肩胸接触优化
- **标签**：`钢丝张力分配` `冗余绳索机器人` `刚度优化` `零空间` `正张成`
- **DOI**：10.1109/tmech.2016.2639053
- **源文件**：`D:\XJTLU\PSP\literature_inbox\41_kinematically_constrained_redundant_cable_driven_parallel_robots_modeling_redundancy_analysis_and_st.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\41_kinematically_constrained_redundant_cable_driven_parallel_robots_modeling_redundancy_an__e69ada.txt`
- **规模**：10 页 / 提取 50890 字符

## 摘要

论文建立带运动学约束的冗余索驱并联机器人柔性模型，把平台刚度写成结构刚度与钢丝张力贡献之和。冗余张力写成特解加 Jacobian 零空间项 tau = J^T+ w + H alpha，在每根钢丝的最小、最大张力约束内优化零空间系数，并证明方向刚度是这些系数的线性函数，因此可用线性规划选择既满足正张力又提高指定方向刚度的解。

## 要点

- 静力方程的通解由最小范数特解和零空间向量组成，冗余主要通过零空间系数处理。
- 钢丝张力不仅平衡外力，还直接改变平台刚度和扰动响应。
- 张力上下界构成零空间系数的多面体可行域，可用线性规划求解。
- 约束钢丝的轴向刚度与所有钢丝预张力共同进入整体刚度矩阵。
- 仅追求最小张力可能牺牲稳定性和刚度，需在可行张力集内选择次级目标。

## 与本项目的关系

V4 当前用有界最小二乘直接求张力，容易只得到低张力可行解。本文提示应保留零空间自由度，用它提高关键方向刚度或稳定裕度，同时把每根钢丝的最小预张力作为硬约束。
