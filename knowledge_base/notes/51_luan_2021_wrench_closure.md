# 51 · luan_2021_wrench_closure

- **中文标题**：索驱动并联机构的力闭合条件
- **原题**：Wrench-Closure Condition of Cable-Driven Parallel Manipulators
- **作者**：Phan Gia Luan, Nguyen Truong Thinh
- **出处**：Applied Sciences, Vol. 11, No. 9, 4228（2021）
- **类型**：期刊论文　**等级**：核心　**分组**：钢丝张力分配与肩胸接触优化
- **标签**：`钢丝张力分配` `正张成` `力闭合` `绳索并联机器人` `接触优化`
- **DOI**：10.3390/app11094228
- **源文件**：`D:\XJTLU\PSP\literature_inbox\51_wrench_closure_condition_of_cable_driven_parallel_manipulators.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\51_wrench_closure_condition_of_cable_driven_parallel_manipulators__e1e16e.txt`
- **规模**：13 页 / 提取 35048 字符

## 摘要

论文研究索驱动并联机构何时能够仅依靠单侧受拉的钢丝平衡任意外力。静力关系写成 Wt = psi_e，其中 W 是 wrench matrix，t 必须逐项非负。作者把 wrench-closure workspace 判定拆成两个条件：W 的列必须张满 R^m，并且列向量的非负组合能够包含原点。随后用锥交集、仿射独立性和 Dines 的 I-minor 方法给出可执行判定，并简化单冗余 n-m=1 的情形。

## 要点

- 钢丝只能受拉，因此所有张力必须满足 t_i >= 0，不能直接使用无约束最小二乘。
- 力闭合的必要条件是 wrench matrix 满秩，充分条件是其列向量的正张成覆盖 R^m。
- 等价判定要求存在严格正的列组合等于零：sum(t_i w_i)=0 且 t_i>0。
- 最少需要 m+1 根仿射独立的钢丝，才可能正张成一个 m 维力空间。
- 当钢丝数只比自由度数多 1 时，剩余向量必须落入前 m 个负列向量所张成的锥中。
- 论文提供基于 I-minor 的递归判定，可用于替代只检查矩阵秩或单个线性规划解的近似方法。

## 与本项目的关系

当前 V4 的 _positive_spanning 使用归一化列向量和线性规划判断零点是否位于凸包内部，方向正确，但没有把满秩、锥交集和钢丝数下限完整分开报告。本文可直接把正张成从软罚分提升为逐姿态硬约束，并为 K2、K3、K7、SLIDE_FWD 给出可解释的失败类型。
