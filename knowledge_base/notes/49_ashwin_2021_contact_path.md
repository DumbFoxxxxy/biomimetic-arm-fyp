# 49 · ashwin_2021_contact_path

- **中文标题**：障碍环境下索驱连续机器人的轮廓与接触力估计
- **原题**：Profile and contact force estimation of cable-driven continuum robots in presence of obstacles
- **作者**：K. P. Ashwin, Soumya Kanti Mahapatra, Ashitava Ghosal
- **出处**：Mechanism and Machine Theory, Vol. 164, 104404（2021）
- **类型**：期刊论文　**等级**：重要　**分组**：钢丝张力分配与肩胸接触优化
- **标签**：`单侧接触` `接触优化` `钢丝路径` `障碍约束` `接触力估计`
- **DOI**：10.1016/j.mechmachtheory.2021.104404
- **源文件**：`D:\XJTLU\PSP\literature_inbox\49_profile_and_contact_force_estimation_of_cable_driven_continuum_robots_in_presence_of_obstacles.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\49_profile_and_contact_force_estimation_of_cable_driven_continuum_robots_in_presence_of_ob__147589.txt`
- **规模**：12 页 / 提取 39787 字符

## 摘要

论文从应变能最小原理建立索驱连续机器人形态优化，把障碍物写成不等式约束，不需要事先知道接触点。拉格朗日乘子非零的位置自动对应接触节点，再用虚功原理估计接触反力。实验验证的形态 RMS 误差约 1.35 mm，接触力误差小于 1 g。

## 要点

- 接触避障可以写成轨迹点始终位于障碍边界之外的不等式约束。
- 约束优化中的拉格朗日乘子可用来识别哪些节点真正发生接触。
- 接触力通过虚功原理从已求得的运动形态反算。
- 180 mm 连续机器人原型的形态 RMS 误差为 1.35 mm。
- 50-350 g 载荷范围内，接触力平均预测误差小于 1 g。
- 方法允许接触点沿障碍表面改变，不要求预先固定接触位置。

## 与本项目的关系

当前 V3/V4 主要用采样点最小距离判断净空。本文提供了更严格的表述：把胸廓和肩胛接触面写成解析不等式，只有乘子非零的点算真实接触；未接触时允许自由脱离。该结构可同时处理 G07 钢丝避障和肩胛接触垫。
