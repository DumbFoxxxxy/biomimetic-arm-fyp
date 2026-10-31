# 13 · sodeyama_2008_kojiro

- **中文标题**：Kojiro 仿人机器人球面胸廓加肩胛骨锁骨的肩部设计与运动
- **原题**：The Designs and Motions of a Shoulder Structure with a Spherical Thorax, Scapulas and Collarbones for Humanoid Kojiro
- **作者**：Yoshinao Sodeyama, Tamaki Nishino, Yuta Namiki, Yuto Nakanishi, Ikuo Mizuuchi, Masayuki Inaba
- **出处**：2008 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)（2008）
- **类型**：会议论文　**等级**：重要　**分组**：肩带机构与机器人肩部实现
- **标签**：`球面胸廓` `肩胛骨` `锁骨` `Kojiro` `滑动面`
- **源文件**：`D:\XJTLU\PSP\The_designs_and_motions_of_a_shoulder_structure_with_a_spherical_thorax_scapulas_and_collarbones_for_humanoid_Kojiro.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\The_designs_and_motions_of_a_shoulder_structure_with_a_spherical_thorax_scapulas_and_colla__a4cf30.txt`
- **规模**：6 页 / 提取 20838 字符

## 摘要

为肌肉驱动机器人 Kojiro 设计新的肩部结构。关键设计是把胸廓做成球面形状，再配合肩胛骨结构，从而让肩胛骨在胸廓上获得期望的运动。作者认为肩胛骨-锁骨结构有四个好处：活动范围大、柔性接近人、与驱动系统隔离更安全、胸腔内可容纳大量部件。论文描述了设计过程与动作测试。

## 要点

- 把胸廓做圆，是让肩胛骨在胸壁上滑动的最直接工程手段。
- 肩胛骨在胸廓上的运动本质是滑动约束，不是铰链。
- 四个优点中，安全性来自驱动系统与人体隔离这一点值得注意。
- 与 2007 年那篇同一系列，是迭代后的版本。

## 与本项目的关系

直接对应你模型里的肩胛胸壁滑动面。你的 thorax 由真实肋骨网格拼成，凹凸不平，接触代理不容易做；这篇的球面胸廓是一个可考虑的简化方向。
