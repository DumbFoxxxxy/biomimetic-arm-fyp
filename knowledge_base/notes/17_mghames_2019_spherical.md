# 17 · mghames_2019_spherical

- **中文标题**：面向仿人机器人与人类的球面主动关节
- **原题**：A Spherical Active Joint for Humanoids and Humans
- **作者**：Sariah Mghames, Manuel Giuseppe Catalano, Antonio Bicchi, Giorgio Grioli
- **出处**：IEEE Robotics and Automation Letters (RA-L), Vol. 4, No. 2, pp. 838-（2019）
- **类型**：期刊论文　**等级**：重要　**分组**：肩带机构与机器人肩部实现
- **标签**：`球关节` `串并联混合` `工作空间` `奇异` `肩关节样机`
- **源文件**：`D:\XJTLU\PSP\A_Spherical_Active_Joint_for_Humanoids_and_Humans.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\A_Spherical_Active_Joint_for_Humanoids_and_Humans__ce0293.txt`
- **规模**：8 页 / 提取 40472 字符

## 摘要

指出球关节有三种运动学实现：串联、并联、混合，各有优劣。本文提出一种混合式主动球关节，结合并联与串联的优点，追求生物关节的四个特征：大工作空间、紧凑尺寸、动态性能和整体球形外形。作者把该关节的工作空间与人体关节对比，发现除了运动学奇异之外几乎可以完全覆盖；并制作了首个样机作为机器人肩关节进行初步测试。

## 要点

- 球关节的三种实现方式（串联、并联、混合）决定了工作空间、刚度与惯量的分布。
- 混合式能同时拿到并联的刚度与串联的大行程。
- 整体球形外形本身就是一个设计目标，因为它贴合生物关节。
- 工作空间的边界通常由运动学奇异而非机械限位决定。

## 与本项目的关系

如果以后要把 MuJoCo 里的简化球铰换成真实机构，这篇的混合式方案是肩关节的候选。
