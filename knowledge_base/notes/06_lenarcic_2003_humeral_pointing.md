# 06 · lenarcic_2003_humeral_pointing

- **中文标题**：仿人肩复合体与肱骨指向运动学
- **原题**：A Humanoid Shoulder Complex and the Humeral Pointing Kinematics
- **作者**：Jadran Lenarcic, Michael Stanisic
- **出处**：IEEE Transactions on Robotics and Automation, Vol. 19, No. 3, pp. 499-（2003）
- **类型**：期刊论文　**等级**：重要　**分组**：肩复合体解剖与运动学建模
- **标签**：`肩带` `并联机构` `盂肱关节` `偏置双指向系统` `冗余自由度` `奇异锥`
- **源文件**：`D:\XJTLU\PSP\A_humanoid_shoulder_complex_and_the_humeral_pointing_kinematics.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\A_humanoid_shoulder_complex_and_the_humeral_pointing_kinematics__d7eae8.txt`
- **规模**：8 页 / 提取 52044 字符

## 摘要

把仿人肩复合体拆成两个子系统：作为内层肩带（肩胸）的并联机构和作为外层球形盂肱关节的串联机构，两者之间用偏置距离和扭转角分开，整体构成一个偏置双指向系统。所谓肱骨指向是肩带位移与肱骨共面、且两个子系统的倾角之比取常数的位形，这和人的肱骨指向一致。肩带保留一个冗余自由度，可用来避开关节的奇异锥。

## 要点

- 肩带与盂肱关节之间存在固定偏置与扭转，不能当成一个单纯的球铰。
- 肱骨指向对应一个固定的两子系统倾角比，这等价于肩肱节律的几何表达。
- 肩带的冗余自由度可以用来优化构型、避开奇异。
- 这是把锁骨与肩胛骨纳入运动学的经典处理方式。

## 与本项目的关系

你在 MuJoCo 里用 SC / AC / ST / GH 四个关节串联表达肩带，并用约 2:1 肩肱节律耦合，理论依据就是这类模型。
