# 20 · kawaharazuka_musashi

- **中文标题**：面向学习控制系统研究的模块化肌肉骨骼仿人平台 Musashi
- **原题**：Component Modularized Design of Musculoskeletal Humanoid Platform Musashi to Investigate Learning Control Systems
- **作者**：Kento Kawaharazuka, Shogo Makino, Kei Tsuzuki, Moritaka Onitsuka, Yuya Nagamatsu, Koki Shinjo, Tasuku Makabe, Yuki Asano, Kei Okada, Koji Kawasaki, Masayuki Inaba
- **出处**：IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)（2019）
- **类型**：会议论文　**等级**：重要　**分组**：肌肉骨骼人形平台（系统级）
- **标签**：`模块化` `关节模块` `肌肉模块` `非线性弹性单元` `通用骨框架` `可重构`
- **源文件**：`D:\XJTLU\PSP\Component_Modularized_Design_of_Musculoskeletal_Humanoid_Platform_Musashi_to_Investigate_Learning_Control_Systems.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Component_Modularized_Design_of_Musculoskeletal_Humanoid_Platform_Musashi_to_Investigate_L__ff758f.txt`
- **规模**：8 页 / 提取 38050 字符

## 摘要

为了做学习控制研究，目标是造出柔性肌肉骨骼结构、传感器冗余、且容易重构的身体。作者开发了可直接测量关节角的关节模块、能实现多种肌肉走线的肌肉模块、以及软结构的非线性弹性单元，随后组装出只由关节模块、肌肉模块、通用骨框架、肌肉钢丝单元和少量附件构成的 MusashiLarm，再扩展到全身形成 Musashi，并做了基础实验与学习控制实验验证其概念。

## 要点

- 关节模块可直测关节角，肌肉模块可实现任意走线路径，两者都是通用件。
- 通用骨框架加少量附件的思路，使改肌肉布置不需要重做骨架。
- 传感器冗余是刻意设计的目标，不是副作用。
- 整机先做单臂（MusashiLarm）再扩展全身，是可直接照搬的推进节奏。

## 与本项目的关系

你现在在骨骼上加夹具来固定钢丝走向，正是这里肌肉模块加通用骨框架思路的轻量版本。若将来要反复改走线，值得把夹具做成标准化模块。
