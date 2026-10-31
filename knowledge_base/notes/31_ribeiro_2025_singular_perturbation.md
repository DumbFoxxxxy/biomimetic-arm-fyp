# 31 · ribeiro_2025_singular_perturbation

- **中文标题**：肌腱驱动软体机器人的奇异摄动控制：理论与实验
- **原题**：Singular-Perturbation Control of a Tendon-Driven Soft Robot: Theory and Experiments
- **作者**：Lucas Novaki Ribeiro, Pablo Borja, Cosimo Della Santina, Bastian Deutschmann
- **出处**：IEEE Transactions on Control Systems Technology, Vol. 33, No. 5, pp. 1929-（2025）
- **类型**：期刊论文　**等级**：参考　**分组**：人工肌肉与软体驱动
- **标签**：`肌腱驱动` `连续体` `奇异摄动` `轨迹跟踪` `驱动器动力学` `低刚度失稳`
- **源文件**：`D:\XJTLU\PSP\Singular-Perturbation_Control_of_a_Tendon-Driven_Soft_Robot_Theory_and_Experiments.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Singular_Perturbation_Control_of_a_Tendon_Driven_Soft_Robot_Theory_and_Experiments__c6096c.txt`
- **规模**：8 页 / 提取 39863 字符

## 摘要

现有肌腱驱动连续体软体机器人的模型控制策略忽略了驱动系统动力学，但驱动动力学对闭环性能有重要影响。本文用奇异摄动理论设计能把驱动动力学纳入的控制器，给出解析形式并做深入实验验证，还用高刚度和低刚度两种肌腱对比所提控制器与传统反馈控制的效果。结果显示奇异摄动控制器能减少振荡、降低误差且不需要高增益；更重要的是，当肌腱刚度较低时，忽略驱动动力学可能导致不稳定。

## 要点

- 肌腱刚度低时，把驱动器当成理想位移源会直接导致闭环失稳。
- 奇异摄动法把快变的驱动动力学与慢变的机械臂动力学分层处理。
- 结论经过高刚度与低刚度肌腱的对照实验验证。
- 不需要高增益，这一点对实际电机很重要。

## 与本项目的关系

你的钢丝属于低刚度长柔性传动，这篇的警告恰好针对你：直接把滑块位置当作钢丝收缩量会失真甚至失稳。做位置-张力联合控制时要从这里取方法。
