# 05 · kawaharazuka_2024_wire_optimization

- **中文标题**：肌腱驱动机器人钢丝布线（含可变绕行点）的数值仿真优化设计
- **原题**：Design Optimization of Wire Arrangement With Variable Relay Points in Numerical Simulation for Tendon-Driven Robots
- **作者**：Kento Kawaharazuka, Shunnosuke Yoshimura, Temma Suzuki, Kei Okada, Masayuki Inaba
- **出处**：IEEE Robotics and Automation Letters (RA-L), Vol. 9, No. 2, pp. 1388-（2024）
- **类型**：期刊论文　**等级**：核心　**分组**：项目直接依据
- **标签**：`绕行点` `力臂` `黑盒优化` `多目标` `可达力空间` `速度空间`
- **源文件**：`D:\XJTLU\PSP\Design_Optimization_of_Wire_Arrangement_With_Variable_Relay_Points_in_Numerical_Simulation_for_Tendon-Driven_Robots.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Design_Optimization_of_Wire_Arrangement_With_Variable_Relay_Points_in_Numerical_Simulation__6a3e91.txt`
- **规模**：8 页 / 提取 34340 字符

## 摘要

指出以前的钢丝布线优化都假设力臂恒定，或者忽略跨多关节、带绕行点的走线。本文把每根钢丝表示为起点加多个绕行点加终点，用黑盒优化求期望的物理性能，并同时考虑可达操作力空间与速度空间的多目标优化。

## 要点

- 绕行点（relay point）是可优化变量，它决定各关节上的力臂，直接影响能发力的大小与方向。
- 力臂随关节角变化，恒力臂假设会让结果失真。
- 目标是力和速度两个空间的折中，属于典型多目标问题。
- 方法可直接用于你后续在骨骼上布置夹具与过线导管的优化。

## 与本项目的关系

你模型里每根钢丝的 joint pulley 是按局部扫描最大化最差力臂自动求出的，思路源于这里。后续在骨骼上加固定钢丝走向的夹具，也应该用这套优化来定位置。
