# 48 · kraus_2015_pulley_friction

- **中文标题**：索驱并联机器人绞盘内置钢丝力测量的滑轮摩擦补偿与验证
- **原题**：Pulley Friction Compensation for Winch-Integrated Cable Force Measurement and Verification on a Cable-Driven Parallel Robot
- **作者**：Werner Kraus, Michael Kessler, Andreas Pott
- **出处**：2015 IEEE International Conference on Robotics and Automation (ICRA)（2015）
- **类型**：会议论文　**等级**：核心　**分组**：钢丝张力分配与肩胸接触优化
- **标签**：`Capstan摩擦` `导轮摩擦` `Dahl模型` `钢丝张力分配` `导轮`
- **DOI**：10.1109/icra.2015.7139406
- **源文件**：`D:\XJTLU\PSP\literature_inbox\48_pulley_friction_compensation_for_winch_integrated_cable_force_measurement_and_verification_on_a_cabl.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\48_pulley_friction_compensation_for_winch_integrated_cable_force_measurement_and_verificat__5ee3f0.txt`
- **规模**：6 页 / 提取 27438 字符

## 摘要

论文建立串联滑轮组的钢丝张力传递模型，使用 Eytelwein/Capstan 关系表示包角引起的张力指数变化，并把库仑、黏性摩擦与 Dahl 预滑移模型组合起来。实验识别 8 根钢丝的摩擦参数，发现滑轮摩擦平均使测力误差达到 8.4%，补偿后误差约降低三分之二。

## 要点

- 单个滑轮的基础张力关系为 F_j = F_(j-1) exp(mu_p alpha)，包角是核心变量。
- 串联多个滑轮时，必须逐级传播张力，不能只对总包角乘一次效率。
- 库仑项负责方向相关的静/动摩擦，黏性项与速度相关，Dahl 项描述预滑移。
- 实验识别的库仑系数约 0.011-0.025，黏性系数约 1.18-3.50 Ns/m。
- 平均钢丝力 6.40 N 时，滑轮摩擦约 0.54 N，未补偿误差约 8.4%。
- 在 50 kg 验证中，摩擦误差可由约 46.5 N 降到 20.4 N。

## 与本项目的关系

V4 当前用 exp(-mu * 总包角) 直接缩放整根钢丝，只能描述静态无回差情况。应改为按每个导轮的有符号包角和运动方向逐级计算，并分别记录导轮两侧张力；准静态无运动时可退化为 Eytelwein 上下界或对称近似。
