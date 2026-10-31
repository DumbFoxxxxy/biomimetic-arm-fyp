# 02 · yoshimura_2026_ais

- **中文标题**：通过参数化设计与 3D 打印柔性结构加速肌肉骨骼机器人开发
- **原题**：Accelerating Musculoskeletal Robotics Through Parametric Design and 3D-Printed Flexible Structures
- **作者**：Shunnosuke Yoshimura, Kento Kawaharazuka, Kei Okada
- **出处**：Advanced Intelligent Systems (Wiley-VCH)（2026）
- **类型**：期刊论文　**等级**：核心　**分组**：项目直接依据
- **标签**：`TPU` `参数化设计` `韧带` `缓冲垫` `肌腱` `下肢` `统一制造流程`
- **源文件**：`D:\XJTLU\PSP\Advanced Intelligent Systems - 2026 - Yoshimura - Accelerating Musculoskeletal Robotics Through Parametric Design and.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Advanced_Intelligent_Systems_2026_Yoshimura_Accelerating_Musculoskeletal_Robotics_Through___426b72.txt`
- **规模**：17 页 / 提取 80702 字符

## 摘要

同一作者把 PSM 思路扩展成统一框架：用同一种 TPU 材料与同一套参数化流程同时生成肌肉、肌腱、韧带和缓冲软垫四类元件，通过调整图案密度和分支厚度在不同身体部位获得不同力学性能。韧带被用作关节约束，配合集成肌腱驱动，实现真人尺度下肢连续协调运动，躯干滚转、俯仰、偏航各可约束在正负 20 度范围。

## 要点

- 肌肉、肌腱、韧带、缓冲垫用一套参数化流程和单一材料完成，避免多材料多工艺。
- 关节活动范围由韧带长度与绷紧角决定，而不是靠硬限位，这是 MuJoCo 韧带条带的建模依据。
- 图案密度与分支厚度是调节局部力学性能的两个主控参数。
- 验证对象是下肢，迁移到右臂需要重新确定每条韧带的布置角。

## 与本项目的关系

你明确说过要参照这篇的方法搭右臂。ligament manifest 里的 81 条韧带条带、缓冲垫思路、以及 lattice_mid / lattice_end 两组晶格参数都来自这里。
