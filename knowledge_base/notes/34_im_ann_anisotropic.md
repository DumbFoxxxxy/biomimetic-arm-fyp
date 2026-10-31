# 34 · im_ann_anisotropic

- **中文标题**：基于计算晶体结构数据的各向异性超弹性材料人工神经网络建模
- **原题**：Artificial neural network modeling for anisotropic hyperelastic materials based on computational crystal structure data
- **作者**：Sunyoung Im, Wonbae Kim, Hyungjun Kim, Maenghyo Cho
- **出处**：Seoul National University, Multiscale Mechanical Design Division（PDF 元数据标题为 AIAA 模板）
- **类型**：期刊或会议论文　**等级**：参考　**分组**：柔性与晶格材料的本构建模（方法储备）
- **标签**：`各向异性` `超弹性` `神经网络` `应变能密度` `本构张量`
- **源文件**：`D:\XJTLU\PSP\Artificial-neural-network-modeling-of-anisotropic-hyperelastic-materials-based-on-computational-crystal-structure-data.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Artificial_neural_network_modeling_of_anisotropic_hyperelastic_materials_based_on_computat__aa69a9.txt`
- **规模**：12 页 / 提取 27095 字符

## 摘要

传统超弹性建模假设应变能密度函数形式，再对应变求导得到应力与本构张量，但对复杂六方、四方晶体结构这类各向异性超弹性材料，很难找到合适的应变能密度函数。本文绕过该假设，直接用计算晶体结构数据里的应变、应力与本构张量数据训练人工神经网络，从而预测各向异性超弹性行为。

## 要点

- 用神经网络直接学习应力-应变映射，跳过猜应变能函数这一步。
- 训练数据来自计算晶体结构模型而不是实验。
- 针对的是强各向异性材料。
- 优点是不需要假设函数形式，代价是可解释性下降。

## 与本项目的关系

X 晶格是强各向异性结构，用解析本构很难写准。这篇提供直接学映射的替代路线，可与第 33 篇组合使用。
