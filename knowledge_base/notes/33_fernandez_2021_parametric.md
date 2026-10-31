# 33 · fernandez_2021_parametric

- **中文标题**：基于机器学习的参数化各向异性有限应变超弹性材料建模及其在超材料优化中的应用
- **原题**：Material modeling for parametric, anisotropic finite strain hyperelasticity based on machine learning with application in optimization of metamaterials
- **作者**：Mauricio Fernandez, Felix Fritzen, Oliver Weeger
- **出处**：International Journal for Numerical Methods in Engineering, DOI 10.1002/nme.6869（2021）
- **类型**：期刊论文　**等级**：参考　**分组**：柔性与晶格材料的本构建模（方法储备）
- **标签**：`超材料` `晶格` `相对密度` `各向异性` `有限应变` `代理模型`
- **源文件**：`D:\XJTLU\PSP\Numerical Meth Engineering - 2021 - Fernández - Material modeling for parametric  anisotropic finite strain hyperelasticity.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Numerical_Meth_Engineering_2021_Ferna_ndez_Material_modeling_for_parametric_anisotropic_fi__1d1ee5.txt`
- **规模**：33 页 / 提取 98261 字符

## 摘要

机械超材料（开闭孔晶格、泡沫、复合材料等）通常可以用微结构参数描述，例如相对密度、长径比、材料、形状或拓扑参数。为了在有限变形范围内建模其等效本构行为并支持高效的多尺度仿真、设计与优化，本文提出一种基于机器学习的本构模型，并在具有立方各向异性的弹性梁晶格上验证：由于微结构失稳与拓扑变化，这类晶格会表现出强非线性等效行为。作者用微结构仿真生成不同参数下的均匀化应力-变形响应作为训练数据，再把参数化超弹性各向异性本构写成神经网络。

## 要点

- 晶格的关键参数就是相对密度、长径比、形状与拓扑，与你 X 晶格的 x/y/z/p 一一对应。
- 训练数据来自单胞微结构仿真加均匀化，不需要做大量样件实验。
- 微结构失稳使等效行为强非线性，所以线性等效模型不适用。
- 训练好的网络可直接用于超材料的快速优化。

## 与本项目的关系

如果你想做晶格参数到刚度与行程的自动优化，这就是标准做法：先单胞仿真，训练代理模型，再在代理模型上做优化。可直接用于你的 (x,y,z,p) 参数搜索。
