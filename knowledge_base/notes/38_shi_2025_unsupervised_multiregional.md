# 38 · shi_2025_unsupervised_multiregional

- **中文标题**：无需应力数据的深度学习用于多区域超弹性属性识别（综述）
- **原题**：Deep learning without stress data on the discovery of multi-regional hyperelastic properties
- **作者**：Ruike Shi, Haitian Yang, Jianxu Chen, Klaus Hackl, Stephane Avril, Yiqian He
- **出处**：Computational Mechanics, Vol. 76, pp. 117-146（REVIEW ARTICLE）（2025）
- **类型**：综述　**等级**：参考　**分组**：柔性与晶格材料的本构建模（方法储备）
- **标签**：`无监督` `多区域` `异质材料` `残余节点力` `自由能` `软组织`
- **DOI**：10.1007/s00466-024-02591-0
- **源文件**：`D:\XJTLU\PSP\Deep learning without stress data on the discovery ofmulti-regional.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Deep_learning_without_stress_data_on_the_discovery_ofmulti_regional__5af28c.txt`
- **规模**：30 页 / 提取 60302 字符

## 摘要

指出 Flaschel 等人提出的无监督框架不需要难以测量的应力标签，但尚未有工作把它应用到异质材料的多区域本构模型上。本文提出无监督多区域本构学习模型：先定义残余节点力指标与自由能指标来识别非均匀界面或材料分布，再通过在扩展应变不变量层加入开关函数实现多区域本构的联合识别，从而让多个区域的本构模型在同一个网络里同时训练。用数值仿真验证，包括承受区域损伤的软生物组织本构建模。

## 要点

- 残余节点力与自由能指标可用来识别材料界面，无需事先知道分区。
- 开关函数让一个大网络同时表示多个区域的不同本构。
- 对异质、损伤演化的软组织也适用。
- 属于综述性质的系统梳理，附完整方法。

## 与本项目的关系

你的整条手臂会同时包含晶格段、致密端块和韧带，是典型的多区域异质结构。若要用一套参数化本构描述不同部位，这篇给出了分区自动识别的方法。
