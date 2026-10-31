# 37 · vlassis_2023_diffusion

- **中文标题**：用去噪扩散算法逆向设计具有精细非线性材料属性的微结构
- **原题**：Denoising diffusion algorithm for inverse design of microstructures with fine-tuned nonlinear material properties
- **作者**：Nikolaos N. Vlassis, WaiChing Sun
- **出处**：Computer Methods in Applied Mechanics and Engineering, Vol. 413, 116126（2023）
- **类型**：期刊论文　**等级**：参考　**分组**：柔性与晶格材料的本构建模（方法储备）
- **标签**：`扩散模型` `生成模型` `逆向设计` `微结构` `CNN代理`
- **源文件**：`D:\XJTLU\PSP\Denoising diffusion algorithm for inverse design of microstructures.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Denoising_diffusion_algorithm_for_inverse_design_of_microstructures__839b06.txt`
- **规模**：20 页 / 提取 75175 字符

## 摘要

用去噪扩散概率模型逆向设计具有精细可调非线性属性的微结构：通过学习马尔可夫扩散过程的逆过程，生成大量拓扑各异、本构响应接近指定非线性行为的候选微结构。为筛掉不满足容许范围的样本，训练一个卷积神经网络代理模型来替代高保真有限元仿真。实验基于开源 mechanical MNIST 数据集，说明该算法不仅能做非线性等效介质的逆向设计，还学到了结构-属性映射。

## 要点

- 生成模型负责提出候选微结构，CNN 代理负责快速筛选，两级分工。
- 目标是给定力学属性反推几何拓扑，与你给定肌肉刚度反推晶格图案是同一类问题。
- 可在潜空间加入额外拓扑或几何约束。
- 训练与验证用公开数据集，便于复现。

## 与本项目的关系

这是你肌肉形状与晶格密度自动设计的终极工具形态。当前不必使用，但若后面要扩展到自动生成晶格图案，这是入口。
