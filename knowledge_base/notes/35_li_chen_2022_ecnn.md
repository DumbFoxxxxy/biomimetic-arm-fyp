# 35 · li_chen_2022_ecnn

- **中文标题**：基于平衡方程的卷积神经网络用于超弹性材料的本构建模
- **原题**：Equilibrium-based convolution neural networks for constitutive modeling of hyperelastic materials
- **作者**：L.F. Li, C.Q. Chen
- **出处**：Journal of the Mechanics and Physics of Solids, Vol. 164, 104931（2022）
- **类型**：期刊论文　**等级**：参考　**分组**：柔性与晶格材料的本构建模（方法储备）
- **标签**：`无监督学习` `等效力平衡` `超弹性` `本构` `应力不可测`
- **DOI**：10.1016/j.jmps.2022.104931
- **源文件**：`D:\XJTLU\PSP\Equilibrium-based convolution neural networks for constitutive.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Equilibrium_based_convolution_neural_networks_for_constitutive__6718ea.txt`
- **规模**：15 页 / 提取 56317 字符

## 摘要

传统机器学习本构方法需要大量多轴应力-应变曲线，但应力只能在应力场简单均匀的有限试件上测出，实验上很难获得；而且多数方法是有监督学习、忽略已知物理定律，导致可解释性与泛化性受限。本文提出基于平衡方程的卷积神经网络 ECNN：训练只输入应变与外载荷，把难以测量的应力当作内部变量，并把平衡方程嵌入网络结构作为约束，因此是无监督的，仅需一个非均匀变形试件就能生成大规模数据。

## 要点

- 把物理定律（平衡方程）写进网络结构，而不是只当作训练数据。
- 训练不需要应力标签，只需要应变和外力。
- 单个非均匀变形试件就能生成大量训练样本。
- 属于物理信息机器学习路线。

## 与本项目的关系

如果你要实测 TPU 打印晶格的本构，应力很难测，这篇的方法正好绕开该难点。与你已决定的等效力耦合路线在思路上一致：宁可改模型内部，也不造无法标定的模型。
