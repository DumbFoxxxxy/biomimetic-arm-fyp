# 36 · vlassis_2020_gdl

- **中文标题**：计算力学中的几何深度学习（一）：各向异性超弹性
- **原题**：Geometric deep learning for computational mechanics Part I: anisotropic hyperelasticity
- **作者**：Nikolaos N. Vlassis, Ran Ma, WaiChing Sun
- **出处**：Computer Methods in Applied Mechanics and Engineering, Vol. 371, 113299（2020）
- **类型**：期刊论文　**等级**：参考　**分组**：柔性与晶格材料的本构建模（方法储备）
- **标签**：`图卷积网络` `加权图描述子` `Sobolev训练` `凸性` `多晶体`
- **源文件**：`D:\XJTLU\PSP\Geometric deep learning for computational mechanics Part I.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Geometric_deep_learning_for_computational_mechanics_Part_I__4c8d13.txt`
- **规模**：39 页 / 提取 130285 字符

## 摘要

把几何深度学习与 Sobolev 训练结合，生成一族能预测训练中未见多晶体等效响应的大变形各向异性超弹性模型。作者指出孔隙率或组分平均取向这类人为设计的微结构度量无法充分表达微结构拓扑，因此引入加权图作为新的高维描述子，用图卷积网络从加权图中提取低维特征并学习其对储能泛函的影响；同时用 Sobolev 训练保证光滑性、避免生成非凸储能函数。

## 要点

- 单纯用孔隙率之类的标量描述微结构不够，拓扑信息需要图来描述。
- Sobolev 训练的目的是防止网络学出非凸的应变能函数。
- 应力通过对能量泛函求方向导数隐式得到，而不是直接拟合。
- 能量凸性是物理可采纳性的要求，不能只看拟合误差。

## 与本项目的关系

当你用代理模型描述 X 晶格时，这条约束很关键：如果代理模型给出的储能非凸，MuJoCo 里就会自发失稳。选代理模型时要检查凸性。
