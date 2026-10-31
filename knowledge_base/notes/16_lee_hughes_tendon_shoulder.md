# 16 · lee_hughes_tendon_shoulder

- **中文标题**：带被动稳定的大工作空间肌腱驱动肩关节
- **原题**：Tendon-Driven Robotic Shoulder Joint With Passive Stabilization for Large Workspace
- **作者**：Sudong Lee, Kieran Gilday, Josie Hughes
- **出处**：IEEE/ASME Transactions on Mechatronics（2024）
- **类型**：期刊论文　**等级**：重要　**分组**：肩带机构与机器人肩部实现
- **标签**：`肌腱布线` `被动螺旋肌腱` `弹簧` `脱位特性` `正负99.87度`
- **源文件**：`D:\XJTLU\PSP\Tendon-Driven_Robotic_Shoulder_Joint_With_Passive_Stabilization_for_Large_Workspace.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Tendon_Driven_Robotic_Shoulder_Joint_With_Passive_Stabilization_for_Large_Workspace__6e3b5d.txt`
- **规模**：11 页 / 提取 64132 字符

## 摘要

提出一种独特的肌腱布线法：把直线驱动肌腱、被动螺旋肌腱和弹簧结合起来，即使在大角度极限工况下也能防止关节脱位，并保证任意位形下可控。改变肌腱数量与布线方式可以调节脱位特性，被动弹簧连接带来固有柔性和恶劣工况下的稳定性。肌腱在这里同时模仿人的韧带与肌肉，以实现结构上与功能上的相似。样机在滚转与俯仰方向达到正负 99.87 度工作空间，并能柔顺承受高冲击。

## 要点

- 直线驱动肌腱管运动，被动螺旋肌腱加弹簧管稳定与防脱位，两类肌腱分工明确。
- 脱位特性是一个可设计量，通过肌腱数量与布线方式调节。
- 被动元件提供固有柔性，冲击下不依赖控制就能保护结构。
- 工作空间正负约 100 度，与人的肩部量级相当。

## 与本项目的关系

你的球铰盂肱关节在仿真里容易绕轴乱转甚至脱位，这篇给的是纯被动的解法（额外加被动肌腱与弹簧），比单靠加大 joint stiffness 更符合生物学。
