# 肘部设计阶段交付与停止边界

本轮停在肘部，不继续肩部／肩胛肌肉设计。交付的是可复现的**肘部原型设计候选与试验准备**，
没有将未验证的TPU驱动登记为完成，也没有冻结制造公差或宣称实物可演示。

## 当前方案

人体肱骨、尺骨和桡骨保留解剖外形，工程肘铰链提供屈伸；尺桡骨作为同一前臂刚体。
FJ1486参考肱肌的屈侧分布，FJ1477参考肱三头肌外侧头的伸侧分布。
两块异形TPU均含内部三维晶格、外表面杆网、通道、端部过渡和夹持端。
远端采用被动转动端盖与前臂支架，两件U形金属压片捕获四枚销轴。

连续通道版入口：
`direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers/continuous_channel_candidate/`。
目录中的 `assembly.json` 是当前127件零位装配清单；只有两块TPU更新，其余125个名义刚性零件保持M10-R1。
两块实际肌肉各自在 `FJ1486/`、`FJ1477/` 下，含OBJ/STL、设计参数和打印摆放候选。
原始分段短导管版保留作对照。完整刚性STEP仍在父目录的 `cad/paired_swivel_arm_review.step`，
125个实体已回读；它不含两块TPU，也不是带原生特征树和配合关系的SolidWorks装配。

## 钢丝、反力和连接关系

每块肌肉是一条内部主动钢丝通路。上游电机拉钢丝，近端刚性盖的导管承靠肩接收导管反力；
钢丝穿过近端盖、TPU通道，到远端盖的止挡／压接坯。止挡靠到远端盖肩部后传递拉力，
两端盖之间的TPU发生压缩变形，近端接口连肱骨、远端接口经转动端盖与支架连前臂。
通道随TPU变形，不作为固定在世界坐标的滑轮；远端转角也不能预设为实际平衡结果。
压接坯目前仅为未压接的名义几何，夹持力、导管安装、摩擦和电机端仍需落实实物测试。
这不是钢丝独立绕关节驱动再套一层装饰TPU。

`elbow_interfaces.json` 单列实际材料路径、骨侧安装变换和端盖归属，优先于早期装配中继承的描述字段。
当前候选钢丝直径0.6mm、通道孔径1.7mm、连续壁外径3.2mm；实际线材／导管规格未冻结。
近端导管候选外径2mm；装配前应以实际采购件复核承靠肩和孔径。

## 已核验与未核验

| 项目 | 数字结果 | 结论边界 |
|---|---|---|
| 刚性肘部规定路径 | 0–150°，266个条件姿态及独立回程 | 各相关规定角步长≤1°；不是实际肌肉驱动范围 |
| 当前127件零位 | 两块新TPU的改变配对无实体相交 | 仅零位；不代表整个运动过程无柔性接触 |
| 连续通道肌肉 | 各为单一连通材料、钢丝名义相交为0、STL回读通过 | 不证明孔道可打印、可变形或摩擦合适 |
| 骨侧防脱结构 | 四枚销轴几何捕获及15条选定子装配路径检查 | 防脱力、防松、工具空间及整机装配未实测 |
| 刚性CAD | 125实体STEP有效 | TPU以实际制造网格单独交付 |
| 打印参数样件 | 四类各一份制造网格、每类计划独立打印3件 | 没有实测行程、力、恢复能力或寿命 |
| 完整驱动范围 | 未知 | 不能将150°填入实际TPU驱动范围 |

屈肌大屈曲区是重点风险：原分段导管假设在140°／150°出现导管相互重叠。
简单表面映射出现局部翻转；表面图形变拟合在132°仍有体积损失、钢丝穿透和金属接触，
这些方法均未被接受为TPU物理变形证据。连续管壁版不能继承它们的运动通过结论。
它们也不能证明真实打印肌肉一定失败，需用实测校正材料和接触模型。

## 依据论文补齐的内容

知识库2026论文§3.1明确使用连续管壁维护钢丝通道，并以实际打印测量建立参数映射。
本候选在原短导管之外加入真实连续管壁；屈／伸肌分别增加约653.882／602.128mm³材料。
外半径1.6mm、孔半径0.85mm是本项目待测参数，没有声称来自尚未取得的附表。
`elbow_test_coupons/README.md` 列出稀疏晶格、短导管、连续壁和参考尺寸四类对照及测量方法。
TPU95A为文献候选；打印机、喷嘴、材料品牌及切片条件仍未确定。

## 查看与复现

在 `D:/XJTLU/PSP` 执行：

```powershell
# 当前实际肌肉的零位MuJoCo查看；由用户手动打开窗口
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_continuous_elbow.py'
# 离屏重新生成真实零位图
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_continuous_elbow.py' --render
# 刚性运动演示，明确不显示TPU变形
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_paired_swivel.py' --keepers --motion
# 重新生成连续通道、对比样件
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/build_continuous_muscle_channels.py' --shifted --extensor-reanchor
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/build_elbow_test_coupons.py'
```

当前图为连续通道版 `elbow_neutral.png`；刚性动画为父目录 `rigid_motion.gif`。
试验文件是下一阶段的准备，不在本轮自动开展采购、打印或肩部设计。
