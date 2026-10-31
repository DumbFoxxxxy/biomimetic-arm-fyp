# S1胸锁连接与新关节候选

本路线独立于旧版本。输出位于：
`mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/`。

原始锁骨与胸骨柄零位实体重叠约319.292mm³，因此不能直接当成可打印的刚体装配。
S1在SC端局部重构球头/颈部和座体，加入左右可拆保持盖及M3紧固件。
原始输入文件保持不变，切除/增加体积、输入哈希和零位检查均有记录。

名义尺寸：球半径10mm、座内半径10.35mm、座外半径13mm、颈半径4.5mm，
保持盖分缝0.2mm、孔径3.4mm。锁骨端22mm范围内作局部重构，其他区域保留原形。
球头和座的间隙意味着真实零件存在微小平移，理想球副模型未模拟这些游隙。
表面简化容差0.01mm；仅清理数值零体积壳，没有删去有效的骨连接来制造通过。

四件聚合物候选：锁骨球头件、胸骨柄座体、左右保持盖。另配两只M3×16、
两只M3螺母和四只垫片；金属件CAD仅为名义包络，螺纹未建模，不是打印螺母图纸。
先装球头，再把两个保持盖从颈部侧面装入，以螺钉连接后座。不能把球头压过保持口。
实际配合、螺钉预紧、防松、机架安装和承载均待验证。

12件零位逐对无检测到的穿透，OBJ/STL回读及12个STEP回读均通过。
`cad/SC_review_assembly.step` 是保持12个独立实体的分面BRep查看装配，
不是原生SolidWorks配合装配或参数特征树。
六方向刚性球头直线移动会在0.2–0.8mm内被阻挡；这是名义几何探测，不能代替防脱载荷试验。

运动检查包含5个移动骨件、32个骨性支撑参考及10个固定SC附件。没有屏蔽相邻移动骨。
7个肋软骨仅作参考，肌肉变形、AC/GH/肘部实际保持器未完成，不作为通过项。
实际四元数间隔已独立检查≤1°（浮点舍入误差除外），首边界细化0.1°。
前屈/外展使用与肱骨长轴正交的工程轴；此前绕全局X轴160°的结果不能作为前屈验收。
本轮名义连续前屈154.56°、外展119.31°、肘150°、联合肱骨方向变化164.05°。
这只验证了指定路径及已包含的实体，不是完整制造装配或TPU驱动通过。
SC单轴前后运动较小，仍须设计SC/AC耦合以避免肩胛/锁骨进入肋骨。

```powershell
& 'C:/Python314/python.exe' 'D:/XJTLU/PSP/mujoco_arm/anatomical_concept/view_reference.py' --SC --supports
```

1肘、2前屈、3外展、4联合抬升、5SC抬升：播放已检查的连续无碰撞前缀；
R停止并回零，T切换支撑参考，M仅在零位查看肌肉参考。运动时自动隐藏未变形肌肉。
加 `--detail` 查看SC局部，加 `--render` 导出静态对照。
joint_review_motion.gif只展示规定刚体姿态，不是电机/肌肉仿真。

生成脚本顺序：repair_inputs.py --supports、fit_joint_candidates.py、build_joint_model.py、
build_SC_connection.py、audit_joint_model.py --SC。CAD导出使用Python3.10运行
export_SC_cad.py及export_SC_assembly.py。最后运行check_SC_assembly.py、
check_SC_capture.py、check_joint_audit.py和write_joint_review.py复核。
