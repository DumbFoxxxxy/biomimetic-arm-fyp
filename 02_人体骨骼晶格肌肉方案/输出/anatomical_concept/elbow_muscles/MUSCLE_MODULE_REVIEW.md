# 肘部异形TPU肌肉制造候选

已生成肱肌参考FJ1486和肱三头肌外侧头参考FJ1477两块候选。它们包含实际三维晶格、参数化外表面菱形网、多个短管状导向单元、穿线孔、TPU端垫、刚性端板和金属附件包络。每块20个组件，TPU肌腹自身为单连通实体；不是长方体占位或纯渲染外壳。当前模型尚未安装到骨架，未验证收缩、驱动或制造性能。

## 与论文方法的对应

依据本地论文《Accelerating Musculoskeletal Robotics Through Parametric Design and 3D‐Printed Flexible Structures》第2.3、4.2节及图4，并对照其PSM前作图4：保留肌腹体积→晶格/外网格→管状导向→端部的组成思路。参考论文p=2mm及15/10/20mm尺寸建立候选；本脚本局部轴依次为横向15、非主要受载方向20、工程收缩方向10mm。

截面来自真实OBJ沿肱骨工程方向的三角面截交线，17个截面拟合周期三次B样条，轴向也用三次B样条；相当于权重为1的曲面参数化。内网格采用规则八分支体单元后裁切，外网格沿参数面生成。尚未逐式复现论文式(3)的体积参数坐标映射，不能直接继承论文的应变、刚度或耐久数值。

钢丝路径通过固定端点、截面内候选点的短步长路径确定，段间在真实导向管处折转。该路径是工程设计，未经解剖肌纤维标定；不采用PCA代替肌纤维。管轴沿相邻直线跨度的角平分线，保留名义管壁。压接后的实际摩擦和受压变形仍需试验。

| 项目 | FJ1486 | FJ1477 |
|---|---:|---:|
| 设计肌腹轴向跨度mm | 180.3393 | 205.9385 |
| TPU材料体积mm³ | 26467.5619 | 32817.2346 |
| 肌腹包络内材料比例 | 0.2102 | 0.1718 |
| 拟合包络位于原肌肉外的体积比例 | 0.0155 | 0.0329 |
| 最大导向折转角° | 7.8704 | 8.9141 |
| 导向管名义壁厚下界mm | 1.0448 | 1.0305 |
| 0.6mm钢丝名义最小间隙mm | 0.1487 | 0.1487 |

拟合前裁去源肌肉轴向两端各10%，用工程端垫和过渡取代。上表“位于原肌肉外的比例”只衡量拟合包络的外溢，不是总体积保留率或完整表面误差证明；局部端垫加宽、表面杆件和连接加强材料另有设计变化。

## 端部与传力

主动功来自电机与钢丝。近端端板固定后提供导管反力接口；主动钢丝经过短导向管，在远端金属压接止挡处终止，拉近两端并压缩晶格。远端端板的横向孔连接被动短钢丝肌腱，再向关节输出。TPU提供可变形体积、弹性和缓冲，当前没有假定其能自行主动收缩。

每块模块：1件TPU肌腹、2件刚性端板、4个名义M3×12螺栓、4个M3螺母、8片垫圈、1个压接止挡，共20件。端垫直径22mm、厚4mm，M3孔距14mm；TPU孔道切削半径0.85mm，检查的保留目标半径0.8mm；0.6mm钢丝仅为首轮候选规格。刚性端板含导管孔/台阶、止挡座及被动肌腱横孔。金属止挡必须实际压接后才可能夹住钢丝，几何包络的通孔不能证明夹持力。

新检查曾发现螺母终态可放置但不能装入；已保留iterations/before_nut_access，并在TPU端部开侧向装入槽。现已检查两端板、四组螺母/内垫圈、四螺栓及止挡共11条1mm采样路径。螺纹啮合、工具可达性、预紧及TPU端部蠕变未验证。

局部分离的生成枝条以记录端点的真实短连接修复，导向孔和螺母装入槽保留为空腔。只对生成过程产生的微小封闭孔腔和松散碎片做记录的局部修整，未删除源解剖模型。各制造STL/OBJ回读均通过；四个刚性端板STEP回读通过。TPU仍以制造网格和可复现生成参数交付。

## 接入骨架仍需完成

解剖候选端点直线牵引仅作为对照筛查，不代表上述弯曲肌腹的实际驱动路线，也未给不存在的滑轮施加约束：

| 参考肌肉 | 直线端点总行程mm | 首次力臂符号不符° | 与骨架相交的姿态数/151 |
|---|---:|---:|---:|
| FJ1486 | 51.893 | 此探针未出现 | 15 |
| FJ1477 | 22.069 | 85 | 133 |

直拉路径不能冻结为最终走线。理想半径8/10/12/14mm卷绕的150°行程分别约20.94/26.18/31.42/36.65mm，只是后续尺寸对照；没有据此生成“已验证的滑轮驱动”结论，也不把行程比例当成TPU可用压缩率。

| 模块 | 原位置刚性干涉对数 | 沿候选外移方向首次采样刚性清空mm |
|---|---:|---:|
| FJ1486 | 10 | 14 |
| FJ1477 | 3 | 6 |

外移结果仅用于设计骨侧支架起点，尚未采用为装配位置；检查只覆盖零位。TPU与骨接触列为变形/摩擦待查，既不按刚性穿透一律判失败，也未当作已通过。近端支架、导管/电机反力固定、远端实际肌腱导向和关节范围内的力臂/行程仍待实现。肩部与肩胛肌群的制造、驱动和真实机架仍属于整臂目标。

## 文件和复现

FJ1486/、FJ1477/内的TPU_muscle.stl为制造几何；proximal_rigid_cap.stl和distal_rigid_cap.stl为刚性端板。其余金属STL是硬件包络，不作为TPU打印件。hard_caps/cad/含端板STEP。print_candidates/为另外旋转/落床的候选文件，变换有记录；这些文件只用于打印准备，装配仍使用父目录公共坐标。切片支撑、材料/孔配合和打印质量未验证。

MuJoCo模块查看器读取与制造STL同源的OBJ，以支持较多三角面；静态钢丝仅示意真实导向之间的名义路径。图片不包含已验证的收缩变形或驱动动画。

```powershell
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/parameterize_elbow_muscles.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/build_shaped_elbow_muscles.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_shaped_muscles.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_muscle_end_assembly.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/screen_elbow_excursion.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_muscle_placement.py
& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_muscle_caps.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/prepare_muscle_prints.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/write_muscle_review.py
& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_shaped_muscles.py --render
& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_shaped_muscles.py --id FJ1486
```

先完成骨侧适配、实际关节走线及打印样件检查，再决定肌腹尺寸和可用范围。测试记录模板中的实测项保持null，未从论文或数字形状填入“实测”结果。
