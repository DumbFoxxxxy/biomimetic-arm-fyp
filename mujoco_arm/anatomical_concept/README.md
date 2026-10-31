# 独立人体骨肌参考入口

当前几何候选入口为输出目录的`CURRENT_CANDIDATE.md`，最新M5B含两块异形TPU
肌肉及可转向远端端头。局部装配/转向检查已记录，实际钢丝导向和整臂驱动未完成。
下文保留各阶段的历史参考导入记录，14肌肉入口不是当前27肌肉来源清单。

查看当前整臂零位（非驱动仿真）：

```powershell
& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_distal_swivels.py --distal-load --whole
```

使用当前目录中已有的5块骨骼和14个人体肌肉，保留共同坐标和原比例。
没有继承旧关节链、锚点或动作指标。此入口是解剖参考，不是驱动仿真。

在 `D:\XJTLU\PSP` 中运行：

```powershell
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/build_reference.py'
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/repair_inputs.py'
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/audit_input_fit.py'
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_reference.py' --render
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_reference.py'
```

M切换肌肉，B切换骨骼，R复位相机。

输出位于 `mujoco_arm/exports/anatomical_concept/`：
INVENTORY.md 包含标签；inventory.json包含坐标、单位、来源哈希和模型清单；
repair_report.json记录清理过程；repaired/保留闭合体积副本。
原始14肌肉中11个包含35个不与主体连通的单三角面，最长边约0.357mm。
副本删除孤立零厚度三角面，不移动其余顶点、不覆盖源文件；重新导入后验证体积拓扑。

input_fit_audit.json记录全部70对骨肌的输入相交和间隙。相交不能直接判为临床错误，
也不能作为合法动态接触直接通过。不得从最近骨骼推断真实起止点，或从PCA推断纤维方向。
现有数据尚未提供完整的肌群、纤维方向或附着区；补充模型可继续加入同坐标入口。

模型来源和许可按OBJ文件头保留：BodyParts3D, (c) The Database Center for Life Science,
CC Attribution-Share Alike2.1 Japan。原许可URL记录在清单中。

## 解剖体积晶格筛查

```powershell
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/lattice_screen.py'
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/lattice_screen.py' --surface-cage
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/analyze_lattice_screen.py'
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_reference.py' --cage
```

两组筛查覆盖全部14肌肉。晶格在真实体积内裁剪，外网格是源表面简化后的三角网格，
尚未采用论文的曲面坐标构造外图案。PCA只定铺网格坐标，不代表纤维方向。
还没有钢丝孔道、肌腱端部、受力或动作验证。材料比例较高的薄片区不应据此冻结。
对照见 LATTICE_SCREEN.md、lattice_comparison.json 和 lattice_detail_comparison.png。

## 27肌肉附着关系与候选区域

当前输入已经扩展为27肌肉＋5骨；原14肌肉的晶格筛查结果不自动推广到新增模型。

```powershell
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/extract_attachment_regions.py'
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_concept/check_attachment_regions.py'
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_reference.py' --attachments --muscle FJ1486
```

附着结构的解剖对应关系及逐项来源在attachment_regions/ATTACHMENT_REVIEW.md。
绿色是起始骨面的接近区域样点，橙色是止端；P切换标记，M/B切换肌肉/骨骼。
显示球体大小仅为可视化符号，不能解释成锚点或钻孔直径。
候选区域以肌肉面片中心到指定骨面的距离≤2mm筛选，5mm仅作敏感性诊断。
这是按照解剖结构约束做几何筛查，不是临床标志自动识别；所有点均未冻结。

目前肩胛、胸廓相关肌群仍缺固定支撑的工程定义；支持参考的源模型可用性见
support_reference_catalog.json。肱二头肌短头止端与桡骨的几何距离约36.45mm，
需要核对共同肌腱/模型分段表达或设计被动肌腱段，不直接挪动骨骼来制造贴合。
三束斜方肌使用全肌附着范围作初筛，不能把每一束都强制连接到所有止端骨骼。

## 支撑参考与新关节候选

用户提供的BodyParts3D目录中已补入39个支撑参考，保存在 `骨骼/支撑参考/`，
包括32个骨参考及7个肋软骨。源文件保持不变，副本来源/哈希在该目录sources.json。
这39件是解剖背景，不意味着必须全部打印，也不是已经设计好的机器人机架。

```powershell
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_reference.py' --supports
```

T切换支撑参考。附着候选查询现已包含这些参考面，项韧带与腹部腱膜仍需明确工程接口。
参考面上的几何接近区域仍不能解释成已确认的骨性标志或机械锚点。

新的关节候选在joint_candidates.json，独立从当前原始骨骼拟合SC/AC/GH/肘中心，
并用真实胸骨面参与SC参考定位。它不是旧链或旧关节点的拷贝。
joint_model/arm.xml试建SC3＋AC3＋GH3＋肘1的球副/铰链约束模型，nv=10，
只是运动学候选，没有制造防脱件或驱动器，不是已完成工程关节。
肘拟合以尺桡骨保持源相对姿态作单项0–150°骨面试探，未决定最终旋前旋后设计。
肩部路径、整体支撑干涉、驱动可控性和柔性肌肉运动尚未验证，不能从模型加载成功
或数学自由度数推断已实现运动。候选模型中的肌肉保持零位，运动展示必须隐藏它们。

## 当前M4可拆力臂候选

`elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15/`
已有111件中立位实际网格、109件刚性STEP、12件新实体及8条安装路径检查。
13条骨架路径已重算，但钢丝/端帽的轴向预测路径存在穿透，驱动尚未完成。
详细结论在`SIDE_LEVER_REVIEW.md`；两件压接止挡仅为包络，尚未固定。

```powershell
& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_side_lever.py --detail
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/review_side_lever.py
```

`rebuild_side_lever.py`记录完整生成顺序，包含耗时的整臂路径复查。
展示保持TPU零位，不用假变形动画或导向控制点冒充实际驱动。
