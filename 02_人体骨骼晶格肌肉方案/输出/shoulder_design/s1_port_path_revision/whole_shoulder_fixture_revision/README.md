# S1夹具修订与首次实际刚性运动检查

2026-11-02。本候选的零位穿线贯通，但运动检查未通过，不能作为可运动或打印定版。

## 实际夹具修订

保持原立杆、横梁、骨侧接口位置及底座40×30×6mm尺寸。外侧底座沿源坐标−X移9mm；在0至20mm的1mm步长候选中，首次满足与20个现有非夹具实体均无超阈相交且底座间隙至少0.5mm。间隙查询上限0.5mm，不表示每一处实际间隙都恰为0.5mm，也未宣称全局最优位移。

原导出脚本将杆、梁、底座分别存STL，没有可拆接头。现明确采用两个一体支架候选：外侧立杆与底座合并，内侧立杆、横梁与底座合并。三处原交叠落实为真实材料连接，不是将碰撞整对加入白名单。两个实际STL均单正连通、水密且回读通过。打印支撑、刚度、夹持反力和安装路径尚未验证。

夹具生成与试选记录见`../fixture_revision/check.json`，脚本为`mujoco_arm/tools/revise_shoulder_fixture.py`。它读取上轮25件装配并核对输入哈希，不修改原文件。

## 22个现有实体的静态复核

`actual_assembly_check.json`包含3骨、4整体TPU、10端盖、3止挡坯和2一体夹具。四组肌肉选型沿用`../whole_shoulder_static/README.md`。

- 六条完整钢丝对全部22件的相交均为0mm³。
- 夹具对全部非夹具实体、两个夹具彼此的相交均为0mm³。
- 超过0.05mm³的实体交叠为0对；**不等于数学零干涉**。最大剩余非零值为冈下肌/肩胛约0.002936mm³，属于仍待定位的贴骨/导出残差。
- TPU壁厚、内部空腔、五金和柔性变形仍未放行。

## 固定GH中心的刚性运动未通过

`rigid_motion_check.json`读取这22件的真实网格，检查肱骨、肱骨侧4个盖及3个止挡的理想随动刚体，分别对固定肩胛、锁骨、其余端盖及两夹具逐对检查。TPU不能作为刚体运动结果来验收，故其变形保留未验证；并非以忽略肌肉碰撞获得运动通过。

SC、AC固定，GH中心强制固定，端盖和止挡假定已正确安装随骨移动；实际开放浅窝、防脱和紧固尚未证明。每方向从零位按1°采样到首次相交，再以0.1°细化。交叠判据为0.001mm³，不是数学全姿态证明。

| 校准轴/符号 | 末个未超阈采样° | 首个超阈采样° | 首个限制 |
|---|---:|---:|---|
| flexion + | 0.1 | 0.2 | 肱骨—FJ1467起端盖 |
| flexion − | 0 | −0.1 | 同上 |
| abduction + | 0 | 0.1 | 同上 |
| abduction − | −0.1 | −0.2 | 同上 |
| longitudinal + | 0.1 | 0.2 | 同上 |
| longitudinal − | −0.3 | −0.4 | 同上 |

符号按源坐标校准轴，完整向量见JSON，不直接等同于已核验的解剖内/外旋标号。以上是当前模型的首次刚性超阈采样，**不是实测可达范围，也不是K1–K7完成结论**。

关键问题：FJ1467中束起端盖此前只按零位肱骨局部修形，几乎没有转动间隙。后续优先按运动包络重设计该端盖及骨侧接口，再同步重建钢丝/TPU配合。不得关闭该骨—盖碰撞检查来获取通过，不降低原动作目标。

## 复现

仓库根目录，项目环境需`.deps`和`mujoco_arm/tools`在Python搜索路径：

1. `revise_shoulder_fixture.py`生成两个新夹具。
2. `check_shoulder_manufactured_assembly.py`重复传入四组`--unified-report`（见上轮README），再加`--fixture-report <fixture_revision/check.json绝对路径>`。
3. `check_shoulder_rigid_motion.py`生成条件刚性运动报告。

静态查看：

```powershell
& 'C:/Python314/python.exe' 'mujoco_arm/tools/view_shoulder_manufactured.py' --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_port_path_revision/whole_shoulder_fixture_revision/actual_assembly_check.json'
```

加`--render`生成下图。场景DOF为0，仅为制造网格静态查看；运动证据来自上述逐姿态实体检查，不来自图片。

![22件真实网格的静态视图](available_parts.png)
