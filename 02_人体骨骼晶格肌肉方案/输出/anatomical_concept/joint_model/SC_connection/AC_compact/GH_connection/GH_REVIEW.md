# G1 盂肱实体连接候选（未定版）

以SC/AC紧凑候选为基础，GH中心保持源拟合位置。将肱骨头局部改为半径16mm的工程球头，颈部半径4.5mm；球窝内半径16.4mm、外半径19.5mm。两条颈部运动让位槽使用半径5.3mm的扫掠圆柱构造，保留两个独立防脱盖，每盖采用名义M3×12螺栓、螺母和两片垫圈。螺纹未建模。

骨侧过渡、球窝、防脱盖及附件共12个连通实体。盖板/附件随肩胛（AC刚体）运动，球头及肱骨随GH刚体运动。全模型25个活动实体对42个固定骨参考/SC附件检查；没有排除相邻骨对。肘仍为理想铰链，未提供实体销轴和防脱。

## 修改与装配修正

初版新增连接与既有锁骨及AC保持盖发生零位相交，已保存于iterations/initial_neutral_conflict/。修订版在真实相交区域开局部间隙，没有关闭碰撞对；之后重新检查零位与全路径。球窝曾有约529mm³的悬空块，加入明确的外侧连接支撑后成为单一实体。局部间隙修正产生0.00406mm³细碎残片，通过记录中心/半径的局部轮廓修整去除；没有删掉实际连接件。具体工具区域、连接端点和几何增减量见assembly.json。

肱骨相对原骨形移除34892.1mm³，肩胛相对AC基线移除675.9mm³并增加6887.2mm³。此方案明显重做了关节近端，不能把它称为未改动的人体关节。原输入均保留。

## 指定路径的连续前缀

每个探针从零位开始，正反向采样≤1°，首次干涉边界细化≤0.1°。未检出是指交叠体积不超过0.001mm³；间隙查询封顶3mm，不是数学上的全姿态证明或打印性能验收。组合运动测量实际肱骨方向，不把关节角相加。表中间隙为已接受前缀内完整1°采样点的最小值，不包含未存入路径数组的细化端点。

| 探针 | 请求输入° | 连续前缀输入° | 肱骨方向变化° | 采样最小间隙mm（封顶3） | 首个限制 |
|---|---:|---:|---:|---:|---|
| elbow_flexion | 150.000 | 150.000 | 0.000 | 0.0000 | 请求路径内未检出 |
| GH_flexion | 160.000 | 145.938 | 145.938 | 0.0000 | FJ3384 / FJ3368 |
| GH_abduction | 160.000 | 123.938 | 123.938 | 0.0000 | FJ3384 / FJ3368 |
| GH_axial_positive | 90.000 | 90.000 | 0.000 | 0.0000 | 请求路径内未检出 |
| GH_axial_negative | 90.000 | 90.000 | 0.000 | 0.0000 | 请求路径内未检出 |
| SC_elevation | 30.000 | 27.500 | 27.497 | 0.0000 | FJ3362 / SC_cap_left |
| SC_protraction | 30.000 | 1.688 | 0.292 | 0.0000 | FJ3384 / FJ3342 |
| SC_retraction | 30.000 | 8.625 | 1.492 | 0.0000 | FJ3362 / FJ3336 |
| SC_rotation_X | 30.000 | 30.000 | 29.540 | 0.0000 | 请求路径内未检出 |
| AC_upward | 50.000 | 35.500 | 35.496 | 0.0000 | FJ3362 / AC_cap_left |
| AC_rotation_X | 30.000 | 30.000 | 29.540 | 0.0000 | 请求路径内未检出 |
| AC_rotation_Z | 30.000 | 30.000 | 5.138 | 0.0000 | 请求路径内未检出 |
| coordinated_raise | 120.000 | 120.000 | 164.994 | 0.0000 | 请求路径内未检出 |

前屈与外展首次完整角度碰撞样本已分区：均发生在新增工程几何，不在保留的原肩胛骨实体中；详见limit_localization.json。后续仍可通过连接外形/支撑路由优化，当前未把范围降低改写为“满足原目标”。

注意：全体间隙最小值为0，零位最近对为工程肩胛/FJ3384与AC_cap_right（同属AC刚体）。虽然未检出超过体积阈值的穿透，这仍不是正的制造间隙证明，不能据此直接冻结公差或打印装配。完整模型仍保留该对的碰撞检查。

## 肌肉附着影响

候选解剖骨面不是冻结的工程锚点。以下区域的原候选点距新骨面超过1mm；其中肩袖近端偏差较大，必须设计新的骨侧连接端、过渡及短肌腱，不能把源肌肉原位叠加当成已经连上球头。

| 肌肉编号 | 骨面 | 超过1mm点数 | 最大偏差mm |
|---|---|---:|---:|
| FJ1467 | scapula | 21 | 2.370 |
| FJ1467 | humerus | 2 | 11.887 |
| FJ1468 | humerus | 2 | 8.096 |
| FJ1500 | humerus | 228 | 11.371 |
| FJ1504 | humerus | 288 | 8.106 |
| FJ1506 | scapula | 2 | 1.683 |
| FJ1506 | humerus | 270 | 12.140 |
| FJ1507 | humerus | 4 | 1.157 |
| FJ1508 | humerus | 96 | 11.271 |
| FJ1554 | scapula | 98 | 2.763 |

## 验证范围

12个独立STEP和12实体复合STEP均通过回读，网格与场景的36项位置检查通过。STEP为三角面BRep，保留独立实体，不含SolidWorks原生特征树/装配配合。球头26方向直拉探针在0.9mm内均遇到球窝/盖；这是有限直线检查，不是所有曲线路径的防脱证明。真实螺栓预紧、打印间隙、工具可达性、强度及寿命尚未验证。

保留更多骨形的GH_preserve_head比较候选有独立文件夹。它保留了一个原肱骨局部，并以弧形连接连续接回骨干，但在很小的前屈/外展范围即相交；肱骨件STEP导出还报Invalid BRep。该候选未被选用，不继承G1的CAD或运动通过结论。

TPU肌腹、钢丝通道、端部载荷路径和驱动尚未建立，本阶段演示隐藏肌肉。7个肋软骨不作为刚性骨验收，尚未验证其柔性接触。下一步是肘的实体连接，以及骨侧肌肉接口和成形肌腹/通道；不是继续把理想约束称为已制造关节。

## 复现

```powershell
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/build_GH_connection.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/audit_joint_model.py --GH
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_SC_assembly.py --GH
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_GH_capture.py
& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_cad.py --GH
& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_assembly.py --GH
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_joint_audit.py --GH
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/diagnose_GH_limits.py
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/write_GH_review.py
& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_reference.py --GH --supports
```

加--detail显示关节细节；按1–5播放已采样前缀。渲染用--render，动作GIF由render_joint_motion.py --GH生成。
