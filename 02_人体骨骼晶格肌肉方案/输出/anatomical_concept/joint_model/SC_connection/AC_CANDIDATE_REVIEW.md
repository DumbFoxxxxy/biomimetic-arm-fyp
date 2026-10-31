# 肩锁实体连接候选对照

两方案均保留原拟合AC中心，不移动GH中心。大方案球半径10mm、球窝外半径13mm；紧凑方案分别6mm、9mm。两者保留完整球头、骨侧连续过渡、球窝、分体保持盖、M3螺栓/螺母/垫圈，均为12个连通实体。紧凑方案螺栓名义M3×12，大方案M3×16；螺纹未建模。

当前后续设计优先使用AC_compact，尚未冻结。它恢复了大方案损失的盂肱抬升路径，代价是孤立AC上回旋探针的范围更小。不是在全部动作上优于大方案，也尚未验证缩小球头的强度。

## 运动采样（数字几何）

下表为从零位连续可达的指定路径前缀；角度为输入旋转向量的模，联合抬升列为实际肱骨方向相对零位的夹角。不是人体临床ROM或已承诺的新验收标准。

| 探针 | 大方案° | 紧凑方案° | 紧凑方案首个限制 |
|---|---:|---:|---|
| elbow_flexion | 150.000 | 150.000 | 请求路径内未检出 |
| GH_flexion | 151.938 | 154.562 | FJ3384 / FJ3368 |
| GH_abduction | 70.938 | 119.312 | FJ3384 / FJ3368 |
| GH_axial_positive | 90.000 | 90.000 | 请求路径内未检出 |
| GH_axial_negative | 90.000 | 90.000 | 请求路径内未检出 |
| SC_elevation | 27.500 | 27.500 | FJ3362 / SC_cap_left |
| SC_protraction | 1.688 | 1.688 | FJ3384 / FJ3342 |
| SC_retraction | 8.625 | 8.625 | FJ3362 / FJ3336 |
| SC_rotation_X | 30.000 | 30.000 | 请求路径内未检出 |
| AC_upward | 47.125 | 35.500 | FJ3362 / AC_cap_left |
| AC_rotation_X | 30.000 | 30.000 | 请求路径内未检出 |
| AC_rotation_Z | 30.000 | 30.000 | 请求路径内未检出 |
| coordinated_raise | 97.536 | 164.049 | FJ3384 / FJ3368 |

## 肌肉接口影响

候选骨面区域不是已冻结锚点。以下区域有采样点距工程修改表面超过1mm，必须重映射；不能直接复用原解剖附着面：

| 肌肉编号 | 骨面 | 超过1mm点数 | 最大偏差mm |
|---|---|---:|---:|
| FJ1467 | scapula | 20 | 2.370 |
| FJ1468 | clavicle | 141 | 4.070 |
| FJ1506 | scapula | 3 | 1.683 |
| FJ1554 | clavicle | 14 | 4.625 |
| FJ1554 | scapula | 98 | 2.763 |

## 验证边界

全部活动零件依其所属刚体运动：肩锁球窝、保持盖和附件随肩胛运动，没有误固定到世界坐标。15个活动实体与42个固定实体检查，邻接骨没有整对排除。角度采样≤1°、边界细化≤0.1°，实体交叠阈值0.001mm³，最小间隙查询封顶3mm。7个肋软骨没有作为刚性实体验收；柔性接触仍未知。回程重新计算，演示仅播放首次碰撞之前的连续前缀。

OBJ/STL回读、各12个独立STEP实体和复合STEP分别检查。STEP为三角面BRep，可导入CAD编辑，但不含原生特征树和装配配合。六方向球体直拉均被固定球窝/保持盖阻挡：仅是有限方向几何检查，不能证明所有逃脱路径或抗拉承载。

SC/AC实体候选已加入模型；GH与肘仍是理想约束，不能据此宣称完整工程骨架已完成。肌肉在运动中隐藏，未模拟TPU、钢丝收缩或力传递。实际机架、装配工具可达性、打印公差、紧固预紧、防脱强度及寿命仍待验证。

## 复现入口

从项目根目录运行（Python环境见已有SC说明）：

```powershell
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/build_AC_connection.py --compact
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/audit_joint_model.py --AC --compact
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_SC_assembly.py --AC --compact
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_AC_capture.py --compact
& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_cad.py --AC --compact
& 'D:/Python3.10.6/python.exe' mujoco_arm/anatomical_concept/export_SC_assembly.py --AC --compact
& 'D:/MujocoEnv/Scripts/python.exe' mujoco_arm/anatomical_concept/check_joint_audit.py --AC --compact
& 'C:/Python314/python.exe' mujoco_arm/anatomical_concept/view_reference.py --AC --compact --supports
```

去掉--compact可复现大方案。视窗按1–5播放已检查路径，M切换肌肉、T切换胸廓参考、R回到零位；运动肌肉隐藏避免误示驱动。
