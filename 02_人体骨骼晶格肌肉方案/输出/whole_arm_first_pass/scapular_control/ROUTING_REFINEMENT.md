# 前锯肌局部走线修整

基于现有整臂做 K6 对照。保留钢丝起点、原胸廓导向、肩胛锚点和原关节轨迹。在原导向与肩胛锚点之间筛选一处带骨面支承的新导向。检查对象为 0.60 mm 名义钢丝、导向安装包络、原整臂刚件和其余局部钢丝。

## 新增导向的候选结果

新增胸廓侧导向落在右第六肋 `FJ3344`，端口源坐标为 **(−92.17, −11.20, 1235.26) mm**，名义支承离骨面约 **2.82 mm**。钢丝顺序为胸廓起点、原胸廓导向、新导向、原肩胛端。新导向在本轮 K6 中的峰值反力约 **6.04 N**。

| 项目 | 检查结果 |
|---|---|
| 几何范围 | 原第六肋钢丝干涉绕开。K6 扫描通过 **27.50°**，下一样点 **27.516°** 的首限为原锁骨与胸锁端盖相交。 |
| 正向力臂 | 两条主动钢丝到 **11.797°** 仍保持正向。**11.813°** 样点的前锯肌等效力臂已变为负值，上斜方肌仍为正。 |
| 主动上回旋 | 在原单轴模型和理想重力预平衡下，钢丝主动拉到 **11.79°**，制动回位至距中立位约 **0.0037°**，往返各检查 49 个姿态。 |
| 装配余量 | 新导向支柱包络在中立位与肩胛的名义间隙为零。上回旋的抽查姿态中无超过 0.001 mm³ 阈值的新相交。 |
| 前伸和后缩 | 后缩抽查未发现新相交。按原前伸控制复查，前锯肌钢丝再次与第六肋相交。 |

这版保留为 **K6 走线候选**。前伸干涉、中立位安装间隙、实际导向眼及摩擦仍需解决，整臂继续保留上一版走线作为对照。[锚点与反力](serratus_guide_points.csv) · [完整复核数据](serratus_guide_refinement.json)

[双视图 MuJoCo 回放](../../../../docs/showcase/media/k6-routing-comparison.mp4) · [峰值截图](../../../../docs/showcase/media/k6-routing-comparison.png)。右侧播放新走线的实际求解轨迹。左侧按同一组姿态显示原走线，便于对照末段与第六肋的位置。彩色叠加线做了视觉加粗。

[同一动作的整臂双视角](../../../../docs/showcase/media/k6-whole-arm.mp4) · [整臂峰值截图](../../../../docs/showcase/media/k6-whole-arm.png)。白色虚线连接中立位的盂肱中心、肘关节中心和前臂末端参考点，白点标出这些点的当前姿态。峰值时盂肱中心上移约 30.9 mm。此轨迹保持盂肱与肘关节相对角度固定。

## 自动绕行对照

把原前锯肌末段当成无摩擦张紧细索，在第六肋局部截面的外侧求最短可见路径，再将结果放回原整臂做三维相交检查。这项试算给出了随姿态变化的绕行点。0–4.0° 的所查点保持直拉。4.5° 起进入绕行分支时，线位仍与第五肋相交，且二维近似后仍留有第六肋三维接触。该局部试算的起止点、绕行侧和障碍范围均已指定。[理想绕行数据](ideal_rib_wrap_probe.json)

## 复现

从项目根目录运行。绿色几何是局部导向的安装包络，用于检查支承位置与空间。进入回放后可拖动旋转。

```powershell
$env:PYTHONPATH = (Resolve-Path .deps).Path
& 'D:/MujocoEnv/Scripts/python.exe' -X utf8 tools/probe_ideal_rib_wrap.py
& 'D:/MujocoEnv/Scripts/python.exe' -X utf8 tools/refine_scapular_routing.py
& 'D:/MujocoEnv/Scripts/python.exe' -X utf8 tools/refine_scapular_routing.py --check
& 'D:/MujocoEnv/Scripts/python.exe' -X utf8 tools/refine_scapular_routing.py --view
python -X utf8 tools/render_scapular_routing.py
python -X utf8 tools/render_scapular_routing.py --whole-arm
```
