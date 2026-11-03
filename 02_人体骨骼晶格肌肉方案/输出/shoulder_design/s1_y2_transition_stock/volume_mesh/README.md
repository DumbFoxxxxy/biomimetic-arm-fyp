# 两块完整TPU肌腹：体积网格前置检查未通过

更新：2026-11-03。当前骨架、肌腹STL和候选指针保持S1-Y2原字节。本轮建立Gmsh体网格输入、日志、数值质量门槛和失败定位流程；**两块完整肌腹均未获得可用四面体网格，尚不能进行可信的TPU连续体变形求解**。

## 结果与源几何约定

[六组尝试与诊断](diagnostics.json)包含直接保留表面、0.005 mm分析副本简化、标准Delaunay及HXT边界恢复试验。全部完整肌腹尝试失败。数值近重合角参数的比较另行记录，没有降低正体积、边界材料体积一致性或单元质量接受门槛。

| 原始STL／导入 | 冈下肌 | 肩胛下肌 |
|---|---:|---:|
| 原始三角面 | 81668 | 85282 |
| 精确坐标去重后的非流形边 | 0 | 6 |
| 现行四位小数分桶导入后的非流形边 | 0 | 0 |
| 导入合并的最大坐标重指派距离 | 0 | 0.00006175 mm |
| 导入删除的退化三角面 | 0 | 8 |

此前通用网格回读使用`mesh_geom.weld(..., decimals=4)`，以四位小数的坐标作合并键，保留桶中首个真实坐标，不把全部顶点吸附到规则网格。此前“零边界／非流形边”结论限定于该导入约定；不能读成原始STL精确坐标拓扑、PLC边界适用性或制造通过。本轮把精确去重和分桶导入结果分别保存，没有修改源STL。

单正材料连通及边计数正常仍不足以证明表面没有非相邻边—面的交点。Gmsh在原表面报告近重合面；采用当前Manifold几何再做独立HXT恢复，仍报告线段与三角面相交，定位在两块肌腹远端接口附近。

## 当前定位结果

| 肌腹 | 报错点源坐标X / Y / Z，mm | 远端接口局部X / Y / Z，mm | 本次局部交点 |
|---|---|---|---:|
| infraspinatus | -179.7301 / -59.2198 / 1321.1309 | -0.1296 / -3.9105 / -2.1090 | 12 |
| subscapularis | -146.6109 / -90.7783 / 1309.6629 | -0.0412 / 4.7787 / -0.9321 | 2 |

报错点来自原生网格器输出的有限精度坐标。附近0.2 mm定位盒选择相交的原三角面，长三角面可能延伸到盒外。[开放诊断STL与交点](diagnostics.json)保留实际边／面坐标。上述交点数只对应此次定位区域；并行边界恢复可以先报不同位置，它不是全肌腹自相交总数。

这些报错点不在S1-Y2四处新增填实空腔的包围盒内。位置更接近既有远端夹持／挂钩过渡，需针对交汇表面修复，不应无依据撤销闭腔处理或移动钢丝锚点。

![局部表面诊断](localized_boundary_failures.png)

投影图展示源三角面及交点，不是应力或位移图。`failure_patch_OPEN.stl`是开放局部面片，只用于定位，不能当完整肌肉或打印件。

## 工具与独立验证

Gmsh 4.15.2已安装到项目忽略的`.deps/gmsh_vendor`，没有打包到仓库。Windows目标目录安装遗漏DLL时，本机从已获批准下载的原wheel缓存恢复，并核对安装RECORD中的SHA256；未修改系统CAD或模型目录。新环境可在自己的几何Python环境安装`gmsh==4.15.2`。

表面导入／体网格构建参考[Gmsh官方离散实体与STL教程](https://gmsh.info/doc/texinfo/gmsh.html#t13)。Python直接读取真实中文路径的STL，再通过API传节点／三角面，不创建兼容目录。

立方体对照得到12个正体积四面体、总体积1000 mm³，独立核对12个边界三角面。这个对照只验证网格API与数据提取路径。[独立验证](verification_checks.json)另用线性方程最小二乘复核局部边—面交点；它不提升完整肌腹失败结果。

当前接受门槛：四面体正体积；相对分析边界的材料体积差小于1e−5；mean-ratio最小质量不低于0.05。它们是本项目数值预检，不代替网格收敛、材料标定或实体验证。完整肌腹尚无四面体结果，门槛保持失败。构建脚本失败退出码为2；诊断脚本可以成功完成失败定位，但不意味着肌腹通过。

[连续体模型待定输入](elasticity_input_status.json)中的材料、打印、边界预紧及接触参数全部未定，未移植论文弹性模量或任意收缩比例。

## 复现

使用具备numpy、scipy、manifold3d、gmsh的几何环境，在项目根目录执行：

```powershell
python mujoco_arm/tools/build_shoulder_tpu_volume_mesh.py --module infraspinatus
python mujoco_arm/tools/build_shoulder_tpu_volume_mesh.py --module subscapularis --mode validated_manifold --algorithm 10 --facet-angle-tol 0.000001
```

上述当前预期失败，详情写入对应`mesh_attempt.json`。以下命令重新捕获两次原生HXT错误、做对照与局部定位：

```powershell
python mujoco_arm/tools/diagnose_shoulder_tpu_volume_mesh.py --replay-hxt
python mujoco_arm/tools/verify_shoulder_volume_mesh_diagnostics.py
python mujoco_arm/tools/render_shoulder_volume_mesh_diagnostics.py
```

各失败目录保存原日志、源哈希、输入边界、生成器快照或其哈希名版本。已有运行结果有历史备份；复算会更新当前诊断，实测记录不在这些命令中修改。

下一步先修复已定位的远端过渡表面，重新验证接口、通道与完整体网格。通过后才进行连续体离散和自由端盖／肌腱接触求解；实际材料和安装数据仍需实测。当前`stage_pass=false`，实际驱动范围未知。
