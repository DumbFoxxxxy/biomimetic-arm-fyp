# 完整肌腹分析边界推进（2026-11-03）

**本轮没有更换S1-Y2，没有合格完整TPU体网格或驱动结果。** 新增显式索引OBJ导入、四个完整表面修整候选，以及一个冈下肌孔道重建候选。所有原模型、锚点、路径及失败保留。

## 新增证据与取舍

| 路线 | 实际结果 | 决定 |
|---|---|---|
| 高精度OBJ输入原全局修复版 | 两块仍未建网；不是只由旧分桶导入造成 | 未采用 |
| 极短边收缩，原顶点簇位移≤0.001 mm | 冈下肌734次、肩胛下肌293次收缩，闭合索引拓扑保留；仍建网失败 | 未采用 |
| 退化面邻域局部修整，位移≤0.02 mm | 退化面减少，但仍有PLC交点；部分三角面法向翻转，保护材料损失另列 | 未采用 |
| 对修整版执行Gmsh t13表面参数化 | 边界拓扑不适于参数化，失败 | 未采用 |
| 冈下肌显式重建较大孔道 | 所用浮点表面扫描不再检出交点；0.75 mm壁厚包络及静态刚体／主动线复核通过 | 仅研究分支，未采用 |
| 新孔道Delaunay完整体网格 | 生成199103个四面体，8个非正体积、73357个平均比质量低于0.05，体积不守恒 | 拒绝用于连续体求解 |
| 新孔道HXT／表面参数化 | 前者原生堆损坏退出；后者重复递归后中断，无最终报告 | 不算完成或通过 |

前六个完整体网格失败均以原生退出码2及日志留档，见[captured_attempts.json](captured_attempts.json)。四个表面候选的独立重放、材料差分、保护区及静态检查见[verification_checks.json](verification_checks.json)。这些检查不是实物验证。

## 为什么重建孔道

[位置分类](crossing_region_classification.json)针对之前报告中已有的浮点检测点：冈下肌107处均距内部路径约0.85065–0.86000 mm，集中在原孔壁；肩胛下肌67处中29处位于保留端部／通道区，另外38处在该区域外。因此两块肌腹不能假定共用同一个修复结论。数量是对应检测方法的记录，不是精确全局自交总数，也不证明唯一成因。

冈下肌[孔道重建候选](../../../02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_tpu_channel_rebuild_revision/infraspinatus/radius_0.9/README.md)明确改变孔径，不以原0.001 mm³判据掩盖材料切除。中心线、锚点与其余190个源不变，孔壁采用正搭接的圆柱／球节点。它没有使制造网格或驱动能力获得通过。

## 工具与复现

`build_shoulder_tpu_volume_mesh.py`现在读取`import_convention=indexed_float64_obj_no_bucket_weld`，对该OBJ保留原索引及双精度坐标，不分桶、不删面。原STL分桶约定保持不变，两块STL的顶点／面数组逐项回归一致；旧工具源在各次尝试目录中归档。

在项目根目录使用已有numpy／scipy／manifold3d／Gmsh环境，令`.deps`可导入。例如PowerShell中先设置`$env:PYTHONPATH = '.deps'`，然后使用该环境的Python：

```powershell
python mujoco_arm/tools/build_shoulder_edge_collapse_surface.py --module infraspinatus
python mujoco_arm/tools/build_shoulder_local_fairing_surface.py --module infraspinatus
python mujoco_arm/tools/build_shoulder_rebuilt_channel_surface.py --module infraspinatus --radius 0.9
```

生成器会写候选目录；需要保留的结果先归档，不能覆盖实测。原生建网统一使用有资源上限的入口，输出目录需新建且不得已有结果：

```powershell
python mujoco_arm/tools/run_shoulder_mesh_worker.py --entry volume --output tmp/shoulder_mesh_run_01 --timeout-seconds 120 --maximum-log-bytes 8388608 -- --module infraspinatus --mode preserve --algorithm 1 --facet-angle-tol 0.000001 --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_tpu_channel_rebuild_revision/infraspinatus/radius_0.9/candidate_assembly.json'
```

资源上限只控制运行和日志量，不改变几何／数值质量门。超限、崩溃、中断或只有初始报告的进程不算完成；监督器保存版本、日志、退出码并回收子进程。已用独立目录的1字节日志上限控制例验证停止与回收，没有把控制例当作肌腹通过。

独立审计命令：

```powershell
python reviews/S1_2026-11-03/indexed_boundary_meshing/verify_candidates.py
python reviews/S1_2026-11-03/indexed_boundary_meshing/verify_rebuilt_channel.py
```

`capture_attempts.py`用于本轮六次原生输出捕获，会重写其尝试目录且拒绝已有捕获备份；它是审查记录，不是日常重复运行入口。新运行应采用上面的监督器和独立输出目录。

## 接下来的设计动作

不继续无界增大补球或放宽质量门。以新孔道为独立对照，从参数化的通道／端部过渡构造与边界三角化重新处理退化薄片，保留真实晶格、端垫、挂钩、锚点和钢丝路径；每个实质修改重新核对接口、壁厚、材料差分与完整体网格。肩胛下肌还需处理保留区外的交点。

只有最终回读表面、数值质量与几何偏差均验收后，才考虑候选替换；之后仍需材料、工艺、预紧、摩擦和单肌肉实测。完整技术文件见[肩部方案](../../../02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)。
