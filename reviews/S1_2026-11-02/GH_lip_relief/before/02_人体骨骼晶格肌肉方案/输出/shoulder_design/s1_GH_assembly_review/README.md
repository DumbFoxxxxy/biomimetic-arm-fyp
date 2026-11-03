# GH装入与捕获审查

更新：2026-11-02。此目录保留问题证据与独立试改，不是已通过的关节定版。

## 结论与文件

| 对象 | 文件 | 结论 |
|---|---|---|
| 187对象前衬垫保持版 | [父级源清单](../s1_front_liner_keepers/assembly_check.json) | 当前集成几何起点；旧GH盖存在轴向脱出反例 |
| 旧盖轴向审查 | [完整审查报告](parent_187/axial_capture_audit.json) | 开盖肱骨、闭盖肱骨及闭盖诊断球沿真实颈轴25 mm受检路径均清 |
| 深度6 mm新盖试改 | [源清单与变化配对检查](retained_depth_6/assembly_check.json) | 替换两盖，保留其余185源；新增与肩胛下肌约6.65813 mm³相交，轴向捕获未复验 |

旧盖检查为裸骨／夹具阶段，固定对象和未装对象逐例列在报告中。不是全部187对象同时参加移动检查；肌肉、锁骨骨段及其他未装模块不能隐含为在场。

## 参数与判据

- 采用原GH构造记录的源坐标颈轴 `[-0.29658719, -0.02805316, -0.95459366]`，不是世界网格近似方向。
- 肱骨使用实际工程网格；独立R16 mm球仅用于诊断，两者分别记录。
- 0–25 mm、步长0.1 mm，每例251姿态；相交体积大于0.001 mm³才报告为受阻。
- 旧盖深3 mm，窝内半径16.4 mm，理想未开槽端口半径约16.1233 mm，大于16 mm球头。该公式辅助解释，不能代替实际网格检查。
- 新盖将球面保持部分延至6 mm，保留原3 mm紧固承靠面，并为实际五金让位；理想未开槽端口半径约15.2630 mm。槽口、五金让位及盖的装入路径仍须验证。

旧26方向有限捕获探针未包含真实颈轴，其结果不能用于否定本反例。新盖与柔性肌肉相交不直接证明肌肉最终不可用，也不能未经变形或实物检查就放行。报告保持 `stage_pass=false`、`actual_drive_range_deg=null`。

## 查看与复算

在项目根目录，使用装有MuJoCo／Pillow的环境静态查看试改网格：

```powershell
python mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_GH_assembly_review/retained_depth_6/assembly_check.json' --check
```

这是加载检查，不进行捕获或TPU变形验收。几何复算需要numpy、scipy和manifold3d；使用项目几何环境，并使根目录`.deps`中的依赖可导入。下列命令写入独立复算目录，保留本次证据：

```powershell
python mujoco_arm/tools/audit_shoulder_GH_assembly.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_front_liner_keepers/assembly_check.json' --out 'tmp/GH_parent_recheck'
python mujoco_arm/tools/audit_shoulder_GH_assembly.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_GH_assembly_review/retained_depth_6/assembly_check.json' --out 'tmp/GH_depth6_recheck'
```

第二条是待执行验证的复现入口，**本次文件整理没有运行，也没有其通过结果**。`build_shoulder_GH_retaining_caps.py`会覆盖此目录下的试改输出；如需重生成，先保存现有文件及哈希。一般阅读无需运行生成器。

## 修订完成条件

1. 协调保持盖与肩胛下肌空间，记录修改对象、局部偏差和实际接触依据。
2. 用重新导出的实际网格分别检查开盖装入、闭盖受阻；扩大逃逸方向与联合运动检查，保留采样边界，不宣称任意姿态证明。
3. 检查盖、螺钉、垫圈、螺母及工具的实际分步装入；重查变化配对与六主动线。
4. 另行验证打印配合、防脱／承载、磨损及实际肌肉驱动。没有实测时保持性能未知。

主技术说明见[SHOULDER_DESIGN第2.10节](../SHOULDER_DESIGN.md)。
