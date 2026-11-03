# GH 独立保持关节对照候选

2026-11-02。目标已恢复为 active，继续按照肩部技术审查稿推进。此处是独立比较候选，尚未替代开放浅窝方案，实际 TPU 驱动范围未知。

## 本轮完成

先核对旧 G1 工程关节与紧凑肌肉候选，发现旧工程肩胛与三角肌重叠 311.909 mm³。其中 GH 新增部分仅贡献约 1.840 mm³，其余主要来自旧 AC 连接上下文。直接套用 G1 会把 AC 冲突误归给 GH。

本候选按 `原解剖肩胛 − (AC 基线 − G1 肩胛) ＋ (G1 肩胛 − AC 基线)` 应用 GH 改动，保持当前原锁骨与夹具，继续使用工程球头、保持盖和名义五金。注意差分也包含 G1 当时的局部上下文让位，不是从头设计的纯 GH 插件；所有改动体积登记在报告中。

两枚原螺栓 OBJ/STL 各有一条非流形边。新螺栓以单一旋转截面重建同名义 M3×12 包络，未建螺纹。肩胛差分导出最初有 81 条非流形边，普通简化后 58 条，重新建立单一几何来源并按 0.01 mm 简化后剩两条四面共边；在两段各不足 4 mm 的接缝上增加 R0.1 局部圆整，最终回读边界／非流形均为 0、单正实体。失败几何保存在 `iterations/unrounded_seam/`，没有通过删除连接部位修复。

更圆的名义螺栓头与支座曾重叠 0.174893 mm³；增加名义 0.15 mm 局部让位，移除约 4.43474 mm³ 支座材料，重新导出后相交为 0。该余量是候选参数，尚未验证打印公差与支座强度。

## 零位与刚性路径结果

| 检查 | 结果 | 适用范围 |
|---|---|---|
| 可用对象 | 34 件，含 3 骨、4 TPU、10 肌肉盖、5 止挡、2 夹具组、10 GH 附件 | 缺肌肉骨侧紧固及完整夹持接触结构 |
| 六路钢丝 | 对 34 件均未检出超过 0.001 mm³ 的相交 | 零位名义路径，钢丝半径 0.3 mm |
| 零位最大刚性相交 | 约 0.00003693 mm³ | 不是正制造间隙证明 |
| 三角肌／新肩胛接触 | 约 1.84044 mm³ | 需柔性接触分析，未判通过或失败 |
| 有限脱出探针 | 26／26 方向在 0.9 mm 内遇到材料 | 0.1 mm 步长、最多 3 mm；仅零位直线平移 |

固定肩胛／锁骨，强制 GH 中心，肱骨侧端盖及止挡理想随骨安装；不把 TPU 当刚体扫掠。1° 常规采样，边界 0.1°，相交阈值 0.001 mm³，间隙查询封顶 2 mm。以下每个末样本最小间隙均记录为 0 mm，不代表有制造余量。

| 探针 | 末个未超阈值样本 | 首个超阈值样本 | 首限 |
|---|---:|---:|---|
| 屈曲＋ | 114.2° | 114.3° | 肩胛下肌远端盖／肩胛 |
| 屈曲－ | −64.2° | −64.3° | 工程肱骨／GH 保持盖 1 |
| 外展＋ | 61.1° | 61.2° | 冈下肌远端盖／肩胛 |
| 外展－ | −42.6° | −42.7° | 工程肱骨／肩胛 |
| 长轴＋ | 68.1° | 68.2° | 肩胛下肌远端盖／GH 螺栓 0 |
| 长轴－ | −47.0° | −47.1° | 肩胛下肌远端盖／肩胛 |

相比原紧凑候选，正屈曲与正轴旋下降，负屈曲、负外展延长，正外展仍受相同端盖限制。不能称为运动覆盖全面改善，也不能把长轴正负直接标为内外旋。

新肩胛相对输入移除约 675.989 mm³、增加约 6855.237 mm³；肱骨仍沿用 G1 工程头，移除约 34892.118 mm³。因此不是未改动的人体关节，也尚未解决解剖附着与工程接口之间的过渡。

当前骨侧名义点的表面距离见父目录 `compatibility.json`；点距只描述该点，不代表整个端垫接触面积，也不能证明已有骨侧安装。完整肌肉固定、三轴耦合、任意姿态防脱、强度与真实驱动都未验证。报告 `stage_pass=false`。

## 复现与查看

在项目根目录、安装 RUNNING.md 所需依赖后运行；本地 manifold3d 位于 `.deps` 时先设置本次 shell 的 Python 搜索路径：

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) '.deps')
python mujoco_arm/tools/audit_shoulder_retained_joint_fit.py
python mujoco_arm/tools/build_shoulder_isolated_retained_joint.py
python mujoco_arm/tools/check_shoulder_rigid_motion.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/GH_only_candidate/assembly_check.json'
python mujoco_arm/tools/check_shoulder_retained_capture.py
python mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/GH_only_candidate/assembly_check.json'
```

查看器是零位静态网格，DOF=0，不是驱动或柔性仿真。添加 `--render` 输出图像。

- [零位装配检查](assembly_check.json)
- [刚性运动检查](rigid_motion_check.json)
- [有限脱出探针](finite_retention_probes.json)
- [静态 MuJoCo 图](available_parts.png)
- [当前肩部技术方案](../../SHOULDER_DESIGN.md)

下一步优先比较肱骨侧端盖／短肌腱接口，重新验证力臂和钢丝传力；并明确晶格工作方向及单肌肉试验。此候选保持独立，不覆盖原骨骼、既有肘部或紧凑端盖候选。
