# 三角肌携支承导入候选：只改两件金属

在原244对象候选上，为整块三角肌安装增加定向导入让位。骨骼、所有TPU、三出线口、孔位和紧固件均保持原文件及变换；只有前束和后束的金属前支承改变，其余242源行不变。

这是条件安装候选，未替换默认191研究模型，S1未定版，线上GitHub模型不变。原244对象及两组安装失败研究均保留，见[整肌安装研究](../s1_deltoid_installation/README.md)。

## 金属改动

沿所选携带方向对支承的骨侧倒扣进行扫掠让位，保留骨／TPU接触及单材料连通。不削减肌腹或扩大其钢丝通道。

| 支承 | 去除金属 mm³ | 原体积 mm³ | 原／新骨侧接触探针响应 mm³ |
|---|---:|---:|---:|
| 前束 | 0.037876 | 675.570 | 0.334745／0.317809 |
| 后束 | 35.915964 | 624.555 | 3.424996／2.292277 |
| 中束 | 0，保留原源文件 | 651.240 | 不变 |
| 共享远端 | 0，保留原源文件 | 721.002 | 不变 |

接触探针为沿局部法向0.02mm的几何相交量，不是承压面积、载荷或强度比。后束去料和接触变化须在后续强度／预紧审查中保留，不能推断机械性能不变。数值清理增减单列在[生成报告](assembly_check.json)，未自动删除孤岛。

首版对无需改动的中束也重新清理导出，产生不必要的数值变化，已留[冻结记录](iterations/numerical_only_exports/frozen_assembly_check.json)。最终中束／共享远端完全沿用原源，只导出两件确有让位的支承。

## 条件装配走廊

[独立复核](verification_checks.json)确认两改件精确OBJ边界、非流形、重复顶点三角形及零面积为零，各一个正材料分量。新485配对无超0.001mm³相交；其余242源与父版完全一致，父版29,646全零位配对证据保留。因此当前零位证据由父版与本次改动复核共同支持，不称重新执行了整套29,646配对。

三角肌与四个前支承作为五件组合，沿源方向`[-0.9428090416,-0.2357022604,0.2357022604]`在0–80mm间移动。0–3mm步长0.1mm，其余1mm，去回共216个整体姿态，无超过判据的实体相交，骨／TPU接触探针仍有响应。无角度扫描，也不构成所有姿态证明。

这一结果依赖明确的阶段库存：

- 191件已在位。
- 5件作为携带组合装入：整块三角肌及四前支承。
- 48件仍在后装清单，最终设计中全部保留。其中新调整为后装的是冈下肌TPU、抱箍上半／上内衬及八个抱箍五金；其他为三角肌端盖、紧固件和止挡套等原后装件，逐件见报告。

尚未验证：五件组合的临时持件／定位结构、冈下肌后装、抱箍上半及五金后装、其后三角肌全部紧固工序、真实穿线／压接、预紧和承载。**216姿态仅证明这个条件下的组合放入走廊，不能判全肩装配通过。** `stage_pass=false`，无实际驱动范围结论。

## 实际几何查看

以下为MuJoCo按源网格及记录位移渲染，TPU未做形状插值或变形，不是驱动演示。阶段场景显示196件，48件后装件暂不显示；完整244件模型仍保留。

80mm起点：

![起点](entry_demo/entry_80mm.png)

0mm落位：

![落位](entry_demo/available_parts.png)

## 复现与交互查看

在项目根目录执行。前两步使用已有几何环境；姿态文件准备仅需Python标准库，不要求SciPy。

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) '.deps'
& D:/MujocoEnv/Scripts/python.exe tools/build_s1_deltoid_support_draft.py
& D:/MujocoEnv/Scripts/python.exe tools/verify_s1_deltoid_support_draft.py
& C:/Python314/python.exe tools/prepare_s1_deltoid_entry_demo.py
```

查看完整244件设计：

```powershell
& C:/Python314/python.exe mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_deltoid_support_draft/assembly_check.json'
```

查看可旋转、缩放的80mm装配阶段场景：

```powershell
& C:/Python314/python.exe mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_deltoid_support_draft/entry_demo/pose_80mm.json'
```

将报告名改为`pose_0mm.json`可看落位阶段。两姿态均为0DOF静态场景，几何放入路径证据来自上述216姿态复核。
