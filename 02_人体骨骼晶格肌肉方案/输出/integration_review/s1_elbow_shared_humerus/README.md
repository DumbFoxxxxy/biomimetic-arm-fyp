# 肩肘共用肱骨候选

后续设计以原版 FJ1486 屈肌和 FJ1477 伸肌为基准。保留其异形晶格、附着点和走线，优先调整共用肱骨上的刚性安装接口。[原版肘部基准及零位冲突](anatomical_elbow_baseline.json)

原版伸肌与当前肩侧三角肌在零位存在约 445.42 mm³ 的几何重叠。该项需要与刚性安装座干涉分开解决。下述转位组合和对应网页展示属于机械集成试验。

当前机械候选将 20 件远端转子、销轴和防脱件恢复到前臂侧，只移动近端安装件。伸肌远端 TPU 过渡段缩至 1 mm，近端采用侧装导向件和后装承托件。[修正版装配与运动记录](extensor_late_reaction_audit.json)

修正版报告记录了按指定关节路径的刚性避让和名义钢丝路径。TPU 变形、带载传力和实际驱动范围仍待验证。下文的零位总成和早期试算保留为设计对照。

## 肘部TPU回读修订

两块肘部TPU的源面索引闭合，合并重合顶点后会产生非流形边。采用半径0.005mm的局部补材，屈肌与伸肌分别增加约0.00001334和0.00000453mm³。补材对各自名义钢丝孔腔的检测相交为零，两件均通过精确OBJ回读。[修复记录](elbow_tpu_readback_repair.json)

将两修复件放回325对象清单，补查647个涉及它们的零位配对，检出三角肌与肘伸肌约445.42mm³相交。结合下文原有12处相交，目前共检出13处跨模块相交。该结果用于碰撞设计，实物接触与变形仍待验证。[补查报告](repaired_elbow_tpu_neutral_audit.json)

两处主要转接座干涉位于三角肌端部区之外，分别与中束和前束钢丝保护区重叠约15.13和24.82mm³。对两座直接试切后，原单实体各分为三个实体。应调整肘侧近端接口位置并重新核对肌腹和走线。[材料映射](interference_material_map.json) · [转接座试切](elbow_adapter_relief_screen.json)

MuJoCo查看时可改用 `assembly_tpu_readback_candidate.json`，该清单已通过325对象加载检查。

这一版共有325个几何对象，肩和肘使用同一根肱骨。肩侧把原抱箍的肌腱座、金属连接形状和过渡材料并入工程肱骨，取代原12件抱箍、内衬和紧固件。零件数包括仍待修复的两块肘部TPU。

## 共用骨接口

肩侧一体化候选在原肱骨上增加约4658.18mm³。精确OBJ回读为单实体，232对象零位中新增骨体未检出超0.001mm³相交。新增材料与六条钢丝走廊无实体相交，三角肌携四支承80mm路径369个姿态未检出超判据相交。承载、材料和完整安装工序待验。详见[肩侧一体化检查](../../shoulder_design/s1_integrated_humerus_anchors/assembly_check.json)。

共用肱骨在局部轴向110mm站位接合。肩侧保留近端，肘侧保留远端，匹配轴段的重叠宽度为2mm。两模块的紧固特征有部分处于同一轴段，单独接合后仍有八处肘侧骨接口干涉。局部转移原肘骨孔腔去除1185.85mm³。精确回读后，为FJ1486一颗螺栓做0.02mm轴向邻域让位，再去除约0.056mm³。共用骨与322个可读回的零件检查通过0.001mm³判据。详见[共用骨审核](shared_humerus_audit.json)。

## 整体零位检查

检查了21,021个肩肘跨模块配对，共检出12处超判据相交。最大的三处均涉及三角肌，分别是FJ1477骨侧转接座约241.68mm³、FJ1486骨侧转接座约168.20mm³、FJ1477近端刚性端盖约92.39mm³。两块肘部TPU因非流形边未进入实体相交计算。详见[跨模块审核](inter_module_neutral_audit.json)。

当前为本地零位候选。优先处理肘侧转接座、近端端盖与三角肌的占位，再修复两块肘部TPU网格并复核全装配及运动。承载与实际驱动仍待验证，发布条件尚未满足。

![肩肘共用肱骨零位](available_parts.png)

在项目根目录查看可旋转的MuJoCo模型。

```powershell
& C:/Python314/python.exe mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus/assembly_check.json'
```

生成步骤依次为 `build_s1_integrated_humerus_anchors.py`、`verify_s1_integrated_humerus_anchors.py`、`study_s1_shared_humerus_seam.py`、`build_s1_shared_humerus.py` 和 `audit_s1_shared_humerus_inter_modules.py`，均在项目 `tools` 目录。
