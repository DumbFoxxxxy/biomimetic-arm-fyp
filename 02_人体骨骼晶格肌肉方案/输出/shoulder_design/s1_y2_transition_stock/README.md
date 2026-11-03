# S1-Y2：保留外形的晶格肌腹过渡区修订

更新：2026-11-02。选为当前几何研究起点，**不是制造定版或实际驱动通过**。唯一选择及报告哈希由[候选指针](../../../../design/S1_CURRENT_CANDIDATE.json)维护。原[S1-Y1](../s1_y_muscle_pair/README.md)、失败搜索及比较候选完整保留。

## 为什么修订

Y1两块肌腹各为单正材料连通，但内部仍有4个封闭空腔。空腔与此前紧凑骨架父级的肌腹空腔体积高度重合，位于内部路径R3.3 mm的保留过渡区，不与名义完整钢丝通道或端接口相交。单材料连通不代表没有内部闭腔。

原开口搜索将刀具向空腔背面伸入一个孔半径，可能越过薄空腔对侧壁。用实际空腔厚度限制后伸后，肩胛下肌0.8／1.2 mm开口能开敞两处，仍保留约0.08114 mm³闭腔；1.6 mm还保留约0.17525 mm³闭腔。原失败仅否定当时的有限搜索，不能据此断言肩胛下肌不能开口。[修订搜索记录](../s1_y_muscle_pair/vent_experiment/outward_start_revision/outward_trials.json)与[独立复核](../s1_y_muscle_pair/transition_stock_comparison/independent_checks.json)保留该区别。

为保留现有外部杆网、通道壁和端部，比较了“封闭过渡空隙填实”和“0.8 mm开口＋最小剩余空隙填实”。本阶段选择前者：没有外部切口，不依赖未确定打印工艺对细孔的还原。该选择只针对几何试制准备，不宣称其力学表现更好。

## 明确的材料变化

每个待填空隙先验证：它是闭合负体积空腔壳、全部位于R3.3 mm保留过渡区、与1.72 mm完整路由孔腔及两端保护区的体积相交不超过0.000001 mm³。仅对通过分类的4处做显式TPU填实；没有把所有负壳一律填掉，也没有交给导出清理或切片器自动修补。

| 模块 | 明确填实的原空隙 | 材料增加 | 相对原肌腹体积 |
|---|---|---:|---:|
| 冈下肌FJ1500 | 1处 | 约13.11811 mm³ | 约0.156% |
| 肩胛下肌FJ1504 | 3处 | 约2.79255 mm³ | 约0.036% |

有符号布尔差分中的约−0.000027 mm³等微小负值作为数值残差保留，不解释为负材料。[生成与逐空腔参数](../s1_y_muscle_pair/transition_stock_comparison/transition_trials.json)记录位置、体积与分类依据。

局部填实会改变该处刚度与折叠方式，影响大小未测。钢丝名义直径0.6 mm、孔径1.72 mm以及原两端接口不变；不能仅凭材料增加比例很小就认定肌肉性能不变。

新增[自由末端／肌腱约束及4端环局部运动](free_terminal_kinematics/README.md)，明确两压接出口间的中央段长度与张力单边约束；理想张紧时仍余7个端部变量待求。局部端环扫描、零位直线探针均不作肩部驱动范围；本轮不改本装配源或指针。

另补[规定端环角的整根肌腱候选](free_strand_routes/README.md)，中央段与两端硬件均在场，含MuJoCo细节查看。此项保留端盖／肱骨零位和未测的名义长度，不作张紧平衡或GH范围。本装配源及选择指针保持不变。

新增[完整TPU体网格前置检查](volume_mesh/README.md)：两肌腹均未通过，已定位远端过渡局部边—面交点；此前回读采用四位小数分桶导入。当前没有完整连续体变形解，源几何未修改，立方体对照不提升肌腹验收。

## 当前数字证据

- [191对象装配回读](assembly_check.json)：只换两件TPU，其余189源及变换不变；所有源单正材料连通、零边界／非流形边。这两块修订肌腹无封闭负壳，其他历史对象并未因此获得无闭腔结论。
- 重新检查全部18145对零位材料及六路完整主动线，无超0.001 mm³判据刚性／主动线相交；既有柔性重叠仍登记。
- [独立单肌肉复核](../s1_y_muscle_pair/transition_stock_comparison/independent_checks.json)检查177刚性对象、六主动线、原R1.8内部通道／端部保护材料、近端0.05 mm接触及远端−0.2／+0.25 mm有限探针。各轴向分层仍有实际材料支持的内部单元与表面杆。
- [裸骨路径](bare_motion_review/bare_joint_motion_checks.json)针对新清单重新执行：正屈曲／外展到120°检查上限；负向−39.6°台面首限／−39.1°夹爪螺钉首限，长轴±120°至上限但未映射为解剖内外旋。常规1°、首界0.1°、独立回程；102对象在场，89未装，部分最小间隙为0。
- [保持与名义放置](readback_keeper_checks.json)重新核对11步4322姿态、裸肱骨0↔10°、防转／侧退与52方向球头有限捕获。真实工具／紧固／承载和完整肌肉安装仍未验收。

零位、裸骨和有限探针都不能代替柔性变形、摩擦、张力与自由末端—肌腱平衡。原K1–K7不降低，实际驱动范围未知，`stage_pass=false`。

## 查看、源文件与复现

默认查看仍用同样命令，读取候选指针核对哈希：

```powershell
python mujoco_arm/tools/view_shoulder_manufactured.py
python mujoco_arm/tools/view_shoulder_bare_GH_motion.py
```

[实际网格静态图](available_parts.png)与[裸骨三姿态图](bare_motion_review/bare_GH_motion.png)不展示TPU变形。两块新STL路径在装配报告的`sources`中；坐标规则与[Y1说明](../s1_y_muscle_pair/README.md)相同，局部STL加其`translation_source_mm`恢复源解剖坐标。

[钢丝／锚点清单](wire_and_anchor_plan.json)只更新模型引用，钢丝和锚点坐标不变。[空白实测记录](single_muscle_observation_template.json)绑定新网格哈希，所有工艺与测量值保持未知；原Y1记录不覆盖或冒充新版本测量。

重算依赖numpy、manifold3d及项目`.deps`，会重写各独立候选目录；已测记录不可随重算覆盖：

```powershell
python mujoco_arm/tools/build_shoulder_Y_outward_vents.py
python mujoco_arm/tools/build_shoulder_Y_transition_stock.py
python mujoco_arm/tools/check_shoulder_Y_transition_stock.py
python mujoco_arm/tools/build_shoulder_Y2_transition_assembly.py
python mujoco_arm/tools/check_shoulder_bare_GH_motion.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/assembly_check.json' --out '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/bare_motion_review'
python mujoco_arm/tools/check_shoulder_GH_nut_keepers.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/assembly_check.json'
```

全部检查后再审阅并同步候选指针、钢丝引用、未测模板及渲染；不能只改哈希绕过检查。Y1晶格生成是本版父级，不应把当前Y2重新输入Y1生成器叠加处理。下一步为实际切片／工艺、穿线与完整肌肉试验，并建立自由末端耦合模型；封闭空隙处理不替代这些工作。
