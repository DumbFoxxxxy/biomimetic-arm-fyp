# 三前衬垫定位与M2捕获候选

2026-11-02。以175对象接收部修订版为父级，改3片前衬垫及2支架，新增12件名义五金，形成187对象独立候选。骨骼、肌腹、钢丝及父级其余170个源文件原样保留。**三前衬垫已有名义定位／止退几何，不再只是在场景中假定预置；实际保持力、打印配合、拧紧和位置精度尚未验证。**

## 1. 每片的结构

衬垫上增加主定位柱（带侧装M2螺母腔）与小定位销。主柱及小销从骨侧穿入接收部，螺母先侧装入柱，再由背侧M2固定。螺母腔处于接收部背部3 mm区域；装入后孔壁包围侧入口，不能把螺母在完全安装后沿原侧口放入作为装配顺序。

| 特征／部件 | 每片数量 | 名义尺寸 |
|---|---:|---|
| TPU主柱 | 1，与衬垫一体 | Ø6.6，接收孔Ø6.8；上端高于接收部背面1 mm |
| TPU小定位销 | 1，与衬垫一体 | Ø2，孔Ø2.2；局部U/V=4.8/4.8 mm |
| 螺母及侧腔 | 1 | M2无螺纹包络，对边4、高1.6；径向0.1、轴向每侧0.05 mm模型间隙 |
| 螺母腔顶TPU屋面 | — | 名义1.65 mm；抗拔／蠕变未知 |
| 背侧隔套 | 1 | OD9／ID6.8×高2.8 mm，给凸出的主柱留空间 |
| 垫圈 | 1 | OD9／ID2.2×厚0.5 mm |
| 螺钉 | 1 | M2×6，无螺纹包络；模型末端超螺母约0.4 mm，不入骨 |

每片新增螺母、螺钉、隔套、垫圈4件，三区12件。不是12个TPU打印件。材料、商品规格、扭矩、螺纹啮合和切片工艺均未确定；Ø6.6柱和1.65屋面也不等于已通过打印或强度验证。

柱根根据原衬垫实际射线交点定位，插入原衬垫约0.4 mm形成实体连接。只修去新增部分可能进入骨面的材料，不移动原骨／肌。新增接收孔、0.01 mm局部法向模型间隙及20 mm衬垫入口连续开敞均记录于[装配报告](assembly_check.json)。模型间隙不是选定的制造公差。

## 2. 为什么需要单独检查衬垫安装

骨可以沿前上方路径装入，不代表衬垫也能沿相同方向随骨进去。三片原衬垫各沿自身法向的1 mm粗采样看似能拆出，但加保持柱后的0.1 mm检查暴露上角／外侧起步干涉约0.03907／0.01331 mm³。

本轮按衬垫自身局部−n、20 mm路径连续开敞接收部，消除该倒扣，真实止退改由定位柱、捕获螺母与M2承担意图，而不借助无法装入的刚性卡接。失败版本及独立骨路径保留于[细查失败记录](iterations/liner_placement_leadin_blocked/readback_retention_checks.json)。初版0.001072 mm³共面读回残差另保存在iterations/initial_contact_residual，不靠提高判据通过。

相对父175对象，三个衬垫增加约191.233／249.888／158.166 mm³；旧衬垫读回移除差值各小于0.0001 mm³。两支架移除约495.801／186.823 mm³，另有约0.001459／0.001222 mm³读回新增差值，不写绝对几何不变或严格零相交。

## 3. 规定的分步安装与证据

[独立回读与保持检查](readback_retention_checks.json)重新加载187对象，检查17391对零位材料，无超0.001 mm³刚性相交；六路内部主动线对全体无超阈值相交。所有新增／替换制造网格均单正材料实体、零边界／非流形边。

按当前实网格重查90对父级接触，另加12对隔套／垫圈／头部／螺母承靠面，共102对有限法向探针均遇材料。螺母到TPU屋面有0.05 mm名义轴向间隙，不能把该检查写成已经预紧。

在无骨、无GH／肌腹、其余定位结构已就位的明确阶段，逐项在0.1 mm下独立计算去／回程，共3396姿态：

1. 衬垫位于接收部前方20 mm时，螺母从侧槽插入／取出6 mm；此时需人工暂持螺母，未验证手部操作。
2. 衬垫与暂持螺母沿该片局部法向装入，反向抽出20 mm；不假定它们粘接。
3. 隔套轴向放入，随后垫圈，再将M2螺钉从背侧放入。受检抽出长度分别12／12／6.1 mm；螺纹解除后的轴向包络清，不等于螺钉已能真实拧紧。
4. 给三个M2轴上方各检查Ø3×长20 mm的候选直线工具入口，阶段内无超阈值相交。该包络不包含真实刀头／内六角窝、完整手／扳手或选定工具。

全部上述受检几何路径清；完整GH、后衬垫、浮动爪、五金、肌肉及整机顺序尚未验收。

## 4. 名义止退与局限

将TPU衬垫作为刚性几何探针，在M2螺母／螺钉名义位置不变时沿−n试拉，三片均在0.1 mm遇到螺母屋面接触；去掉螺钉、沿侧槽方向试移螺母，安装后的框架在约1.1–1.2 mm遇到材料。各方向的有限探针见检查报告。

这些只证明所选方向存在几何捕获，不证明任意脱出路径、真实TPU拉出力、保持位置精度、疲劳或螺纹握持。径向0.1 mm、螺母0.05 mm间隙也不是实际装配游隙；不以固定金属组的位置假设替代真实力学验证。

## 5. 保持件参与的肩胛装入复查

[当前装入细查](insertion_review/selected_path_dense_check.json)保留全部12保持五金：步骤场景57对象（含肩胛），187源对象中130个尚未装入，包括GH／肱骨／锁骨骨段／肌腹／后爪等。肩胛前上方20 mm去／回各201姿态，整段0.1 mm，均未超原0.001 mm³判据；最大残差约0.00014782 mm³、最小间隙含0。不是严格零穿透或制造间隙验收。

[MuJoCo单步骤图](insertion_review/scapula_insertion_demo.png)和交互查看包含57对象，骨由mocap平移，nq=0／nmocap=1，无动力学步进。保持件是名义装配，实际预紧与TPU变形未建模，不作为肌肉驱动演示。

[完整187对象静态图](available_parts.png)、[仅衬垫打印几何细节](part_detail.png)分别显示实际网格；后者隐藏所有金属／框架，不能理解成装配已完成。

## 6. 复现

从项目根目录，在装有numpy／Manifold的环境生成和检查；命令写本候选，先保存需保留的结果：

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) '.deps')
python mujoco_arm/tools/build_shoulder_front_liner_keepers.py
python mujoco_arm/tools/check_shoulder_front_liner_keepers.py
python mujoco_arm/tools/audit_shoulder_scapula_insertion.py
python mujoco_arm/tools/verify_shoulder_scapula_insertion.py
```

审查和查看默认已指向本187对象候选；历史175版用显式--assembly／--out、--audit-dir、--review-dir隔离输出。使用装有MuJoCo／Pillow的环境查看：

```powershell
python mujoco_arm/tools/view_shoulder_scapula_insertion.py
python mujoco_arm/tools/view_shoulder_scapula_insertion.py --render
python mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_front_liner_keepers/assembly_check.json'
```

单独重绘衬垫背侧定位柱与螺母入口：

```powershell
python mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_front_liner_keepers/assembly_check.json' --part 'clamp:medial_mid/front_liner' --render --azimuth 270 --elevation 25 --distance 0.045
```

空格播放／暂停，0就位，1取出20 mm。参数可编辑于脚本／JSON，没有新的SolidWorks原生特征树、公差或受力放行。

下一步需完整GH／肱骨、后衬垫／浮动爪和五金／工具安装，并打印验证定位柱／嵌入螺母承载及保持；实际驱动范围仍未知。原肌肉软接触、通道壁与自由末端问题保留，stage_pass=false，不继承旧动作通过。
