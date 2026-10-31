# 右臂仿生项目文献知识库

- 源目录：`D:\XJTLU\PSP`（含 `literature_inbox/` 新文献投放目录）
- PDF 文件：**50** 个 → 去重后 **47** 篇独立文献
- 全文提取：`fulltext/`（每篇一个纯文本文件，按页分隔，可全文检索）
- 单篇笔记：`notes/`（含中文摘要、要点、与本项目的关系）
- 机器可读：`manifest.json`（原始提取）、`manifest_enriched.json`（含注解）、`documents.csv`
- 关键词检索：`python search.py 关键词`

等级含义：**核心** = 方案直接依据；**重要** = 肩肘机构与平台设计主要参考；**参考** = 方法储备与对照路线；**背景** = 相关但非设计文献；**离题** = 与本项目无关。

## 项目直接依据

X 晶格人工肌肉（PSM）、TPU 柔性结构参数化设计、本项目毕设基线、带韧带的肩部肌肉骨骼仿真、钢丝布线优化。这五篇决定你手臂方案的基本做法。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 01 | 2020 | [图案化结构肌肉（PSM）：利用各向异性柔性结构的任意形状钢丝驱动机器人肌肉](notes/01_psm_yoshimura_iros2020.md) | IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) | 核心 |
| 02 | 2026 | [通过参数化设计与 3D 打印柔性结构加速肌肉骨骼机器人开发](notes/02_yoshimura_2026_ais.md) | Advanced Intelligent Systems (Wiley-VCH) | 核心 |
| 03 | 2025 | [由仿生肌肉结构驱动的高仿人机械臂设计（西交利物浦大学本科毕业设计）](notes/03_wang_2025_fyp.md) | XJTLU SAT301 Final Year Project, Bachelor of Engineering | 核心 |
| 04 | 2021 | [含韧带的肩复合体肌肉骨骼仿真构建及其基于模型预测控制的验证](notes/04_sahara_shoulder_ligaments.md) | IEEE 会议论文（PDF 未标注具体会议） | 核心 |
| 05 | 2024 | [肌腱驱动机器人钢丝布线（含可变绕行点）的数值仿真优化设计](notes/05_kawaharazuka_2024_wire_optimization.md) | IEEE Robotics and Automation Letters (RA-L), Vol. 9, No. 2, pp. 1388- | 核心 |

## 肩复合体解剖与运动学建模

把肩带（胸锁、肩锁、肩胛胸壁、盂肱）抽象成运动学模型的工作，是关节中心、轴与自由度取法的依据。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 06 | 2003 | [仿人肩复合体与肱骨指向运动学](notes/06_lenarcic_2003_humeral_pointing.md) | IEEE Transactions on Robotics and Automation, Vol. 19, No. 3, pp. 499- | 重要 |
| 07 | 2005 | [肩复合体建模及其在仿人机器人上肢设计中的应用](notes/07_tondu_2005_shoulder_complex.md) | 2005 5th IEEE-RAS International Conference on Humanoid Robots, pp. 313- | 重要 |
| 08 | 2026 | [人体关节的设计特征及对仿人机器人设计的建议（综述）](notes/08_burgess_2026_human_joints.md) | Bioinspiration and Biomimetics, Vol. 21, 041005（Topical Review, Open Access） | 重要 |
| 09 | 2012 | [用于肌肉骨骼机械臂的仿人肩复合体](notes/09_ikemoto_2012_iros.md) | 2012 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), pp. 4892- | 重要 |
| 10 | 2014 | [仿人肌肉骨骼机械臂的肩复合体连杆机构](notes/10_ikemoto_2014_linkage_preprint.md) | Bioinspiration and Biomimetics 投稿预印本（2014 年 2 月） | 重要 |

## 肩带机构与机器人肩部实现

别人实际做出来的机器人肩部：球窝关节、肩胛骨滑移、肌腱被动稳定。对应你在 SolidWorks 与 MuJoCo 里要复现的结构。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 11 | — | [带肩胛骨的肌肉驱动机器人的肩部结构](notes/11_sodeyama_blade.md) | IEEE 会议论文（PDF 未保留出处信息，正文引用止于 2005 年） | 重要 |
| 12 | 2007 | [采用肩胛骨-锁骨结构的宽活动范围肩部设计与运动](notes/12_sodeyama_2007_bladebone.md) | Proceedings of the 2007 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), pp. 3629- | 重要 |
| 13 | 2008 | [Kojiro 仿人机器人球面胸廓加肩胛骨锁骨的肩部设计与运动](notes/13_sodeyama_2008_kojiro.md) | 2008 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) | 重要 |
| 14 | 2024 | [受肌肉骨骼机械智能启发的高仿生机器人肩部研制与特性](notes/14_yang_2024_tor_shoulder.md) | IEEE Transactions on Robotics, Vol. 40 | 重要 |
| 15 | 2023 | [具备类人生物力学能力与性能的仿生拟人机械臂设计与实证研究（博士论文）](notes/15_yang_2023_phd_thesis.md) | University of Manchester, PhD Thesis, Department of Mechanical, Aerospace and Civil Engineering | 重要 |
| 16 | 2024 | [带被动稳定的大工作空间肌腱驱动肩关节](notes/16_lee_hughes_tendon_shoulder.md) | IEEE/ASME Transactions on Mechatronics | 重要 |
| 17 | 2019 | [面向仿人机器人与人类的球面主动关节](notes/17_mghames_2019_spherical.md) | IEEE Robotics and Automation Letters (RA-L), Vol. 4, No. 2, pp. 838- | 重要 |
| 18 | 2019 | [面向仿人机器人应用的主动球关节](notes/18_folgheraiter_2019_3pss1s.md) | 2019 IEEE Conference on Cybernetics and Intelligent Systems (CIS) and Robotics, Automation and Mechatronics (RAM), pp. 571- | 参考 |

## 肌肉骨骼人形平台（系统级）

东京大学 JSK 一系（Kojiro/Kenzoh/Kengoro/Musashi）与 ECCE 等整机平台，提供肌腱驱动、冗余自由度与模块化的工程经验。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 19 | 2011 | [带非线性弹簧单元可调柔顺肌腱的肌肉骨骼仿人机器人 Kenzoh 研制](notes/19_nakanishi_2011_kenzoh.md) | Proceedings of the 2011 IEEE International Conference on Robotics and Biomimetics (ROBIO), Phuket, Thailand | 重要 |
| 20 | 2019 | [面向学习控制系统研究的模块化肌肉骨骼仿人平台 Musashi](notes/20_kawaharazuka_musashi.md) | IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) | 重要 |
| 21 | 2016 | [面向真实世界物理交互的仿人肌肉骨骼机器人 Kengoro](notes/21_asano_2016_kengoro.md) | 2016 IEEE-RAS 16th International Conference on Humanoid Robots (Humanoids), pp. 876- | 重要 |
| 22 | — | [ECCE1：系列拟人肌肉骨骼上躯干中的第一台](notes/22_marques_ecce1.md) | IEEE 会议论文（PDF 未保留出处信息） | 参考 |
| 23 | 2020 | [含桡尺关节与柔性机加工弹簧手指的仿人前臂与手部设计](notes/23_kawaharazuka_2020_radioulnar.md) | Journal of Robotics and Mechatronics（收稿 2018-09-18，录用 2020-01-19） | 重要 |
| 24 | 2016 | [由复丝肌肉驱动的肌肉骨骼下肢机器人](notes/24_kurumaya_2016_multifilament.md) | Robotics and Mechatronics (Robomech Journal), 3:18 | 参考 |

## 仿真、优化与联合建模

多体仿真、优化与仿真器选型的参考，含主动-被动外骨骼联合优化与人形机器人协同设计。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 25 | 2022 | [主动-被动肩部外骨骼的优化设计：人机交互的计算建模](notes/25_nasr_2022_exoskeleton.md) | Multibody System Dynamics, Vol. 57, pp. 73-106 | 重要 |
| 26 | 2024 | [面向机器人与强化学习的 PyBullet 与 Isaac Sim 评测](notes/26_berrocal_pybullet_isaacsim.md) | IEEE 会议论文（TECNALIA / UPV-EHU） | 参考 |
| 27 | 2025 | [基于代理模型引导优化的仿人机器人设计与运动控制协同适应](notes/27_du_2025_coadaptation.md) | Biomimetic Intelligence and Robotics, Vol. 5, 100255 | 参考 |

## 钢丝张力分配与肩胸接触优化

面向 V4 钢丝路径、正张成、导轮摩擦、单侧接触和肩胛节律建模的同行评审文献。新增论文先进入候选清册，取得可检索全文后再写入正式注解。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 41 | 2017 | [运动学约束冗余索驱并联机器人：建模、冗余分析与刚度优化](notes/41_jamshidifar_2017_stiffness.md) | IEEE/ASME Transactions on Mechatronics, Vol. 22, No. 2, pp. 921-930 | 重要 |
| 42 | 2021 | [面向康复任务的索驱并联机器人多目标优化设计](notes/42_benhamida_2021_nsga2.md) | Mechanism and Machine Theory, Vol. 156, 104141 | 重要 |
| 43 | 2012 | [用动力学仿真确定特定运动的最优多关节肌肉布置](notes/43_asaoka_2012_muscle_arrangement.md) | Proceedings of the 12th IEEE-RAS International Conference on Humanoid Robots, pp. 216-221 | 重要 |
| 47 | 2018 | [完成日常生活活动所需的肩肘活动范围：系统综述](notes/47_oosterwijk_2018_adl_rom.md) | Physiotherapy Theory and Practice, Vol. 34, No. 7, pp. 505-528 | 重要 |
| 48 | 2015 | [索驱并联机器人绞盘内置钢丝力测量的滑轮摩擦补偿与验证](notes/48_kraus_2015_pulley_friction.md) | 2015 IEEE International Conference on Robotics and Automation (ICRA) | 核心 |
| 49 | 2021 | [障碍环境下索驱连续机器人的轮廓与接触力估计](notes/49_ashwin_2021_contact_path.md) | Mechanism and Machine Theory, Vol. 164, 104404 | 重要 |
| 50 | 2024 | [椭球面是否适合建模肩胸滑动平面？](notes/50_blache_2024_scapula_ellipsoid.md) | Journal of Biomechanics, Vol. 164, 111989 | 核心 |
| 51 | 2021 | [索驱动并联机构的力闭合条件](notes/51_luan_2021_wrench_closure.md) | Applied Sciences, Vol. 11, No. 9, 4228 | 核心 |

## 人工肌肉与软体驱动

气动、织物、连续体等其它人工肌肉路线，以及肌腱驱动软体机器人的控制理论。用于对比与借用方法。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 28 | 2012 | [弹性体折纸：可编程纸-弹性体复合气动驱动器](notes/28_martinez_2012_origami.md) | Advanced Functional Materials, Vol. 22, pp. 1376-1384 | 参考 |
| 29 | 2022 | [气动软体机器人：挑战与优势（综述）](notes/29_su_2022_pneumatic_review.md) | Actuators (MDPI), Vol. 11, 92 | 参考 |
| 30 | 2024 | [编码缝纫软体织物机器人](notes/30_guo_2024_encoded_sewing.md) | Science Advances, Vol. 10, eadk3855 | 参考 |
| 31 | 2025 | [肌腱驱动软体机器人的奇异摄动控制：理论与实验](notes/31_ribeiro_2025_singular_perturbation.md) | IEEE Transactions on Control Systems Technology, Vol. 33, No. 5, pp. 1929- | 参考 |
| 32 | — | [用仿人机器人对培养于软生物反应器中的人体细胞施加力学刺激](notes/32_mouthuy_humanoid_bioreactor.md) | 期刊论文（含软生物反应器与仿人机器人肩关节） | 背景 |

## 柔性与晶格材料的本构建模（方法储备）

超弹性与各向异性本构建模、微结构逆向设计的机器学习方法。若后面要把 X 晶格失效应变、非线性刚度写进仿真，或做晶格自动优化，从这里取方法。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 33 | 2021 | [基于机器学习的参数化各向异性有限应变超弹性材料建模及其在超材料优化中的应用](notes/33_fernandez_2021_parametric.md) | International Journal for Numerical Methods in Engineering, DOI 10.1002/nme.6869 | 参考 |
| 34 | — | [基于计算晶体结构数据的各向异性超弹性材料人工神经网络建模](notes/34_im_ann_anisotropic.md) | Seoul National University, Multiscale Mechanical Design Division（PDF 元数据标题为 AIAA 模板） | 参考 |
| 35 | 2022 | [基于平衡方程的卷积神经网络用于超弹性材料的本构建模](notes/35_li_chen_2022_ecnn.md) | Journal of the Mechanics and Physics of Solids, Vol. 164, 104931 | 参考 |
| 36 | 2020 | [计算力学中的几何深度学习（一）：各向异性超弹性](notes/36_vlassis_2020_gdl.md) | Computer Methods in Applied Mechanics and Engineering, Vol. 371, 113299 | 参考 |
| 37 | 2023 | [用去噪扩散算法逆向设计具有精细非线性材料属性的微结构](notes/37_vlassis_2023_diffusion.md) | Computer Methods in Applied Mechanics and Engineering, Vol. 413, 116126 | 参考 |
| 38 | 2025 | [无需应力数据的深度学习用于多区域超弹性属性识别（综述）](notes/38_shi_2025_unsupervised_multiregional.md) | Computational Mechanics, Vol. 76, pp. 117-146（REVIEW ARTICLE） | 参考 |

## 无关文献

与手臂项目无关，仅因同目录被收录。

| 编号 | 年 | 文献 | 出处 | 等级 |
| --- | --- | --- | --- | --- |
| 39 | 2019 | [非线性薛定谔方程的高阶微扰分析](notes/39_bonetti_2019_nls.md) | Communications in Nonlinear Science and Numerical Simulation, Vol. 72, pp. 152-161, DOI 10.1016/j.cnsns.2018.12.010 | 离题 |

## 标签索引

- `肩胛骨`：04 11 12 13 50
- `钢丝张力分配`：41 42 43 48 51
- `锁骨`：04 07 12 13
- `TPU`：01 02 03
- `各向异性`：01 33 34
- `接触优化`：49 50 51
- `正张成`：41 42 51
- `X晶格`：01 03
- `代理模型`：27 33
- `仿真`：04 12
- `冗余自由度`：06 19
- `冗余驱动`：08 22
- `单侧接触`：49 50
- `并联机构`：06 18
- `气动人工肌肉`：09 10
- `盂肱关节`：06 09
- `绳索并联机器人`：42 51
- `肌腱驱动`：11 31
- `肩复合体`：04 07
- `超弹性`：34 35
- `连杆机构`：08 10
- `韧带`：02 04
- `2D到3D`：30
- `3PSS-1S`：18
- `ADL`：47
- `Benchmark`：26
- `Blade`：11
- `CNN代理`：37
- `Capstan摩擦`：48
- `Dahl模型`：48
- `ECA`：15
- `ECCEROBOT`：22
- `EMG控制`：10
- `FDM打印`：01
- `ICA`：15
- `K7`：43
- `Kengoro`：21
- `Kenzoh`：19
- `Kojiro`：13
- `Kojiro后续`：19
- `MCA`：15
- `MISA`：15
- `MPC`：04
- `McKibben肌肉`：24
- `NSGA-II`：42
- `PPO`：26
- `PSM`：01
- `Pareto`：42
- `Sobolev训练`：36
- `Stable-Baselines3`：26
- `V-REP`：18
- `下肢`：02
- `不完整球窝`：14
- `丝杠滑台`：03
- `串并联混合`：17
- `九自由度`：07
- `代谢能量`：25
- `仿真器选型`：26
- `优化效率`：27
- `优化设计`：25
- `传感器驱动一体肌肉模块`：21
- `低刚度失稳`：31
- `保真度`：26
- `偏置双指向系统`：06
- `光通信`：39
- `关节主动力矩`：25
- `关节指标`：47
- `关节模块`：20
- `关节类型`：08
- `冗余系统`：24
- `冗余绳索机器人`：41
- `几何约束`：12
- `凸性`：36
- `分节`：08
- `刚度优化`：41
- `前臂`：23
- `前臂肘部`：15
- `力学刺激`：32
- `力臂`：05
- `力闭合`：51
- `功能活动范围`：47
- `加权图描述子`：36
- `动力学优化`：43
- `协同设计`：27
- `单压力源`：30
- `参数化设计`：02
- `变刚度肌腱`：19
- `可变绕行点`：42
- `可编程变形`：28
- `可达力空间`：05
- `可重构`：20
- `图卷积网络`：36
- `基线`：03
- `复丝肌肉`：24
- `多体动力学`：25
- `多关节肌肉`：43
- `多功能骨骼结构`：21
- `多区域`：38
- `多晶体`：36
- `多目标`：05
- `多目标优化`：42
- `多项式拟合`：03
- `大活动范围`：11
- `奇异`：17
- `奇异摄动`：31
- `奇异锥`：06
- `密集附着`：24
- `导轮`：48
- `导轮摩擦`：48
- `工作空间`：17
- `工作空间参数`：18
- `应力不可测`：35
- `应变能密度`：34
- `康复基准`：47
- `异质材料`：38
- `弹簧`：16
- `强化学习`：27
- `微型骨肌模块`：23
- `微扰法`：39
- `微结构`：37
- `快速成型`：29
- `悬浮系统`：08
- `扩散模型`：37
- `折纸结构`：28
- `抛球实验`：09
- `拟人化`：22
- `接触力估计`：49
- `控制系统`：29
- `旋前旋后`：23
- `无监督`：38
- `无监督学习`：35
- `晶格`：33
- `最大肌力`：04
- `有限应变`：33
- `本构`：35
- `本构张量`：34
- `机加工弹簧手指`：23
- `机械智能`：14
- `柔性躯干`：22
- `柔性驱动器`：15
- `桡尺关节`：23
- `椭球接触面`：50
- `模块化`：20
- `正负99.87度`：16
- `步进电机`：03
- `残余节点力`：38
- `气动`：28
- `气动软体驱动器`：29
- `滑动面`：13
- `球关节`：17
- `球面胸廓`：13
- `生成模型`：37
- `生物反应器`：32
- `相对密度`：33
- `神经网络`：34
- `等效力平衡`：35
- `等效运动学`：07
- `纸弹性体复合`：28
- `线性驱动`：03
- `组织工程`：32
- `织物机器人`：30
- `绕行点`：05
- `统一制造流程`：02
- `综述`：29
- `缓冲垫`：02
- `缝纫约束`：30
- `肌肉布置`：43
- `肌肉模块`：20
- `肌肉骨骼仿人`：32
- `肌肉骨骼模型`：25
- `肌腱`：02
- `肌腱布线`：16
- `肘奇异`：07
- `肩关节`：32
- `肩关节样机`：17
- `肩带`：06
- `肩肘ROM`：47
- `肩胛胸壁关节`：09
- `肩胛运动`：10
- `肩胸滑动`：50
- `肩部`：15
- `肱桡盂肱耦合`：14
- `肱骨头`：14
- `背驱动性`：24
- `胸腔腾空`：11
- `脱位特性`：16
- `自由度`：08
- `自由能`：38
- `自锁`：14
- `表面肌电`：10
- `被动螺旋肌腱`：16
- `计算效率`：26
- `认知`：22
- `调制不稳定性`：39
- `超材料`：33
- `轨迹跟踪`：31
- `转录组`：32
- `软体机器人`：28
- `软组织`：38
- `软部件`：29
- `运动生成`：12
- `连续体`：31
- `逆向设计`：37
- `逆运动学`：18
- `通用骨框架`：20
- `速度空间`：05
- `钢丝路径`：49
- `钢丝驱动`：01
- `障碍约束`：49
- `零空间`：41
- `非线性弹性单元`：20
- `非线性弹簧`：19
- `非线性波动`：39
- `韧带关节囊软骨`：15
- `韧带软骨肌腱`：14
- `顺序驱动`：30
- `驱动器动力学`：31
- `高维设计空间`：27
- `黑盒优化`：05

## 重复与同源文件

- **10 ikemoto_2014_linkage_preprint** 与：Humanlike_shoulder_complex_for_musculoskeletal_robot_arms.pdf（同一工作的会议版）
- **15 yang_2023_phd_thesis** 与：Development_and_Characteristics_of_a_Highly_Biomimetic_Robotic_Shoulder...（同作者的期刊版）
- **22 marques_ecce1** 与：ECCE1_The_first_of_a_series_of_anthropom.pdf（完全相同的文件）
- **24 kurumaya_2016_multifilament** 与：Musculoskeletal lower‑limb robot driven.pdf（完全相同的文件）
- **32 mouthuy_humanoid_bioreactor** 与：Humanoid robots to mechanically stress human.pdf（完全相同的文件）

说明：3 组文件是逐字节相同的重复副本（ECCE1、Humanoid bioreactor、Multifilament lower-limb），已在上面标出；另外 Ikemoto 的肩复合体与 Haosen Yang 的肩设计各有一对同源的不同版本，内容互补而非重复。
