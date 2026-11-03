# 项目当前状态

更新：2026-11-02。当前推进肩部S1；**最新集成几何候选为175对象接收部装入修正版（原锁骨支承版及父级保留），尚未定版，预测和实测驱动范围均未知。**

## 路线与阶段

采用真实解剖外形＋工程关节＋钢丝／异形TPU晶格肌肉，允许短肌腱。S1理想模型固定锁骨、肩胛骨，先研究GH肩关节，首轮手拉输入。肩胛／锁骨固定已有几何候选，但真实夹紧、整体安装与工作台适配未验收；固定骨骼仍是理想模型假设。

| 部分 | 当前交付 | 尚未完成 |
|---|---|---|
| 肘部 | 既有候选及试验准备独立保留 | 实际驱动能力；本轮不重新验收 |
| 肩部S1 | GH保持关节对照、60 mm抱箍、两块肩袖肌自由末端及近端支承 | 完整固定、柔性运动、样件及性能验证；新S1未选定 |
| 晶格方向样件 | 8类20对象、24份空白测量记录 | 打印工艺与实际测量；未自动替换整肌腹 |
| 整肩／整臂 | 长期目标保留 | 肩胛运动、多肌肉集成及整臂驱动 |

## 本轮整骨装入审查

[接收部修订与步骤审查](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_fixture_assembly_review/README.md)修正了1 mm采样漏过的起步干涉：两支架局部移除约0.242／0.062 mm³，其余173个源文件不变。175对象零位15225对无超阈值刚性相交，90对父级接触和六肩胛M3孔路复核通过，六主动线保持净空。

在三前衬垫预置、GH／肱骨／锁骨骨段／肌腹等尚未装入的明确条件下，肩胛前上方20 mm路径已独立完成去／回程各201姿态、0.1 mm检查；最大残差0.00014782 mm³，最小间隙含0，不是制造间隙合格。前衬垫实际保持、余下骨架／五金和工具步骤未验收。

[MuJoCo三姿态图](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_fixture_assembly_review/scapula_insertion_demo.png)与交互脚本只显示45对象部分步骤，单骨mocap平移，不是肩部驱动。当前完整175对象清单见[修订报告](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_fixture_assembly_review/receiver_relief/assembly_check.json)。

## 父175对象锁骨支承证据

[175对象锁骨抱夹／支承候选](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/bench_mount_candidate/scapula_clamp_candidate/clavicle_mount_candidate/README.md)新增两半夹体、双衬垫、骨外M3、支承及两组M4背板接口，原骨／肌腹／钢丝保持不变。175对象回读、15225对零位检查无超阈值刚性相交，六主动线相交为0，21对有限接触和四紧固孔已检查。

所需螺钉轴向清出及前夹体20 mm抽出已检查；未装肌肉时，锁骨20 mm候选取出路径清，完整固定未变形肌腹场景0.6 mm接触三角肌。当前需先骨／夹具、后肌腹；完整装配、工具、松衬垫保持及真实夹紧/滑移/支承强度仍未验证。[独立检查](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/bench_mount_candidate/scapula_clamp_candidate/clavicle_mount_candidate/readback_contact_checks.json)不代表驱动或真实固定通过。

## 父155对象肩胛夹持证据

[155对象肩胛夹持候选](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/bench_mount_candidate/scapula_clamp_candidate/README.md)新增三区双面TPU衬垫、浮动爪、骨外M3，并将固定爪连接进两支架。骨骼、肌腹和钢丝源文件不变。155对象回读、11935对零位检查无超阈值刚性相交，六主动线相交为0，36对接触与六骨外M3孔路已检查。

三浮动爪在取下本区五金后，可沿受检方向抽出20 mm；六M3在未装肌肉场景可按所需轴向行程清出，但外侧两枚装肌肉后16.0／22.3 mm碰三角肌，当前需先装夹具再装肌腹。完整骨骼装入、螺母／垫圈、工具和衬垫保持未验收。锁骨固定、预紧／摩擦／强度和驱动仍未知，不能视为完成真实固定。

## 父122对象台架接口证据

[122对象台架接口候选](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/bench_mount_candidate/README.md)新增公共背板／支脚、隔套、穿透螺栓及局部台面接口，替换两支架基座孔座，保留骨侧位置。122对象回读、7381对零位检查无超阈值刚性相交，六路主动线相交为0。骨侧夹爪和锁骨固定仍缺，台面片只是外部固定支承假设。

六根螺钉在仅台架场景可沿轴抽出；完整未变形场景中外侧第二M4抽出7.7 mm碰肩胛。应先装台架再装骨骼，工具及完整装配顺序尚未验证。[独立检查](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/bench_mount_candidate/readback_access_checks.json)不代表运动、强度或制造通过。

## 父92对象肌肉／近端支承证据

- [92对象配合报告](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/assembly_check.json)。应读取fit_revision，父目录初版保留已知占位。
- [近端支承结构、修改量及复现](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/README.md)。
- [独立孔路／接触检查](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/mount_contact_checks.json)、[核验汇总](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/verification_checks.json)。

已回读92对象；四M3孔轴通畅，8处有限法向接触探针遇到材料，六路内部主动钢丝对全部对象未超0.001 mm³相交阈值。新增近端前支承、背板及紧固件共20对象，不代表20件TPU打印件或实际安装完成。

原三角肌／肩胛约1.84044 mm³柔性重叠尚未判定，通道壁布尔／导出差值也未解决。螺钉预紧、薄骨承压、压接、卡扣和台架夹持均未验证。报告保持stage_pass=false；MuJoCo图为实际名义网格静态查看，DOF=0。

## 历史结果如何使用

24对象紧凑端盖、34对象GH保持关节、46对象抱箍、72对象肌腱连接均保留为独立对照。24对象的屈曲＋120°（检查上限）、外展＋61.1°及长轴＋79.3°／−47.0°是旧条件刚性探针；34／46对象有另一组范围。

72／92对象的两个肌腹末端改为自由末端，不能套用此前刚性角度或24对象的零位直线力矩覆盖结论。旧刚性脚本会拒绝这种错误约束。完整比较见[肩部技术文件第4节](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)。

## 下一阶段优先级

1. 完善三前衬垫实际保持，检查后衬垫／浮动爪／GH／五金和工具的完整安装顺序，补其余模块固定与实际工作台接口。
2. 完成通道壁与终端检查，确定工艺并开展方向样件、压接／卡扣及单肌肉试验。
3. 建立自由末端—肌腱—柔性肌腹运动与受力模型，保留未知参数，与实测对照。
4. 在保持GH上先验证一对拮抗肌，再扩大布局；重算往返及联合动作，随后进入整肩／整臂。

机械保持关节与新接口已形成候选，用户已授权必要时重设计；尚未作新S1最终取舍。K1–K7目标与原底线保留，不用降低指标或转移旧通过结果制造新通过结论。

## 文档入口

- [肩部技术方案与证据边界](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)
- [有效规格与验收边界](design/CURRENT_BASELINE.md)
- [设计任务书及历史](design/DESIGN_BRIEF.md)
- [方向／通道样件](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_actuation_calibration/axis_coupons/README.md)
- [技术主导规则](PROJECT_CONTROL.md)、[工作计划](task_plan.md)、[工作日志](progress.md)
- [本次整理前快照](reviews/S1_2026-11-02/document_optimization_before/SHA256.json)
