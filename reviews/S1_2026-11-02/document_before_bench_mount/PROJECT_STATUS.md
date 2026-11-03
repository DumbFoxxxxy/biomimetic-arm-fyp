# 项目当前状态

更新：2026-11-02。当前推进肩部S1；**最新集成几何候选为92对象近端支承配合修正版，尚未定版，预测和实测驱动范围均未知。**

## 路线与阶段

采用真实解剖外形＋工程关节＋钢丝／异形TPU晶格肌肉，允许短肌腱。S1理想模型固定锁骨、肩胛骨，先研究GH肩关节，首轮手拉输入。现有台架夹持与工作台连接不完整，固定骨骼仍是模型假设。

| 部分 | 当前交付 | 尚未完成 |
|---|---|---|
| 肘部 | 既有候选及试验准备独立保留 | 实际驱动能力；本轮不重新验收 |
| 肩部S1 | GH保持关节对照、60 mm抱箍、两块肩袖肌自由末端及近端支承 | 完整固定、柔性运动、样件及性能验证；新S1未选定 |
| 晶格方向样件 | 8类20对象、24份空白测量记录 | 打印工艺与实际测量；未自动替换整肌腹 |
| 整肩／整臂 | 长期目标保留 | 肩胛运动、多肌肉集成及整臂驱动 |

## 最新几何与证据入口

- [92对象配合报告](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/assembly_check.json)。应读取fit_revision，父目录初版保留已知占位。
- [近端支承结构、修改量及复现](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/README.md)。
- [独立孔路／接触检查](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/mount_contact_checks.json)、[核验汇总](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/verification_checks.json)。

已回读92对象；四M3孔轴通畅，8处有限法向接触探针遇到材料，六路内部主动钢丝对全部对象未超0.001 mm³相交阈值。新增近端前支承、背板及紧固件共20对象，不代表20件TPU打印件或实际安装完成。

原三角肌／肩胛约1.84044 mm³柔性重叠尚未判定，通道壁布尔／导出差值也未解决。螺钉预紧、薄骨承压、压接、卡扣和台架夹持均未验证。报告保持stage_pass=false；MuJoCo图为实际名义网格静态查看，DOF=0。

## 历史结果如何使用

24对象紧凑端盖、34对象GH保持关节、46对象抱箍、72对象肌腱连接均保留为独立对照。24对象的屈曲＋120°（检查上限）、外展＋61.1°及长轴＋79.3°／−47.0°是旧条件刚性探针；34／46对象有另一组范围。

72／92对象的两个肌腹末端改为自由末端，不能套用此前刚性角度或24对象的零位直线力矩覆盖结论。旧刚性脚本会拒绝这种错误约束。完整比较见[肩部技术文件第4节](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)。

## 下一阶段优先级

1. 补完整肩胛／锁骨夹持、工作台连接、其余模块骨侧固定及工具进入空间。
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
