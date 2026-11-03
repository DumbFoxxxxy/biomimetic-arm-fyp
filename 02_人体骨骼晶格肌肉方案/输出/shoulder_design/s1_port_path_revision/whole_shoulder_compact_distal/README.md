# 紧凑远端盖候选：暂停前几何对照

> 状态说明：此处保留候选输入及检查结果，不代表已采纳的重设计方案。设计工作暂停；当前决策状态以 [肩部技术审查稿](../../SHOULDER_DESIGN.md) 为准。

2026-11-02。当前已建模刚性组件的前屈检查到120°；不代表整肩定版或实际驱动范围。

## 当前组合

| 模块 | 检查报告 |
|---|---|
| 三角肌 | `../unified_deltoid_compact_master/check.json` |
| 冈下肌 | `../unified_infraspinatus_compact_distal_fitted/check.json` |
| 肩胛下肌 | `../unified_subscapularis_compact_distal_fitted/check.json` |
| 大圆肌 | `../unified_teres_major_compact_master/check.json` |
| 两夹具 | `../fixture_revision/check.json` |
| 钢丝主路径 | `../routes_compact_distal_candidate.json` |

未改两块远端接口的骨侧名义点、骨骼或关节中心。改变端部高度和形状，同步TPU端垫与通道。

## 结构与修订依据

FJ1500/FJ1504远端均采用两圆端相连的窄盖：长11.4、宽6.4、厚3mm；钢丝孔R0.86与M3候选孔R1.7相距5mm，止挡腔R1.4深2.6mm，保留0.4mm肩部。各配一枚未压接金属止挡坯，外径2.5、内径0.9、长2.5mm。与原大圆盖相比，主要减少关节附近的横向刚性占位。

冈下肌堆叠高度减少4mm，肩胛下肌减少5.5mm；新TPU端垫匹配两圆端轮廓，R4.2、中心距5mm，按骨面修形，补贯通的候选M3孔。骨骼安装孔尚未加工/建模，不能把这些孔当成已完成的骨侧固定。

首次仅按盖体避碰的筛选曾选择7/5.97mm，后续发现冈下肌钢丝入口挤入肱骨且碰止挡。该组完整失败数据保留于`../compact_distal_unfiltered_revision`，未采用。当前筛选增加了1mm轴向入口与R1.61通道壁守护区，重建从端垫侧进入的连接段，实际钢丝/骨和钢丝/止挡相交均0。

原部件、旧主路径及上侧起端候选均保留。三角肌与大圆肌只有在各自路线记录严格相等、网格哈希未变时才复用，整装重新计算。

## 当前证据

24个现有实体包括3骨、4TPU、10盖、5枚止挡坯和2夹具。14个TPU/盖及新增2枚止挡网格均单正连通、水密且回读通过。内部空腔保留并记录，打印支持未验收。

- 6条完整钢丝对全部24件相交均为0；无超过0.05mm³的实体交叠，最大非零贴骨残差约0.002013mm³。
- 新两止挡零位与盖不相交，向内0.04mm不碰肩部、0.06mm被挡；这不是压接握持力或强度试验。
- 新增2项几何回归通过：钢丝/M3轴向通路、止挡间隙与肩部捕获。实物紧固件头部未包括在该通路测试中。

固定GH中心、SC/AC固定、端盖与止挡理想安装。步长1°，首次限制细化0.1°；交叠阈值0.001mm³：

| 校准轴 | 上轮末个未超阈° | 当前末个未超阈° | 首限 |
|---|---:|---:|---|
| flexion + | 24.8 | 120（采样上限） | 受检区间未发现首限，不是最大可达值 |
| flexion − | −51.6 | −51.6 | 肱骨—肩胛骨，−51.7°首次超阈 |
| abduction + | 54.1 | 61.1 | 冈下肌远端盖—肩胛骨，61.2°首次超阈 |
| abduction − | −25.0 | −25.0 | 肱骨—肩胛骨，−25.1°首次超阈 |
| longitudinal + | 80.3 | 79.3 | 冈下肌远端盖—肩胛骨，79.4°首次超阈；回退1°已登记 |
| longitudinal − | −38.4 | −47.0 | 肩胛下肌远端盖—肩胛骨，−47.1°首次超阈 |

轴符号按源坐标，尚未替代解剖内/外旋标定。报告保存末采样最小间隙；接触或未超阈微交叠时可为0，不能当作正制造间隙。原指标保留，不因该结果修改目标或宣布K1–K7通过。

壁厚守护区仍有数值残差，冈下肌约0.008316、肩胛下肌约0.002327mm³；骨侧固定、实际螺钉/头部空间、防转、压接、导管反力、防脱、安装、柔性变形、摩擦、强度和材料行程均未验收。**本方案不能按整套可打印或驱动定版使用。**

## 复现与查看

仓库根目录，以项目环境运行；本机需`.deps`与`mujoco_arm/tools`加入Python搜索路径：

1. `search_shoulder_distal_stacks.py`保留原圆盖高度比较。
2. `search_shoulder_compact_distal.py`筛选紧凑盖，并检查入口壁厚空间。
3. `design_shoulder_compact_distal.py`生成新路线、局部盖及止挡。
4. 冈下肌、肩胛下肌分别运行`build_shoulder_unified_muscle.py --module <英文名> --fit-terminals --routes <routes_compact_distal_candidate.json绝对路径>`。
5. 用`rebind_shoulder_unchanged_routes.py --tag compact_master`重绑定三角肌上侧edge14版及大圆肌connected版；任一模块自身路径改变即拒绝。
6. `check_shoulder_manufactured_assembly.py`传入上表4个`--unified-report`、夹具`--fixture-report`与本目录`--output-dir`；随后运行`check_shoulder_rigid_motion.py --assembly <本目录actual_assembly_check.json>`。

```powershell
& 'C:/Python314/python.exe' 'mujoco_arm/tools/view_shoulder_manufactured.py' --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_port_path_revision/whole_shoulder_compact_distal/actual_assembly_check.json'
```

加`--render`渲染整装，再加`--part 'cap:FJ1500/distal'`只查看该盖。场景DOF为0，不是实际拉线或TPU变形仿真。

![24件静态查看](available_parts.png)

![独立紧凑端盖，周边已隐藏](part_detail.png)
