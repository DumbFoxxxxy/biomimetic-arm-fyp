# GH局部入口与裸肱骨安装姿态候选

更新：2026-11-02。187名义对象，独立研究分支；尚未替代集成起点，`stage_pass=false`，真实驱动范围未知。

## 修改与证据

在父级[保持盖让位版](../retained_depth_6_relief/assembly_check.json)基础上，仅修改工程肩胛GH附近入口；其余186源及变换不变。相对原187对象前衬垫保持版，累计改变两盖和工程肩胛三源，其他184源不变。没有改肌腹、肱骨、关节中心或骨长。

[当前源清单及全体零位检查](assembly_check.json)记录工程肩胛移除约71.19736 mm³、导出新增差值约0.000465 mm³。四组连续平移包络名义横向／端部让位0.1 mm；源GH中心55 mm邻域外移除为0。肩胛仍为单一正材料实体，回读边界与非流形边均为0。全187对象17391对无超阈值刚性相交，六主动线净空，原柔性相交保留。

初版曾因整骨额外简化及更换STL原点引入邻域外0.01216 mm³差值而被拒绝。最终沿用父级局部坐标原点，不再额外全骨简化；失败网格及脚本保存在[初版记录](../../../../../reviews/S1_2026-11-02/GH_entry_initial_rejected/failure.txt)。材料去除不等于剩余壁厚、骨侧承压或强度通过。

## 六步名义放置

前置条件：台架、夹具和肩胛已处于报告指定阶段；裸肱骨绕既定GH屈曲轴到10°，由外部暂持。未装肌肉、骨干抱箍、肌腱和锁骨骨段等对象逐步列明。0°和10°是检查模型姿态，非实物驱动结果。

| 顺序 | 子组件 | 从最终位置向外的受检方向／行程 |
|---|---|---|
| 1 | 螺母0＋下垫圈0 | 源方向[0.484693, −0.865685, −0.125151]，25 mm |
| 2 | 螺母1＋下垫圈1 | 源方向[−0.484693, 0.865685, 0.125151]，25 mm |
| 3 | 保持盖0＋上垫圈0 | 源方向normalize[−1,−1,−1]，25 mm |
| 4 | 保持盖1＋上垫圈1 | 真实颈轴方向，25 mm |
| 5 | 螺钉0 | 真实颈轴方向，15 mm |
| 6 | 螺钉1 | 真实颈轴方向，15 mm |

安装沿反方向进行；拆卸按6→1逆序。每步固定上下文包括之前已装对象，两片盖／垫圈和侧入螺母需暂持，不预设已具备保持结构。

[10°姿态顺序报告](ordered_posed_placement_checks.json)独立重算各步去／回程，0.1 mm共2612个子组件姿态均未超0.001 mm³判据。[裸肱骨姿态报告](humerus_assembly_pose_search.json)在全部GH五金与固定上下文在场时检查0↔10°往返，每1°采样，共22姿态清。它只验证裸肱骨这一段规定路径，不是最新肩部肌肉驱动范围。

原零位直线顺序、26方向粗搜索及648条两段平移搜索的失败均保留；不能把搜索失败称为数学不可装配。本版通过改变裸肱骨安装姿态解决第二片盖的远端肱骨障碍。

## 球头保持与工具包络

[真实颈轴复验](axial_audit/axial_capture_audit.json)：开盖实际肱骨25 mm清；闭盖肱骨与独立R16诊断球均在1.1 mm首个超阈值采样处受阻。[52方向球头探针](finite_capture_direction_checks.json)包含世界格点、GH基底格点及真实颈轴，全部在0.4–1.1 mm采样处遇材料。仅为有限直线捕获证据，不是任意平移／转动逃逸证明或抗拔强度。

直线长柄工具0包络碰肱骨，不能据此宣称工具可进入。[L形工具包络](L_key_envelope_checks.json)采用2.5 mm六角钥匙候选的R1.5圆形杆包络、短臂8 mm、长臂50 mm、弯部R1.6包络，两螺钉均找到名义静态及±20°、1°步长转动净空。当前螺钉模型无真实驱动孔，包络未验证批头啮合、工具进入、手部空间、螺母反力扳手、材料或扭矩；不能写为已能拧紧。

## 复现与查看

在根目录使用已配置numpy／scipy／manifold3d及`.deps`的几何环境，按以下顺序复算。生成器会覆盖本分支输出，先保存需要保留的结果。

```powershell
python mujoco_arm/tools/open_shoulder_GH_assembly_entries.py
python mujoco_arm/tools/search_shoulder_GH_assembly_pose.py
python mujoco_arm/tools/check_shoulder_GH_entry_sequence.py
python mujoco_arm/tools/check_shoulder_GH_L_key_access.py
python mujoco_arm/tools/check_shoulder_GH_entry_capture.py
python mujoco_arm/tools/audit_shoulder_GH_assembly.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_GH_assembly_review/assembly_entries/assembly_check.json' --out 'tmp/GH_entries_axial_recheck'
```

使用MuJoCo／Pillow环境查看已有网格，不重生成设计：

```powershell
python mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_GH_assembly_review/assembly_entries/assembly_check.json'
```

[187对象零位静态图](available_parts.png)与[工程肩胛实际网格](part_detail.png)为DOF0显示，不表现安装动作、肌肉变形或实际驱动。

## 仍需完成

落实侧入螺母／垫圈暂持、防转与扳手反力；核对实际工具及螺钉驱动孔、进入与预紧；验证新槽剩余壁厚与强度、打印配合；将裸骨步骤接入完整肩部安装，并验证钢丝／TPU柔性驱动。保持原K指标，不转移历史刚性角度。
