# 三角肌共享远端夹持候选：244对象

在234对象中束候选上补齐共享远端的骨侧夹持。新增10件、替换肱骨／三角肌TPU／共享端盖3件，其余231父源行不变。骨骼定位、三线出口和原附着中心保持；S1尚未定版，191默认研究模型与GitHub线上资产不变。

## 连接与改动

端盖双孔从原14mm对称布局改为局部`[4,-4]`、`[12,0]`mm，间距约8.944mm，向肘侧避开三条出线结构。采用两套名义M3×40、前后垫圈、背侧螺母和贴骨背板；前支承位于TPU与肱骨之间。实际肱骨双孔Ø3.4mm。来源筛选见[布局对照](../s1_deltoid_distal_study/README.md)。

端盖及TPU补真实承压耳，连接翼名义法向−5.5至−0.03mm；前支承名义−7.5至−5.5mm，按实际骨面让位，TPU侧留0.01mm名义接口间隙。不是只画孔或把端盖理想绑定到骨上。

力的名义传递路径为：三条钢丝各自终端→共享端盖／TPU端部→前支承→肱骨；双螺钉、背板和螺母闭合夹持。三原套筒仍是未验证止挡参考，不能依据此图确认握持力。真实螺纹、预紧、保持力及强度均未验证。

骨去除约376.876mm³；TPU增591.373、去699.279mm³。修改限制在声明端部R18、法向−18至0.2mm，端部区外去除约7.88e−6mm³；三线R2保护区损失分别约4.25e−7／8.19e−7／7.19e−7mm³，未有意削减通道保护壁。改动包含端部重构，不称原肌腹完全不变。

首版加承压耳时填入后束套筒沉孔，产生2.85976mm³相交，冻结在[失败版本](iterations/filled_stop_seat/frozen_assembly_check.json)。最终按原三处R1.4、深2.6mm沉孔重新开出，保留套筒承压台阶；不是删去止挡套或忽略该碰撞。

## 检查结果

- [生成／回读报告](assembly_check.json)：3081个涉及修改件的配对，无新增或增加超过0.001mm³判据的相交；六主动线最终位置畅通。
- [独立复核](verification_checks.json)：13件精确OBJ回读边界、非流形、重复顶点三角形、零面积均为零，每件一个正材料分量。骨／前支承／TPU／端盖接触链、前后螺钉承压及四切向保持探针通过局部几何判据。
- 三个原套筒零位均不撞端盖，向内0.06mm会遇承压台阶。仅证明沉孔几何保留，不证明套筒压接后能握住钢丝。
- 前侧细工具杆和背侧OD8、内AF5.7的名义套筒包络未发现相交；真实工具啮合、操作及扭矩未验证。
- 244件MuJoCo静态加载通过，0DOF；该查看场景不计算TPU变形或驱动动作。

额外执行[整肩零位审计](whole_neutral_audit.json)：244源哈希一致，29,646个源对象配对全部检查，**无整对排除骨骼或接触**，未发现超过0.001mm³判据的实体相交；六主动线对所有其他源检查通过，仅跳过与自身重合的同一钢丝参考几何。所有当前源均有体归属，无遗留退役键。旧TPU—肱骨约0.002013mm³重叠在本次端部修订后减至约1.55e−7mm³。

上述结论只针对规定数字零位与导入精度，不是最小间隙保证、全动作范围、完整安装或实物验收。

## 条件装配路径

7个已采样阶段，共1034个去回姿态无实体干涉：

1. 未装三角肌和共享盖时，前支承沿法向20mm进入，0.5mm步长。
2. 共享盖沿源方向`[-0.0931338034,-0.4963254975,-0.8631263495]`侧向16mm装入，0.1mm步长。
3. 背板从背侧35mm进入，0.5mm步长。
4. 0号螺钉／前垫圈沿法向45mm进入，0.5mm步长。
5. 1号螺钉／前垫圈同样装入。
6. 0号背垫圈／螺母从背侧15mm进入，0.5mm步长。
7. 1号背垫圈／螺母同样装入。

步骤1与2之间三角肌需放到最终位置，**该放入过渡尚未验证**。端盖在144条15°方向粗筛中只有一条16mm走廊通过，选定方向又以0.1mm往返检查；局部走廊通过不证明操作者可从任意外部位置取放。

侧向装盖会跨过三条规划钢丝，因此必须**先装盖，再穿三条钢丝及装终端**。三个旧止挡套保留在最终源清单中，但作为未验证后装参考归位，没有计作已通过的安装阶段。实际穿线／压接、三角肌整体安装及全肩工序联动仍待解决，`stage_pass=false`，无驱动范围结论。

前侧实际几何：

![前侧](detail/available_parts.png)

背侧实际几何：

![背侧](detail/rear_bearing.png)

## 复现

项目根目录PowerShell，使用已有几何环境及`.deps`：

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) '.deps'
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_deltoid_distal_fasteners.py
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_deltoid_distal_offset.py
& D:/MujocoEnv/Scripts/python.exe tools/build_s1_deltoid_distal_mount.py
& D:/MujocoEnv/Scripts/python.exe tools/verify_s1_deltoid_distal_mount.py
& D:/MujocoEnv/Scripts/python.exe tools/audit_s1_mechanical_candidate.py
& C:/Python314/python.exe mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_deltoid_distal_mount/assembly_check.json'
```

下一步检查整块三角肌／其余肌肉安装和钢丝终端，再统一肩肘共用肱骨。打印及完整有限元不新增为机械阶段或GitHub展示的前置门槛。
