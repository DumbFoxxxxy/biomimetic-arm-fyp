# 三角肌中束侧装螺母夹持候选：234对象

在[前束227对象工作候选](../s1_anterior_mount/README.md)上补中束近端固定，保留原骨骼位置、肌肉中心和出线口。新增7件、修改3件，其余224父源行不变。默认191研究模型和GitHub线上展示不变，S1尚未定版。

## 机械结构

中束肩胛背侧紧邻肱骨，未采用穿骨后加螺母。名义Ø5×6盲孔嵌件也因侧壁空间不足而未采用，筛选失败保留在[方案对照](../s1_middle_blind_study/README.md)。

改用两块侧向滑入骨槽的六角M3螺母，前侧螺钉穿过端盖、TPU连接翼及贴骨前支承后拧入螺母。螺母槽上方承压面闭合反力链，下方和侧面限制位移／转动；螺钉装入前需人工暂时保持螺母，不宣称永久卡扣。螺纹只采用名义包络，实际啮合、预紧和承载尚未验证。

| 项目 | 当前几何 |
|---|---|
| 两孔局部坐标 | `[0,-8]`、`[8,-2]` mm，原端盖坐标基 |
| 孔间距 | 10mm；端盖中心不动 |
| 螺钉 | 名义M3×14、M3×12，各配前垫圈 |
| 螺母／槽 | AF5.5×2.4mm；AF5.7×2.7mm槽 |
| 两侧装走廊 | 同一局部切向，13mm／10mm |
| 名义局部骨包络 | 槽顶1.5mm、槽底1mm、Ø8.6mm；切槽后另查材料见证区 |
| 前支承名义层 | 法向−3.5至−1.5mm，按真实骨形让位 |
| TPU连接翼 | 法向−1.5至−0.03mm，补两个孔的承压耳 |
| 支承—TPU名义让位 | 0.01mm，未改变判据 |

新增：两螺母、两螺钉、两垫圈和一个前支承。修改肩胛骨、中束端盖及一体三角肌TPU。骨去除约411.300mm³；TPU增加68.151、去除623.938mm³，属于声明的端部R12、法向−18至0.2mm内的重做。端部区外损失约2.19e−6mm³，FJ1467钢丝R2保护区损失约2.94e−7mm³。此为有意的端部重构，不称原TPU完全不变。

首版清理回读产生0.001589mm³支承—TPU重叠，超0.001mm³判据，已保存[首版冻结报告](iterations/contact_roundtrip/frozen_assembly_check.json)及修改网格；未放松判据。最终增加实际0.01mm名义接口让位，数值清理增减量仍在报告中单列。

## 回读与局部装配

[assembly_check.json](assembly_check.json)：2285个涉及修改件的配对无新增／增加超过判据的相交，六主动线最终位置清楚。继承的TPU—肱骨约0.002013mm³相交仍存在，未据本次检查判为合法接触；原TPU—肩胛约1.840436mm³减少至约5.75e−6mm³。报告是相对父版的局部改动检查，不是全肩零碰撞结论。

[verification_checks.json](verification_checks.json)重新回读10件修改／新增精确OBJ：边界、非流形、重复顶点三角形和零面积均为零，每件一个正材料分量。234件MuJoCo静态加载通过，0DOF，不表示TPU变形模拟。

切完侧槽后，两螺母各自指定的槽顶、槽底材料环缺失体积均为零；螺母上下移动0.2mm会遇骨，±5°转动会遇限转面；螺钉装入后沿侧槽移动0.3mm会遇螺杆。上述探针证明名义几何限制存在，不证明保持力或强度。

六个条件装配阶段，共1012个去回采样姿态无实体干涉：

1. 三角肌及中束盖尚未装入时，0号螺母沿13mm侧槽装入，步长0.1mm。
2. 1号螺母沿10mm侧槽装入，步长0.1mm。
3. 前支承沿法向15mm进入，步长0.5mm。
4. 中束盖沿源方向`[0.2485632982,-0.9552239751,0.1605099506]`侧向12mm装入，步长0.1mm。
5. 0号螺钉／垫圈沿法向30mm进入，步长0.5mm。
6. 1号螺钉／垫圈同样装入。

步骤3与4之间需要把三角肌放到最终位置，**该放入过渡没有验证**。端盖方向经144条15°粗采样筛选，只有纯切向找到11条粗采样候选，选定路径再以0.1mm往返检查。该路径跨过FJ1467钢丝，所以必须先装盖、后穿中束钢丝；实际喂线仍未验证。

前支承—骨、支承—TPU和端盖—TPU接触探针均有响应；前侧细工具杆包络未发现相交。实际工具啮合、手部操作、人工保持螺母、预紧、TPU加载及全肩安装仍未验，`stage_pass=false`，不报告驱动范围。

![中束实际接口局部](detail/available_parts.png)

## 复现

项目根目录PowerShell，沿用已配置几何环境和`.deps`：

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) '.deps'
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_middle_blind_mount.py
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_middle_offset_inserts.py
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_middle_side_nuts.py
& D:/MujocoEnv/Scripts/python.exe tools/build_s1_middle_side_nut_mount.py
& D:/MujocoEnv/Scripts/python.exe tools/verify_s1_middle_side_nut_mount.py
& C:/Python314/python.exe mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_middle_side_nut_mount/assembly_check.json'
```

下一步为共享远端接口、整块肌肉安装及全肩机械复核，再统一肩肘共用肱骨。本候选不触发GitHub定版发布。
