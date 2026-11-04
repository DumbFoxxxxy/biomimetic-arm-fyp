# 三角肌前束夹持候选：227对象

本候选在后束217对象工作结构上补齐前束近端的名义机械连接。原191默认模型及全部父版本保留；S1仍未定版，不替换GitHub展示。

## 结构与改动

- 保持原端盖中心、钢丝出口、骨骼定位及肌肉主体；14mm双孔线转80°，不是旋转整块肌肉。
- 锁骨实际双孔Ø3.4mm，名义M3×45穿骨紧固；中心射线骨厚约22.20／22.49mm，不能据此推强度或孔缘裕量。
- 新增贴骨前支承、背板与两套螺钉／前后垫圈／螺母，共10件；修改锁骨、三角肌TPU及前束端盖3件，其余214父源行不变。
- 前支承补足骨到TPU反力路径；TPU端部增加连接翼，端盖补实际双孔承压耳。钢丝通道及端面名义0.03mm间隙保留。
- 锁骨去除397.0634mm³，TPU增加75.1242、去除72.1678mm³；声明端部区外去除约8.70e−8mm³，原钢丝R2保护区损失约1.27e−7mm³。数值清理增减另见报告。

## 检查结果与安装顺序

[assembly_check.json](assembly_check.json)记录2860个涉及修改件的配对，无新增或增加超过0.001mm³判据的相交；六主动线最终位置均畅通。13件精确OBJ回读的边界、非流形、重复顶点三角形和零面积检查通过，各有一个正材料分量；227件MuJoCo加载，0DOF静态查看。

前支承—骨、前支承—TPU、背板—骨和端盖—TPU的规定微位移接触探针均有响应。探针证明名义几何接触，不证明接触压力、夹紧力或TPU保持力。

原整盖沿法向20mm进入失败，26个去回姿态撞自身TPU，最大约4.7564mm³；原记录和校验器保存在[iterations/normal_cap_entry](iterations/normal_cap_entry/verification_checks.json)，不隐去失败。

[cap_entry_study.json](cap_entry_study.json)对144条15°采样方向检查零位连续12mm走廊，6条纯切向候选通过粗采样。选定源方向`[0,-0.0424696178,0.9990977588]`，按0.1mm复核往返242姿态通过源对象碰撞检查。该路径会跨过尚未安装的FJ1468主动线，因此必须**先装端盖，再穿前束钢丝**；其他五主动线路径未在本次局部移动中发生超过判据的相交。

[side_entry_checks.json](side_entry_checks.json)按以下条件顺序重放7阶段970姿态，未发现实体干涉：

1. 三角肌尚未装入时，前支承沿法向20mm进入。
2. 端盖沿选定侧向12mm走廊装入；本阶段前束钢丝未装。
3. 背板从背侧35mm进入。
4. 0号前侧螺钉／垫圈沿法向50mm进入。
5. 1号前侧螺钉／垫圈沿法向50mm进入。
6. 0号背侧垫圈／螺母从背侧12mm进入。
7. 1号背侧垫圈／螺母从背侧12mm进入。

步骤1与2之间需将三角肌放到最终位置，**此放入过渡未验证，不能计为装配通过**。之后的实际穿线与预紧也未验证，不计入上述7个已采样阶段。其余对象均保留在阶段账目中，未整对忽略相邻骨骼。

前侧细工具杆名义包络清楚；真实内六角啮合、后螺母反夹持工具、手部空间、螺纹／预紧、TPU加载变形和整肩安装均未验证。本候选只完成局部名义连接与条件装配检查，`stage_pass=false`，无驱动范围结论。

![局部实际几何](detail/available_parts.png)

## 复现与查看

在项目根目录PowerShell执行，几何检查使用现有含SciPy环境和`.deps`；可替换为自行配置的等效Python环境。

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) '.deps'
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_anterior_mount.py
& D:/MujocoEnv/Scripts/python.exe tools/build_s1_anterior_mount.py
& D:/MujocoEnv/Scripts/python.exe tools/verify_s1_anterior_mount.py
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_anterior_cap_entry.py
& D:/MujocoEnv/Scripts/python.exe tools/verify_s1_anterior_mount.py --side-entry
& C:/Python314/python.exe mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_anterior_mount/assembly_check.json'
```

前束支持的名义力链为：端盖／TPU连接翼→前支承→锁骨，双螺钉与背板闭合夹持。该路径尚需实体材料、预紧和防脱验证。本轮不选择打印参数。
