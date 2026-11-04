# 整块三角肌安装研究

本轮以244对象共享远端候选为输入，检查真实完整TPU网格的安装，不以端盖局部路径代替整肌安装。

## 两种未通过的刚体顺序

[只移动肌肉](muscle_entry_study.json)：保留四个前支承及中束骨内螺母，端盖及对应紧固件列为后装。12条满足四端法向外移条件的直线方向，均在首个0.1mm采样处出现实体相交，主要是中束／共享远端支承。各阶段逐件清单以报告为准，有限方向采样不能证明柔性装入不可能。

[肌肉携四支承](support_bundle_entry_study.json)：五件保持相对位置一起移动，相同12方向仍在首个0.1mm采样处失败，主要转为锁骨／肩胛与贴骨支承的导入干涉。临时持件方式尚未设计，本研究不假定操作者能无约束保持所有部件。

两报告中，最终模型的部件未删除；每件均明确为当前固定、移动或后装。失败仅适用于所测刚性直线路径，**不证明实际TPU不可弯曲装入**，也不能未经变形验证称柔性安装已通过。

## 诊断与结构处理

[完整走廊干涉分类](entry_relief_analysis.json)检查源方向`[-0.9428090416,-0.2357022604,0.2357022604]`的0–80mm走廊；0–2mm用0.1mm，其余用1mm采样。越过首次碰撞后的姿态只用于诊断，明确不算可达。

该走廊相交涉及：三处前支承、抱箍上半及1号五金、冈下肌。前支承及冈下肌相交位于声明端部区，抱箍相交超出端部区；直接从TPU扣除整条干涉扫掠不符合当前局部改动边界。所查干涉未进入R2保护区，但这不构成削减任意工作结构的理由。

实际处理采用[金属支承导入让位候选](../s1_deltoid_support_draft/README.md)，只改两件金属，不削TPU、不改骨。其条件路径要求冈下肌、抱箍上半与全部抱箍紧固件后装；相关后装步骤仍需验证，不能据此称整肩安装完成。

## 复现

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) '.deps'
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_deltoid_muscle_entry.py
& D:/MujocoEnv/Scripts/python.exe tools/study_s1_deltoid_bundle_entry.py
& D:/MujocoEnv/Scripts/python.exe tools/analyze_s1_deltoid_entry_relief.py
```

上述研究不要求打印或完整有限元；当前结果不是驱动、强度或实物装配验收。
