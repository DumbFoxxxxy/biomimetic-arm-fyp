# 方案二：人体解剖外形＋工程关节＋异形TPU晶格肌肉

以人体骨骼和肌肉模型作为独立几何参考，用工程关节明确运动，再用钢丝牵引TPU晶格肌腹。当前在肘部设计候选与试验准备阶段停止，不继续肩部肌肉设计。

| 目录 | 内容 |
|---|---|
| `输入/骨骼/` | 原人体骨骼及支撑参考，原骨形与修订几何明确区分 |
| `输入/muscle/` | 已筛选的人体肌肉参考模型 |
| `输出/anatomical_concept/` | 当前独立人体骨骼路线，含肘部候选、接口、试验件与诊断 |
| `输出/anatomical_arm/` | 早期A1探索，属于历史，不作为当前定版 |
| `脚本/anatomical_concept/` | 当前方案脚本入口，源文件保留原运行位置 |
| `脚本/anatomical_arm/` | 早期A1脚本入口 |

- [当前肘部交付与验证边界](输出/anatomical_concept/ELBOW_HANDOFF.md)
- [肘部设计思路](输出/anatomical_concept/肘部设计思路.md)
- [当前候选入口](输出/anatomical_concept/CURRENT_CANDIDATE.md)
- [打印对比样件与记录方法](输出/anatomical_concept/elbow_test_coupons/README.md)

当前肘部方案含两块真实异形TPU肌肉、连续内部钢丝通道、工程肘铰链和可转动远端接口。0–150°只表示规定刚性路径检查，实际TPU驱动范围未知，未制造定版。

在 `D:/XJTLU/PSP` 中可继续使用原查看命令：

```powershell
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_continuous_elbow.py'
```

也可通过本目录的脚本入口调用同一文件。完整BP3D源库仍位于外部 `D:/XJTLU/isa_BP3D_4.0_obj_99`，本次未移动外部资料。

注意：当前MuJoCo原生接口直接打开含中文的真实XML路径会失败。请使用上述查看脚本或原ASCII兼容路径；新目录中的脚本入口已验证可用。
