# A1 人体骨骼与异形肌肉原型

这是独立的研究候选，不能替代 V2-S1。原始骨骼没有改写。

## 复现

在 `D:\XJTLU\PSP` 的 PowerShell 中运行：

```powershell
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_arm/fit_elbow.py'
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_arm/build.py'
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_arm/audit.py'
& 'D:/MujocoEnv/Scripts/python.exe' 'mujoco_arm/anatomical_arm/audit_scapula_tasks.py'
& 'D:/Python3.10.6/python.exe' 'mujoco_arm/anatomical_arm/export_cad.py'
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_arm/view_arm.py' --render
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_arm/render_motion.py'
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_arm/view_arm.py'
```

查看器：M 切换真实肘部晶格；L 切换全臂候选体积；H 切换工程关节参考件；R 回零。
修改关节坐标后自动隐藏未变形肌肉，避免把刚体肌肉误认为柔性运动仿真。
默认没有重力或驱动控制，使用理想关节查看。MuJoCo 不用骨骼凸包碰撞做验收；
验收脚本使用原始三角网格的实体布尔运算。

修改 `muscle_recipe.json` 可重生成肘部肌腹位置、宽度和厚度。
`build.py` 中的 muscle() 接收 A/C/B、横截面尺寸，生成弯曲变截面 BCC 晶格、
外部交叉网格、连续孔道及带孔端板；不是长方体占位。
表面简化容差0.01毫米。STEP 是分面 BRep，可编辑但没有还原原生参数特征树。

## 不能据此得出的结论

- 骨骼动作通过，不意味着关节参考件已经装得上或具有防脱能力。
- A1 中的销轴和轴套仍是参考几何；缺少配套孔腔、连接桥和完整 GH 保持器。
- 全臂半透明体积只用于布局，未完成骨面投影、肌肉间避让和驱动设计。
- 肘部晶格是可导出的几何样件，安装端、导管反力支架和骨骼夹具仍需集成。
- A/C/B 是候选接口，不是可以直接钻孔的坐标。
- 固定零位 C 到运动 B 的距离变化仅是路径参考，不能等同于 TPU 肌腹缩短量。
- 肌肉力、恢复能力、实际可用行程和疲劳尚未测量。
- SC2/AC1 机构的指定滑移/内外旋姿态不一定可实现；另有跟踪残差报告，
  不用 SC 转角冒充毫米滑移，也不新增隐藏自由度。

输出统一位于 `mujoco_arm/exports/anatomical_arm/A1`。
