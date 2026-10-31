# 仿生机械手臂 FYP

GitHub：[DumbFoxxxxy/biomimetic-arm-fyp](https://github.com/DumbFoxxxxy/biomimetic-arm-fyp)（私有）。

| 入口 | 内容 |
|---|---|
| [方案一](01_工程骨架方案/README.md) | 工程骨架、V2-S1及旧钢丝/PSM路线 |
| [方案二](02_人体骨骼晶格肌肉方案/README.md) | 人体骨骼、工程关节、异形TPU晶格肌肉；当前停在肘部候选 |
| [共用资料](00_共用资料/README.md) | 文献、历史任务书、依赖位置 |
| [查看与复现](RUNNING.md) | 无兼容目录的运行命令、克隆及验证方法 |

模型和输出实际存放在各方案目录，源码实际存放在mujoco_arm。无需目录联接或兼容文件夹。
历史记录中的旧地址由project_runtime集中解析到实际目录；新生成的清单写入真实相对路径。
当前肘部仍未制造定版，实际TPU驱动范围待实测。文件与路径改造不改变设计验证结论。
