# 查看与复现：无需兼容目录

本项目直接使用两套方案的实际目录，不需要目录联接、符号链接或恢复旧文件夹。源码实际保留在 `mujoco_arm/`。

## 当前肘部

在项目根目录执行（使用已安装所需依赖的Python环境）：

```powershell
python mujoco_arm/anatomical_concept/view_continuous_elbow.py
python mujoco_arm/anatomical_concept/view_paired_swivel.py --keepers --motion
```

本机也可以继续指定原解释器：

```powershell
& 'C:/Python314/python.exe' 'mujoco_arm/anatomical_concept/view_continuous_elbow.py'
```

第一条显示实际肌肉制造网格的零位，第二条仅展示规定刚性运动。二者都不证明TPU实际驱动范围。
脚本可通过绝对路径从其他工作目录调用，项目根目录从源码位置确定。

仅检查模型加载、或离屏渲染：

```powershell
python mujoco_arm/anatomical_concept/view_continuous_elbow.py --check
python mujoco_arm/anatomical_concept/view_continuous_elbow.py --render
python mujoco_arm/anatomical_concept/view_paired_swivel.py --keepers --motion --check
```

## 生成候选与试验样件

这些命令会写入真实方案输出目录；仅在需要重新生成时执行。

```powershell
python mujoco_arm/anatomical_concept/build_elbow_test_coupons.py
python mujoco_arm/anatomical_concept/build_continuous_muscle_channels.py --shifted --extensor-reanchor
```

本机原 `D:/MujocoEnv/Scripts/python.exe` 环境可以继续使用。源码保持原有数值算法和输出类型；
依赖如numpy、scipy、manifold3d、mujoco、Pillow，以及CAD导出的OpenCascade仍需另外安装。
本次没有将CAD软件、Python环境或外部完整BP3D库打包进仓库。

## 克隆到新目录

```powershell
git clone https://github.com/DumbFoxxxxy/biomimetic-arm-fyp.git
cd biomimetic-arm-fyp
git lfs install
git lfs pull
```

配置所需Python环境后直接运行上述脚本，**不要再运行目录恢复操作**。
原 `tools/restore_workspace_links.ps1` 已退役为提示入口，不创建任何联接。
已验证带中文和空格的替代项目目录。MuJoCo资源由Python读取并通过内存VFS提供给原生接口，
请用查看脚本加载；不能据此保证其他软件原生的“打开XML文件”入口支持中文路径。

## 路径和历史记录

- 原始人体骨骼：`02_人体骨骼晶格肌肉方案/输入/骨骼/`。
- 肌肉参考：`02_人体骨骼晶格肌肉方案/输入/muscle/`。
- 当前输出：`02_人体骨骼晶格肌肉方案/输出/anatomical_concept/`。
- 旧方案CAD与输出：`01_工程骨架方案/SolidWorks/`、`01_工程骨架方案/输出/`。

`project_runtime/` 集中解析历史清单里的旧地址，包括历史工作站根目录。
它把地址转换为真实路径，不创建虚拟目录，也不修改第三方pathlib或MuJoCo。
历史证据JSON及制造网格保持原字节；新生成的清单记录实际相对路径。
旧文档可能仍描述当时的兼容目录做法，以本说明为准。

## 路径验证

```powershell
python -m unittest discover -s tests -v
python tools/verify_link_free_runtime.py
```

第二条会在忽略上传的 `tmp/` 中建立小型普通目录副本，验证中文根目录和不同工作目录，
不改生产制造网格。SolidWorks COM连接、GUI装配外部引用和所有历史优化任务没有在本次逐一重跑。
