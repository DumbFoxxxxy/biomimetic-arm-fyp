# S1四组肌肉统一零位候选

2026-11-02。范围为肩部S1，不是整臂或实物放行。

| 模块 | 选用报告 | TPU正材料连通分量 |
|---|---|---|
| 三角肌 | `../unified_deltoid_outlet_socket/check.json` | 1 |
| 冈下肌 | `../unified_infraspinatus_fitted/check.json` | 1 |
| 肩胛下肌 | `../unified_subscapularis_fitted/check.json` | 1 |
| 大圆肌 | `../unified_teres_major_connected/check.json` | 1 |

六条钢丝共用`../routes_outlet_candidate.json`。组合检查强制四份报告使用同一主路径哈希，并逐项核对实际STL。14个肌肉/端盖网格均水密、正材料单连通且回读通过。封闭空腔保留并单独登记，尚未验证打印支持。

## 修改及选择理由

- 按主路径生成连续通道及端部；原端垫和短腱袖并入TPU制造件，不是删除连接件。
- 冈下肌旧源已有断连晶格。平口终止TPU袖口，消除薄盖外不必要的圆头悬浮体；另外裁掉19.852300、9.461000、0.148871mm³三个已分类的独立晶格尾端。它们全部来自旧晶格，与锚点端垫、有效通道壁及过渡区相交均为0。不是以体积阈值删除；原文件保留。
- `unified_infraspinatus_connected`保留了2肋实验，但仍有2片未连接且增加结构，未选入本候选。
- 肩胛下肌首次导出有5条非流形边；以同一0.001mm清理精度再次回读后通过，未改变水密标准。
- 大圆肌保留约2105.08mm³完整端部，以R1、路径长6.3643mm的真实拟合肋连接，新增材料16.3114mm³。肋对骨、盖及相邻肌肉相交为0；强度、柔顺性及新壁厚未验收。

## 联合检查

`actual_assembly_check.json`包含25件：3骨、4整体TPU、10端盖、3枚未压接止挡坯、5夹具件。六条完整钢丝探针对全部25件的相交均为0。实例数减少来自端垫和腱袖合并，不能单凭数量与旧41/36件报告作性能比较。

仍有4对超过0.05mm³的相交：

| 对象 | mm³ | 待处理 |
|---|---:|---|
| 外侧底座—冈下肌 | 57.445934 | 夹具避让或真实变形接触验证；未通过 |
| 外侧底座—外侧立杆 | 36.325895 | 核对合法连接与实际装配 |
| 内侧底座—内侧立杆 | 32.194543 | 同上 |
| 内侧横梁—内侧立杆 | 288.114169 | 同上 |

所有非零交叠也单独保留。骨侧TPU残差约0.00003–0.00294mm³，需区分拟合接触、导出误差与真实穿透；六路净空不是骨架运动通过。

自由壁厚守护区仍有约0.00188–0.00798mm³残差；大圆肌加肋后尚未重算，旧值明确标为父版结果。打印支持、导管反力固定、紧固/压接、防脱、安装路径、摩擦、连续运动及实际驱动范围均未验收。

## 复现

在仓库根目录的项目Python环境运行，本机需将`.deps`及`mujoco_arm/tools`加入Python搜索路径：

1. 三角肌：`build_shoulder_unified_muscle.py --routes <主路径绝对路径>`。
2. 其余三组：`build_shoulder_unified_muscle.py --module <英文模块名> --fit-terminals --routes <主路径绝对路径>`。
3. 大圆肌再运行`connect_shoulder_lattice_fragments.py --module teres_major --parent-report <unified_teres_major_fitted/check.json绝对路径>`。
4. `check_shoulder_manufactured_assembly.py`重复传入4个`--unified-report <上表报告绝对路径>`，生成本目录报告。

查看命令（本机MuJoCo/Pillow环境）：

```powershell
& 'C:/Python314/python.exe' 'mujoco_arm/tools/view_shoulder_manufactured.py' --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_port_path_revision/whole_shoulder_static/actual_assembly_check.json'
```

加`--render`生成下图。场景DOF为0，只显示实际制造网格，不是拉线或柔性变形仿真。

![25件MuJoCo静态网格](available_parts.png)
