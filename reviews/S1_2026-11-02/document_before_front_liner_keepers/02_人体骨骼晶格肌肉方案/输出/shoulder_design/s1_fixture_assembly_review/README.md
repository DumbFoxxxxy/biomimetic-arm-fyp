# 肩胛整骨装入步骤：接收部微调与路径审查

2026-11-02。不是肩关节运动验收。本轮从175对象锁骨支承候选出发，只修改两件接收支架，骨、肌腹、钢丝及其他173个源文件保持不变。完整候选仍为175对象，未定版，驱动范围未知。

## 1. 粗采样漏检与修复

原175对象的9个受检前向平移方向中，归一化[0,−1,0.3]在1 mm采样下能到20 mm；但独立0.1 mm去／回程检查发现起步0.1–0.7 mm会穿过两刚性接收部，最大约0.040747 mm³，到1 mm后才消失。这不算从零位可达，不能以粗通过得出装配通过。原审查、细查和脚本保留在[失败记录](iterations/coarse_hidden_barrier/selected_path_dense_check.json)。

当前以两个夹区的实际骨块，沿20 mm路径作连续Minkowski开敞，横向及端部另留0.02 mm模型间隙；不是已选制造公差。未移动骨骼、关节中心或锚点，也未裁剪肌腹：

| 修订对象 | 相对父级移除 | 读回新增差值 |
|---|---:|---:|
| 内侧缘中段接收部 | 0.242427 mm³ | 0.00000418 mm³ |
| 外侧缘上部接收部 | 0.061515 mm³ | 0.00001052 mm³ |

上角接收部不主动改型。两支架仍是单正材料实体；没有以丢弃断片获得通过。精确局部裁剪范围、原点和扫掠方向见[新装配报告](receiver_relief/assembly_check.json)。当前实体为制造网格和参数化脚本，未交付新的SolidWorks特征树或壁厚／强度放行。

## 2. 所选步骤的对象和前置条件

这是先装台架／接收部、再放入肩胛的阶段；不是将完整机械手臂整体平移：

- 移动：仅工程肩胛骨，骨长和朝向不变。
- 固定：44个实际对象，包括两支架、背板／台面及五金、三片前衬垫、空的锁骨后支承及其M4连接。
- 尚未装入：130对象，包括肩胛浮动爪、后衬垫及M3，GH捕获五金、肱骨、锁骨骨段、肌腹及相关接口。完整逐项清单见[步骤场景清单](selected_stage_scene.json)，不以“合法接触”或整对关闭碰撞掩盖它们。
- 三片前衬垫需预先留在名义位置，实际临时保持、压装／定位和手／工具空间尚未设计验证，不能假定它们粘在夹具上。

同时比较了衬垫随骨临时协同移动和骨无衬垫的情形。随骨带前衬垫的9个受检方向都起步受阻；不能用这种操作替代未完成的衬垫保持设计。“骨无衬垫”即使路径清，也不能证明以后衬垫装得进去。

## 3. 路径与复核结果

取出方向源坐标为[0,−0.9578263,0.2873479]。20 mm路径对应源Y前移约19.157 mm、Z上移约5.747 mm；装入沿反方向从20 mm回到零位。它不是肩胛运动范围、GH角度或任何K指标结果。

[独立细查](selected_path_dense_check.json)在去程和回程各重新变换并计算201姿态，共402姿态，整段步长0.1 mm。所有检查姿态相交体积均不超0.001 mm³；最大残差约0.00014782 mm³。最小间隙记录包含0，末端查询达到2 mm封顶，**没有证明严格零相交、正制造间隙或数学意义的连续无碰撞**。

[当前175对象独立回读](receiver_relief/readback_contact_checks.json)检查15225对零位材料，无超阈值刚性相交；六条主动线对全部当前材料不超阈值。按新支架几何重新核对90对父级名义接触声明，0.05 mm法向探针均遇到材料，六条肩胛M3孔轴仍清。接触声明不是预紧、摩擦、承载或骨夹紧证明。

原三角肌／肩胛柔性重叠、通道壁残差、两自由肌腹末端、材料与驱动未知等均保留。旧动作角度、旧原位维护受阻和其他骨架装配步骤不被本检查覆盖。

## 4. MuJoCo查看

[三姿态实网格图](scapula_insertion_demo.png)显示20／10／0 mm三姿态，蓝色为预置前衬垫。场景仅45对象，来源是175对象清单；无肱骨、锁骨骨段、GH五金或肌腹。骨由mocap平移，nq=0、nmocap=1；没有动力学步进，显示场景关闭碰撞，不用于几何验收，验收来自上述实际实体检查。

从项目根目录、装有MuJoCo和Pillow的Python运行：

```powershell
python mujoco_arm/tools/view_shoulder_scapula_insertion.py
```

空格播放／暂停，0回名义装好姿态，1回20 mm取出姿态。鼠标控制视角。它是单步装配查看，不是驱动演示；实物临时夹持与手部操作尚未验证。

```powershell
python mujoco_arm/tools/view_shoulder_scapula_insertion.py --render
python mujoco_arm/tools/view_shoulder_scapula_insertion.py --check
```

完整175对象静态查看另用：

```powershell
python mujoco_arm/tools/view_shoulder_manufactured.py --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_fixture_assembly_review/receiver_relief/assembly_check.json'
```

[完整175对象静态图](receiver_relief/available_parts.png)不代表全部已装配验证。

## 5. 可复现生成顺序

生成命令会写本审查目录；旧175父级和失败记录保留，先保存需要保留的后续结果。使用装有numpy／Manifold等依赖的环境：

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) '.deps')
python mujoco_arm/tools/relieve_shoulder_fixture_for_scapula_insertion.py
python mujoco_arm/tools/audit_shoulder_scapula_insertion.py --assembly '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_fixture_assembly_review/receiver_relief/assembly_check.json'
python mujoco_arm/tools/verify_shoulder_scapula_insertion.py
python mujoco_arm/tools/check_shoulder_fixture_receiver_revision.py
```

审查脚本默认检查当前receiver_relief修正版；上方显式输入也可复现。若用--assembly检查其他历史候选，会覆盖本审查输出，必须重新细查，不能将粗采样标签沿用。原失败记录及父级哈希是生成依据，若父级变更，应重新审查，不能绕过哈希。

输出使用本项目实际短目录，不建兼容目录或联接。后续先补三前衬垫保持、后衬垫／浮动爪／五金与GH组件完整安装、工具路径，再继续实际肌肉驱动约束和样件标定。所有报告stage_pass=false，整机目标未完成。
