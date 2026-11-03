# 上侧起端安装区候选：当前S1集成入口

2026-11-02。解决原FJ1467起端盖几乎无转动间隙的问题；整肩运动和实物驱动仍未通过。

## 设计依据及真实改变

下侧端盖缩短凸台1.5mm后，局部外展仍约32°受弯管限制；90°侧出口研究仅改善到约46°，未选入本方案。更深埋管碰肩胛骨，不能采用。

在同一原始肩胛骨网格上发现邻近上侧骨面，选择焊接后顶点505，面积加权外法向朝上。这里是工程安装区域候选，不代表已证明生物肌腱的精确附着点。

| 参数 | 原候选 | 当前候选 |
|---|---|---|
| 源坐标安装点mm | −161.26, −63.69, 1344.01 | −162.823, −64.9978, 1348.45 |
| 位移 | — | 4.885377mm |
| 法向 | −0.281906, −0.269910, −0.920694 | −0.278712, 0.088169, 0.956319 |
| TPU端垫 | 旧下侧端垫及3mm凸台 | R11，名义高2.5mm，骨侧轮廓拟合；不加旧凸台 |
| 刚性盖 | R6，厚1mm及R3弯管 | R9.5，厚1.5mm；保留两M3候选孔，移除旧弯管 |
| 钢丝入口 | 下方绕出 | 向上出线；TPU过渡绕向外缘连接原有效肌腹 |

骨骼、骨长、关节中心及有效肌腹核心路径未移动。新连接段长度和实际应变/力矩需要重算，不能继承旧25%假设。

初版过渡外缘参数13.5mm在与原晶格交接处产生一条非流形边。R0.2/R0.5局部填补试验未解决，保留未选。将实际过渡外缘改为14.0mm（外调0.5mm）后，真实几何消除狭缝，网格通过；未放宽边界边/非流形边必须为0的标准。

## 唯一当前集成组合

| 模块 | 报告 |
|---|---|
| 三角肌 | `../unified_deltoid_upper_origin_edge14p0/check.json` |
| 冈下肌 | `../unified_infraspinatus_upper_edge14_master/check.json` |
| 肩胛下肌 | `../unified_subscapularis_upper_edge14_master/check.json` |
| 大圆肌 | `../unified_teres_major_upper_edge14_master/check.json` |
| 两个夹具 | `../fixture_revision/check.json` |
| 钢丝主路径 | `../routes_upper_origin_edge14p0_candidate.json` |

其余三组只在逐项确认自己的路径记录完全相等、网格哈希未变后重绑定主路径版本，不继承改变邻件之后的整装结论。完整22件及刚性运动均重新计算。

## 检查结果与限制

- 14个肌肉/端盖网格通过水密、正材料单连通和回读检查；两个夹具也通过。
- 三角肌TPU为1个正材料分量，14个封闭空腔单独记录，没有为了通过而填腔。制造支持仍未验收。
- 六条完整钢丝探针对22个现有实体相交均为0mm³。新端盖/TPU、夹具/肌肉相交也为0。
- 无超过0.05mm³的实体交叠，最大未分类贴骨残差约0.002936mm³；不宣称全模型数学零干涉。
- 当前自由TPU壁厚守护区仍有前/后/中束约0.00804/0.00325/0.00664mm³残差，未通过壁厚验收。

刚性运动以固定GH中心、固定SC/AC、端盖及止挡理想安装为条件；不把TPU当刚体或把静态钢丝强行固定为运动约束。角度步长1°，首次限制细化0.1°，超阈标准0.001mm³。

| 校准轴及符号 | 末个未超阈采样° | 首个超阈采样° | 首限 |
|---|---:|---:|---|
| flexion + | 24.8 | 24.9 | 肩胛下肌远端盖FJ1504—肩胛骨 |
| flexion − | −51.6 | −51.7 | 肱骨—肩胛骨 |
| abduction + | 54.1 | 54.2 | 冈下肌远端盖FJ1500—肩胛骨 |
| abduction − | −25.0 | −25.1 | 肱骨—肩胛骨 |
| longitudinal + | 80.3 | 80.4 | 冈下肌远端盖FJ1500—肩胛骨 |
| longitudinal − | −38.4 | −38.5 | 肩胛下肌远端盖FJ1504—肩胛骨 |

各末个采样的最小表面间隙为0，表示已接触或存在未超过报告阈值的交叠，不能当作有制造间隙。JSON同时保存首限体积、配对和源轴向量；符号不直接替代已核验的解剖内/外旋标签。

这只是当前理想刚性模型的连续采样边界，**不是TPU实际驱动范围，也不是K1–K7通过**。原指标保留比较，下一步优先缩减/重布置FJ1500、FJ1504远端刚性接口占位，不能通过关闭骨—盖检查通过。

实际骨侧紧固、夹持及压接、导管反力、安装工具路径、防脱、TPU变形接触、摩擦、材料行程、强度和打印工艺均未验收。新骨面上的孔位只是候选，不代表骨骼已开孔或紧固件已设计。

## 复现

项目环境需将`.deps`及`mujoco_arm/tools`加入Python搜索路径，在仓库根目录运行：

1. `search_shoulder_origin_stack.py`：保留凸台缩短比较与裸骨采样基准。
2. `design_shoulder_upper_origin.py --edge-offset 14.0`：生成骨面候选及主路径。
3. `build_shoulder_unified_muscle.py --routes <routes_upper_origin_edge14p0_candidate.json绝对路径>`。
4. `rebind_shoulder_unchanged_routes.py --routes <上述路径> --tag upper_edge14_master`，重复传入三份旧模块`--report`：冈下肌/肩胛下肌fitted、大圆肌connected。模块路径不同即拒绝复用。
5. `check_shoulder_manufactured_assembly.py`传入上表4份`--unified-report`、`--fixture-report`及`--output-dir <本目录绝对路径>`。
6. `check_shoulder_rigid_motion.py --assembly <本目录actual_assembly_check.json绝对路径>`。

静态查看：

```powershell
& 'C:/Python314/python.exe' 'mujoco_arm/tools/view_shoulder_manufactured.py' --report '02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_port_path_revision/whole_shoulder_upper_edge14/actual_assembly_check.json'
```

加`--render`生成下图。DOF为0，仅查看实际网格；运动结果来自逐姿态实体报告。

![当前22件真实网格](available_parts.png)
