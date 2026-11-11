# 肘部供线端反力接口候选

2026-11-08。以208个刚性件的 `assembly.json` 为当前局部候选清单，保留原骨骼、原肘肌腹外形及内部钢丝走线。主骨架指针未切换；本目录不是已完成的胸廓供线方案。最终装配无工程背板。

## 增加了什么

原端盖有半径0.8mm的钢丝通孔，但没有外导管的专用承压／保持接口。在原端盖上增加阶梯座和向空余侧伸出的两固定耳，另配金属导管端头、可拆压板及两套名义M2x14紧固件。每块端盖新增约1050.56mm³，原端盖材料仅有数值噪声级求差，没有移动原安装点或削掉原连接。

| 项目 | 当前候选尺寸及作用 |
|---|---|
| 钢丝 | 直径0.6mm候选；穿过导管、端头及原端盖，不固定在近端端盖 |
| 导管参考截面 | 外半径1.5mm、孔半径0.6mm；外导管材料、轴向及弯曲刚度未定 |
| 金属端头 | 主体外半径2mm、孔半径0.65mm；上段容纳导管的孔半径1.55mm、有效插入4mm |
| 端头承压座 | 半径2.1mm容纳主体；下方环形肩面有效承压面积约10.2966mm² |
| 防退出结构 | 端头法兰半径3mm，压板孔半径1.7mm，最大轴向浮动0.2mm |
| 双固定耳 | 位于原支架空余侧，侧插M2螺母，4.2mm宽槽约束4mm螺母对边；螺母直接承压，不埋松散垫片 |
| 尚未确定 | 导管与端头的压接／粘接防脱、真实螺纹、公差、预紧、强度、寿命及实际工具 |

端头法兰的机械保持已设计；导管自身的防拔出仍待验证，不能把二者混为一项。图中绿管只是名义短导管参考，不证明其材料能承受预测反力。

名义压缩反力传递：导管端面 → 金属端头内台阶 → 端头底面 → 一体端盖承压座 → 原端盖M3连接 → 原骨侧支架 → 肱骨安装座。外导管远端和电机侧还未连接，钢丝远端终端本轮未验证。

## 真实供线边界及坐标

[route_endpoints.json](route_endpoints.json)是后续走线应读取的端口角色记录。坐标为BP3D源坐标，单位mm：

| 路径 | 导管实际承压台阶中心 X/Y/Z | 向外单位轴线 |
|---|---|---|
| FJ1486，肘屈肌 | −185.71493 / −97.82312 / 1224.77969 | 0.172585 / −0.014178 / 0.984893 |
| FJ1477，肘伸肌 | −179.07486 / −45.22341 / 1294.07952 | 同上 |

端头出孔位于承压台阶沿轴线向外4mm处。压缩端面支承存在，但以这4mm插入段作为理想切向夹持仍需保持／接触验证。其余外导管节点必须由自身弹性、实际导向或骨面接触决定运动。

`assembly.interfaces.supply_exit_source_mm`沿用生成器字段，实际含义是名义绿管末端；它不是新的刚性支承或固定路径点。`route_endpoints.json`已明确将它登记为`reference_only_stub_tip`且`reference_tip_fixed=false`。当前模型没有把这个末端固定为钢丝导向点。

内部等效梁仍把原首个TPU导向节点夹持为近端边界，尚未解析完整晶格—端盖过渡。新增端头没有把这个近似自动升级为外供线受力验证。

## 几何及安装证据

- `fit_screen.json`：新增／修改部件配对检查、两条进线及两个名义导管短段均无>0.001mm³相交；新增材料各为单一连通体。
- `export_screen.json`：208个刚体精确重载后复核，配合仍清；18个新增／修改网格的导出体积误差均为0，无自动删除孤岛。
- `installation_screen.json`：上方骨紧固件全装好时，屈肌压板直线下降失败；绕导管轴旋转的54个分段候选也失败，均保留。
- `installation_late_bone0.json`：屈肌上方bone0螺栓、前后垫片及螺母四件明确后装。先装端头、侧插M2螺母、压板及M2，再将原bone0前组沿原轴45→0mm、后组35→0mm装回。两套接口共18阶段、562个位置，0.5mm采样均无超阈值相交；没有排除整对邻接部件。

压板须在导管端头固定前先穿到导管上并临时停放；停放姿态及进入中的端头／参考短管也做了配对检查。另一骨螺栓及手动持件保持原支架姿态，临时夹持工装、螺纹旋入、扭矩工具、完整长导管及原整块肌腹装入没有验证。因此这只是局部工序几何候选，不是整机可装配结论。

旧对称固定耳碰原骨侧支架，早期复杂阶梯螺母孔导出未通过。现版改为单一矩形侧插槽，严格导出复核通过；失败记录保存在`fit_history_*`、`export_failure_diagnostic.json`及诊断OBJ，均不属于选定装配。

## 接入局部驱动模型后的结果

`model.xml`为294速度自由度、215可视几何。208刚体包含真实新端头／压板／紧固件，新增质量按打印件1200、钢件7800kg/m³估计；肱骨刚体质量约0.61513kg，基础惯量仍是假设。两绿管网格为零质量、无接触响应的参考显示，不是动态导管或刚性路径约束。模型仍缺少胸廓电机、外导管及其反力。

`local_cycle.json`：2s自由受力，屈肌卷线2mm再放回，其余26路保持名义线长；肘角度未给定、软件肘角度限位关闭，实体止挡采用已记录的0.005mm数值接触缓冲与1e6软接触参数。

- 预测肘范围0..1.62235°，最终0.02864°；无警告、无60N饱和、无记录接触压入。
- 屈／伸峰值张力约18.0864／4.6735N；它们是内部无摩擦局部驱动结果，不是已计算的外导管反力。若仅以18.0864N估算端头肩面平均压应力，约1.76MPa，尚无强度通过结论。
- 最大端部平移等式分量残差0.09298mm；1012保存姿态，最大旋转采样增量0.52042°。
- `local_cycle_solid_screen.json`：两钢丝对208刚体、16364跨体刚体对的1012姿态均无>0.001mm³相交。完整TPU、外导管和孔腔变形未检查，不是数学连续运动证明。

当前约1.62°只是这套单路卷线、其它线长保持的预测；两条肩胛—前臂长头仍为直线代理，也在牵制肘。不能将它作为骨架／肌肉能力上限，原K指标不变。

## 复现及整臂查看

在项目根目录PowerShell中使用当前环境：

```powershell
$env:PYTHONPATH='tools;.deps'
D:/MujocoEnv/Scripts/python.exe tools/build_channel_first_elbow_supply_interface.py
D:/MujocoEnv/Scripts/python.exe tools/describe_channel_first_elbow_ports.py
D:/MujocoEnv/Scripts/python.exe tools/build_channel_first_elbow_supply_model.py
D:/MujocoEnv/Scripts/python.exe tools/audit_channel_first_elbow_supply_installation.py --late-flexor-bone0
D:/MujocoEnv/Scripts/python.exe tools/simulate_channel_first_elbow_channel_beams.py --folder '02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_supply_interface' --cycle --duration 2 --contact-stiffness-scale .01 --contact-margin-mm .005 --result-name local_cycle
D:/MujocoEnv/Scripts/python.exe tools/audit_channel_first_elbow_channel_beams.py --folder '02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_supply_interface' --result-name local_cycle
D:/MujocoEnv/Scripts/python.exe tools/view_channel_first_elbow_cycle.py --folder '02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_supply_interface' --result-name local_cycle
```

回放支持鼠标旋转／平移／缩放，显示整胸廓和手臂；两原TPU肌腹没有作为已验证变形显示。`--check`已重建1012姿态，GUI本轮未启动。单独启动XML不会施加Python中的肩胛支承和卷线控制，不能重现这次动力学。

`render_channel_first_elbow_supply.py`可重新生成下图。橙色为新一体端盖与压板，金属件灰色、名义导管短段绿色；完整肌腹未显示。

![伸肌供线接口零位](FJ1477_inlet_55.png)

后续优先确定胸廓／桌面底座上的驱动与反力连接，设计实际外供线段并让它自由弯曲、接触和传力，再协调跨肩肘的各路收放线。此局部候选不能替代这些全局工作。
