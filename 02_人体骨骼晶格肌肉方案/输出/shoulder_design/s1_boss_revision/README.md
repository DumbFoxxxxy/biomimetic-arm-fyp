# S1 端垫凸台连接修订

日期：2026-11-02。

已修复原10个端垫中凸台与主体相隔1.5mm的问题。根因是将Manifold默认从0到h的圆柱当成中心对称圆柱，额外加了h/2位移。
凸台现从主体外表面下方0.01mm开始构造并与主体取并，顶面仍为主体高度+名义3mm。
该嵌入用于避免旋转后共面浮点误差造成分离，不增加名义外部高度。

## 验证与交付

- 平端垫、减骨随形端垫、任意旋转端垫3项回归测试均通过。
- 20个实际导出STL（10端垫+10端盖）均为单一连通材料，水密、边界边/非流形边为0，回读体积误差不超过1e-5。
- 已逐文件核对哈希后复制到生产shoulder_parts目录，记录见deployment.json。
- 修复前20件STL及2份生产报告保存在before/，不与修复后件混用。
- 当前s1_parts_export.json增加内存实体与STL回读连通件数量，并在检查失败时返回非零状态。

## 零位子装配检查

check_shoulder_assembly.py新增--manufactured-parts，读取真实STL，不以参数重建代替制造件。
生产文件复核：27件（3骨、4肌肉包络、10端垫、10端盖），20处超过0.05mm3的接触，18处在已有白名单内。
其余2处：肱骨/三角肌中束起点端盖0.102mm3；冈下肌远端盖/三角肌包络0.242mm3。
均小于原设计采用的0.5mm3预算，但仍保留记录，未称为无干涉或可直接制造。
旧0.484mm3报告使用参数重建，本次读实际STL；二者不作为相同对象的受控前后对比。

范围不含腱段、夹具、五金、钢丝、软衬及运动过程；白名单内的夹持接触也尚未通过实物装配确认。
本修订只解决端垫断连及其检查缺口，不代表S1驱动、完整装配或强度已通过。

复核命令（需numpy、manifold3d；本机.deps可提供manifold3d）：

```powershell
python -m unittest discover -s tests -p test_shoulder_pad_connections.py -v
python mujoco_arm/tools/check_shoulder_assembly.py --manufactured-parts "02_人体骨骼晶格肌肉方案/输出/anatomical_concept/shoulder_parts"
```

原始审查报告保留修复前事实；本文件与当前生产JSON记录修复后的状态。
