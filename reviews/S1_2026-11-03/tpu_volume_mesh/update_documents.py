"""Record full-muscle tetrahedralization failures and clarify import precision."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'
def edit(p,pairs):
    s=p.read_text(encoding='utf8')
    for a,b in pairs:
        if a not in s:raise ValueError(f'Missing expected text: {p}: {a[:60]}')
        s=s.replace(a,b,1)
    p.write_text(s,encoding='utf8')

section='''### 2.19 完整TPU体积网格前置检查：未通过

[体网格尝试与失败定位](s1_y2_transition_stock/volume_mesh/README.md)将当前两块完整肌腹输入Gmsh 4.15.2，完成6组表面保留／分析副本／边界恢复算法试验；均未获得通过数值质量门槛的完整肌腹四面体。对照立方体有12个正体积四面体、总体积1000 mm³，独立验证仅支持API与数据提取流程，不提升肌腹结果。

此前回读默认使用四位小数坐标分桶合并，保留首个真实坐标；其零边界／非流形边结果限定于该导入约定。本轮原始STL精确去重：冈下肌非流形边0，肩胛下肌6；后者导入最大重指派约0.00006175 mm、删除8个退化面后边计数为0。即使此时为单正材料连通且边计数正常，仍不足以证明PLC表面不相交或适用于体网格。

标准边界恢复报告近重合面；调整数值恢复参数及独立HXT试验仍出现线段与三角面相交。已保存原生输出坐标、实际开放局部面片和[独立交点验证](s1_y2_transition_stock/volume_mesh/verification_checks.json)。本次定位在两肌腹远端接口附近，不在四处新增填实空腔包围盒内；定位数量随首先报告的位置变化，只是局部示例，不能当全模型自交总数。

分析副本0.005 mm简化没有替换制造网格；体积、表面及导入偏差另列。0.001 mm³实体碰撞判据、四面体正体积和数值质量门槛均未据失败放宽。网格构建失败退出码2，诊断脚本成功定位不等于肌腹通过。

当前源STL、骨骼、锚点与候选指针均不改。先修复被定位的远端过渡表面，重新核对接口／通道、表面与完整体网格；通过后再进入连续体离散与自由端盖／肌腱接触求解。材料本构、实际打印、预紧与摩擦仍未测，本轮没有产生TPU应力、位移或驱动范围。`stage_pass=false`。

'''
edit(D/'SHOULDER_DESIGN.md',[
 ('版本：v0.9.15 · 2026-11-02','版本：v0.9.16 · 2026-11-03'),
 ('## 3. 肌肉／近端支承父级已验证内容与继承边界',section+'## 3. 肌肉／近端支承父级已验证内容与继承边界'),
 ('各源均为单正材料连通、零边界及非流形边。','各源按现行分桶导入检查为单正材料连通、零边界及非流形边；原始STL精确拓扑及体网格适用性另见2.19节。'),
 ('| 4 | 完善已建立的自由末端／肌腱端点与张力单边约束，接入TPU平衡与变姿态路由 |','| 4 | 先修复完整肌腹体网格输入，再将已建立的端点／张力约束接入TPU平衡与变姿态路由 |'),
])
with (D/'SHOULDER_DESIGN.md').open('a',encoding='utf8') as f:f.write('\n2026-11-03 v0.9.16：新增完整肌腹体网格前置门、6组实际失败、对照API与局部边—面交点独立诊断；澄清四位小数分桶导入的拓扑结论边界。未修改源几何或输出虚构变形。[更新前快照](../../../reviews/S1_2026-11-03/tpu_volume_mesh/before_SHA256.json)。\n')
edit(ROOT/'PROJECT_STATUS.md',[
 ('更新：2026-11-02 · 技术文件v0.9.15。','更新：2026-11-03 · 技术文件v0.9.16。'),
 ('[肩部技术文件v0.9.15]','[肩部技术文件v0.9.16]'),
 ('## 最新进展与限制\n','## 最新进展与限制\n\n新增[完整TPU体网格前置检查](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/volume_mesh/README.md)：6组完整肌腹尝试均未通过，已定位远端过渡附近的边—三角面交点。立方体API对照通过，仅验证工具流程；当前无可用完整TPU四面体或变形解。源几何与指针不变，先修复分析边界。此前网格零边／非流形结论采用四位小数分桶导入，不是原始STL精确拓扑或体网格通过证明。\n'),
 ('## 当前优先问题\n','## 当前优先问题\n\n**新增前置门：完整肌腹分析表面尚未通过体网格恢复。** 需对已定位的远端夹持／挂钩过渡做局部修复并复核接口／通道；保持既有失败记录。此门未通过前，不计算或展示TPU应力、位移与实际驱动结果。\n'),
 ('先核对S1-Y2骨架紧固与两块Y肌腹的切片／穿线可行性；','先修复两块完整肌腹的体网格输入并复核接口／通道，再核对S1-Y2骨架紧固与切片／穿线可行性；'),
])
edit(ROOT/'design/CURRENT_BASELINE.md',[
 ('更新：2026-11-02','更新：2026-11-03'),
 ('技术文件2.14–2.18节','技术文件2.14–2.19节'),
 ('## 文档分工','当前通用STL导入采用四位小数坐标分桶合并；零边界／非流形边只对应该导入约定，不保证原始精确拓扑或体网格适用性。完整肌腹6次体网格恢复仍失败，须修复已定位远端过渡并复验后再进行连续体求解。数值网格预检与材料标定分别验收，立方体对照不作肌腹通过。\n\n## 文档分工'),
])
brief=ROOT/'design/DESIGN_BRIEF.md';raw=brief.read_bytes();marker='## 历史记录阅读说明'.encode('utf8');head,tail=raw.split(marker,1)
s=head.decode('utf8').replace('## 当前入口（2026-11-02）','## 当前入口（2026-11-03）').replace('技术文件v0.9.15','技术文件v0.9.16')
s+='[完整肌腹体网格前置检查](../02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/volume_mesh/README.md)未通过，已定位远端过渡边—面交点；此前回读拓扑采用四位小数分桶导入。先修复并验证完整分析边界，当前源STL及指针不改，没有TPU变形或驱动通过。\n\n'
brief.write_bytes(s.encode('utf8')+marker+tail)
for p in [ROOT/'README.md',ROOT/'02_人体骨骼晶格肌肉方案/README.md']:
    s=p.read_text(encoding='utf8').replace('当前状态（2026-11-02）','当前状态（2026-11-03）').replace('技术文件2.14–2.18节','技术文件2.14–2.19节')
    if p.parent.name=='02_人体骨骼晶格肌肉方案':s=s.replace('v0.9.15','v0.9.16')
    p.write_text(s,encoding='utf8')
edit(D/'s1_y2_transition_stock/README.md',[
 ('## 当前数字证据','新增[完整TPU体网格前置检查](volume_mesh/README.md)：两肌腹均未通过，已定位远端过渡局部边—面交点；此前回读采用四位小数分桶导入。当前没有完整连续体变形解，源几何未修改，立方体对照不提升肌腹验收。\n\n## 当前数字证据'),
])
print('Technical v0.9.16 records failed complete-muscle meshing and the import-precision boundary.')
