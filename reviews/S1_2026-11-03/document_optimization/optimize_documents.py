"""Reorganize documentation only; preserve selected geometry and old technical text."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / '02_人体骨骼晶格肌肉方案/输出/shoulder_design'
HERE = Path(__file__).resolve().parent
old = (HERE / 'before/02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md').read_text(encoding='utf-8')

def section(start, end):
    return old.split(start, 1)[1].split(end, 1)[0].strip()

def paragraph(text):
    return text.strip() + '\n\n'

parts = [paragraph('''
# 肩部 S1 技术方案与验证依据

版本：v0.9.17 · 2026-11-03 · 技术审查稿，未定版。

**当前采用S1-Y2作为研究起点：真实解剖骨形、GH工程球铰、钢丝与异形TPU晶格肌肉，共191个名义对象。** 已有零位配合、有限裸骨运动和局部安装证据，尚不能确认可制造、可靠装配或稳定驱动。预测及实测驱动范围均未知，`stage_pass=false`。

本次只重组技术文件，不修改几何、锚点、钢丝、判据或候选选择。原稿逐字节保存在[v0.9.16完整归档](SHOULDER_DESIGN_v0.9.16_ARCHIVE.md)，保留全部迭代参数、失败与版本记录。

| 阅读目的 | 位置 |
|---|---|
| 了解目标、当前结构及输入 | 第1–2节 |
| 理解晶格、钢丝与传力 | 第3–4节 |
| 判断验证结果与动作范围 | 第5–6节 |
| 判断可行性、后续验收与复现 | 第7–9节 |

## 1. 目标与研究边界
''')]
parts.append(paragraph(section('## 1. 目标、范围与当前决策', '## 2. 最新候选：结构与实际传力关系')))
parts.append(paragraph('''
“几何结果”仅对应指定输入与采样；“预测”需要材料、约束与模型；“实测”需要实物记录。对象数包含五金、钢丝和台面，不是打印件数或完成度。

## 2. 当前骨架、支承与装配

### 2.1 唯一输入及坐标

| 文件 | 用途 |
|---|---|
| [候选选择记录](../../../design/S1_CURRENT_CANDIDATE.json) | 当前报告路径与SHA256；默认查看核对它 |
| [S1-Y2说明](s1_y2_transition_stock/README.md) | 当前模型、几何修改及生成顺序 |
| [191对象源清单](s1_y2_transition_stock/assembly_check.json) | 网格、哈希、装配变换、归属及配对结果 |
| [钢丝／锚点清单](s1_y2_transition_stock/wire_and_anchor_plan.json) | 当前线位、端点及上游路由引用 |

单位为mm。局部STL顶点加源清单中各自的`translation_source_mm`恢复源解剖坐标，不能把不同局部坐标直接混用。研究分支若声明不同导入约定，须读取其记录。

可编辑参数位于生成脚本和JSON；目前没有覆盖全部对象的新SolidWorks原生特征树及装配配合交付。重新导出须复核文件、变换和适用检查，不能只替换单个STL或修改指针哈希。

### 2.2 GH骨架及台架连接

GH采用机械保持球铰，球头名义半径16 mm、窝内半径16.4 mm；保留原始骨形对照，工程肩胛有局部差分修订。骨干抱箍位于相对GH定义的60 mm站位。名义间隙不是实物配合公差，有限捕获不是承载或任意姿态防脱证明。

紧凑GH骨架将M2头／螺母埋入新增凸台、取消四片外置垫圈，第0侧横销内收0.3 mm，并对新增凸台让位约4.68953 mm³。骨长、GH中心与肱骨未因此移动。实际M2头型、头下9 mm长度、批头、预紧及座壁／横销承载尚未验证；详细材料差分见[骨架父级记录](s1_GH_assembly_review/compact_recessed_keepers/README.md)。

两块肩袖肌经近端支承、穿骨紧固件及背板连到肩胛骨，骨外夹爪与锁骨抱夹再接公共台架。主体背板225×148×6 mm；18 mm厚局部台面片代表外部支承假设，尚未选定实际工作台。

### 2.3 安装步骤的适用范围

已有11步名义放置、裸肱骨0↔10°回零及有限防转／侧退证据。先放部分螺母与保持五金，再装裸肱骨、保持盖及螺钉；**步骤不含完整肌腹、抱箍、肌腱安装或实际拧紧**。后衬垫、浮动夹爪、工具进入及全套顺序仍待整合。

历史台架步骤发现完整场景内M4抽出受肩胛阻挡，安装顺序不能只按孔路推定。零间隙接触探针只说明存在材料，不证明夹紧力、抗滑或薄骨承压安全。局部装入、紧固和保持的完整数值见[旧稿第2.5–2.14节](SHOULDER_DESIGN_v0.9.16_ARCHIVE.md)。

## 3. 异形TPU肌肉

### 3.1 布局与附着
'''))
layout = section('### 2.2 骨骼、关节与肌肉布局', '### 2.3 两块肩袖肌的名义力链')
layout = layout[layout.index('| 模块 |'):]
parts.append(paragraph(layout))
parts.append(paragraph('### 3.2 共形晶格与连续通道'))
y = section('### 2.15 S1-Y1父级：两块共形Y晶格肌腹', '### 2.16 当前S1-Y2：封闭过渡空隙的明确填实')
# Keep detailed actual parameters and wall-envelope qualifications; move historical trials to archive.
y = y[y.index('沿内部通道中心线'):y.index('骨骼未改，但针对Y1清单')]
y = y.replace('TPU体积 |', 'Y1父级TPU体积 |').replace('另见2.19节', '另见第5节')
parts.append(paragraph(y))
parts.append(paragraph('### 3.3 当前Y2过渡区修订'))
y2 = section('### 2.16 当前S1-Y2：封闭过渡空隙的明确填实', '### 2.17 自由末端／肌腱约束与局部端环转动')
y2 = y2[:y2.index('[新裸骨路径]')]
parts.append(paragraph(y2))
parts.append(paragraph('### 3.4 论文方法与标定'))
paper = section('## 5. 论文方法、晶格方向与标定', '## 6. 目标／底线／最新实际范围')
paper = paper.replace('见2.15节', '见3.2节')
parts.append(paragraph(paper))
parts.append(paragraph('## 4. 钢丝、短肌腱与传力模型\n\n### 4.1 名义力链及反力'))
parts.append(paragraph(section('### 2.3 两块肩袖肌的名义力链', '### 2.4 关键名义参数')))
parts.append(paragraph('### 4.2 近端接口与名义长度'))
params = section('### 2.4 关键名义参数', '### 2.5 支架到台面的名义连接')
parts.append(paragraph(params.replace('见2.17节', '见4.3节')))
parts.append(paragraph('### 4.3 自由端约束与局部端环检查'))
ends = section('### 2.17 自由末端／肌腱约束与局部端环转动', '### 2.18 规定端环角的完整肌腱路由候选')
parts.append(paragraph(ends))
parts.append(paragraph('### 4.4 完整肌腱候选与斜向出线'))
routes = section('### 2.18 规定端环角的完整肌腱路由候选', '### 2.19 完整TPU体积网格前置检查：未通过')
parts.append(paragraph(routes.replace('补上2.17节', '补上4.3节')))
parts.append(paragraph('''
## 5. 验证证据与阻塞项

### 5.1 当前S1-Y2的检查范围

| 检查及输入 | 已有结果 | 不能据此推出 |
|---|---|---|
| [191对象零位](s1_y2_transition_stock/assembly_check.json) | 18145配对；六主动线无超0.001 mm³判据刚性相交 | 全姿态、严格零穿透、实际穿线或驱动 |
| [102对象裸骨](s1_y2_transition_stock/bare_motion_review/bare_joint_motion_checks.json) | 101固定＋单肱骨；连续往返和首限见第6节 | 89对象未装，不含肌腹／抱箍／肌腱约束 |
| [保持与名义放置](s1_y2_transition_stock/readback_keeper_checks.json) | 11步0.1 mm往返共4322姿态；裸肱骨0↔10°共22姿态；52方向球头有限捕获 | 实际拧紧、工具、强度及任意联合逃逸 |
| [两肌腹独立复核](s1_y_muscle_pair/transition_stock_comparison/independent_checks.json) | 覆盖、保护区、有限接触、177刚体与六主动线 | 实物刚度、受压通道、工艺与肌肉性能 |
| 端环／整根肌腱候选 | 条件几何和端点约束已建立，见第4节 | 自由末端—肌腱—TPU耦合平衡 |
| 完整肌腹体网格 | 原S1-Y2六组尝试均失败，见下节 | 尚无合格完整TPU四面体或变形解 |

0.001 mm³是相交报告判据，不是制造公差。接触探针、单件孔路和有限捕获均不能替代真实装配、夹紧或承载验证；部分名义路径末端间隙为0。

### 5.2 完整肌腹分析边界：未通过
'''))
parts.append(paragraph(section('### 2.19 完整TPU体积网格前置检查：未通过', '## 3. 肌肉／近端支承父级已验证内容与继承边界')))
parts.append(paragraph('''
### 5.3 修复实验分支：未采用

另有[局部增料修复试验](s1_tpu_distal_surface_revision/junction_repair_trials.json)与[全表面修复试验](s1_tpu_surface_global_revision/global_repair_trials.json)，均记录`adopted_as_current=false`，不替代S1-Y2，也不新增通过结论。导出及回读仍报告交点，内部某轮扫描为0不能当作最终文件通过；不同扫描和导入方法的数量不直接比较。

分支须独立核对来源、导入精度、材料差分、接口及通道，并对最终回读文件重新验证。体网格构建后还需数值质量和收敛检查、材料标定及耦合求解。本次文档整理不重跑或验收这些实验。

### 5.4 其他必须保留的问题
'''))
unresolved = section('**未解决项必须保留：**', '## 4. 候选谱系与证据适用范围')
parts.append(paragraph(unresolved))
parts.append(paragraph('## 6. 动作目标与当前范围\n\n### 6.1 裸骨几何结果'))
bone = section('### 2.14 S1-Y1沿用的紧凑191对象骨架父级及裸骨范围', '### 2.15 S1-Y1父级：两块共形Y晶格肌腹')
bone = bone[bone.index('| 裸骨探针 |'):bone.index('11步名义放置')]
parts.append(paragraph('以下引用当前[S1-Y2裸骨复核](s1_y2_transition_stock/bare_motion_review/bare_joint_motion_checks.json)，不是直接继承父级结果。只移动裸肱骨，101对象固定、89对象未装。'))
parts.append(paragraph(bone))
parts.append(paragraph('### 6.2 原指标与驱动状态'))
metrics = section('## 6. 目标／底线／最新实际范围', '## 7. 下一阶段工作与完成判据')
metrics = metrics[metrics.index('| 编号／动作 |'):]
parts.append(paragraph(metrics))
parts.append(paragraph('''
## 7. 可行性判断与后续验收

**方案已有可继续研究的几何基础，尚不足以确认稳定、可重复的三轴肩部驱动。** 数字分析先处理完整肌腹边界；实物准备可并行核对紧固与工艺。体网格失败不等于已经证明不能打印，表面闭合也不等于已经证明可打印。

| 工作 | 交付及通过条件 | 依赖 |
|---|---|---|
| 修复完整分析输入 | 独立候选、材料差分、接口／通道复核、最终导出回读、体网格质量与收敛依据；失败保留 | 连续体求解前完成 |
| 落实GH与骨骼固定 | 实际紧固件／工具、进入与拧紧顺序、座壁／横销承载、预紧和抗滑依据 | 带载实物演示前完成 |
| 收敛完整安装 | 肌腹、抱箍、肌腱、其余模块及工具均在场的安装过程与复核 | 局部步骤不能拼作整机通过 |
| 确定打印及穿线工艺 | 设备、喷嘴、材料、方向、切片、通道与端部固定、可重复穿线记录 | 可与数字边界修复并行 |
| 单肌肉及端部标定 | 输入位移、缩短、张力、恢复、滞后、卡线／损伤、重复性；安装长度、预紧、压接保持 | 不填假定收缩比例 |
| 求解自由末端耦合 | TPU本构、肌腱／接触模型、边界条件、力矩、行程及实测对照 | 合格离散与标定数据 |
| 验证拮抗驱动 | 单轴、往返和联合路径的范围、首限、间隙与失败，再评估扩大布局 | 实际驱动试验 |

[单肌肉模板](s1_y2_transition_stock/single_muscle_observation_template.json)与[端部测量模板](s1_y2_transition_stock/free_terminal_kinematics/measurement_template.json)均未测。未取得数据时只记录假设或预测，不关闭性能验收项。

## 8. 查看、复现与维护

从根目录、已装MuJoCo／Pillow的Python环境运行，仅读取现有模型：

```powershell
python mujoco_arm/tools/view_shoulder_manufactured.py
python mujoco_arm/tools/view_shoulder_bare_GH_motion.py
```

第一条是191对象零位静态模型，第二条是102对象裸骨规定运动；均非TPU变形或驱动仿真。裸骨使用单mocap骨、DOF=0，无动力学步进。环境与操作见[RUNNING](../../../RUNNING.md)。

生成顺序见[S1-Y2复现](s1_y2_transition_stock/README.md)。重新生成会写入候选目录，应先保留需要追溯的结果和实测记录；源清单、变换、钢丝与检查报告须一致，不能只换图片或改哈希绕过复验。

仅核对文档及源一致性可运行：

```powershell
python reviews/S1_2026-11-02/Y_muscle_pair/verify_documents.py
```

该检查核对链接、源哈希、原K指标与任务书历史字节，不重新执行几何或实物验证。

## 9. 历史与文档职责

| 入口 | 维护内容 |
|---|---|
| [PROJECT_STATUS](../../../PROJECT_STATUS.md) | 当前进度、阻塞项与下一步 |
| [CURRENT_BASELINE](../../../design/CURRENT_BASELINE.md) | 有效要求及验收边界 |
| 本文 | 当前结构、传力、参数及证据含义 |
| [DESIGN_BRIEF](../../../design/DESIGN_BRIEF.md) | 当前入口和原样保留的历史任务书 |
| [v0.9.16归档](SHOULDER_DESIGN_v0.9.16_ARCHIVE.md) | 全部谱系、失败、材料修改及版本记录 |
| [本次修改前快照](../../../reviews/S1_2026-11-03/document_optimization/before_SHA256.json) | 文档与指针的修改前哈希 |

历史谱系包括24／34／46／72／92／122／155／175／187／195及191对象版、Y1与Y2。195版裸骨屈曲首限96.2°、首次191试改99.1°只适用于对应输入；旧刚性角度与力矩代理不转移到自由末端结构。历史步骤及失败完整保留，未删除连接部位或降低指标获得通过。

v0.9.17：按当前结构、传力、证据和验收重组正文，历史全文归档；同步入口，明确修复实验未采用，区分体网格、可打印性与驱动验证。没有新增几何或性能通过结论。
'''))
new = ''.join(parts)
# Replace residual old-section references in retained technical paragraphs.
new = new.replace('见2.14节', '见第2节').replace('另见2.19节', '另见5.2节')
new = new.replace('本轮原始STL精确去重', '该次原始STL精确去重')
new = new.replace('不足以证明PLC表面不相交', '不足以证明分片线性边界（PLC）不相交')
new = new.replace('本轮没有产生TPU应力、位移或驱动范围', '截至所列证据，没有TPU应力、位移或驱动范围结果')
(D / 'SHOULDER_DESIGN.md').write_text(new, encoding='utf-8', newline='\n')

for rel in ['README.md', '02_人体骨骼晶格肌肉方案/README.md', 'design/CURRENT_BASELINE.md', 'design/DESIGN_BRIEF.md']:
    p = ROOT / rel
    raw = (HERE / 'before' / rel).read_bytes()
    marker = '## 历史记录阅读说明'.encode('utf-8')
    if rel.endswith('DESIGN_BRIEF.md'):
        active, history = raw.split(marker, 1)
        active = active.decode('utf-8').replace('v0.9.16', 'v0.9.17')
        active = active.replace('下一步核对真实紧固件／工具与座壁强度', '数字分析先修复完整肌腹边界；并行核对真实紧固件／工具与座壁强度')
        p.write_bytes(active.encode('utf-8') + marker + history)
        continue
    text = raw.decode('utf-8').replace('\r\n', '\n').replace('v0.9.16', 'v0.9.17')
    text = text.replace('技术文件2.14–2.19节', '技术文件第2–6节')
    p.write_text(text, encoding='utf-8', newline='\n')

print('Technical document:', len(old.splitlines()), '->', len(new.splitlines()), 'lines')

status = '''# 项目当前状态

更新：2026-11-03 · 技术文件v0.9.17。

**当前推进肩部S1，采用S1-Y2作为研究起点，尚未定版。** 路线为真实解剖外形＋工程关节＋钢丝／异形TPU晶格肌肉，允许短肌腱。S1以固定肩胛／锁骨、手拉钢丝为研究条件；肘部既有候选独立保留。

预测与实测驱动范围均未知，`stage_pass=false`。本次仅优化技术文件，没有更改模型、锚点、钢丝或验收指标。

## 当前采用什么

| 入口 | 用途 |
|---|---|
| [S1-Y2说明](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/README.md) | 紧凑GH骨架＋两块共形Y肌腹，共191名义对象；非制造冻结 |
| [候选指针](design/S1_CURRENT_CANDIDATE.json) | 当前源清单、裸骨及保持报告路径与哈希 |
| [肩部技术文件v0.9.17](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md) | 结构、参数、传力、验证证据及后续验收 |
| [有效规格](design/CURRENT_BASELINE.md) | 原K1–K7目标／底线与结论边界 |
| [查看／复现](RUNNING.md) | 当前静态、裸骨规定运动和历史安装步骤分别说明 |

191是含骨、肌腹、五金、钢丝及台面的对象数，不是打印件数或完成度。两块肩袖肌为共形Y格，三角肌与大圆肌保留旧Z格。Y2对4处已分类的非通道／非端接口闭腔显式填实，局部刚度影响未测；其余189源及变换、骨骼、钢丝和锚点不变。

## 已有证据及边界

| 项目 | 当前证据 | 尚未验证 |
|---|---|---|
| 191对象零位 | 18145配对、六主动线无超0.001 mm³判据刚性相交 | 全姿态、严格零穿透、实际穿线与驱动 |
| 裸骨运动 | 正屈曲／外展到120°检查上限；负向−39.6°台面首限／−39.1°夹爪螺钉首限；长轴±120°到检查上限 | 89对象未装，不能转为带肌肉范围；内外旋正负映射未核定 |
| 名义放置及保持 | 11步、0.1 mm往返4322姿态；0↔10°裸肱骨22姿态；52方向球头有限捕获 | 真实紧固／工具、强度、任意联合逃逸及完整肌肉安装 |
| 自由端／肌腱 | 端点与张力单边约束、局部端环检查、规定端环角的完整肌腱候选 | 安装长度、预紧、压接与自由端盖—TPU耦合平衡 |
| 肌腹与工艺 | 两块实际共形晶格、连续通道；方向样件8类20对象 | 设备／材料／切片、穿线、实际收缩／张力；24份测量记录仍为空 |
| 肘部 | 既有候选和试验准备独立保留 | 本轮未复验，不能据此宣称实际驱动通过 |

裸骨常规1°、首界0.1°，从零位连续采样并独立回程。部分最小间隙为0；120°不是最大范围。0.001 mm³是报告阈值，不是制造公差。详细输入及报告见[技术文件第5–6节](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)。

## 当前阻塞项

1. **完整TPU分析边界未通过。** [原模型六组体网格尝试](02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/volume_mesh/README.md)均失败，远端过渡已定位边—面交点。立方体对照仅验证工具流程。当前没有合格完整肌腹四面体或TPU变形解。
2. **固定及完整安装未验收。** 实际紧固件、工具、预紧、座壁／横销承载、后衬垫与夹爪、肌腹及五金的完整安装仍待确认。三角肌／肩胛约1.84044 mm³未变形柔性重叠继续保留。
3. **真实驱动缺少数据。** 安装长度、压接、摩擦、材料与预紧均未标定；肩胛下肌约16.98°斜向出线需弯折／侧向接触试验。不能固定自由端盖或借用裸骨角度计算通过。

现有“零边界／非流形边”基于四位小数分桶导入，不保证原始STL精确拓扑或体网格适用性。另有局部和全表面修复实验，均未采用；导出回读仍有交点，见[技术文件5.3节](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)。默认选择与源几何保持不变。

## 下一步及验收

数字分析先修复完整肌腹边界，核对接口／通道、最终文件与体网格质量；并行落实实际紧固、完整安装及打印／穿线工艺。网格失败不等于已经证明不能打印，可打印性也不能只凭表面闭合判断。

随后标定单肌肉、压接与预紧，求解自由端盖—肌腱—TPU变形和接触，再验证拮抗驱动、往返与联合动作。具体交付条件见[技术文件第7节](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN.md)。原K1–K7不降低，旧刚性范围和论文性能不填作当前驱动能力。

## 文档及历史

- [设计任务书](design/DESIGN_BRIEF.md)：当前入口以下历史原样保留，不同时生效。
- [v0.9.16技术全文归档](02_人体骨骼晶格肌肉方案/输出/shoulder_design/SHOULDER_DESIGN_v0.9.16_ARCHIVE.md)：详细谱系、修改量、失败和版本记录。
- [本次整理前快照](reviews/S1_2026-11-03/document_optimization/before_SHA256.json)：文档与候选指针哈希。
- [工作规则](PROJECT_CONTROL.md)、[计划](task_plan.md)、[日志](progress.md)。
'''
(ROOT / 'PROJECT_STATUS.md').write_text(status, encoding='utf-8', newline='\n')
