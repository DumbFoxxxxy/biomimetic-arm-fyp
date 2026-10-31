# 文献知识库使用说明

把 `D:\XJTLU\PSP` 根目录和 `literature_inbox/` 下的 PDF 整理成了一个可检索、
可扩展的知识库。当前已提取 **50 个 PDF、47 篇独立文献**。

## 目录结构

| 路径 | 内容 |
| --- | --- |
| `INDEX.md` | 主索引。按 8 个主题分组，每篇给出一行摘要与等级 |
| `notes/` | 47 篇单篇笔记，每篇含中文摘要、要点、与本项目的关系 |
| `fulltext/` | 47 篇提取出的纯文本全文，按页分隔，供检索 |
| `documents.csv` | 扁平总表，可直接用 Excel 排序筛选 |
| `manifest.json` | 原始提取结果（页数、字符数、PDF 元数据、去重指纹） |
| `manifest_enriched.json` | 上面的结果加上全部注解，供脚本或以后接 RAG |
| `annotations.json` | **人工撰写的注解源文件**，修改后重新生成即可 |
| `search.py` | 关键词检索工具 |
| `tmp/` | 中间产物（摘要清单 `digest.txt` 等），可随时删 |

当前新增分组：

| 分组 | 内容 |
| --- | --- |
| `cable_contact` | 钢丝张力分配、正张成、导轮摩擦、肩胸接触与接触优化 |

## 怎么用

查文献：

```
D:\MujocoEnv\Scripts\python.exe search.py 肩胛骨 节律
D:\MujocoEnv\Scripts\python.exe search.py --tag X晶格
D:\MujocoEnv\Scripts\python.exe search.py --tier 核心
D:\MujocoEnv\Scripts\python.exe search.py --doc 14
D:\MujocoEnv\Scripts\python.exe search.py 韧带 --full
```

不加 `--full` 只搜标题、作者、标签、摘要与要点，很快。加 `--full` 会把每篇 PDF 的
正文全文一起扫一遍，用于找具体数值、公式附近的原话。

用编辑器看：直接打开 `INDEX.md`，点表格里的标题就进单篇笔记。

## 等级怎么定的

- **核心**：方案直接依据。改这几篇的做法就等于改方案。
- **重要**：肩肘机构与整机平台设计的主要参考。
- **参考**：方法储备或对照路线（软体驱动、本构建模、优化方法）。
- **背景**：相关但不是设计文献。
- **离题**：与本项目无关。

## 重复与同源关系

50 个文件里有 3 组是逐字节相同的重复副本，已在索引中标注：

- `ECCE1 the first of a series...` = `ECCE1_The_first_of_a_series_of_anthropom.pdf`
- `Humanoid robots to mechanically stress human cells...` = `Humanoid robots to mechanically stress human.pdf`
- `Musculoskeletal lower-limb robot driven by multifilament muscles.pdf` = `Musculoskeletal lower‑limb robot driven.pdf`

另外两对**不是**重复，而是同一工作的不同版本，内容互补，两篇都保留了：

- 09（Ikemoto, IROS 2012 会议版）与 10（同作者 2014 期刊预印本，多了机构推证与肌电实验）
- 14（Haosen Yang, IEEE T-RO 2024）与 15（同作者曼彻斯特博士论文，多了驱动器 ICA/ECA/MCA 与完整术语表）

## 加入新 PDF

1. 把新 PDF 按候选编号前缀放进 `D:\XJTLU\PSP\literature_inbox\`。
   候选编号和 DOI 见 `literature_inbox\candidates.csv`。
2. 运行 `rebuild_knowledge_base.bat`，它会重新提取全文并刷新索引与总表。
3. 打开 `annotations.json`，按现有格式给新文献补一条中文注解（`id` 接着编，`file` 写文件名，
   `slug` 从刷新后的 `manifest.json` 里对应该文件的记录抄过来）。
4. 再跑一次 `generate_knowledge_base.py`，笔记与索引就会带上新条目。

正式写入 `annotations.json` 前先检查 PDF：

```powershell
D:\MujocoEnv\Scripts\python.exe knowledge_base\check_literature_intake.py
```

没有文本层的扫描件会标为 `OCR_REQUIRED`，不会进入正式文献笔记。

### 当前扩充状态

`literature_inbox/candidates.csv` 固定了候选编号 40–51。12 篇文献的题录和 DOI
已经核验；其中 8 篇取得全文并已入库。文献 40、44、45、46 经用户确认无法取得
学校权限，保留题录并标记为 `unavailable`。

## 重新生成

```
rebuild_knowledge_base.bat
```

或分步：

```
D:\MujocoEnv\Scripts\python.exe build_knowledge_base.py   # 提取全文与元数据
D:\MujocoEnv\Scripts\python.exe make_digest.py            # 生成 tmp/digest.txt
D:\MujocoEnv\Scripts\python.exe generate_knowledge_base.py # 生成索引、笔记、总表
```

注意：`generate_knowledge_base.py` 会覆盖 `notes/` 与 `INDEX.md`，手工改过的内容要写进
`annotations.json` 而不是直接改生成结果。
