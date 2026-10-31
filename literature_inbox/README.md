# 文献投放与入库目录

本目录用于接收“肩部系统与钢丝/接触优化知识库扩充”计划中的新论文。

## 文件命名

每个 PDF 使用候选编号作为前缀：

```text
40_kobayashi_1998_redundant_tendons.pdf
41_jamshidifar_2017_constrained_cable_robot.pdf
...
```

候选编号、题录和检索状态以 `candidates.csv` 为准。不要修改编号。

## 获取规则

1. 优先使用出版社、Crossref、OpenAlex 或作者公开的正式版本。
2. 受限论文由用户使用学校账号手动登录并获得全文，然后保存到本目录。
3. 不记录、不保存、不请求学校账号密码。
4. 普通预印本不进入正式知识库。
5. 扫描件必须有文本层；否则先 OCR，再入库。

## 检查

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\knowledge_base\check_literature_intake.py
```

检查报告写入：

```text
D:\XJTLU\PSP\knowledge_base\literature_intake_report.json
```

达到 12 篇后，再把每篇文献写入 `knowledge_base/annotations.json`，然后运行：

```powershell
D:\XJTLU\PSP\knowledge_base\rebuild_knowledge_base.bat
```

当前状态：12 篇候选的 DOI 已通过 Crossref 核验。8 篇取得全文并已入库；
文献 40、44、45、46 经用户确认无法通过学校权限获取，保留为题录候选。
