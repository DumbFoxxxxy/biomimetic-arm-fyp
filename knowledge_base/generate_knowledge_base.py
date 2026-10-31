"""Build the browsable knowledge base from annotations.json + manifest.json.

Outputs:
  INDEX.md              grouped master index (Chinese)
  notes/NN_<short>.md   one note per document
  documents.csv         flat table with tags and summaries
  manifest_enriched.json  manifest + annotation fields, for scripts/RAG
  search.py             regenerated keyword search helper
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(r"D:\XJTLU\PSP")
KB = ROOT / "knowledge_base"
NOTES = KB / "notes"
TIER_ORDER = ["核心", "重要", "参考", "背景", "离题"]


def main() -> int:
    NOTES.mkdir(parents=True, exist_ok=True)
    ann = json.loads((KB / "annotations.json").read_text(encoding="utf-8"))
    manifest = json.loads((KB / "manifest.json").read_text(encoding="utf-8"))
    by_file = {m["file"]: m for m in manifest}
    groups = {g["key"]: g for g in ann["groups"]}
    docs = sorted(ann["documents"], key=lambda d: d["id"])

    enriched = []
    for d in docs:
        m = dict(by_file.get(d["file"], {}))
        m.update(
            id=d["id"], group=d["group"], tier=d["tier"], short=d["short"],
            title=d["title"], title_zh=d["title_zh"],
            authors=d["authors"], venue=d["venue"],
            year=d["year"], kind=d["kind"], tags=d["tags"],
            doi=d.get("doi", m.get("meta", {}).get("doi", "")),
            summary_zh=d["summary_zh"], key_points=d["key_points"],
            relevance_zh=d["relevance_zh"], duplicates=d.get("duplicates", []),
            note=f"notes/{d['id']:02d}_{d['short']}.md",
        )
        enriched.append(m)
    (KB / "manifest_enriched.json").write_text(
        json.dumps(enriched, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # ---- per-document notes -------------------------------------------------
    for d in docs:
        m = by_file.get(d["file"], {})
        lines = [
            f"# {d['id']:02d} · {d['short']}",
            "",
            f"- **中文标题**：{d['title_zh']}",
            f"- **原题**：{d['title']}",
            f"- **作者**：{d['authors']}",
            f"- **出处**：{d['venue']}"
            + (f"（{d['year']}）" if d["year"] else ""),
            f"- **类型**：{d['kind']}　**等级**：{d['tier']}"
            f"　**分组**：{groups[d['group']]['name_zh']}",
            "- **标签**：" + " ".join(f"`{t}`" for t in d["tags"]),
        ]
        if d.get("doi") or m.get("meta", {}).get("doi"):
            lines.append(
                f"- **DOI**：{d.get('doi') or m.get('meta', {}).get('doi')}"
            )
        lines += [
            f"- **源文件**：`{ROOT / d['file']}`",
            f"- **提取全文**：`{KB / 'fulltext' / (d['slug'] + '.txt')}`",
            f"- **规模**：{m.get('pages', '?')} 页 / 提取 {m.get('chars', '?')} 字符",
            "",
            "## 摘要",
            "",
            d["summary_zh"],
            "",
            "## 要点",
            "",
        ]
        lines += [f"- {p}" for p in d["key_points"]]
        lines += ["", "## 与本项目的关系", "", d["relevance_zh"]]
        if d.get("duplicates"):
            lines += ["", "## 重复或同源文件", ""]
            lines += [f"- {x}" for x in d["duplicates"]]
        (NOTES / f"{d['id']:02d}_{d['short']}.md").write_text(
            "\n".join(lines) + "\n", encoding="utf-8"
        )

    # ---- flat table ---------------------------------------------------------
    with (KB / "documents.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["编号", "等级", "分组", "年份", "原题", "中文标题", "作者",
                    "出处", "类型", "DOI", "标签", "页数", "源文件", "笔记", "重复"])
        for d in docs:
            m = by_file.get(d["file"], {})
            w.writerow([
                f"{d['id']:02d}", d["tier"], groups[d["group"]]["name_zh"],
                d["year"] or "", d["title"], d["title_zh"], d["authors"],
                d["venue"], d["kind"],
                d.get("doi", m.get("meta", {}).get("doi", "")),
                "; ".join(d["tags"]), m.get("pages", ""),
                d["file"], f"notes/{d['id']:02d}_{d['short']}.md",
                " / ".join(d.get("duplicates", [])),
            ])

    # ---- master index -------------------------------------------------------
    out = [
        "# 右臂仿生项目文献知识库",
        "",
        f"- 源目录：`{ROOT}`（含 `literature_inbox/` 新文献投放目录）",
        f"- PDF 文件：**{len(manifest)}** 个 → 去重后 **{len(docs)}** 篇独立文献",
        "- 全文提取：`fulltext/`（每篇一个纯文本文件，按页分隔，可全文检索）",
        "- 单篇笔记：`notes/`（含中文摘要、要点、与本项目的关系）",
        "- 机器可读：`manifest.json`（原始提取）、`manifest_enriched.json`（含注解）、`documents.csv`",
        "- 关键词检索：`python search.py 关键词`",
        "",
        "等级含义：**核心** = 方案直接依据；**重要** = 肩肘机构与平台设计主要参考；"
        "**参考** = 方法储备与对照路线；**背景** = 相关但非设计文献；**离题** = 与本项目无关。",
        "",
    ]

    for key, g in groups.items():
        sel = [d for d in docs if d["group"] == key]
        if not sel:
            continue
        out += [f"## {g['name_zh']}", "", g["desc"], "",
                "| 编号 | 年 | 文献 | 出处 | 等级 |", "| --- | --- | --- | --- | --- |"]
        for d in sel:
            link = f"[{d['title_zh']}](notes/{d['id']:02d}_{d['short']}.md)"
            if d["year"]:
                yr = str(d["year"])
            else:
                yr = "—"
            out.append(f"| {d['id']:02d} | {yr} | {link} | {d['venue']} | {d['tier']} |")
        out.append("")

    # tag index
    tags: dict[str, list[int]] = {}
    for d in docs:
        for t in d["tags"]:
            tags.setdefault(t, []).append(d["id"])
    out += ["## 标签索引", ""]
    for t in sorted(tags, key=lambda x: (-len(tags[x]), x)):
        out.append(f"- `{t}`：" + " ".join(f"{i:02d}" for i in tags[t]))
    out.append("")

    # duplicate map
    out += ["## 重复与同源文件", ""]
    for d in docs:
        if d.get("duplicates"):
            out.append(f"- **{d['id']:02d} {d['short']}** 与："
                       + "；".join(d["duplicates"]))
    out += [
        "",
        "说明：3 组文件是逐字节相同的重复副本（ECCE1、Humanoid bioreactor、"
        "Multifilament lower-limb），已在上面标出；"
        "另外 Ikemoto 的肩复合体与 Haosen Yang 的肩设计各有一对同源的不同版本，"
        "内容互补而非重复。",
        "",
    ]
    (KB / "INDEX.md").write_text("\n".join(out), encoding="utf-8")

    print(f"wrote INDEX.md, {len(docs)} notes, documents.csv, manifest_enriched.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
