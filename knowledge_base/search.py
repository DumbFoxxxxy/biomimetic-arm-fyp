"""Keyword search over the literature knowledge base.

Examples:
    python search.py scapula
    python search.py 肩胛骨 节律
    python search.py --tag X晶格
    python search.py --tier 核心
    python search.py --doc 14
    python search.py --list

Without --full it searches notes, titles, tags and Chinese summaries only.
With --full it also scans every word of the extracted PDF text.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

KB = Path(r"D:\XJTLU\PSP\knowledge_base")


def load() -> list[dict]:
    return json.loads((KB / "manifest_enriched.json").read_text(encoding="utf-8"))


def haystack(doc: dict, full: bool) -> str:
    parts = [
        doc.get("file", ""), doc.get("title", ""), doc.get("title_zh", ""),
        doc.get("authors", ""), doc.get("venue", ""), doc.get("kind", ""),
        doc.get("doi", ""),
        " ".join(doc.get("tags", [])), doc.get("summary_zh", ""),
        " ".join(doc.get("key_points", [])), doc.get("relevance_zh", ""),
    ]
    if full:
        ft = KB / "fulltext" / f"{doc.get('slug', '')}.txt"
        if ft.exists():
            parts.append(ft.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(p for p in parts if p)


def snippets(text: str, terms: list[str], width: int = 170) -> list[str]:
    out = []
    low = text.lower()
    for term in terms:
        for m in re.finditer(re.escape(term.lower()), low):
            a = max(0, m.start() - width // 2)
            b = min(len(text), m.end() + width // 2)
            frag = re.sub(r"\s+", " ", text[a:b]).strip()
            out.append(("..." if a else "") + frag + ("..." if b < len(text) else ""))
            if len(out) >= 3:
                return out
    return out


def show(doc: dict, terms: list[str], full: bool, verbose: bool) -> None:
    print(f"\n[{doc['id']:02d}] {doc['title_zh']}  ({doc['tier']})")
    print(f"     {doc['title']}")
    print(f"     {doc['venue']}  {doc.get('year') or ''}")
    print(f"     源文件: D:\\XJTLU\\PSP\\{doc['file']}")
    print(f"     笔记:   {KB / doc['note']}")
    if verbose:
        print(f"     标签:   {', '.join(doc['tags'])}")
        print(f"     摘要:   {doc['summary_zh'][:220]}...")
    if terms:
        hits = snippets(haystack(doc, full), terms)
        for h in hits:
            print(f"       · {h}")


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("terms", nargs="*", help="关键词，可多个（任一匹配即可）")
    ap.add_argument("--full", action="store_true", help="同时检索 PDF 全文")
    ap.add_argument("--tag", help="按标签精确筛选")
    ap.add_argument("--tier", help="按等级筛选：核心/重要/参考/背景/离题")
    ap.add_argument("--group", help="按分组 key 筛选，如 core/shoulder/platform")
    ap.add_argument("--doc", type=int, help="显示指定编号的文献")
    ap.add_argument("--tags", action="store_true", help="列出所有标签及数量")
    ap.add_argument("--list", action="store_true", help="列出全部文献")
    args = ap.parse_args()

    docs = load()

    if args.tags:
        counts: dict[str, list[int]] = {}
        for d in docs:
            for t in d["tags"]:
                counts.setdefault(t, []).append(d["id"])
        for t in sorted(counts, key=lambda x: (-len(counts[x]), x)):
            print(f"{len(counts[t]):2d}  {t:<32} {' '.join(f'{i:02d}' for i in counts[t])}")
        return 0

    if args.doc:
        sel = [d for d in docs if d["id"] == args.doc]
        if not sel:
            print(f"没有编号 {args.doc} 的文献")
            return 1
        d = sel[0]
        show(d, [], False, True)
        print("\n要点:")
        for p in d["key_points"]:
            print(f"  - {p}")
        print(f"\n与本项目的关系: {d['relevance_zh']}")
        note = KB / d["note"]
        if note.exists():
            print(f"\n完整笔记: {note}")
        return 0

    sel = docs
    if args.tag:
        sel = [d for d in sel if args.tag in d["tags"]]
    if args.tier:
        sel = [d for d in sel if d["tier"] == args.tier]
    if args.group:
        sel = [d for d in sel if d["group"] == args.group]

    terms = [t for t in args.terms if t.strip()]
    if terms:
        keep = []
        for d in sel:
            text = haystack(d, args.full).lower()
            if any(t.lower() in text for t in terms):
                keep.append(d)
        sel = keep

    if args.list or (not terms and not args.tag and not args.tier and not args.group):
        for d in sorted(sel, key=lambda x: x["id"]):
            print(f"{d['id']:02d}  {d['tier']:<3} {d['title_zh'][:46]:<48} "
                  f"{str(d.get('year') or '—'):<5} {d['venue'][:40]}")
        print(f"\n共 {len(sel)} 篇")
        return 0

    for d in sorted(sel, key=lambda x: x["id"]):
        show(d, terms, args.full, False)
    print(f"\n共匹配 {len(sel)} 篇；单篇详情用 --doc 编号")
    return 0


if __name__ == "__main__":
    sys.exit(main())
