"""Extract text and metadata from every PDF under D:\\XJTLU\\PSP.

Outputs (under knowledge_base/):
  fulltext/<slug>.txt   per-document plain text, page-delimited
  manifest.json         machine-readable record for every document
  metadata.csv          flat table
  tmp/abstracts.txt     first-page text of every document, for human review
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
import sys
import unicodedata
from pathlib import Path

from pypdf import PdfReader

# the console on this machine is GBK; keep printing safe
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
logging.getLogger("pypdf").setLevel(logging.ERROR)
logging.getLogger("pypdf._reader").setLevel(logging.ERROR)

ROOT = Path(r"D:\XJTLU\PSP")
KB = ROOT / "knowledge_base"
FULLTEXT = KB / "fulltext"
TMP = KB / "tmp"


def slugify(name: str) -> str:
    stem = Path(name).stem
    stem = unicodedata.normalize("NFKD", stem)
    stem = stem.replace("\u2010", "-").replace("\u2011", "-").replace("\u2013", "-")
    stem = re.sub(r"[^0-9A-Za-z]+", "_", stem)
    stem = re.sub(r"_+", "_", stem).strip("_")
    stem = stem[:90]
    digest = hashlib.md5(name.encode("utf-8")).hexdigest()[:6]
    return f"{stem}__{digest}"


def clean(text: str) -> str:
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def norm_for_fingerprint(text: str) -> str:
    return re.sub(r"[^0-9a-z]+", "", text.lower())[:1500]


def main() -> int:
    FULLTEXT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)

    pdfs = sorted(
        p for p in ROOT.rglob("*")
        if p.suffix.lower() == ".pdf" and KB not in p.parents
    )
    print(f"found {len(pdfs)} pdf files")

    records = []
    abstracts = []

    for path in pdfs:
        rel = str(path.relative_to(ROOT))
        rec = {
            "file": rel,
            "name": path.name,
            "size_mb": round(path.stat().st_size / 1e6, 2),
            "slug": slugify(path.name),
        }
        try:
            reader = PdfReader(str(path))
            rec["pages"] = len(reader.pages)
            rec["encrypted"] = bool(reader.is_encrypted)
            meta = reader.metadata or {}
            rec["meta"] = {
                k.lstrip("/"): str(v) for k, v in meta.items() if v is not None
            }
            pages = []
            for i, page in enumerate(reader.pages):
                try:
                    pages.append(page.extract_text() or "")
                except Exception as exc:  # noqa: BLE001
                    pages.append("")
                    rec.setdefault("page_errors", {})[i + 1] = repr(exc)
        except Exception as exc:  # noqa: BLE001
            rec["error"] = repr(exc)
            rec["pages"] = 0
            pages = []

        body = clean("\n\n".join(
            f"\n===== PAGE {i + 1} =====\n{t}" for i, t in enumerate(pages)
        ))
        rec["chars"] = len(body)
        rec["chars_per_page"] = round(rec["chars"] / rec["pages"], 1) if rec["pages"] else 0
        rec["has_text_layer"] = rec["chars_per_page"] >= 200

        out = FULLTEXT / f"{rec['slug']}.txt"
        out.write_text(body, encoding="utf-8")

        first = clean(pages[0]) if pages else ""
        rec["first_page_chars"] = len(first)
        rec["fingerprint"] = hashlib.md5(
            norm_for_fingerprint(first or body).encode("utf-8")
        ).hexdigest()[:12]
        rec["sha1"] = hashlib.sha1(path.read_bytes()).hexdigest()[:12]

        abstracts.append(
            f"\n{'=' * 78}\nFILE: {rel}\nPAGES: {rec['pages']}  "
            f"CHARS: {rec['chars']}  TEXTLAYER: {rec['has_text_layer']}\n{'=' * 78}\n"
            + (first[:3000] if first else "[no extractable text on page 1]")
        )
        records.append(rec)
        print(f"  {rel}  pages={rec['pages']} chars={rec['chars']}")

    # duplicate detection: identical bytes, or identical opening text
    by_sha1: dict[str, list[str]] = {}
    by_fp: dict[str, list[str]] = {}
    for rec in records:
        if rec["chars"] == 0:
            continue
        by_sha1.setdefault(rec["sha1"], []).append(rec["file"])
        by_fp.setdefault(rec["fingerprint"], []).append(rec["file"])
    exact = [v for v in by_sha1.values() if len(v) > 1]
    same_open = [v for v in by_fp.values() if len(v) > 1]
    for rec in records:
        rec["exact_duplicate_of"] = next(
            (g for g in exact if rec["file"] in g and len(g) > 1), None
        )
        rec["shares_opening_with"] = next(
            (g for g in same_open if rec["file"] in g and len(g) > 1), None
        )

    (KB / "manifest.json").write_text(
        json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (TMP / "abstracts.txt").write_text("\n".join(abstracts), encoding="utf-8")

    with (KB / "metadata.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "pages", "chars", "chars_per_page", "has_text_layer",
                    "title", "author", "subject", "doi", "size_mb", "slug"])
        for rec in records:
            m = rec.get("meta", {})
            w.writerow([rec["file"], rec["pages"], rec["chars"],
                        rec["chars_per_page"], rec["has_text_layer"],
                        m.get("Title", ""), m.get("Author", ""), m.get("Subject", ""),
                        m.get("doi", "") or m.get("DOI", ""), rec["size_mb"], rec["slug"]])

    missing = [r["file"] for r in records if not r["has_text_layer"]]
    print(f"\nno usable text layer ({len(missing)}):")
    for f in missing:
        print("  -", f)
    print(f"\nexact duplicate groups ({len(exact)}):")
    for g in exact:
        print("  -", " == ".join(g))
    print(f"\nsame-opening groups ({len(same_open)}):")
    for g in same_open:
        print("  -", " ~ ".join(g))
    return 0


if __name__ == "__main__":
    sys.exit(main())
