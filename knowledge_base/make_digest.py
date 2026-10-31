"""Condense each extracted full text into a compact title/authors/abstract digest."""

from __future__ import annotations

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml

import json
import re
import sys
from project_runtime import WorkspacePath as Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

KB = Path(r"D:\XJTLU\PSP\knowledge_base")

NOISE_LINE = re.compile(r"^[\s\-_=.*|]{20,}$")
NOISE_FRAG = (
    "creative commons", "provided the original work is properly cited",
    "this is an open access article", "terms of the creative commons",
    "all rights reserved", "downloaded from", "see the terms and conditions",
)


def tidy_block(text: str) -> str:
    lines = []
    for raw in text.split("\n"):
        line = raw.strip()
        if not line or NOISE_LINE.match(line):
            continue
        if any(frag in line.lower() for frag in NOISE_FRAG):
            continue
        lines.append(line)
    return " ".join(lines)


def main() -> int:
    records = json.loads((KB / "manifest.json").read_text(encoding="utf-8"))
    chunks = []
    for rec in records:
        slug = rec["slug"]
        fpath = KB / "fulltext" / f"{slug}.txt"
        text = fpath.read_text(encoding="utf-8")
        body = re.sub(r"===== PAGE \d+ =====", "", text)

        head = tidy_block(body[:2500])[:700]

        abstract = ""
        for pat in (r"\bABSTRACT\b", r"\bAbstract\b", r"\bA B S T R A C T\b"):
            m = re.search(pat, body)
            if m:
                abstract = tidy_block(body[m.end():m.end() + 2600])[:1500]
                break
        if not abstract:
            abstract = tidy_block(body[:4000])[:1200]

        chunks.append(
            f"\n{'#' * 76}\n## {rec['file']}\n"
            f"pages={rec['pages']}  chars={rec['chars']}  slug={slug}\n"
            f"--HEAD--\n{head}\n--ABSTRACT--\n{abstract}\n"
        )
    out = KB / "tmp" / "digest.txt"
    out.write_text("\n".join(chunks), encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size} bytes, {len(records)} docs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
