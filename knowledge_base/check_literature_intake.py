"""Validate newly supplied literature before adding it to annotations.json."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
from difflib import SequenceMatcher
from pathlib import Path

from pypdf import PdfReader


logging.getLogger("pypdf").setLevel(logging.ERROR)


ROOT = Path(__file__).resolve().parents[1]
KB = ROOT / "knowledge_base"
INBOX = ROOT / "literature_inbox"
CANDIDATES = INBOX / "candidates.csv"
REPORT = KB / "literature_intake_report.json"


def normalize(text: str) -> str:
    return re.sub(r"[^0-9a-z]+", "", text.lower())[:1500]


def load_existing() -> tuple[set[str], set[str]]:
    records = json.loads((KB / "manifest.json").read_text(encoding="utf-8"))
    records = [
        item
        for item in records
        if not item.get("file", "").replace("/", "\\").startswith(
            "literature_inbox\\"
        )
    ]
    return (
        {item["sha1"] for item in records if item.get("sha1")},
        {item["fingerprint"] for item in records if item.get("fingerprint")},
    )


def inspect_pdf(path: Path) -> dict:
    reader = PdfReader(str(path))
    pages = len(reader.pages)
    extracted = []
    for page in reader.pages[: min(3, pages)]:
        extracted.append(page.extract_text() or "")
    text = "\n".join(extracted)
    chars_per_page = len(text) / max(1, min(3, pages))
    return {
        "file": path.name,
        "meta_title": str((reader.metadata or {}).get("/Title", "") or ""),
        "pages": pages,
        "sample_chars": len(text),
        "sample_chars_per_page": round(chars_per_page, 1),
        "has_text_layer": chars_per_page >= 200,
        "sha1": hashlib.sha1(path.read_bytes()).hexdigest()[:12],
        "fingerprint": hashlib.md5(normalize(text).encode("utf-8")).hexdigest()[:12],
        "sample_text": re.sub(r"\s+", " ", text)[:500],
        "match_text": re.sub(r"\s+", " ", text)[:12000],
    }


def title_similarity(left: str, right: str) -> float:
    def words(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()

    return SequenceMatcher(None, words(left), words(right)).ratio()


def match_candidate(pdf: Path, info: dict, candidates: list[dict]) -> dict | None:
    prefix = pdf.name[:2]
    by_prefix = next(
        (item for item in candidates if item["candidate_id"] == prefix),
        None,
    )
    if by_prefix is not None:
        return by_prefix

    haystack = f"{info['meta_title']} {info['match_text']}"
    haystack_lower = haystack.lower()
    for item in candidates:
        doi = item.get("doi", "").strip().lower()
        if doi and doi in haystack_lower:
            return item

    scored = []
    for item in candidates:
        title = item.get("title", "").strip()
        if not title:
            continue
        score = max(
            title_similarity(title, info["meta_title"]),
            title_similarity(title, info["match_text"][:1000]),
        )
        scored.append((score, item))
    if not scored:
        return None
    score, candidate = max(scored, key=lambda pair: pair[0])
    return candidate if score >= 0.72 else None


def candidate_filename(candidate: dict) -> str:
    short = re.sub(r"[^A-Za-z0-9]+", "_", candidate["title"]).strip("_").lower()
    return f"{candidate['candidate_id']}_{short[:100]}.pdf"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-count", type=int, default=None)
    parser.add_argument("--rename-matched", action="store_true")
    args = parser.parse_args()
    with CANDIDATES.open(newline="", encoding="utf-8-sig") as handle:
        candidates = list(csv.DictReader(handle))

    existing_sha1, existing_fingerprints = load_existing()
    pdfs = sorted(INBOX.glob("*.pdf"))
    report_rows = []
    errors = []
    warnings = []
    matched_ids: set[str] = set()
    seen_sha1: dict[str, str] = {}
    seen_fingerprints: dict[str, str] = {}

    for pdf in pdfs:
        try:
            info = inspect_pdf(pdf)
        except Exception as exc:
            errors.append(f"{pdf.name}: cannot read PDF: {exc}")
            report_rows.append({"file": pdf.name, "status": "ERROR", "error": str(exc)})
            continue

        candidate = match_candidate(pdf, info, candidates)
        if candidate is None:
            info["status"] = "IGNORED_UNRELATED"
            warnings.append(f"{pdf.name}: no candidate title or DOI match; ignored")
            report_rows.append(info)
            continue

        candidate_id = candidate["candidate_id"]
        if args.rename_matched and not pdf.name.startswith(candidate_id):
            target = INBOX / candidate_filename(candidate)
            if target.exists():
                errors.append(f"{pdf.name}: cannot rename, target exists: {target.name}")
                info["status"] = "RENAME_CONFLICT"
                report_rows.append(info)
                continue
            info["original_file"] = pdf.name
            pdf.rename(target)
            pdf = target
            info["file"] = target.name
        matched_ids.add(candidate_id)
        info["candidate_id"] = candidate_id
        info["title"] = candidate["title"]
        info["status"] = "READY"

        if info["pages"] <= 3:
            info["status"] = "REJECTED"
            errors.append(f"{pdf.name}: page count must exceed 3")
        if not info["has_text_layer"]:
            info["status"] = "OCR_REQUIRED"
            errors.append(f"{pdf.name}: OCR_REQUIRED, no usable text layer")
        if info["sha1"] in existing_sha1:
            info["status"] = "DUPLICATE_EXISTING"
            errors.append(f"{pdf.name}: duplicate of an existing knowledge-base PDF")
        if info["sha1"] in seen_sha1:
            info["status"] = "DUPLICATE_INBATCH"
            errors.append(f"{pdf.name}: duplicate of {seen_sha1[info['sha1']]}")
        if info["fingerprint"] in existing_fingerprints:
            info["status"] = "DUPLICATE_EXISTING"
            errors.append(f"{pdf.name}: opening text duplicates an existing PDF")
        if info["fingerprint"] in seen_fingerprints:
            info["status"] = "DUPLICATE_INBATCH"
            errors.append(
                f"{pdf.name}: opening text duplicates {seen_fingerprints[info['fingerprint']]}"
            )

        seen_sha1[info["sha1"]] = pdf.name
        seen_fingerprints[info["fingerprint"]] = pdf.name
        report_rows.append(info)

    nonempty_dois = [
        item["doi"].strip().lower()
        for item in candidates
        if item.get("doi", "").strip()
    ]
    if len(nonempty_dois) != len(set(nonempty_dois)):
        errors.append("candidate list contains duplicate DOI values")

    pending = [
        item["candidate_id"]
        for item in candidates
        if item["status"] != "unavailable"
        and item["candidate_id"] not in matched_ids
    ]
    unavailable = [
        item["candidate_id"] for item in candidates if item["status"] == "unavailable"
    ]
    accessible_count = len(candidates) - len(unavailable)
    required_count = (
        accessible_count if args.require_count is None else args.require_count
    )
    if len(matched_ids) < required_count:
        errors.append(
            f"matched {len(matched_ids)} PDFs; {required_count} required"
        )

    report = {
        "candidate_file": str(CANDIDATES),
        "inbox": str(INBOX),
        "required_count": required_count,
        "matched_candidate_ids": sorted(matched_ids),
        "pending_candidate_ids": pending,
        "unavailable_candidate_ids": unavailable,
        "pdf_count": len(pdfs),
        "warnings": warnings,
        "errors": errors,
        "files": report_rows,
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        f"matched={len(matched_ids)} pending={len(pending)} "
        f"unavailable={len(unavailable)} "
        f"pdfs={len(pdfs)} errors={len(errors)}"
    )
    print(f"report: {REPORT}")
    if pending:
        print("pending: " + ", ".join(pending))
    if unavailable:
        print("unavailable: " + ", ".join(unavailable))
    for warning in warnings:
        print(f"WARNING: {warning}")
    for error in errors:
        print(f"ERROR: {error}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
