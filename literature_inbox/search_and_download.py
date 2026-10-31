"""Search DOI metadata and retrieve legal open-access PDFs for the intake list."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INBOX = ROOT / "literature_inbox"
CANDIDATES = INBOX / "candidates.csv"
RESULTS_JSON = INBOX / "candidate_search_results.json"
RESULTS_CSV = INBOX / "candidate_search_results.csv"

EMAIL = "Ziqian.Wang23@student.xjtlu.edu.cn"
USER_AGENT = f"FYPKnowledgeBase/1.0 (mailto:{EMAIL})"

SEARCH_QUERIES = {
    "48": "capstan friction tendon sheath cable-driven robot",
    "49": "cable routing unilateral contact tendon-driven mechanism optimization",
    "50": "scapulothoracic contact sliding simulation shoulder",
    "51": "redundant cable-driven parallel robot positive tension distribution wrench closure",
}

ACCEPTED_TYPES = {
    "article",
    "journal-article",
    "proceedings-article",
    "review",
}

DIRECT_OA_OVERRIDES = {
    "10.1080/09593985.2017.1422206": (
        "https://pure.rug.nl/ws/files/54074116/2_8_2018_Shoulder_a.pdf"
    ),
    "10.3390/app11094228": (
        "https://mdpi-res.com/d_attachment/applsci/applsci-11-04228/"
        "article_deploy/applsci-11-04228.pdf"
    ),
}


def request_json(url: str) -> dict:
    last_error = None
    for attempt in range(4):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=40) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code != 429 or attempt == 3:
                raise
            time.sleep(2.0 * (attempt + 1))
    raise last_error


def request_bytes(url: str) -> tuple[bytes, str]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read(), response.headers.get("Content-Type", "")


def crossref_search(query: str, rows: int = 5) -> list[dict]:
    encoded = urllib.parse.quote(query)
    data = request_json(
        "https://api.crossref.org/works"
        f"?query.bibliographic={encoded}&rows={rows}"
        "&select=DOI,title,author,issued,container-title,type,URL"
    )
    output = []
    for item in data["message"]["items"]:
        title = (item.get("title") or [""])[0]
        date_parts = item.get("issued", {}).get("date-parts", [[None]])
        year = date_parts[0][0] if date_parts and date_parts[0] else None
        authors = []
        for author in item.get("author", []):
            name = " ".join(
                part for part in (author.get("given"), author.get("family")) if part
            )
            if name:
                authors.append(name)
        output.append(
            {
                "doi": item.get("DOI", ""),
                "title": title,
                "authors": ", ".join(authors),
                "year": year,
                "venue": (item.get("container-title") or [""])[0],
                "type": item.get("type", ""),
                "landing_url": item.get("URL", ""),
            }
        )
    return output


def openalex_by_doi(doi: str) -> dict:
    if not doi:
        return {}
    return request_json(
        "https://api.openalex.org/works/https://doi.org/"
        + urllib.parse.quote(doi, safe="/:;()._-")
        + f"?mailto={urllib.parse.quote(EMAIL)}"
    )


def unpaywall_by_doi(doi: str) -> dict:
    if not doi:
        return {}
    return request_json(
        f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi, safe='/')}"
        f"?email={urllib.parse.quote(EMAIL)}"
    )


def openalex_search(query: str, rows: int = 8) -> list[dict]:
    data = request_json(
        "https://api.openalex.org/works"
        f"?search={urllib.parse.quote(query)}&per-page={rows}"
        f"&mailto={urllib.parse.quote(EMAIL)}"
    )
    output = []
    for item in data.get("results", []):
        open_access = item.get("open_access") or {}
        best = item.get("best_oa_location") or {}
        authors = [
            author.get("author", {}).get("display_name", "")
            for author in item.get("authorships", [])
        ]
        output.append(
            {
                "doi": (item.get("doi") or "").replace("https://doi.org/", ""),
                "title": item.get("display_name", ""),
                "authors": ", ".join(name for name in authors if name),
                "year": item.get("publication_year"),
                "venue": (
                    (item.get("primary_location") or {}).get("source") or {}
                ).get("display_name", ""),
                "type": item.get("type", ""),
                "is_oa": bool(open_access.get("is_oa")),
                "oa_url": best.get("pdf_url") or open_access.get("oa_url") or "",
                "landing_url": best.get("landing_page_url") or item.get("id", ""),
            }
        )
    return output


def oa_pdf_urls(doi: str) -> list[str]:
    urls: list[str] = []
    if doi.lower() in DIRECT_OA_OVERRIDES:
        urls.append(DIRECT_OA_OVERRIDES[doi.lower()])

    try:
        semantic = request_json(
            "https://api.semanticscholar.org/graph/v1/paper/DOI:"
            + urllib.parse.quote(doi, safe="/:;()._-")
            + "?fields=openAccessPdf"
        )
        semantic_url = (semantic.get("openAccessPdf") or {}).get("url")
        if semantic_url:
            urls.append(semantic_url)
    except Exception as exc:
        print(f"  Semantic Scholar lookup failed for {doi}: {exc}", file=sys.stderr)

    try:
        work = openalex_by_doi(doi)
        candidates = [work.get("best_oa_location") or {}]
        candidates.extend(work.get("locations") or [])
        for location in candidates:
            if location.get("pdf_url"):
                urls.append(location["pdf_url"])
    except Exception as exc:
        print(f"  OpenAlex lookup failed for {doi}: {exc}", file=sys.stderr)

    try:
        data = unpaywall_by_doi(doi)
        candidates = [data.get("best_oa_location") or {}]
        candidates.extend(data.get("oa_locations") or [])
        for location in candidates:
            if location.get("url_for_pdf"):
                urls.append(location["url_for_pdf"])
    except Exception as exc:
        print(f"  Unpaywall lookup failed for {doi}: {exc}", file=sys.stderr)

    unique = []
    for url in urls:
        if not url:
            continue
        for candidate in resolve_pdf_urls(url):
            if candidate not in unique:
                unique.append(candidate)
    return unique


def resolve_pdf_urls(url: str) -> list[str]:
    if (
        url.lower().endswith(".pdf")
        or "/pdf" in url.lower()
        or "article_deploy" in url.lower()
    ):
        return [url]
    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "Mozilla/5.0 (compatible; FYPKnowledgeBase/1.0)"},
        )
        with urllib.request.urlopen(request, timeout=40) as response:
            content_type = response.headers.get("Content-Type", "")
            if "pdf" in content_type.lower():
                return [url]
            body = response.read().decode("utf-8", errors="replace")
            final_url = response.geturl()
    except Exception as exc:
        print(f"  landing-page lookup failed for {url}: {exc}", file=sys.stderr)
        return [url]

    found = []
    citation_match = re.search(
        r'<meta[^>]+name=["\']citation_pdf_url["\'][^>]+content=["\']([^"\']+)',
        body,
        re.I,
    )
    if citation_match:
        found.append(html.unescape(citation_match.group(1)))
    for match in re.findall(r'href=["\']([^"\']+)["\']', body, re.I):
        lowered = match.lower()
        if any(token in lowered for token in ("pdf", "download", "bitstream", "blob")):
            found.append(urllib.parse.urljoin(final_url, html.unescape(match)))
    return found or [url]


def safe_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    return text[:72] or "paper"


def download_pdf(candidate_id: str, item: dict) -> str:
    for url in oa_pdf_urls(item.get("doi", "")):
        try:
            content, content_type = request_bytes(url)
        except Exception as exc:
            print(f"  download failed {url}: {exc}", file=sys.stderr)
            continue
        if not content.startswith(b"%PDF"):
            print(f"  skipped non-PDF response from {url} ({content_type})")
            continue
        suffix = f"{candidate_id}_{safe_name(item.get('title', 'paper'))}.pdf"
        path = INBOX / suffix
        path.write_bytes(content)
        return path.name
    return ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--top", type=int, default=5)
    args = parser.parse_args()

    with CANDIDATES.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    output = []
    for row in rows:
        candidate_id = row["candidate_id"]
        print(f"[{candidate_id}]")
        if row.get("title", "").strip():
            results = crossref_search(row["title"], rows=args.top)
            results.sort(
                key=lambda item: -_title_similarity(row["title"], item["title"])
            )
            selected = results[0] if results else {}
        else:
            query = SEARCH_QUERIES[candidate_id]
            try:
                results = [
                    item
                    for item in openalex_search(query, rows=max(args.top * 2, 10))
                    if item.get("doi") and item.get("type") in ACCEPTED_TYPES
                ]
            except Exception as exc:
                print(f"  OpenAlex search failed: {exc}", file=sys.stderr)
                results = [
                    item
                    for item in crossref_search(query, rows=max(args.top * 2, 10))
                    if item.get("doi") and item.get("type") in ACCEPTED_TYPES
                ]
            oa_results = [item for item in results if item.get("is_oa")]
            selected = (oa_results or results)[0] if results else {}

        entry = {
            "candidate_id": candidate_id,
            "priority": row["priority"],
            "category": row["category"],
            "method_target": row["method_target"],
            "query_title": row.get("title", ""),
            "selected": selected,
            "top_results": results,
            "downloaded_filename": (
                row.get("pdf_filename", "")
                if row.get("status", "") == "downloaded"
                else ""
            ),
        }
        if selected:
            print(
                f"  {selected.get('doi')} | {selected.get('year')} | "
                f"{selected.get('title')}"
            )
            if args.download and selected.get("doi"):
                entry["downloaded_filename"] = download_pdf(candidate_id, selected)
                if entry["downloaded_filename"]:
                    print(f"  saved {entry['downloaded_filename']}")
                else:
                    print("  no legal open-access PDF found")
                time.sleep(0.25)
        output.append(entry)

    RESULTS_JSON.write_text(
        json.dumps(output, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    fields = [
        "candidate_id",
        "priority",
        "category",
        "doi",
        "year",
        "title",
        "authors",
        "venue",
        "type",
        "landing_url",
        "downloaded_filename",
    ]
    with RESULTS_CSV.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for entry in output:
            selected = entry["selected"]
            writer.writerow(
                {
                    "candidate_id": entry["candidate_id"],
                    "priority": entry["priority"],
                    "category": entry["category"],
                    "doi": selected.get("doi", ""),
                    "year": selected.get("year", ""),
                    "title": selected.get("title", ""),
                    "authors": selected.get("authors", ""),
                    "venue": selected.get("venue", ""),
                    "type": selected.get("type", ""),
                    "landing_url": selected.get("landing_url", ""),
                    "downloaded_filename": entry["downloaded_filename"],
                }
            )
    print(f"wrote {RESULTS_JSON}")
    print(f"wrote {RESULTS_CSV}")
    return 0


def _title_similarity(left: str, right: str) -> float:
    from difflib import SequenceMatcher

    def normalize(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()

    return SequenceMatcher(None, normalize(left), normalize(right)).ratio()


if __name__ == "__main__":
    raise SystemExit(main())
