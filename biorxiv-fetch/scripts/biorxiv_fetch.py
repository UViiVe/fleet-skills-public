#!/usr/bin/env python3
"""
biorxiv_fetch.py — Fetch metadata and PDF from bioRxiv or medRxiv.

Inputs:
    --doi       DOI of the paper (e.g. 10.1101/2024.04.19.590240)
    --search    Search query (e.g. "TP53 cancer")
    --start     Start date for range query (YYYY-MM-DD)
    --end       End date for range query (YYYY-MM-DD)
    --server    biorxiv or medrxiv (default: biorxiv)
    --outdir    Output directory for PDFs (default: /tmp)
    --no-pdf    Skip PDF download, metadata only

Outputs:
    PDF file saved to outdir
    Metadata printed to stdout

Example invocations:
    python3 biorxiv_fetch.py --doi 10.1101/2024.04.19.590240
    python3 biorxiv_fetch.py --search "TP53 cancer" --days 60
    python3 biorxiv_fetch.py --start 2026-05-01 --end 2026-06-11 --no-pdf
"""

import argparse
import os
import re
import sys
import time
from datetime import date, timedelta

try:
    import requests
except ImportError:
    sys.exit("Install with: pip install requests")


PDF_MIN_BYTES = 1024
DEFAULT_USER_AGENT = "Mozilla/5.0 (compatible; biorxiv-fetch/1.0)"


def clean_doi(doi: str) -> str:
    """Strip version suffix and whitespace from DOI."""
    return re.sub(r'v\d+$', '', doi.strip())


def extract_doi_from_url(url: str) -> str:
    """Extract DOI from a biorxiv.org URL."""
    match = re.search(r'(10\.1101/[^\sv?]+)', url)
    if not match:
        raise ValueError(f"Could not extract DOI from URL: {url}")
    return clean_doi(match.group(1))


def get_metadata(doi: str, server: str = "biorxiv") -> dict:
    """Fetch paper metadata from bioRxiv API."""
    doi = clean_doi(doi)
    servers = [server, "medrxiv" if server == "biorxiv" else "biorxiv"]
    for srv in servers:
        url = f"https://api.biorxiv.org/details/{srv}/{doi}/na/json"
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        collection = r.json().get("collection", [])
        if collection:
            paper = collection[-1]
            paper["_server"] = srv
            return paper
    raise ValueError(f"No paper found for DOI: {doi}")


def download_pdf(doi: str, paper: dict, outdir: str) -> str:
    """Download PDF and return local path."""
    doi = clean_doi(doi)
    version = paper.get("version", "1")
    server = paper.get("_server", "biorxiv")
    pdf_url = f"https://www.biorxiv.org/content/{doi}v{version}.full.pdf"

    os.makedirs(outdir, exist_ok=True)
    safe_name = doi.replace("/", "_").replace(".", "-")
    pdf_path = os.path.join(outdir, f"{safe_name}.pdf")

    headers = {"User-Agent": os.environ.get("BIORXIV_FETCH_USER_AGENT", DEFAULT_USER_AGENT)}
    r = requests.get(pdf_url, timeout=60, stream=True, headers=headers)
    r.raise_for_status()
    content_type = r.headers.get("Content-Type", "").lower()
    if content_type and "pdf" not in content_type:
        raise ValueError(f"Expected PDF response, got Content-Type: {content_type}")

    tmp_path = f"{pdf_path}.part"
    with open(tmp_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)

    with open(tmp_path, "rb") as f:
        magic = f.read(5)
    size = os.path.getsize(tmp_path)
    if magic != b"%PDF-" or size < PDF_MIN_BYTES:
        os.remove(tmp_path)
        raise ValueError(
            f"Downloaded content from {pdf_url} did not validate as a PDF "
            f"(magic={magic!r}, size={size} bytes). Try --no-pdf or browser mode."
        )

    os.replace(tmp_path, pdf_path)

    return pdf_path


def search_papers(query: str, days: int = 30, server: str = "biorxiv", max_results: int = 10) -> list:
    """Search for papers matching query in recent days."""
    end = date.today().isoformat()
    start = (date.today() - timedelta(days=days)).isoformat()
    return fetch_range(start, end, server=server, query=query, max_results=max_results)


def fetch_range(start: str, end: str, server: str = "biorxiv", query: str = None, max_results: int = 100) -> list:
    """Fetch all papers in a date range, optionally filtered by query."""
    results = []
    cursor = 0
    query_terms = query.lower().split() if query else []

    while len(results) < max_results:
        url = f"https://api.biorxiv.org/details/{server}/{start}/{end}/{cursor}/json"
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        collection = r.json().get("collection", [])
        if not collection:
            break

        for paper in collection:
            if query_terms:
                text = (paper.get("title", "") + " " + paper.get("abstract", "")).lower()
                if not all(t in text for t in query_terms):
                    continue
            results.append(paper)
            if len(results) >= max_results:
                break

        cursor += 100
        if len(collection) < 100:
            break
        time.sleep(0.2)  # polite rate limiting

    return results


def print_paper(paper: dict, pdf_path: str = None) -> None:
    """Print paper metadata in readable format."""
    print(f"\n{'='*60}")
    print(f"Title:   {paper.get('title', 'N/A')}")
    print(f"Authors: {paper.get('authors', 'N/A')[:120]}")
    print(f"Date:    {paper.get('date', 'N/A')}")
    print(f"DOI:     {paper.get('doi', 'N/A')}")
    print(f"Category:{paper.get('category', 'N/A')}")
    print(f"Version: {paper.get('version', 'N/A')}")
    if pdf_path:
        print(f"PDF:     {pdf_path} ({os.path.getsize(pdf_path)//1024} KB)")
    print(f"\nAbstract:\n{paper.get('abstract', 'N/A')[:600]}")
    print(f"{'='*60}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch bioRxiv/medRxiv papers")
    parser.add_argument("--doi", help="DOI or biorxiv URL")
    parser.add_argument("--search", help="Search query")
    parser.add_argument("--start", help="Start date YYYY-MM-DD (for range)")
    parser.add_argument("--end", help="End date YYYY-MM-DD (for range)")
    parser.add_argument("--days", type=int, default=30, help="Days back for search (default 30)")
    parser.add_argument("--server", default="biorxiv", choices=["biorxiv", "medrxiv"])
    parser.add_argument("--outdir", default="/tmp", help="PDF output directory")
    parser.add_argument("--no-pdf", action="store_true", help="Skip PDF download")
    parser.add_argument("--max", type=int, default=10, help="Max results for search/range")
    args = parser.parse_args()

    if args.doi:
        # Mode A: fetch by DOI
        doi = extract_doi_from_url(args.doi) if "biorxiv.org" in args.doi else args.doi
        paper = get_metadata(doi, server=args.server)
        pdf_path = None
        if not args.no_pdf:
            pdf_path = download_pdf(doi, paper, args.outdir)
        print_paper(paper, pdf_path)
        print(f"✓ biorxiv_fetch completed: 1 paper fetched")

    elif args.search:
        # Mode B: search
        papers = search_papers(args.search, days=args.days, server=args.server, max_results=args.max)
        print(f"✓ Found {len(papers)} papers matching '{args.search}'")
        for i, paper in enumerate(papers, 1):
            print(f"\n[{i}] {paper.get('title','N/A')[:80]}")
            print(f"     {paper.get('doi','N/A')} | {paper.get('date','N/A')}")
            if not args.no_pdf:
                try:
                    pdf_path = download_pdf(paper["doi"], paper, args.outdir)
                    print(f"     PDF: {pdf_path}")
                except Exception as e:
                    print(f"     PDF download failed: {e}")
        print(f"\n✓ biorxiv_fetch completed: {len(papers)} papers")

    elif args.start and args.end:
        # Mode C: date range
        papers = fetch_range(args.start, args.end, server=args.server, max_results=args.max)
        print(f"✓ Found {len(papers)} papers from {args.start} to {args.end}")
        for p in papers:
            print(f"  {p.get('doi','?')} | {p.get('date','?')} | {p.get('title','?')[:70]}")
        print(f"\n✓ biorxiv_fetch completed: {len(papers)} papers fetched")

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
