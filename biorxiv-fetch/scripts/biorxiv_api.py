#!/usr/bin/env python3
"""
biorxiv_api.py — bioRxiv/medRxiv METADATA helpers (Mode A/B/C).
For PDF download see Mode D in SKILL.md (browser automation — NOT in this file,
because PDFs require an in-browser fetch to clear Cloudflare; pure-python cannot
get the PDF — it always returns ~5KB Cloudflare HTML).
"""
import requests, re
from datetime import date, timedelta

API = "https://api.biorxiv.org/details"


def fetch_biorxiv_meta(doi: str) -> dict:
    """Mode A — metadata + abstract by DOI. Always works (API host, not Cloudflare-walled)."""
    doi_clean = re.sub(r'v\d+$', '', doi.strip())
    for server in ("biorxiv", "medrxiv"):
        r = requests.get(f"{API}/{server}/{doi_clean}/na/json", timeout=15)
        collection = r.json().get("collection", [])
        if collection:
            p = collection[-1]
            v = p.get("version", "1")
            return {
                "title": p.get("title"), "authors": p.get("authors"),
                "abstract": p.get("abstract"), "doi": doi_clean, "date": p.get("date"),
                "category": p.get("category"), "version": v, "server": server,
                "html_url": f"https://www.biorxiv.org/content/{doi_clean}v{v}.full",
                "pdf_url": f"https://www.biorxiv.org/content/{doi_clean}v{v}.full.pdf",
            }
    raise ValueError(f"No paper found for DOI: {doi_clean}")


def search_biorxiv(query: str, days: int = 30, server: str = "biorxiv", max_results: int = 10) -> list:
    """Mode B — keyword search over a recent date window (client-side AND filter)."""
    end = date.today().isoformat()
    start = (date.today() - timedelta(days=days)).isoformat()
    terms = [t.lower() for t in query.split()]
    out, cursor = [], 0
    while len(out) < max_results:
        r = requests.get(f"{API}/{server}/{start}/{end}/{cursor}/json", timeout=15)
        r.raise_for_status()
        coll = r.json().get("collection", [])
        if not coll:
            break
        for p in coll:
            text = (p.get("title", "") + " " + p.get("abstract", "")).lower()
            if all(t in text for t in terms):
                out.append(p)
                if len(out) >= max_results:
                    break
        cursor += 100
        if len(coll) < 100:
            break
    return out


def fetch_biorxiv_range(start: str, end: str, server: str = "biorxiv") -> list:
    """Mode C — all metadata records in a YYYY-MM-DD/YYYY-MM-DD window."""
    out, cursor = [], 0
    while True:
        r = requests.get(f"{API}/{server}/{start}/{end}/{cursor}/json", timeout=15)
        r.raise_for_status()
        coll = r.json().get("collection", [])
        if not coll:
            break
        out.extend(coll)
        cursor += 100
        if len(coll) < 100:
            break
    return out


if __name__ == "__main__":
    import sys, json
    if len(sys.argv) > 1:
        print(json.dumps(fetch_biorxiv_meta(sys.argv[1]), indent=2)[:1500])
    else:
        print("usage: biorxiv_api.py <DOI>   # e.g. 10.1101/2025.08.14.670059")
