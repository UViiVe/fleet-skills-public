#!/usr/bin/env python3
"""
bioRxiv/medRxiv search via the official biorxiv.org API.
No API key required. Rate limit: ~3 req/sec.

Usage:
  python3 biorxiv_search.py "VEGF-A peptide" --server biorxiv --days 30 --max 10
  python3 biorxiv_search.py "CRISPR zebrafish" --server medrxiv --days 14
"""

import argparse
import json
import sys
import time
from datetime import datetime, timedelta

import requests

BASE_URL = "https://api.biorxiv.org"


def search_by_interval(server: str, start: str, end: str, max_results: int = 20) -> list:
    """Fetch preprints from biorxiv/medrxiv in a date range."""
    results = []
    cursor = 0
    while len(results) < max_results:
        url = f"{BASE_URL}/details/{server}/{start}/{end}/{cursor}/json"
        resp = requests.get(url, timeout=15)
        if resp.status_code == 429:
            time.sleep(2)
            continue
        resp.raise_for_status()
        data = resp.json()
        collection = data.get("collection", [])
        if not collection:
            break
        results.extend(collection)
        if len(collection) < 100:
            break
        cursor += 100
        time.sleep(0.4)
    return results[:max_results]


def filter_by_query(papers: list, query: str) -> list:
    """Client-side filter: keep papers whose title/abstract contain query terms."""
    terms = [t.lower() for t in query.split()]
    out = []
    for p in papers:
        text = (p.get("title", "") + " " + p.get("abstract", "")).lower()
        if all(t in text for t in terms):
            out.append(p)
    return out


def format_paper(p: dict) -> dict:
    return {
        "doi": p.get("doi", ""),
        "title": p.get("title", ""),
        "authors": p.get("authors", ""),
        "date": p.get("date", ""),
        "category": p.get("category", ""),
        "abstract": p.get("abstract", "")[:400] + ("..." if len(p.get("abstract", "")) > 400 else ""),
        "url": f"https://www.biorxiv.org/content/{p.get('doi', '')}v{p.get('version', '1')}",
        "server": p.get("server", ""),
    }


def main():
    parser = argparse.ArgumentParser(description="Search bioRxiv/medRxiv preprints")
    parser.add_argument("query", help="Search terms (space-separated, all must match)")
    parser.add_argument("--server", choices=["biorxiv", "medrxiv"], default="biorxiv")
    parser.add_argument("--days", type=int, default=30, help="How many days back to search (default: 30)")
    parser.add_argument("--max", type=int, default=10, help="Max results to return (default: 10)")
    parser.add_argument("--json", action="store_true", help="Output raw JSON")
    args = parser.parse_args()

    end_date = datetime.today().strftime("%Y-%m-%d")
    start_date = (datetime.today() - timedelta(days=args.days)).strftime("%Y-%m-%d")

    try:
        papers = search_by_interval(args.server, start_date, end_date, max_results=500)
        filtered = filter_by_query(papers, args.query)
        formatted = [format_paper(p) for p in filtered[:args.max]]
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)

    if args.json:
        print(json.dumps(formatted, indent=2))
    else:
        print(f"Found {len(formatted)} results for '{args.query}' on {args.server} (last {args.days} days)\n")
        for i, p in enumerate(formatted, 1):
            print(f"{i}. {p['title']}")
            print(f"   Authors: {p['authors'][:80]}{'...' if len(p['authors']) > 80 else ''}")
            print(f"   Date: {p['date']} | Category: {p['category']}")
            print(f"   URL: {p['url']}")
            print(f"   Abstract: {p['abstract']}\n")


if __name__ == "__main__":
    main()
