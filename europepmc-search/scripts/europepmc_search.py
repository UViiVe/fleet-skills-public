#!/usr/bin/env python3
"""
Europe PMC search — covers 32 preprint servers including bioRxiv, medRxiv, SSRN,
ChemRxiv, OSF, Research Square, and more. Also indexes PubMed/MEDLINE.

Requires: pip3 install pyeuropepmc --break-system-packages

Usage:
  python3 europepmc_search.py "VEGF-A skin redness" --max 10
  python3 europepmc_search.py "zebrafish mitochondria" --preprints-only --max 5
  python3 europepmc_search.py "AI peptide design" --days 30 --json
"""

import argparse
import json
import sys
from datetime import datetime, timedelta


def search(query: str, max_results: int = 10, preprints_only: bool = False,
           days: int = None) -> list:
    try:
        from pyeuropepmc import SearchClient
    except ImportError:
        return [{"error": "pyeuropepmc not installed. Run: pip3 install pyeuropepmc --break-system-packages"}]

    client = SearchClient()

    full_query = query
    if preprints_only:
        full_query += " AND src:PPR"
    if days:
        cutoff = (datetime.today() - timedelta(days=days)).strftime("%Y-%m-%d")
        full_query += f" AND (FIRST_PDATE:[{cutoff} TO *])"

    try:
        data = client.search(full_query, page_size=min(max_results, 100))
        results = data.get("resultList", {}).get("result", [])
    except Exception as e:
        return [{"error": str(e)}]

    out = []
    for r in results[:max_results]:
        abstract = r.get("abstractText", "") or ""
        out.append({
            "pmid": r.get("pmid", ""),
            "pmcid": r.get("pmcid", ""),
            "doi": r.get("doi", ""),
            "title": r.get("title", ""),
            "authors": r.get("authorString", ""),
            "journal": r.get("journalTitle", ""),
            "pub_year": r.get("pubYear", ""),
            "pub_date": r.get("firstPublicationDate", ""),
            "source": r.get("source", ""),
            "is_preprint": r.get("source", "") == "PPR",
            "abstract": abstract[:400] + ("..." if len(abstract) > 400 else ""),
            "url": f"https://europepmc.org/article/{r.get('source','')}/{r.get('id','')}",
        })
    return out


def main():
    parser = argparse.ArgumentParser(description="Search Europe PMC (32 preprint servers + PubMed)")
    parser.add_argument("query", help="Search query (Europe PMC query syntax supported)")
    parser.add_argument("--max", type=int, default=10, help="Max results (default: 10)")
    parser.add_argument("--preprints-only", action="store_true", help="Restrict to preprints only")
    parser.add_argument("--days", type=int, default=None, help="Restrict to papers from last N days")
    parser.add_argument("--json", action="store_true", help="Output raw JSON")
    args = parser.parse_args()

    results = search(args.query, max_results=args.max,
                     preprints_only=args.preprints_only, days=args.days)

    if args.json:
        print(json.dumps(results, indent=2))
        return

    if results and "error" in results[0]:
        print(f"Error: {results[0]['error']}")
        sys.exit(1)

    label = "preprints" if args.preprints_only else "results"
    days_label = f" (last {args.days} days)" if args.days else ""
    print(f"Found {len(results)} {label} for '{args.query}'{days_label}\n")
    for i, r in enumerate(results, 1):
        src_tag = f"[{r['source']}]" if r['source'] else ""
        print(f"{i}. {r['title']} {src_tag}")
        print(f"   Authors: {r['authors'][:80]}{'...' if len(r['authors']) > 80 else ''}")
        print(f"   Journal: {r['journal']} | Year: {r['pub_year']}")
        if r['doi']:
            print(f"   DOI: https://doi.org/{r['doi']}")
        print(f"   URL: {r['url']}")
        if r['abstract']:
            print(f"   Abstract: {r['abstract']}")
        print()


if __name__ == "__main__":
    main()
