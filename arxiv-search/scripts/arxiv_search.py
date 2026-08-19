#!/usr/bin/env python3
"""
arXiv search script for the fleet.
Usage: python3 arxiv_search.py <query> [--max N] [--sort relevance|date|updated] [--category cs.AI]

Returns JSON array of paper results.
"""

from __future__ import annotations

import argparse
import json
import sys

def search(query: str, max_results: int = 10, sort_by: str = "relevance", category: str | None = None) -> list[dict]:
    try:
        import arxiv
    except ImportError:
        print(json.dumps({"error": "arxiv package not installed. Run: pip3 install arxiv --break-system-packages"}))
        sys.exit(1)

    sort_map = {
        "relevance": arxiv.SortCriterion.Relevance,
        "date": arxiv.SortCriterion.SubmittedDate,
        "updated": arxiv.SortCriterion.LastUpdatedDate,
    }

    # Optionally scope to a category (e.g., cs.AI, q-bio.GN, physics.bio-ph)
    full_query = f"cat:{category} AND ({query})" if category else query

    client = arxiv.Client()
    search_obj = arxiv.Search(
        query=full_query,
        max_results=min(max_results, 50),
        sort_by=sort_map.get(sort_by, arxiv.SortCriterion.Relevance),
    )

    results = []
    for r in client.results(search_obj):
        results.append({
            "arxiv_id": r.entry_id.split("/abs/")[-1],
            "title": r.title,
            "authors": [a.name for a in r.authors[:6]],
            "abstract": r.summary[:500] + ("..." if len(r.summary) > 500 else ""),
            "published": r.published.isoformat() if r.published else "",
            "updated": r.updated.isoformat() if r.updated else "",
            "pdf_url": r.pdf_url or "",
            "categories": r.categories[:5],
            "journal_ref": r.journal_ref or "",
            "doi": r.doi or "",
        })
    return results


def main():
    parser = argparse.ArgumentParser(description="Search arXiv")
    parser.add_argument("query", help="Search query")
    parser.add_argument("--max", type=int, default=10, help="Max results (default 10, max 50)")
    parser.add_argument("--sort", choices=["relevance", "date", "updated"], default="relevance")
    parser.add_argument("--category", default=None, help="arXiv category filter (e.g. cs.AI, q-bio.GN)")
    args = parser.parse_args()

    results = search(args.query, max_results=args.max, sort_by=args.sort, category=args.category)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
