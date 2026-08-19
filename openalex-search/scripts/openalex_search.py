#!/usr/bin/env python3
"""
openalex_search.py — Search OpenAlex for scholarly works including preprints.
No API key required (polite pool). 450M+ works. Semantic + keyword search.
Covers arXiv, bioRxiv, medRxiv, Zenodo, institutional repos, journals, and more.

Usage:
  python3 openalex_search.py "cystic fibrosis CFTR" --preprints --days 30
  python3 openalex_search.py "VEGF-A peptide binding" --days 14 --max 20
  python3 openalex_search.py "zebrafish mitochondria CRISPR" --from 2026-01-01 --max 10
  python3 openalex_search.py "UV skin redness VEGF" --open-access --summary
  python3 openalex_search.py --doi 10.1101/2024.01.23.576694

API docs: https://docs.openalex.org
"""

import argparse
import json
import sys
import urllib.request
import urllib.parse
from datetime import datetime, timedelta

BASE_URL = "https://api.openalex.org/works"
MAILTO = "atlas@uviive.com"  # polite pool — increases rate limits

# OpenAlex work types
WORK_TYPES = [
    "article", "preprint", "book-chapter", "book", "dataset",
    "dissertation", "editorial", "erratum", "grant", "letter",
    "paratext", "peer-review", "reference-entry", "report",
    "review", "standard", "supplementary-materials",
]

def fetch_json(url):
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": f"atlas-fleet-preprint-scout/1.0 (mailto:{MAILTO})"
        })
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        print(json.dumps({"error": str(e), "url": url}))
        sys.exit(1)

def build_filter(date_from=None, date_to=None, work_type=None,
                 preprints=False, open_access=False, concepts=None):
    filters = []
    if preprints:
        filters.append("type:preprint")
    elif work_type:
        filters.append(f"type:{work_type}")
    if open_access:
        filters.append("open_access.is_oa:true")
    if date_from:
        filters.append(f"from_publication_date:{date_from}")
    if date_to:
        filters.append(f"to_publication_date:{date_to}")
    if concepts:
        # Concept IDs or display names
        filters.append(f"concepts.display_name.search:{concepts}")
    return ",".join(filters) if filters else None

def build_url(query, filter_str, page, per_page, sort):
    params = {
        "per-page": per_page,
        "page": page,
        "mailto": MAILTO,
    }
    if query:
        params["search"] = query
    if filter_str:
        params["filter"] = filter_str
    if sort:
        params["sort"] = sort
    return BASE_URL + "?" + urllib.parse.urlencode(params)

def search(query, filter_str, max_results=25, sort="publication_date:desc"):
    results = []
    page = 1
    total_count = 0

    while len(results) < max_results:
        per_page = min(50, max_results - len(results))
        url = build_url(query, filter_str, page, per_page, sort)
        data = fetch_json(url)
        meta = data.get("meta", {})
        total_count = meta.get("count", 0)
        batch = data.get("results", [])
        if not batch:
            break
        results.extend(batch)
        # Check if we've exhausted results
        if len(results) >= total_count:
            break
        page += 1

    return results[:max_results], total_count

def get_best_oa_url(work):
    """Extract the best available open access URL."""
    oa = work.get("open_access", {})
    if oa.get("oa_url"):
        return oa["oa_url"]
    # Try primary location
    primary = work.get("primary_location", {})
    if primary.get("pdf_url"):
        return primary["pdf_url"]
    if primary.get("landing_page_url"):
        return primary["landing_page_url"]
    return None

def get_source_name(work):
    """Get the journal/server/repo name."""
    primary = work.get("primary_location", {})
    source = primary.get("source") or {}
    return source.get("display_name", "")

def format_paper(w, index=None):
    doi = w.get("doi", "")
    if doi:
        doi = doi.replace("https://doi.org/", "")
    
    # Authors (first 5)
    authorships = w.get("authorships", [])
    authors = "; ".join(
        a.get("author", {}).get("display_name", "")
        for a in authorships[:5]
    )
    if len(authorships) > 5:
        authors += f" +{len(authorships) - 5} more"

    # Concepts (top 5)
    concepts = [c.get("display_name", "") for c in w.get("concepts", [])[:5]]

    # Abstract (inverted index → reconstruct)
    abstract = ""
    inv_idx = w.get("abstract_inverted_index")
    if inv_idx:
        try:
            positions = {}
            for word, pos_list in inv_idx.items():
                for pos in pos_list:
                    positions[pos] = word
            abstract = " ".join(positions[i] for i in sorted(positions.keys()))
            if len(abstract) > 500:
                abstract = abstract[:500] + "..."
        except Exception:
            abstract = ""

    oa_url = get_best_oa_url(w)
    work_url = w.get("id", "").replace("https://openalex.org/", "https://openalex.org/works/")

    return {
        "index": index,
        "title": w.get("title", ""),
        "authors": authors,
        "date": w.get("publication_date", ""),
        "year": w.get("publication_year", ""),
        "type": w.get("type", ""),
        "source": get_source_name(w),
        "doi": doi,
        "url": oa_url or (f"https://doi.org/{doi}" if doi else work_url),
        "openalex_id": w.get("id", ""),
        "abstract": abstract,
        "cited_by": w.get("cited_by_count", 0),
        "open_access": w.get("open_access", {}).get("is_oa", False),
        "oa_status": w.get("open_access", {}).get("oa_status", ""),
        "concepts": concepts,
        "language": w.get("language", ""),
    }

def main():
    parser = argparse.ArgumentParser(
        description="Search OpenAlex — 450M+ scholarly works including preprints"
    )
    parser.add_argument("keywords", nargs="*",
                        help="Search terms (semantic + keyword search)")
    parser.add_argument("--days", type=int, default=None,
                        help="Search last N days")
    parser.add_argument("--from", dest="date_from", help="Start date YYYY-MM-DD")
    parser.add_argument("--to", dest="date_to", help="End date YYYY-MM-DD")
    parser.add_argument("--max", type=int, default=25,
                        help="Max results (default: 25)")
    parser.add_argument("--sort", default="publication_date:desc",
                        help="Sort field:dir (default: publication_date:desc)")
    parser.add_argument("--preprints", action="store_true",
                        help="Preprints only (type:preprint)")
    parser.add_argument("--type", dest="work_type",
                        help=f"Work type filter. Options: {', '.join(WORK_TYPES)}")
    parser.add_argument("--open-access", action="store_true",
                        help="Open access only")
    parser.add_argument("--doi", help="Fetch specific work by DOI")
    parser.add_argument("--summary", action="store_true",
                        help="Print summary stats only")
    parser.add_argument("--list-types", action="store_true",
                        help="List work type options")
    args = parser.parse_args()

    if args.list_types:
        print(json.dumps({"work_types": WORK_TYPES}))
        return

    if args.doi:
        doi_clean = args.doi.replace("https://doi.org/", "")
        url = f"https://api.openalex.org/works/https://doi.org/{doi_clean}?mailto={MAILTO}"
        data = fetch_json(url)
        formatted = format_paper(data, 1)
        print(json.dumps({"count": 1, "results": [formatted]}, indent=2))
        return

    date_to = args.date_to
    date_from = args.date_from
    if args.days and not date_from:
        date_from = str((datetime.now() - timedelta(days=args.days)).date())
        date_to = date_to or str(datetime.now().date())

    filter_str = build_filter(
        date_from=date_from,
        date_to=date_to,
        work_type=args.work_type,
        preprints=args.preprints,
        open_access=args.open_access,
    )

    query = " ".join(args.keywords) if args.keywords else None
    results, total_count = search(query, filter_str, max_results=args.max, sort=args.sort)
    formatted = [format_paper(w, i+1) for i, w in enumerate(results)]

    # Count by type and source
    type_counts = {}
    source_counts = {}
    for r in formatted:
        t = r["type"]
        type_counts[t] = type_counts.get(t, 0) + 1
        s = r["source"]
        if s:
            source_counts[s] = source_counts.get(s, 0) + 1

    if args.summary:
        print(json.dumps({
            "query": query,
            "filter": filter_str,
            "date_range": f"{date_from or 'any'} to {date_to or 'today'}",
            "total_in_openalex": total_count,
            "returned": len(formatted),
            "by_type": type_counts,
            "top_sources": dict(sorted(source_counts.items(), key=lambda x: -x[1])[:10]),
        }, indent=2))
    else:
        print(json.dumps({
            "query": query,
            "filter": filter_str,
            "date_range": f"{date_from or 'any'} to {date_to or 'today'}",
            "total_in_openalex": total_count,
            "returned": len(formatted),
            "by_type": type_counts,
            "results": formatted,
        }, indent=2))

if __name__ == "__main__":
    main()
