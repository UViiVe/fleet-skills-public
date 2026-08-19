#!/usr/bin/env python3
"""
Waterfall abstract fetcher: PubMed → Europe PMC → Semantic Scholar → Crossref

Maximizes abstract coverage for biomedical literature by trying multiple sources.
"""

import json
import sys
import time
import urllib.request
import urllib.parse
from typing import Optional, Dict, Any, List

class AbstractFetcher:
    def __init__(self, enrich_pmids=False):
        self.enrich_pmids = enrich_pmids
        self.stats = {
            'total': 0,
            'pubmed': 0,
            'europepmc': 0,
            'semantic_scholar': 0,
            'crossref': 0,
            'failed': 0,
            'pmids_enriched': 0
        }
    
    def enrich_pmid(self, paper: Dict[str, Any]) -> Optional[str]:
        """Try to find PMID from DOI or title using PubMed search"""
        if not self.enrich_pmids:
            return None
        
        try:
            query = None
            doi = (paper.get('doi') or '').replace('https://doi.org/', '')
            title = paper.get('title', '')
            
            if doi:
                query = f'"{doi}"[DOI]'
            elif title:
                # Clean title for search
                query = f'"{title}"[Title]'
            else:
                return None
            
            encoded_query = urllib.parse.quote(query)
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={encoded_query}&retmode=json"
            with urllib.request.urlopen(url, timeout=10) as response:
                data = json.loads(response.read().decode('utf-8'))
                pmids = data.get('esearchresult', {}).get('idlist', [])
                if pmids:
                    self.stats['pmids_enriched'] += 1
                    return pmids[0]
        except Exception:
            pass
        return None
    
    def fetch_from_pubmed(self, pmid: str) -> Optional[str]:
        """Fetch abstract from PubMed E-utilities using PMID"""
        try:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml&rettype=abstract"
            with urllib.request.urlopen(url, timeout=10) as response:
                xml_data = response.read().decode('utf-8')
                # Simple XML parsing for AbstractText
                if '<AbstractText>' in xml_data:
                    start = xml_data.find('<AbstractText>') + 14
                    end = xml_data.find('</AbstractText>', start)
                    abstract = xml_data[start:end]
                    # Clean up XML entities
                    abstract = abstract.replace('&lt;', '<').replace('&gt;', '>').replace('&amp;', '&')
                    if len(abstract) > 50:
                        self.stats['pubmed'] += 1
                        return abstract
        except Exception:
            pass
        return None
    
    def fetch_from_europepmc(self, query: str) -> Optional[str]:
        """Fetch abstract from Europe PMC"""
        try:
            # Try by PMID or DOI
            encoded_query = urllib.parse.quote(query)
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={encoded_query}&format=json&resulttype=core"
            with urllib.request.urlopen(url, timeout=10) as response:
                data = json.loads(response.read().decode('utf-8'))
                if data.get('resultList', {}).get('result'):
                    result = data['resultList']['result'][0]
                    abstract = result.get('abstractText', '')
                    if len(abstract) > 50:
                        self.stats['europepmc'] += 1
                        return abstract
        except Exception:
            pass
        return None
    
    def fetch_from_semantic_scholar(self, doi: Optional[str] = None, title: Optional[str] = None) -> Optional[str]:
        """Fetch abstract from Semantic Scholar API"""
        try:
            if doi:
                url = f"https://api.semanticscholar.org/graph/v1/paper/DOI:{doi}?fields=abstract"
            elif title:
                # Search by title
                encoded_title = urllib.parse.quote(title)
                url = f"https://api.semanticscholar.org/graph/v1/paper/search?query={encoded_title}&fields=abstract&limit=1"
            else:
                return None
                
            with urllib.request.urlopen(url, timeout=10) as response:
                data = json.loads(response.read().decode('utf-8'))
                if doi:
                    abstract = data.get('abstract', '')
                else:
                    # Search result
                    papers = data.get('data', [])
                    abstract = papers[0].get('abstract', '') if papers else ''
                
                if abstract and len(abstract) > 50:
                    self.stats['semantic_scholar'] += 1
                    return abstract
        except Exception:
            pass
        return None
    
    def fetch_from_crossref(self, doi: str) -> Optional[str]:
        """Fetch abstract from Crossref API"""
        try:
            url = f"https://api.crossref.org/works/{doi}"
            with urllib.request.urlopen(url, timeout=10) as response:
                data = json.loads(response.read().decode('utf-8'))
                abstract = data.get('message', {}).get('abstract', '')
                if abstract and len(abstract) > 50:
                    self.stats['crossref'] += 1
                    return abstract
        except Exception:
            pass
        return None
    
    def waterfall_fetch(self, paper: Dict[str, Any]) -> Optional[str]:
        """
        Try all sources in waterfall order until we find an abstract.
        
        Args:
            paper: dict with keys: title, doi (optional), pmid (optional), year, authors
        
        Returns:
            Abstract text or None
        """
        self.stats['total'] += 1
        doi = (paper.get('doi') or '').replace('https://doi.org/', '')
        pmid = paper.get('pmid')
        title = paper.get('title', '')
        
        # 0. Try to enrich PMID if missing
        if not pmid and self.enrich_pmids:
            pmid = self.enrich_pmid(paper)
            if pmid:
                paper['pmid'] = pmid  # Cache for future use
            time.sleep(0.1)  # Rate limiting for PubMed search
        
        # 1. Try PubMed (fastest, best for biomedical)
        if pmid:
            abstract = self.fetch_from_pubmed(pmid)
            if abstract:
                return abstract
            time.sleep(0.1)  # Rate limiting
        
        # 2. Try Europe PMC (broader coverage)
        if pmid:
            abstract = self.fetch_from_europepmc(f"ext_id:{pmid} src:med")
            if abstract:
                return abstract
        elif doi:
            abstract = self.fetch_from_europepmc(f"doi:{doi}")
            if abstract:
                return abstract
        time.sleep(0.1)
        
        # 3. Try Semantic Scholar (cross-disciplinary)
        if doi:
            abstract = self.fetch_from_semantic_scholar(doi=doi)
            if abstract:
                return abstract
        else:
            abstract = self.fetch_from_semantic_scholar(title=title)
            if abstract:
                return abstract
        time.sleep(0.2)  # Semantic Scholar rate limit: 100 req/5min
        
        # 4. Last resort: Crossref
        if doi:
            abstract = self.fetch_from_crossref(doi)
            if abstract:
                return abstract
        
        self.stats['failed'] += 1
        return None
    
    def fetch_all(self, papers: List[Dict[str, Any]], progress: bool = True) -> List[Dict[str, Any]]:
        """
        Fetch abstracts for a list of papers.
        
        Returns:
            Updated papers list with 'abstract' field added
        """
        results = []
        for i, paper in enumerate(papers):
            if progress and i % 10 == 0:
                print(f"Progress: {i}/{len(papers)}", file=sys.stderr)
            
            abstract = paper.get('abstract')
            if not abstract or len(abstract) < 50:
                abstract = self.waterfall_fetch(paper)
            
            paper_copy = paper.copy()
            paper_copy['abstract'] = abstract if abstract else None
            results.append(paper_copy)
        
        return results
    
    def print_stats(self):
        """Print success rate statistics"""
        print("\n=== Abstract Retrieval Stats ===", file=sys.stderr)
        print(f"Total papers: {self.stats['total']}", file=sys.stderr)
        if self.enrich_pmids:
            print(f"  PMIDs enriched: {self.stats['pmids_enriched']}", file=sys.stderr)
        print(f"  PubMed: {self.stats['pubmed']}", file=sys.stderr)
        print(f"  Europe PMC: {self.stats['europepmc']}", file=sys.stderr)
        print(f"  Semantic Scholar: {self.stats['semantic_scholar']}", file=sys.stderr)
        print(f"  Crossref: {self.stats['crossref']}", file=sys.stderr)
        success = self.stats['total'] - self.stats['failed']
        rate = 100.0 * success / self.stats['total'] if self.stats['total'] > 0 else 0
        print(f"\nSuccess rate: {success}/{self.stats['total']} ({rate:.1f}%)", file=sys.stderr)


def main():
    if len(sys.argv) < 2:
        print("Usage: fetch_abstracts.py <papers.json> [output.json] [--enrich-pmids]", file=sys.stderr)
        print("\nExpected input format:", file=sys.stderr)
        print('[{"title": "...", "doi": "...", "pmid": "...", ...}, ...]', file=sys.stderr)
        print("\nOptions:", file=sys.stderr)
        print("  --enrich-pmids    Search PubMed for missing PMIDs (slower but better coverage)", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith('--') else None
    enrich_pmids = '--enrich-pmids' in sys.argv
    
    with open(input_file, 'r') as f:
        papers = json.load(f)
    
    fetcher = AbstractFetcher(enrich_pmids=enrich_pmids)
    results = fetcher.fetch_all(papers, progress=True)
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults written to: {output_file}", file=sys.stderr)
    else:
        print(json.dumps(results, indent=2))
    
    fetcher.print_stats()


if __name__ == '__main__':
    main()
