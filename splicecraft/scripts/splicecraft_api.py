"""
splicecraft_api.py - Helper for driving SpliceCraft via its local agent API.

Source: SpliceCraft docs/agent-api.md + docs/cli.md
        https://github.com/Binomica-Labs/SpliceCraft

Security note: Uses the HTTP API exclusively - never imports splicecraft directly,
never calls _save_* functions. This is the sanctioned path for programmatic access.

Data-dir safety: The API server is one of only four sanctioned callers for all
save operations. This client is read-only by default; write operations require
explicit opt-in at client construction time.

Usage:
    from splicecraft_api import SpliceCraftClient
    sc = SpliceCraftClient()
    sc.fetch_ncbi("L09137")
    features = sc.list_features()
"""

import json
import os
import platform
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Token file location (platform-specific)
# ---------------------------------------------------------------------------

def _default_data_dir() -> Path:
    """Return the default SpliceCraft data directory for the current platform."""
    override = os.environ.get("SPLICECRAFT_DATA_DIR")
    if override:
        return Path(override)
    system = platform.system()
    if system == "Darwin":
        return Path.home() / "Library" / "Application Support" / "splicecraft"
    elif system == "Windows":
        appdata = os.environ.get("APPDATA", "")
        return Path(appdata) / "splicecraft"
    else:  # Linux, WSL2, etc.
        xdg = os.environ.get("XDG_DATA_HOME", str(Path.home() / ".local" / "share"))
        return Path(xdg) / "splicecraft"


def _read_token_file(data_dir: Path) -> tuple[int, str]:
    """
    Read port and bearer token from the agent_token file.
    Returns (port, token).
    Raises FileNotFoundError if the server hasn't been started yet.
    Raises ValueError if the token file is malformed.
    """
    token_path = data_dir / "agent_token"
    if not token_path.exists():
        raise FileNotFoundError(
            f"SpliceCraft agent token file not found at {token_path}.\n"
            "Start the agent server first: splicecraft --agent"
        )
    content = token_path.read_text().strip().splitlines()
    if len(content) < 2:
        raise ValueError(
            f"Malformed token file at {token_path}: expected 2 lines (port, token)."
        )
    try:
        port = int(content[0])
    except ValueError:
        raise ValueError(f"Non-integer port in token file: {content[0]!r}")
    token = content[1].strip()
    if not token:
        raise ValueError("Empty bearer token in token file.")
    return port, token


# ---------------------------------------------------------------------------
# Base client
# ---------------------------------------------------------------------------

class SpliceCraftClient:
    """
    HTTP client for the SpliceCraft local agent API.

    Automatically discovers port + bearer token from the platform token file.
    Re-reads the token file on each request (token rotates each server launch).

    Parameters
    ----------
    data_dir : Path or None
        Override the data directory. Default: platform-standard location or
        $SPLICECRAFT_DATA_DIR.
    """

    def __init__(self, data_dir: Path | None = None, allow_write: bool = False):
        self._data_dir = Path(data_dir) if data_dir else _default_data_dir()
        self._allow_write = allow_write

    def _require_write_enabled(self, action: str) -> None:
        if not self._allow_write:
            raise PermissionError(
                f"{action} requires SpliceCraftClient(allow_write=True). "
                "Review the active record and library backup before enabling writes."
            )

    # ------------------------------------------------------------------
    # Low-level HTTP
    # ------------------------------------------------------------------

    def _auth_headers(self, token: str) -> dict[str, str]:
        return {
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
        }

    def get(self, endpoint: str, **params) -> Any:
        """
        GET request to the agent API.

        Parameters
        ----------
        endpoint : str
            e.g. "/status", "/list-features"
        params : dict
            Query string parameters (unused by most GET endpoints; included for completeness).

        Returns
        -------
        Parsed JSON response (dict, list, or primitive).
        """
        port, token = _read_token_file(self._data_dir)
        url = f"http://127.0.0.1:{port}{endpoint}"
        req = urllib.request.Request(url, headers=self._auth_headers(token))
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            body_text = e.read().decode(errors="replace")
            raise RuntimeError(
                f"SpliceCraft API error {e.code} on {endpoint}: {body_text}"
            ) from e
        except urllib.error.URLError as e:
            if "Connection refused" in str(e):
                raise ConnectionRefusedError(
                    "Cannot connect to SpliceCraft agent server. "
                    "Start it with: splicecraft --agent"
                ) from e
            raise

    def post(self, endpoint: str, body: dict | None = None) -> Any:
        """
        POST request to the agent API.

        Parameters
        ----------
        endpoint : str
            e.g. "/fetch", "/design-primers"
        body : dict or None
            JSON request body.

        Returns
        -------
        Parsed JSON response.
        """
        port, token = _read_token_file(self._data_dir)
        url = f"http://127.0.0.1:{port}{endpoint}"
        data = json.dumps(body or {}).encode()
        req = urllib.request.Request(
            url,
            data=data,
            headers=self._auth_headers(token),
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            body_text = e.read().decode(errors="replace")
            raise RuntimeError(
                f"SpliceCraft API error {e.code} on {endpoint}: {body_text}"
            ) from e
        except urllib.error.URLError as e:
            if "Connection refused" in str(e):
                raise ConnectionRefusedError(
                    "Cannot connect to SpliceCraft agent server. "
                    "Start it with: splicecraft --agent"
                ) from e
            raise

    # ------------------------------------------------------------------
    # Discovery
    # ------------------------------------------------------------------

    def list_tools(self) -> list[dict]:
        """Return the live endpoint inventory from /tools."""
        return self.get("/tools")

    def status(self) -> dict:
        """Return the current record snapshot (name, length, accession, etc.)."""
        return self.get("/status")

    # ------------------------------------------------------------------
    # Records: load / fetch
    # ------------------------------------------------------------------

    def fetch_ncbi(self, accession: str) -> dict:
        """
        Fetch a plasmid from NCBI by accession and load it into the active session.

        Parameters
        ----------
        accession : str
            e.g. "L09137" (pUC19), "MW463917.1" (pACYC184)
        """
        return self.post("/fetch", {"accession": accession})

    def load_file(self, path: str | Path) -> dict:
        """
        Load a local plasmid file (.gb, .gbk, .dna, .embl, FASTA, .ab1).

        The path must be accessible from the machine running the SpliceCraft server.
        """
        return self.post("/load", {"path": str(path)})

    def load_library_entry(self, name: str, force: bool = False) -> dict:
        """Load a library entry by name into the active session."""
        return self.post("/load-entry", {"name": name, "force": force})

    # ------------------------------------------------------------------
    # Records: sequence + features
    # ------------------------------------------------------------------

    def get_sequence(self, start: int, end: int, bottom_strand: bool = False) -> dict:
        """
        Extract a DNA sequence range from the loaded record.

        Parameters
        ----------
        start, end : int
            1-based bp range (inclusive).
        bottom_strand : bool
            If True, return the reverse complement.
        """
        return self.post("/get-sequence", {
            "start": start,
            "end": end,
            "bottom": bottom_strand,
        })

    def list_features(self) -> list[dict]:
        """Return all features on the loaded record."""
        return self.get("/list-features")

    def add_feature(
        self,
        start: int,
        end: int,
        label: str,
        feature_type: str = "misc_feature",
        strand: int = 1,
    ) -> dict:
        """Add a feature annotation to the loaded record."""
        self._require_write_enabled("add_feature")
        return self.post("/add-feature", {
            "start": start,
            "end": end,
            "label": label,
            "type": feature_type,
            "strand": strand,
        })

    def delete_feature(self, index: int) -> dict:
        """Delete the feature at the given 0-based index."""
        self._require_write_enabled("delete_feature")
        return self.post("/delete-feature", {"index": index})

    def find_orfs(self, min_length: int = 300) -> list[dict]:
        """Find open reading frames in the loaded record."""
        return self.post("/find-orfs", {"min_length": min_length})

    # ------------------------------------------------------------------
    # Records: save / export
    # ------------------------------------------------------------------

    def save(self) -> dict:
        """Save the loaded record to the library."""
        self._require_write_enabled("save")
        return self.post("/save")

    def export_genbank(self, output_path: str | Path) -> dict:
        """
        Export the loaded record as GenBank format.

        Parameters
        ----------
        output_path : str or Path
            Server-side path to write the .gb file.
        """
        return self.post("/export-genbank", {"path": str(output_path)})

    def export_fasta(self, output_path: str | Path) -> dict:
        """Export the loaded record sequence as FASTA."""
        return self.post("/export-fasta", {"path": str(output_path)})

    def export_gff3(self, output_path: str | Path) -> dict:
        """Export the loaded record annotations as GFF3."""
        return self.post("/export-gff3", {"path": str(output_path)})

    # ------------------------------------------------------------------
    # Library
    # ------------------------------------------------------------------

    def list_library(self) -> list[dict]:
        """List all saved plasmid library entries."""
        return self.get("/list-library")

    def search_library(self, query: str) -> list[dict]:
        """Search across all collections for plasmids matching the query."""
        return self.post("/search-library", {"query": query})

    def list_collections(self) -> list[dict]:
        """List collections and the active one."""
        return self.get("/list-collections")

    def set_active_collection(self, name: str) -> dict:
        """Switch the active collection. Required before cross-collection lookups."""
        return self.post("/set-active-collection", {"name": name})

    # ------------------------------------------------------------------
    # Restriction enzyme analysis
    # ------------------------------------------------------------------

    def restriction_scan(
        self,
        filter_mode: str = "all",  # "all" | "unique" | "6bp" | "4bp" | "goldenbraid"
        collection: str | None = None,
    ) -> list[dict]:
        """
        Scan the loaded record for restriction enzyme cut sites.

        Parameters
        ----------
        filter_mode : str
            "all" — all enzymes
            "unique" — single-cutters only
            "6bp" — 6+ bp recognition sequences
            "4bp" — 4+ bp recognition sequences
            "goldenbraid" — Golden Braid connectors only
        collection : str or None
            Use a named enzyme collection (if configured).

        Returns
        -------
        List of restriction site dicts with position, enzyme name, cut positions.

        [VERIFY: exact parameter names on /list-restriction-sites endpoint —
         confirm via `splicecraft-cli tools` on a running instance]
        """
        body: dict = {"filter": filter_mode}
        if collection:
            body["collection"] = collection
        return self.post("/list-restriction-sites", body)

    # ------------------------------------------------------------------
    # Primer design
    # ------------------------------------------------------------------

    def design_primers(
        self,
        task: str = "detection",  # "detection" | "cloning" | "goldenbraid" | "generic"
        target_start: int | None = None,
        target_end: int | None = None,
        **kwargs,
    ) -> dict:
        """
        Design primers using Primer3.

        Parameters
        ----------
        task : str
            Primer type. One of: "detection", "cloning", "goldenbraid", "generic"
        target_start, target_end : int
            Target region (1-based bp).
        **kwargs
            Additional Primer3 parameters passed to the API.
            [VERIFY: exact parameter names — confirm via `splicecraft-cli tools`]

        Returns
        -------
        dict with designed primer sequences, Tm, GC%, etc.
        """
        body: dict = {"task": task}
        if target_start is not None:
            body["target_start"] = target_start
        if target_end is not None:
            body["target_end"] = target_end
        body.update(kwargs)
        return self.post("/design-primers", body)

    def design_mutagenesis(
        self,
        feature_index: int,
        mutation: str,
    ) -> dict:
        """
        Design SOE-PCR mutagenesis primers for a site-directed mutation.

        Parameters
        ----------
        feature_index : int
            Index of the CDS feature to mutate (from list_features()).
        mutation : str
            Mutation descriptor, e.g. "L54A" (amino acid change).

        Returns
        -------
        dict with primer sequences for the SOE-PCR strategy.

        [VERIFY: exact parameter names on /design-mutagenesis — confirm via `splicecraft-cli tools`]
        """
        return self.post("/design-mutagenesis", {
            "feature_index": feature_index,
            "mutation": mutation,
        })

    # ------------------------------------------------------------------
    # Assembly simulation
    # ------------------------------------------------------------------

    def simulate_assembly(
        self,
        grammar: str,  # "goldenbraid" | "moclo" | custom grammar name
        parts: list[str],  # part names from the Parts Bin
        level: int | None = None,
        **kwargs,
    ) -> dict:
        """
        Simulate a Golden Braid or MoClo assembly.

        Parameters
        ----------
        grammar : str
            Assembly grammar: "goldenbraid", "moclo", or a custom grammar name.
        parts : list[str]
            Part names from the active Parts Bin.
        level : int or None
            Assembly level (L0, L1, L2, L3). Grammar-dependent.
        **kwargs
            Additional assembly parameters.

        Returns
        -------
        dict with the assembled construct sequence, features, and overhangs.

        [VERIFY: confirm endpoint name ("/simulate-assembly" vs "/design-gb-part")
         and parameter schema via `splicecraft-cli tools` on a running instance]
        """
        body: dict = {"grammar": grammar, "parts": parts}
        if level is not None:
            body["level"] = level
        body.update(kwargs)
        return self.post("/simulate-assembly", body)

    def design_gb_part(
        self,
        grammar: str,
        sequence: str,
        part_type: str,
        **kwargs,
    ) -> dict:
        """
        Design a Golden Braid / MoClo L0 part with correct overhangs.

        Parameters
        ----------
        grammar : str
            "goldenbraid" or "moclo"
        sequence : str
            DNA sequence of the part (without overhangs).
        part_type : str
            Part type / position in the grammar (e.g. "promoter", "CDS", "terminator").
            [VERIFY: valid values depend on the active grammar definition]
        """
        body = {"grammar": grammar, "sequence": sequence, "part_type": part_type}
        body.update(kwargs)
        return self.post("/design-gb-part", body)

    def simulate_gibson(
        self,
        fragment_names: list[str],
        **kwargs,
    ) -> dict:
        """
        Simulate a Gibson assembly from named library fragments.

        [VERIFY: exact parameter schema via `splicecraft-cli tools`]
        """
        return self.post("/simulate-gibson", {"fragments": fragment_names, **kwargs})

    # ------------------------------------------------------------------
    # BLAST / HMM search
    # ------------------------------------------------------------------

    def blast(
        self,
        query: str,
        program: str = "blastn",  # "blastn" | "blastp" | "blastx" | "tblastn" | "tblastx"
        db: str = "local",  # "local" | "ncbi" | "pfam"
        **kwargs,
    ) -> dict:
        """
        Run BLAST against the local library or NCBI.

        Parameters
        ----------
        query : str
            Query sequence (DNA or protein).
        program : str
            BLAST program.
        db : str
            "local" — in-process pyhmmer against your library
            "ncbi" — NCBI BLAST online
            "pfam" — Pfam/EMBL-EBI

        Returns
        -------
        dict with hits table.
        """
        return self.post("/blast", {"query": query, "program": program, "db": db, **kwargs})

    # ------------------------------------------------------------------
    # Codon optimization
    # ------------------------------------------------------------------

    def optimize_protein(self, amino_acid_sequence: str, codon_table: str | None = None) -> dict:
        """
        Codon-optimize an amino acid sequence to DNA using the active codon table.

        Parameters
        ----------
        amino_acid_sequence : str
            Single-letter AA sequence.
        codon_table : str or None
            Named codon table to use. If None, uses the active table.
        """
        body: dict = {"sequence": amino_acid_sequence}
        if codon_table:
            body["codon_table"] = codon_table
        return self.post("/optimize-protein", body)

    # ------------------------------------------------------------------
    # Sanger verification
    # ------------------------------------------------------------------

    def align_plasmidsaurus_zip(
        self,
        zip_path: str | Path,
        target_name: str | None = None,
    ) -> dict:
        """
        Align Plasmidsaurus sequencing results to the loaded (or named) plasmid.

        Parameters
        ----------
        zip_path : str or Path
            Path to the Plasmidsaurus results .zip file.
        target_name : str or None
            Name of the target library entry. If None, aligns to the loaded record.

        [VERIFY: exact parameter schema via `splicecraft-cli tools`]
        """
        body: dict = {"path": str(zip_path)}
        if target_name:
            body["target"] = target_name
        return self.post("/align-plasmidsaurus-zip", body)

    # ------------------------------------------------------------------
    # History / provenance
    # ------------------------------------------------------------------

    def get_history(self) -> dict:
        """Return the construct lineage as nested JSON (fragment-by-fragment provenance)."""
        return self.get("/get-history")

    # ------------------------------------------------------------------
    # Data safety
    # ------------------------------------------------------------------

    def list_backups(self) -> list[dict]:
        """List available data backups."""
        return self.get("/list-backups")


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    print("SpliceCraft API Client — Example Usage")
    print("=" * 50)

    try:
        sc = SpliceCraftClient()
    except Exception as e:
        print(f"[ERROR] Failed to initialize client: {e}")
        sys.exit(1)

    # 1. Check server status
    try:
        stat = sc.status()
        print(f"[OK] Server running. Loaded record: {stat.get('name', 'none')}")
    except ConnectionRefusedError as e:
        print(f"[ERROR] {e}")
        print("Run: splicecraft --agent")
        sys.exit(1)

    # 2. Fetch pUC19 from NCBI
    print("\nFetching pUC19 (L09137) from NCBI...")
    result = sc.fetch_ncbi("L09137")
    print(f"[OK] Loaded: {result}")

    # 3. List features
    print("\nListing features...")
    features = sc.list_features()
    print(f"[OK] {len(features)} features found:")
    for i, f in enumerate(features[:5]):
        print(f"  [{i}] {f.get('label', '?')} ({f.get('type', '?')}) "
              f"{f.get('start', '?')}..{f.get('end', '?')}")
    if len(features) > 5:
        print(f"  ... ({len(features) - 5} more)")

    # 4. Restriction scan
    print("\nRunning restriction enzyme scan...")
    sites = sc.restriction_scan(filter_mode="unique")
    print(f"[OK] {len(sites)} unique restriction sites found")

    # 5. Design detection primers for a region
    print("\nDesigning detection primers for bp 100-300...")
    try:
        primers = sc.design_primers(task="detection", target_start=100, target_end=300)
        print(f"[OK] Primers designed: {primers}")
    except Exception as e:
        print(f"[WARN] Primer design: {e}")

    # 6. Export GenBank
    output_path = "/tmp/pUC19_export.gb"
    print(f"\nExporting GenBank to {output_path}...")
    export = sc.export_genbank(output_path)
    print(f"[OK] Exported: {export}")

    print("\nDone. SpliceCraft API client is working correctly.")
