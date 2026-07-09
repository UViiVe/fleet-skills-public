---
name: splicecraft
id: "skill_79f66d3d81b148e297a80d2a931883f0"
display-name: "SpliceCraft"
short-description: "Programmatic plasmid design, cloning simulation, and primer design via the SpliceCraft local API. Use when construct work should stay local."
category: wet_lab
keywords: [plasmid, cloning, primer design, Golden Gate, SpliceCraft]
last-updated: "2026-07-09"
version: "1.0.1"
author: "Fleet Skills contributors"
created: 2026-05-29
description: "Programmatic plasmid design, cloning simulation, and primer design via the SpliceCraft local API. Supports Golden Braid/MoClo, restriction analysis, BLAST, and SOE-PCR mutagenesis. Local-only: no data leaves the machine."
status: active
promoted: true
license: MIT
changelog:
  - version: "1.0.1"
    date: "2026-07-09"
    event: public-cleanup
    summary: "Removed internal machine notes and documented read-only default API helper behavior."
---

# SpliceCraft

SpliceCraft is a terminal-based local plasmid workbench. Its local REST API can
be used for scripted construct workflows without sending plasmid data to a cloud
service.

**Project:** https://github.com/Binomica-Labs/SpliceCraft

---

## When to Use

| Task | SpliceCraft action |
|---|---|
| Plasmid design or construct planning | Load `.gb` file, inspect map, edit sequence |
| Golden Gate or MoClo simulation | `simulate-assembly` |
| Restriction enzyme analysis | `list-restriction-sites` |
| Primer design | `design-primers` or `design-mutagenesis` |
| BLAST search | `blast` endpoint, local engine, or NCBI online |
| Sanger read verification | `align-plasmidsaurus-zip` |
| NCBI accession fetch | `splicecraft L09137` or `/fetch` API endpoint |
| Export | `export-genbank`, `export-fasta`, or `export-gff3` |

Do not use this skill to bypass SpliceCraft's normal save paths or to write to a
shared plasmid library without an explicit backup and review step.

---

## Installation

```bash
pipx install splicecraft
splicecraft --version
```

Alternative install:

```bash
uv tool install splicecraft
```

Requirements: Python 3.10+ and a terminal at least 100x30.

---

## Start the Local API Server

```bash
splicecraft --agent
splicecraft --agent --agent-port 6800
```

The server writes a token file at startup:

- macOS: `~/Library/Application Support/splicecraft/agent_token`
- Linux: `~/.local/share/splicecraft/agent_token`
- Windows: `%APPDATA%\splicecraft\agent_token`

Token file format:

```text
<port>
<bearer_token>
```

The token rotates each launch. Re-read it each session.

---

## Python Helper

Use `scripts/splicecraft_api.py` for token discovery, authentication, and common
operations.

```python
from splicecraft_api import SpliceCraftClient

sc = SpliceCraftClient()  # read-only by default
sc.fetch_ncbi("L09137")
features = sc.list_features()
sites = sc.restriction_scan()
primers = sc.design_primers(task="detection", target_start=100, target_end=300)
sc.export_genbank("/path/to/output.gb")
```

Write operations require an explicit opt-in:

```python
sc = SpliceCraftClient(allow_write=True)
sc.add_feature(start=100, end=200, label="example")
sc.save()
```

Enable writes only after confirming the active record and backup state.

---

## CLI Sidecar

```bash
splicecraft-cli status
splicecraft-cli tools
splicecraft-cli features
splicecraft-cli fetch L09137
splicecraft-cli load-file /path/to/plasmid.gb
splicecraft-cli list-restriction-sites
splicecraft-cli export-genbank /path/output.gb
```

Most subcommands support `--json` for machine-readable output.

---

## Safety Rules

- Never import private SpliceCraft modules from ad-hoc scripts to write user data.
- Use the HTTP API for programmatic access.
- Keep plasmid libraries backed up before enabling write operations.
- Confirm endpoint signatures with `splicecraft-cli tools` on the running version.
- Inspect exported records before using them in wet-lab work.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `Connection refused` | Start the server with `splicecraft --agent` |
| `Token file not found` | Start the server or check `SPLICECRAFT_DATA_DIR` |
| Port conflict | Use `splicecraft --agent --agent-port 6800` |
| Command not found | Reinstall with `pipx install splicecraft` and run `pipx ensurepath` |
| Python version error | Use Python 3.10+ |
| Need to verify reads | Use `splicecraft-cli align-plasmidsaurus-zip` or the GUI |

---

## References

- Agent API docs: https://github.com/Binomica-Labs/SpliceCraft/blob/master/docs/agent-api.md
- CLI sidecar docs: https://github.com/Binomica-Labs/SpliceCraft/blob/master/docs/cli.md
- Install docs: https://github.com/Binomica-Labs/SpliceCraft/blob/master/docs/install.md
- Full feature reference: https://github.com/Binomica-Labs/SpliceCraft/blob/master/docs/features.md
