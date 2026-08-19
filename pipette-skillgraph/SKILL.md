---
id: skill_bd8b652a9a43790131ba1095296cf1d9
name: pipette-skillgraph
display-name: Pipette SkillGraph
short-description: Query the Pipette.bio SkillGraph MCP server to discover bioinformatics
  workflows, pipeline transitions, and tool documentation across 78 skills and 483
  pipeline edges. Use when you need to find, chain, or document a bioinformatics workflow
  that may already be catalogued.
category: methodology
keywords: pipette, skillgraph, bioinformatics workflows, workflow discovery, MCP
version: 1.0.0
last-updated: '2026-08-19'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
changelog:
- version: 1.0.0
  date: '2026-07-07'
  event: schema-compliance
  summary: 'Backfilled public schema. De-agentified: removed fleet agent callouts,
    internal paths, Chester SCP commands. Vienna public subset.'
description: Query the Pipette.bio SkillGraph MCP server to discover bioinformatics
  workflows, pipeline transitions, and tool documentation across 78 skills and 483
  pipeline edges. Use when you need to find, chain, or document a bioinformatics workflow
  that may already be catalogued.
---

# Pipette SkillGraph — Bioinformatics Workflow Discovery

A reference skill for querying the [Pipette.bio SkillGraph](https://skillgraph.pipette.bio) — a knowledge graph of 78 bioinformatics workflows, 483 pipeline transitions, and ~300 tools, built from open-access literature and the EDAM ontology.

**Primary use cases:**
- **Manuscript writing:** Find and cite canonical pipeline steps for a given omics data type
- **Workflow design:** Tool selection for multi-omics integration projects
- Discover upstream/downstream skills before designing a new analysis pipeline

---

## When to Use This Skill

✅ **Use when:**
- You need to know what tools/steps are standard for a given bioinformatics workflow (e.g., "what comes after variant-calling?")
- You're designing a multi-step pipeline and want to confirm the canonical transitions
- You want SKILL.md-style documentation for a workflow not available locally (ATAC-seq, metagenomics, ADMET, molecular docking, etc.)
- You're writing a methods section and need to cite the standard tool stack
- You want to check whether a skill you'd build already exists in the graph before building it

❌ **Do not use when:**
- You need to *run* a workflow — this skill is for discovery and documentation only
- The workflow is already well-documented in your local skills (check your local skills directory first)
- You need real-time data (this graph is static; last-known update reflects the repo state)

---

## Quick Start

The SkillGraph is a hosted MCP server. Query it via HTTP POST or directly via the MCP protocol.

### Direct HTTP (no MCP client needed)

```bash
# List all skills
curl -s -X POST https://skillgraph.pipette.bio/mcp \
  -H "Content-Type: application/json" \
  -d '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"list_skills","arguments":{}}}' \
  | python3 -m json.tool

# Get full SKILL.md for a specific workflow
curl -s -X POST https://skillgraph.pipette.bio/mcp \
  -H "Content-Type: application/json" \
  -d '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"get_skill","arguments":{"skill_id":"atacseq"}}}' \
  | python3 -m json.tool

# Search by keyword
curl -s -X POST https://skillgraph.pipette.bio/mcp \
  -H "Content-Type: application/json" \
  -d '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"search_skills","arguments":{"query":"CRISPR"}}}' \
  | python3 -m json.tool

# Find a pipeline path between two skills
curl -s -X POST https://skillgraph.pipette.bio/mcp \
  -H "Content-Type: application/json" \
  -d '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"find_path","arguments":{"source":"wgs-alignment","target":"pathway-enrichment"}}}' \
  | python3 -m json.tool
```

### MCP Protocol (Claude Code clients)

```bash
claude mcp add --transport http --scope user skillgraph https://skillgraph.pipette.bio/mcp
```

---

## Inputs

**Required:**
- A workflow name, keyword, or data type to query (e.g., `"atacseq"`, `"CRISPR"`, `"VCF"`)

**Optional:**
- Source + target skill IDs for pipeline path queries
- Domain keyword filter for `list_skills` (e.g., `"single-cell"`, `"drug-discovery"`)

---

## Outputs

Depending on the tool called:

| Tool | Output |
|---|---|
| `list_skills` | JSON list of skill IDs + brief descriptions, optionally filtered by domain |
| `get_skill` | Full SKILL.md content — triggers, tools, I/O formats, workflow steps, critical rules, output artifacts |
| `search_skills` | Ranked matches across skill IDs, keywords, and tool names |
| `get_transitions` | Upstream + downstream skills with edge evidence (paper counts, data types, EDAM validation) |
| `find_path` | Shortest pipeline path between two skills, annotated with data flow types (BAM→VCF→etc.) |
| `get_graph_stats` | Summary: total skills, edges, domain breakdown |

---

## ⚠️ CRITICAL: ASK THIS FIRST

**Does a local skill already cover this?**

Before querying the SkillGraph, check your local skills directory. Common locally available skills include:
- scRNA-seq (scanpy, seurat, trajectory inference, disease/drug discovery)
- Bulk RNA-seq (DESeq2, differential expression, counts-to-DE)
- ChIP-Atlas workflows
- GWAS, Mendelian randomization, polygenic risk scores
- Pathway/functional enrichment
- Multi-omics integration, proteomics

If a local skill exists, use it — it's already tuned to your conventions. Only query SkillGraph for workflows **not covered locally** or to cross-check a pipeline design.

---

## Available Skill Domains (as of 2026-05-27)

78 skills across these domains (from SkillGraph `get_graph_stats`):
- **Single-cell:** scanpy, seurat, trajectory inference, cell-cell communication
- **Bulk RNA-seq:** differential expression, co-expression networks
- **Epigenomics:** ATAC-seq, ChIP-seq, methylation
- **Genetics/Genomics:** WGS alignment, variant calling, GWAS analysis, polygenic risk scores
- **Drug discovery:** ADMET prediction, molecular docking, target identification
- **Metagenomics:** metagenome binning, taxonomic profiling
- **Proteomics:** mass-spec workflows
- **Multi-omics:** integration pipelines

Notable skills **not commonly available locally** (highest value targets to pull):
- `atacseq` — ATAC-seq peak calling, chromatin accessibility
- `metagenome-binning` — metagenome assembly + binning
- `admet-prediction` — drug ADMET property prediction
- `molecular-docking` — AutoDock / docking workflows
- `wgs-alignment` — WGS read alignment (bwa/samtools stack)

---

## Standard Workflow

### Step 1 — Identify the gap
Confirm the workflow you need isn't already available in your local skills directory.

### Step 2 — Query SkillGraph
Use `search_skills` or `list_skills` with a domain filter to find the relevant skill ID.

### Step 3 — Pull the SKILL.md
Call `get_skill` with the skill ID to retrieve the full documentation — tools, I/O, workflow steps, critical rules.

### Step 4 — Check transitions
Call `get_transitions` to see what feeds into and out of this skill. Essential for pipeline design.

### Step 5 — Optionally adapt for local use
If the skill is high-value and will be reused, adapt the retrieved SKILL.md to your local skill format and save it to your local skills directory for reuse.

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| MCP POST returns 404 | Wrong endpoint path | Use `https://skillgraph.pipette.bio/mcp` (no trailing slash) |
| `get_skill` returns empty | Skill ID wrong | Use `search_skills` first to get exact ID; IDs are kebab-case |
| Skill doc is terse / missing steps | Graph auto-generated from literature | Cross-check with primary tool docs; use as scaffold, not gospel |
| Pipeline path not found | Disconnected nodes in graph | Try intermediate steps manually; not all transitions are in graph |
| `list_skills` returns 78 but you need more detail | List-only view | Follow up with `get_skill` on the IDs of interest |
| Tool names differ from local conventions | Graph uses canonical names; local may have wrappers | Map graph tool name → local equivalent before assuming unavailable |
| EDAM ontology terms unfamiliar | Domain-specific vocabulary | See [references/edam-primer.md](references/edam-primer.md) |
| MCP not available (no Claude Code client) | Using a different agent framework | Use direct HTTP POST (see Quick Start) — same API |
| Adapted skill diverges from graph original | Local customization | Version-stamp the adaptation; note source + date in frontmatter |
| Graph data stale for fast-moving tools | Built from literature; lags new releases | Use graph for pipeline architecture; verify tool versions independently |
| `find_path` returns unexpectedly long path | Graph routes through indirect nodes | Check intermediate nodes — sometimes a shortcut exists outside the graph |
| Local skill overlaps with SkillGraph skill | Independently authored | Keep local skill canonical; use SkillGraph for cross-check only |

---

## Decision Points

| Decision | Default | Override when |
|---|---|---|
| Local vs. SkillGraph | Use local skill | No local skill exists for the workflow |
| Pull full SKILL.md vs. summary | Pull full doc | Only need tool names → `list_skills` or `search_skills` |
| Adapt to local format vs. use as-is | Use as-is for one-off | Will be reused → adapt and save to your local skills directory |
| MCP vs. HTTP POST | HTTP POST (works everywhere) | Claude Code client available → native MCP |

---

## Suggested Next Steps

After finding a workflow in SkillGraph:

1. **Adapt + install:** If the skill is worth reusing, use the `skill-creator` skill to package the retrieved SKILL.md into your local format, then save it to your local skills directory.
2. **Manuscript writing use case:** Paste the `get_skill` output directly into a manuscript methods section as a drafting scaffold — it includes tools, steps, and output artifacts.
3. **Pipeline design use case:** Use `find_path` to design multi-step pipeline architectures before committing to tool selection.

---

## References

**External:**
- SkillGraph repo: https://github.com/variomeanalytics/bioinformatics-agent-skills
- Live graph explorer: https://skillgraph.pipette.bio
- MCP endpoint: https://skillgraph.pipette.bio/mcp
- Pipette.bio: https://pipette.bio

**Bundled:**
- [references/edam-primer.md](references/edam-primer.md) — EDAM ontology quick reference for bioinformatics data types
- [references/available-skills.md](references/available-skills.md) — snapshot of SkillGraph skill IDs as of 2026-05-27
