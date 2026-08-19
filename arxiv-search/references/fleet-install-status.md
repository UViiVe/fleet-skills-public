# Fleet arXiv Install Status

Install command: `pip3 install arxiv --break-system-packages`
Package version: arxiv 2.4.1

| Agent | Machine | Status | Date |
|-------|---------|--------|------|
| Atlas 🏛️ | Mac Studio M3 Ultra | ✅ Installed | 2026-03-23 |
| Buster 🐕 | Linux/Ubuntu (RTX 3080 Ti) | ⏳ Pending | — |
| Zevo 🦓 | Mac Mini | ⏳ Pending | — |
| Uvy 🦾 | — | ⏳ Pending | — |
| S.A.R.A.H. 🏠 | M1 Mac Mini | ⏳ Pending | — |
| Pip 🐣 | MacBook | ⏳ Pending | — |
| Mendel 🧬 | MacBook Neo | ⏳ Pending | — |

## Deploy SOP (per machine)

1. SSH to machine (requires Steve named approval per fleet SOP)
2. Run: `pip3 install arxiv --break-system-packages`
3. Verify: `python3 -c "import arxiv; print('OK')"`
4. Copy skill to `~/.openclaw/workspace/skills/arxiv-search/`
5. Update this file with ✅ + date

## Skill Distribution

Option A — Manual copy via SSH/SCP  
Option B — Drive-based upgrade via 1000 Acres Skills Library manifest (preferred once manifest is live)
