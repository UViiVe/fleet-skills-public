# EDAM Ontology — Quick Reference for SkillGraph Data Types

EDAM (EMBRACE Data And Methods) is a bioinformatics ontology used to annotate tool inputs/outputs. The SkillGraph uses EDAM terms to validate pipeline edge compatibility.

## Most Common Data Types You'll Encounter

| EDAM Term | What It Is | Typical File Format |
|---|---|---|
| `data:2044` | Sequence reads | FASTQ |
| `data:1385` | Sequence alignment | BAM, SAM, CRAM |
| `data:3498` | Sequence variations (variants) | VCF, BCF |
| `data:3112` | Gene expression matrix | TSV, H5AD, RDS |
| `data:3917` | Count matrix (raw) | MTX, H5, CSV |
| `data:3002` | Annotation track (peaks, intervals) | BED, narrowPeak, broadPeak |
| `data:3108` | Phenotype data | TSV, CSV |
| `data:1772` | Pathway/ontology annotation | GMT, GMX |
| `data:2976` | Protein structure | PDB |
| `data:2536` | Pathway | SBML, GPML |

## Pipeline Compatibility Rule (from SkillGraph)

An edge A→B is valid when:
- The *output* EDAM type of skill A matches the *input* EDAM type of skill B
- The edge has ≥1 supporting paper OR is in the curated ground-truth graph

## Why This Matters

When `get_transitions` returns an edge, the data types listed tell you what format to pass between steps. If your intermediate file is in a different format than the graph expects, you need a conversion step not shown in the graph.

## Reference

EDAM ontology: https://edamontology.org
OWL file: https://github.com/edamontology/edamontology
