# Cell Type Marker Genes (Seurat-flavored)

A curated set of canonical markers for cluster annotation. Use 3–5 markers per cell type — single-marker calls are routinely wrong.

The lists below apply equally to scanpy and Seurat workflows. Where naming conventions differ (Seurat capitalizes gene symbols), both forms are valid for human / mouse data depending on rownames in your object.

## Immune cells

### T cells

- **Pan-T:** `CD3D`, `CD3E`, `CD3G`, `CD2`
- **CD4+ helper:** `CD4`, `IL7R` (CD127), `LDHB`
- **CD8+ cytotoxic:** `CD8A`, `CD8B`, `GZMK`, `GZMB`
- **Naive:** `CCR7`, `SELL` (CD62L), `LEF1`, `TCF7`
- **Memory:** `IL7R` (central), `S100A4`, `GZMK` (effector memory)
- **Regulatory:** `FOXP3`, `IL2RA` (CD25), `CTLA4`, `IKZF2`
- **Activated / exhausted:** `TNFRSF9` (CD137), `CD69`, `IFNG`, `PDCD1`, `HAVCR2`
- **gamma-delta:** `TRDV1`, `TRDV2`, `TRGV9`
- **MAIT:** `SLC4A10`, `KLRB1`, `RORC`

### B cells

- **Pan-B:** `CD79A`, `CD79B`, `MS4A1` (CD20), `CD19`
- **Naive:** `IGHD`, `IGHM`, `TCL1A`, `FCER2`
- **Memory:** `CD27`, `TNFRSF13B`, `AIM2`
- **Plasma cells:** `MZB1`, `SSR4`, `JCHAIN`, `XBP1`, `IGHA1`, `IGHG1`, `PRDM1`

### Monocytes / macrophages

- **Classical (CD14+):** `CD14`, `S100A8`, `S100A9`, `FCN1`, `VCAN`
- **Non-classical (CD16+):** `FCGR3A` (CD16), `CDKN1C`, `MS4A7`, `LST1`
- **Intermediate:** `CD14` + `FCGR3A` co-expression
- **Tissue macrophage:** `CD68`, `C1QA`, `C1QB`, `APOE`, `MARCO`
- **Alveolar:** `MARCO`, `MCEMP1`, `MSR1`

### Dendritic cells

- **cDC1:** `CLEC9A`, `XCR1`, `IDO1`
- **cDC2:** `FCER1A`, `CD1C`, `CLEC10A`
- **pDC:** `LILRA4`, `CLEC4C`, `IL3RA`, `IRF7`, `IRF8`, `GZMB`

### NK and ILC

- **NK:** `GNLY`, `NKG7`, `GZMB`, `NCAM1` (CD56), `KLRB1` (CD161), `KLRD1`, `PRF1`
- **CD56-bright (regulatory):** `NCAM1` high, `GZMK`, low `FCGR3A`
- **CD56-dim (cytotoxic):** `NCAM1` low, `FCGR3A` high, `GNLY`
- **ILC1 / ILC2 / ILC3:** `IL7R`, no CD3, then `TBX21`, `GATA3`, `RORC` respectively. SingleR often mis-calls ILCs as NK.

### Other

- **Megakaryocytes / platelets:** `PPBP`, `PF4`, `GP9`, `GP1BB`, `ITGA2B`
- **Mast cells:** `TPSAB1`, `TPSB2`, `CPA3`, `KIT` (CD117)
- **Eosinophils:** `CCR3`, `IL5RA`, `EPX`, `PRG2`
- **Neutrophils** (rarely captured in 10X due to fragility): `FCGR3B`, `CSF3R`, `S100A8/9`

## Brain and neural tissue

### Neurons

- **Pan-neuronal:** `RBFOX3` (NeuN), `SNAP25`, `SYT1`, `MAP2`
- **Excitatory:** `SLC17A7` (VGLUT1), `SLC17A6` (VGLUT2), `CAMK2A`, `SATB2`
- **Inhibitory:** `GAD1`, `GAD2`, `SLC32A1` (VGAT), `DLX1`, `DLX2`
- **Interneuron subtypes:** `SST`, `PVALB`, `VIP`, `LAMP5`
- **Dopaminergic:** `TH`, `SLC6A3` (DAT), `DRD2`, `NR4A2`
- **Serotonergic:** `TPH2`, `SLC6A4` (SERT), `FEV`
- **Cholinergic:** `CHAT`, `SLC18A3` (VAChT)

### Glia

- **Astrocytes:** `AQP4`, `GFAP`, `SLC1A2` (GLT1), `SLC1A3` (GLAST), `ALDH1L1`
- **Oligodendrocytes (mature):** `MBP`, `MOG`, `PLP1`, `MAG`
- **OPCs:** `PDGFRA`, `CSPG4` (NG2), `SOX10`
- **Microglia:** `CX3CR1`, `P2RY12`, `TMEM119`, `AIF1` (IBA1), `CSF1R`
- **Ependymal:** `FOXJ1`, `DNAH11`, `CFAP44`

## Epithelial — lung

- **AT1:** `AGER`, `PDPN`, `CLIC5`, `HOPX`
- **AT2:** `SFTPC`, `SFTPB`, `SFTPA1`, `ABCA3`
- **Ciliated:** `FOXJ1`, `DNAH5`, `DNAI2`, `TPPP3`
- **Club (Clara):** `SCGB1A1`, `SCGB3A2`, `CYP2F2`
- **Goblet:** `MUC5AC`, `MUC5B`, `TFF3`
- **Basal:** `KRT5`, `KRT14`, `TP63`
- **Tuft / brush:** `POU2F3`, `TRPM5`, `DCLK1`

## Epithelial — intestine

- **Enterocytes:** `FABP1`, `APOA1`, `SI`, `KRT20`
- **Goblet:** `MUC2`, `TFF3`, `CLCA1`
- **Paneth:** `LYZ`, `DEFA5`, `DEFA6`
- **Enteroendocrine:** `CHGA`, `TPH1`, `NEUROD1`
- **Tuft:** `DCLK1`, `TRPM5`, `POU2F3`
- **Stem cells (LGR5+):** `LGR5`, `OLFM4`, `ASCL2`
- **Transit-amplifying:** `MKI67`, `TOP2A` + enterocyte markers

## Epithelial — pancreas

- **Acinar:** `PRSS1`, `CTRB1`, `CPA1`, `CPB1`
- **Ductal:** `KRT19`, `MUC1`, `CFTR`
- **Beta cells:** `INS`, `IAPP`, `MAFA`, `NKX6-1`
- **Alpha cells:** `GCG`, `IRX2`, `MAFB`
- **Delta cells:** `SST`, `HHEX`
- **PP cells:** `PPY`
- **Epsilon cells:** `GHRL`

## Stromal and vascular

### Fibroblasts

- **Pan-fibroblast:** `COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `PDGFRA`
- **Myofibroblasts:** `ACTA2` (αSMA), `TAGLN`, `MYH11`
- **Inflammatory CAFs:** `IL6`, `CXCL12`, `HAS1`
- **Matrix CAFs:** `COL1A1`, `COL3A1`, `POSTN`

### Endothelial

- **Pan-endothelial:** `PECAM1` (CD31), `VWF`, `CDH5` (VE-cadherin), `CLDN5`
- **Arterial:** `GJA5`, `EFNB2`, `DLL4`, `HEY1`
- **Venous:** `NR2F2`, `EPHB4`, `ACKR1`
- **Lymphatic:** `PROX1`, `FLT4` (VEGFR3), `LYVE1`, `PDPN`
- **Capillary:** `CA4`, `RGCC`, `KDR` (VEGFR2)

### Mural

- **Pericytes:** `RGS5`, `PDGFRB`, `CSPG4` (NG2)
- **Smooth muscle:** `ACTA2`, `TAGLN`, `MYH11`, `CNN1`, `DES`

## Stem and progenitor

- **HSC:** `CD34`, `KIT` (CD117), `THY1` (CD90), `HLF`, `MLLT3`
- **CMP:** `MPO`, `CSF1R`, `MEF2C`
- **CLP:** `IL7R`, `DNTT`, `RAG1/2`
- **MSC:** `THY1` (CD90), `ENG` (CD105), `NT5E` (CD73), `PDGFRB`
- **Neural stem/progenitor:** `NES`, `SOX2`, `FABP7` (BLBP), `HES1`, `HES5`

## Visualization in Seurat

```r
# Feature plots
FeaturePlot(seurat_obj, features = c("CD3D", "CD4", "CD8A"), cols = c("lightgrey", "red"))

# Dot plot (best for many markers across clusters)
DotPlot(seurat_obj, features = c("CD3D", "CD4", "CD8A", "CD14", "MS4A1", "NKG7")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Violin plots
VlnPlot(seurat_obj, features = c("CD3D", "CD8A"), pt.size = 0)

# Heatmap of top markers
top10 <- markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 10)
DoHeatmap(seurat_obj, features = top10$gene)
```

## Annotation workflow

1. Run `find_all_cluster_markers()` (Wilcoxon).
2. For each cluster, look at the top 10–20 genes by `avg_log2FC`.
3. Cross-reference against this table and tissue-specific atlases.
4. Visualize with dot plots and feature plots to confirm specificity.
5. Assign labels via `annotate_clusters_manual()`.
6. Optionally validate with `annotate_with_singler()` or `annotate_with_azimuth()`.
7. Re-inspect: clusters without clear identity may be doublets, transitional states, or genuinely novel.

## Pitfalls

- **Single-marker calls.** `CD4` is expressed on monocytes; `CD8A` shows up on NK cells at low levels. Confirm with multiple markers.
- **Ambient RNA contamination.** Hemoglobin (`HBA1`, `HBB`) in non-erythroid clusters is almost always soup — SoupX removes it.
- **State vs identity.** Activated CD8 looks different from resting CD8 but is the same lineage. Distinguish type from state.
- **Mixed clusters.** A cluster with T-cell *and* B-cell markers is almost certainly doublets.
- **Ribosomal / mitochondrial dominance.** If a cluster's top markers are all `RPS*`, `RPL*`, `MT-*`, it's a technical artifact.

## SingleR and Azimuth references

- **HPCA** (Human Primary Cell Atlas): broad human reference, good for PBMC and most blood/immune.
- **BlueprintEncode**: human tissue, useful for hematopoietic + stromal.
- **ImmGen**: mouse immune.
- **MonacoImmuneData**: high-resolution human immune.
- **DICE**: human PBMC subtypes.
- **Azimuth references**: tissue-specific, pre-trained (pbmc, lung, heart, etc.).

Pick the reference that matches the tissue. Mixed references give worse results than the right specific one.

## Citations

1. Domínguez Conde C, Xu C, Jarvis LB, et al. *Cross-tissue immune cell analysis reveals tissue-specific features in humans (CellTypist).* Science 376:eabl5197 (2022).
2. Aran D, Looney AP, Liu L, et al. *Reference-based analysis of lung single-cell sequencing (SingleR).* Nat Immunol 20:163–172 (2019).
3. Hao Y, et al. *Integrated analysis of multimodal single-cell data (Azimuth).* Cell 184:3573–3587.e29 (2021).
4. Franzén O, Gan LM, Björkegren JLM. *PanglaoDB.* Database 2019:baz046 (2019).
5. Tabula Sapiens Consortium. *The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans.* Science 376:eabl4896 (2022).
