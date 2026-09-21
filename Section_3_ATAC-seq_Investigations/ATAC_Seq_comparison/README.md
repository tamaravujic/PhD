# HDLEC vs imLEC ATAC-seq comparison

ATAC-seq chromatin accessibility comparison between primary lymphatic endothelial cells (HDLEC) and the immortalised imLEC line, plus an ERG cis-regulatory element (CRE) check, to assess whether imLEC is a valid model system for CRISPRi work targeting ERG binding sites.

Pipeline created by Dr Hannah Maude in the Cebola lab. Adapted by Tamara Vujic.

## Contents

| File | Description |
|---|---|
| `ATAC_HDLEC_imLEC_Comparison.ipynb` | Main analysis notebook (Bash + R). Start here — it walks through QC, the genome-wide comparison, and the ERG-specific comparison, and tells you when to run each script below. |
| `corr_replicates.sh` | PBS script: deepTools `multiBigwigSummary` + `plotCorrelation` across individual replicate bigwigs (HDLEC1, HDLEC2, imLEC1, imLEC2), 10 kb bins. |
| `corr_pooled.sh` | PBS script: same correlation QC, but on pooled (replicate-merged) HDLEC vs imLEC bigwigs. |
| `computeMatrix_HDLEC_imLEC_identity.sh` | PBS script: deepTools `computeMatrix` over the genome-wide HDLEC/imLEC-specific and shared CRE sets (Part 1 heatmap input). |
| `computeMatrix_ERG_functional_HDLEC_imLEC.sh` | PBS script: deepTools `computeMatrix` over the retained/lost functional ERG CRE sets (Part 2 heatmap input). |

## Requirements

**Conda environments** (see the notebook's Requirements cell for creation commands):
- `diffbind` — R 4.3.3, built from `diffbind.yml` (DiffBind, GenomicRanges, tidyverse, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, clusterProfiler, etc.). Used for all R / VS Code sections.
- `ATAC` — built from `ATAC.yml` (deepTools, BEDTools). Used for all terminal sections (`computeMatrix`, `plotHeatmap`, `multiBigwigSummary`, `plotCorrelation`, `intersectBed`).

**Input files** (see the notebook's "Files required" cell for the full list), broadly:
- Per-replicate BAMs (blacklist-filtered, indexed) and MACS2 peak calls for HDLEC and imLEC
- Per-replicate and pooled BigWigs (`*.BPM.bw`)
- HDLEC and imLEC consensus peak/CRE BED files
- ERG ChIP-seq narrowPeak calls and a putative tiered ERG CRE BED file
- `hg38-blacklist.v2.bed` (ENCODE hg38 blacklist)

**Note on paths:** the PBS scripts and sample sheet use absolute paths specific to the original HPC project directory (`/rds/general/user/tv722/...`). Update `DIR`, `REFDIR`, and the sample sheet paths before running on a different system.

## Workflow

1. **Setup** — create the `diffbind` and `ATAC` conda environments (notebook Requirements cell).
2. **QC** — run `corr_replicates.sh` and `corr_pooled.sh` for a genome-wide sanity check of replicate and cell-type agreement.
3. **Part 1: genome-wide comparison** — build a consensus CRE set, count reads with DiffBind, classify HDLEC-specific / imLEC-specific / shared regions, then visualise with `computeMatrix_HDLEC_imLEC_identity.sh` and `plotHeatmap`.
4. **Part 2: ERG-specific comparison** — intersect ERG ChIP peaks with HDLEC ATAC to define functional ERG CREs, check which are retained/lost in imLEC, then visualise with `computeMatrix_ERG_functional_HDLEC_imLEC.sh` and `plotHeatmap`.
