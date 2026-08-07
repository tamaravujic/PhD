# Primary Lymphedema Genomic Analysis Pipeline

## Overview

Collection of scripts used to identifying non-coding genetic variants in the ERG gene, associated with Primary Lymphedema (PL). Integrates variant annotation, structural variant analysis, chromatin accessibility prediction, and transcription factor binding predictions across multiple genome assemblies.

This work is for a PhD project funded by the British Heart Foundation at Imperial College London
Supervised by Dr Graeme Birsdey and Dr Inês Cebola

---

## Scripts by Analysis Section

---

## SECTION 1: Creating a Cohort and Identifying Variants

**GEL-specific scripts** - Scripts prefixed with `GEL_` require access to Genomics England secure research environment. They will not run outside this environment.

| Script | Purpose |
|--------|---------|
| `GEL_COHORT_VARIANT_SEARCH.sh` | Search for specific variants in GEL cohort |
| `GEL_Variant_search_batch.sh` + `GEL_Variant_search_batch.conf` | Batch search for variants in multiple genomic regions|
| `GEL_vcf_variant_checker.sh` | Check presence/genotypes of variants in VCF files |
| `GEL_denovo_variant_filtering.sh` | Filter de novo variants in ERG regions |
| `GEL_SV_SEARCH.sh` | Search for structural variants in GEL cohort |
| `GEL_find_SV_breakpoints.sh` | Find SV breakpoints in GEL data |
| `GEL_VEP.sh` | VEP annotation inside GEL environment |
| `GEL_relatedness.txt` | Documentation: checking kinship in GEL data |

**General variant discovery pipeline**

---

## SECTION 2: Variant Pathogenicity Analysis
| Script | Purpose |
|--------|---------|
| `Variant_Pathogenicity_threshold_filtering.R` | Filter variants by pathogenicity scores (CADD, REVEL, DeepSEA) |

### Subsection 2.1: MotifbreakR - Transcription Factor Binding Prediction

| Script | Purpose |
|--------|---------|
| `MotifbreakR_SNV_pipeline.R` | Predict TF binding disruption for SNVs (4 databases: HOCOMOCO, HOMER, ENCODE, FactorBook) |
| `MotifbreakR_INDEL_pipeline.R` | Predict TF binding disruption for INDELs |
| `motifbreakR_output_sort.R` | Consolidate motifbreakR results with pathogenicity scores |

### Subsection 2.2: ChromBPNet - Deep Learning Chromatin Accessibility

| Script | Purpose |
|--------|---------|
| `chrombpnet.ipynb` | ChromBPNet analysis: predictions, SHAP scores, visualisation |
| `chrombpnet_pipeline_sif.sh` | Main ChromBPNet prediction pipeline (Singularity container) |
| `ChrombpNet_train_bias_model_sif.sh` | Train bias models across 5-fold cross-validation |
| `ChromBPNet_variant_prediction.sh` | Score variants with ChromBPNet |
| `chrombpnet_variant_scores_analysis.R` | Generate volcano plots and group comparisons |

### Subsection 2.3: DeepSEA - Non-coding Variant Impact

| Script | Purpose |
|--------|---------|
| `deepsea_column_extraction.sh` | Extract DeepSEA predictions from annotation files |

### Subsection 2.4: RegVar - Tissue-specific Regulatory Predictions

| Script | Purpose |
|--------|---------|
| `regvar_postprocessing.sh` | Filter RegVar output by variant-gene pair list |

### Subsection 2.5: AlphaGenome - Non-coding Variant Prediction

| Script | Purpose |
|--------|---------|
| `Alphagenome_github.ipynb` | Predict effects of non-coding variants using AlphaGenome API |

### Subsection 2.6: VEP - Variant Effect Predictor

| Script | Purpose |
|--------|---------|
| `VEP_annotation.sh` | Annotate variants with VEP 111 (CADD, REVEL, gnomAD, SpliceAI plugins) |


---

## SECTION 3: ATAC-seq Investigations

### Subsection 3.1: imLEC ATAC-seq analysis

| Script | Purpose |
|--------|---------|
| `imlec_atac_seq.ipynb` | ATAC-seq processing pipeline for immortalised lymphatic endothelial cells |
| `imLEC_ATAC_venn_diagrams.R` | Generate Venn diagrams for peak overlaps between cell types |

### Subsection 3.2: Enhancer Cell type specifcity 

| Script | Purpose |
|--------|---------|
| `Enhancer_specificity_heatmaps-2.ipynb` | DiffBind analysis and enhancer specificity heatmap generation |
| `Submit_dbaCount-2.sh` | PBS job for DiffBind peak counting (dba.count step) |
| `Binary_enhancer_specificity_heatmaps_github.R` | Generate binary enhancer specificity heatmaps |
| `ENCODE_ATAC_data_processing_for_heatmap.sh` | Process ENCODE ATAC-seq data for comparison with in-house data |


---

## Citation

If using this pipeline in research, please cite:

Vujic T (2026). Primary Lymphedema genomic analysis pipeline. GitHub.

And cite relevant tool papers (see ACKNOWLEDGMENTS.md).

---

## Contact

- **Author:** Tamara Vujic, Imperial College London
- tv722@ic.ac.uk

---

**Last Updated:** August 2026
