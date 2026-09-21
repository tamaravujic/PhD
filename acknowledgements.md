# Acknowledgments and Attribution

This repository builds upon work from many researchers, bioinformaticians, and tool developers. We gratefully acknowledge the following contributions:

---

### Dr Hannah Maude (Cebola Lab, Imperial College London)

Dr Maude contributed significantly to the development of analysis workflows, particularly:
- **ChromBPNet Analysis Pipeline:** Complete analysis workflow adaptation of ChromBPNet deep learning model for variant effect prediction
- **ATAC-seq Analysis:** Complete ATAC-seq analysis workflow from raw sequencing to peak calling and enhancer characterisation
- **HDLEC vs imLEC ATAC-seq Comparison Pipeline:** Genome-wide and ERG-specific chromatin accessibility comparison pipeline, adapted by Tamara Vujic
- **Guidance for all computational analysis**


----
## External Tools and Pipelines

### ChromBPNet
**Official ChromBPNet Pipeline:** https://github.com/kundajelab/chrombpnet?tab=readme-ov-file

---
### ATAC-seq Analysis Pipeline

**CebolaLab ATAC-seq Pipeline:** https://github.com/CebolaLab/ATAC-seq

This pipeline provides best-practice workflows for processing ATAC-seq data, including:
- Quality control (FastQC)
- Read alignment (Bowtie2)
- Peak calling (MACS2)
- Enhancer characterisation

---

## Bioinformatics Tools (Alphabetical Order)

We acknowledge the developers of the following essential bioinformatics tools:


| Tool | Reference |
|------|-----------|
| **AlphaGenome** | Avsec, Ž. et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 649, 1206–1218 (2026) |
| **BEDTools** | Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26, 841–842 (2010) |
| **CADD (Combined Annotation Dependent Depletion)** | Tenywa, J. F. et al. Genome region aware CADD thresholds for noncoding variant prioritization. *NAR Genomics and Bioinformatics* 7 (2025) |
| **ChIPseeker** | Yu, G., Wang, L.-G. & He, Q.-Y. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. *Bioinformatics* 31, 2382–2383 (2015) |
| **ChromBPNet** | Pampari, A. et al. ChromBPNet: bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants. *bioRxiv* (2025) doi:10.1101/2024.12.25.630221 |
| **DeepSEA** | Zhou, J. & Troyanskaya, O. G. Predicting effects of noncoding variants with deep learning-based sequence model. *Nature Methods* 12, 931–934 (2015) |
| **deepTools** | Ramírez, F., Ryan, D. P., Grüning, B. et al. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Research* 44, W160–W165 (2016) |
| **DiffBind** | Stark, R. & Brown, G. DiffBind: differential binding analysis of ChIP-Seq peak data. Bioconductor vignette (2011). https://www.bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf |
| **Fathmm** | Shihab, H. A. et al. An integrative approach to predicting the functional effects of non-coding and coding sequence variation. *Bioinformatics* 31, 1536–1543 (2015) |
| **motifbreakR** | Coetzee, S. G., Coetzee, G. A. & Hazelett, D. J. motifbreakR: an R/Bioconductor package for predicting variant effects at transcription factor binding sites. *Bioinformatics* btv470 (2015) |
| **RegVar** | Lu, H. et al. RegVar: Tissue-specific Prioritization of Non-coding Regulatory Variants. *Genomics Proteomics Bioinformatics* 21, 385–395 (2023) |
| **SpliceAI** | Jaganathan, K. et al. Predicting Splicing from Primary Sequence with Deep Learning. *Cell* 176, 535–548.e24 (2019) |
| **VEP (Variant Effect Predictor)** | Dyer, S. C. et al. Ensembl 2025. *Nucleic Acids Research* 53, D948–D957 (2024) |

---

## Data Resources

### Genomics England (GEL)

We acknowledge Genomics England for providing access to whole genome sequencing data from the 100,000 Genomes Project
**Citation:** https://www.genomicsengland.co.uk/

### ENCODE Project

ATAC-seq data for lymphatic endothelial cell types was downloaded from:
- ENCODE Consortium portal: https://www.encodeproject.org/
- Used for enhancer specificity comparisons


**Last Updated:** September 2026
