# Acknowledgments and Attribution

This repository builds upon work from many researchers, bioinformaticians, and tool developers. We gratefully acknowledge the following contributions:

---

### Dr Hannah Maude (Cebola Lab, Imperial College London)

Dr Maude contributed significantly to the development of analysis workflows, particularly:
- **ChromBPNet Analysis Pipeline:** Complete analysis workflow adaptation of ChromBPNet deep learning model for variant effect prediction
- **ATAC-seq Analysis:** Complete ATAC-seq analysis workflow from raw sequencing to peak calling and enhancer characterisation
- **Guidance for all computational analysis 


----
## External Tools and Pipelines

### ChromBPNet
**ChromBPNet:** The official github for ChromBPNet can be found here https://github.com/kundajelab/chrombpnet?tab=readme-ov-file
The Pre print https://www.biorxiv.org/content/10.1101/2024.12.25.630221v2

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
| **VEP (Variant Effect Predictor)** | Dyer, S. C. et al. Ensembl 2025. *Nucleic Acids Research* 53, D948–D957 (2024) |
| **CADD (Combined Annotation Dependent Depletion)** | Tenywa, J. F. et al. Genome region aware CADD thresholds for noncoding variant prioritization. *NAR Genomics and Bioinformatics* 7 (2025) |
| **Fathmm** | Shihab, H. A. et al. An integrative approach to predicting the functional effects of non-coding and coding sequence variation. *Bioinformatics* 31, 1536–1543 (2015) |
| **RegVar** | Lu, H. et al. RegVar: Tissue-specific Prioritization of Non-coding Regulatory Variants. *Genomics Proteomics Bioinformatics* 21, 385–395 (2023) |
| **DeepSEA** | Zhou, J. & Troyanskaya, O. G. Predicting effects of noncoding variants with deep learning-based sequence model. *Nature Methods* 12, 931–934 (2015) |
| **SpliceAI** | Jaganathan, K. et al. Predicting Splicing from Primary Sequence with Deep Learning. *Cell* 176, 535–548.e24 (2019) |
| **motifbreakR** | Coetzee, S. G., Coetzee, G. A. & Hazelett, D. J. motifbreakR: an R/Bioconductor package for predicting variant effects at transcription factor binding sites. *Bioinformatics* btv470 (2015) |
| **ChromBPNet** | Pampari, A. et al. ChromBPNet: bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants. *bioRxiv* (2025) doi:10.1101/2024.12.25.630221 |
| **AlphaGenome** | Avsec, Ž. et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 649, 1206–1218 (2026) |

---

## Data Resources

### Genomics England (GEL)

We acknowledge Genomics England for providing access to whole genome sequencing data from the 100,000 Genomes Project
**Citation:** https://www.genomicsengland.co.uk/

### ENCODE Project

ATAC-seq data for lymphatic endothelial cell types was downloaded from:
- ENCODE Consortium portal: https://www.encodeproject.org/
- Used for enhancer specificity comparisons


**Last Updated:** August 2026