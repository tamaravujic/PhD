# ==============================================================================
# Script: imLEC_ATAC_venn_diagrams.R
#
# Purpose:
#   Generate Venn diagrams comparing ATAC-seq peak accessibility between
#   immortalised lymphatic endothelial cells (imLEC) and primary human dermal
#   lymphatic endothelial cells (HDLEC). Creates three comparative visualisations:
#   genome-wide peaks, ERG-associated peaks, and shared ERG cis-regulatory elements.
#
# Workflow:
#   1. Load peak files for imLEC and HDLEC samples
#   2. Calculate overlap statistics (shared peaks)
#   3. Generate Venn diagrams for each dataset comparison
#   4. Visualise cell-type-specific and shared chromatin accessibility
#
# Input Files:
#   imlec.replicated_broadPeak.sorted.bed - imLEC ATAC-seq peaks
#   LEC_CREs.sorted.bed - HDLEC ATAC-seq peaks
#   shared_peaks.bed - peaks overlapping in both cell types
#   imlec_ERG.bed, LEC_CREs_ERG.bed - ERG region peaks
#   shared_ERG.bed - shared ERG peaks
#   (And corresponding CRE-specific files)
#
# Output:
#   Three Venn diagram visualisations (displayed in R graphics window)
#
# Notes:
#   - Venn diagrams use pairwise comparison (2-circle format)
#   - Visualisations show both cell-type-specific and shared peaks
#   - Colour scheme: blue (imLEC), salmon (HDLEC)
# ==============================================================================

# ==============================================================================
# SECTION 1: LOAD REQUIRED PACKAGES
# ==============================================================================
install.packages("VennDiagram")
library(VennDiagram)
library(grid)
library(GenomicRanges)

# ==============================================================================
# SECTION 2: VENN DIAGRAM 1 - GENOME-WIDE ATAC-SEQ PEAKS
# ==============================================================================
# Compare all ATAC-seq peaks detected in imLEC versus HDLEC

# Load peak files for both cell types
imlec <- read.table("/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/imlec.replicated_broadPeak.sorted.bed", 
                    header=FALSE)
hdlec <- read.table("/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/LEC_CREs.sorted.bed", 
                    header=FALSE)
shared <- read.table('/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/shared_peaks.bed', 
                     header=FALSE)

# Count peaks in each dataset
imlec_n <- nrow(imlec)
hdlec_n <- nrow(hdlec)
shared_n <- nrow(shared)

# Create Venn diagram with pairwise comparison
# area1 = imLEC peaks, area2 = HDLEC peaks, cross.area = shared peaks
VENN_1 <- draw.pairwise.venn(
  area1 = imlec_n,
  area2 = hdlec_n,
  cross.area = shared_n,
  category = c("imLEC", "HDLEC"),
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 3,
  cat.cex = 2.5,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)

# ==============================================================================
# SECTION 3: VENN DIAGRAM 2 - ERG-ASSOCIATED PEAKS
# ==============================================================================
# Compare ATAC-seq peaks specifically within ERG genomic regions

# Load ERG-specific peak files
imlec_ERG <- read.table("/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/imlec_ERG.bed", 
                        header=FALSE)
hdlec_ERG <- read.table('/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/LEC_CREs_ERG.bed', 
                        header=FALSE)
shared_ERG <- read.table('/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/shared_ERG.bed', 
                         header=FALSE)

# Count ERG-associated peaks
imlec_ERG_n <- nrow(imlec_ERG)
hdlec_ERG_n <- nrow(hdlec_ERG)
shared_ERG_n <- nrow(shared_ERG)

# Create Venn diagram for ERG peaks
VENN_2 <- draw.pairwise.venn(
  area1 = imlec_ERG_n,
  area2 = hdlec_ERG_n,
  cross.area = shared_ERG_n,
  category = c("imLEC", "HDLEC"),
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 3,
  cat.cex = 2.5,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)

# ==============================================================================
# SECTION 4: VENN DIAGRAM 3 - SHARED ERG CIS-REGULATORY ELEMENTS
# ==============================================================================
# Compare peaks that are both ERG-associated AND shared between cell types

# Load shared ERG cis-regulatory element (CRE) files
imlec_ERG_cre <- read.table('/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/shared_ERG_cres_imlec.bed', 
                            header=FALSE)
hdlec_ERG_cre <- read.table("/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/shared_ERG_cres_lec.bed", 
                            header=FALSE)
shared_ERG_cre <- read.table('/Users/tamaravujic/Desktop/BHF_1+3/Laboratory /CRISPR/CRISPRi/imLEC ATAC/ERG_CRES_COMMON_IMLEC_HDLEC.bed', 
                             header=FALSE)

# Count shared ERG CRE peaks
imlec_ERG_cre_n <- nrow(imlec_ERG_cre)
hdlec_ERG_cre_n <- nrow(hdlec_ERG_cre)
shared_ERG_cre_n <- nrow(shared_ERG_cre)

# Create Venn diagram for shared ERG CREs
VENN_3 <- draw.pairwise.venn(
  area1 = imlec_ERG_cre_n,
  area2 = hdlec_ERG_cre_n,
  cross.area = shared_ERG_cre_n,
  category = c("imLEC", "HDLEC"),
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 3,
  cat.cex = 2.5,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)
