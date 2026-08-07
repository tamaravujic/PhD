# ==============================================================================
# Script: SNV_motifbreakR_pipeline.R
#
# Purpose:
#   Complete pipeline for analyzing single nucleotide variants (SNVs) using
#   motifbreakR. Predicts transcription factor binding disruption across 4
#   motif databases for SNV datasets.
#
# Workflow:
#   1. Load motifbreakR and motif databases (HOCOMOCO, HOMER, ENCODE, FactorBook)
#   2. Read SNVs from  BED file
#   3. Run motifbreakR analysis across all 4 databases
#   4. Calculate p-values for motif disruption
#   5. Export results to CSV files
#   6. Optional: visualise specific variants
#
# Input:
#  patients_SNVs.bed - BED file with SNVs (already in VCF-compliant format)
#     Format: chr, start, end, variant_ID, score, strand
#     Example: chr21   39313477   39313478   chr21:39313478:G:A   0   +
#
# Output:
#   HG38_SNV_motifbreakR_HOCOMOCO_with_pvalues.csv
#   HG38_SNV_motifbreakR_HOMER_with_pvalues.csv
#   HG38_SNV_motifbreakR_ENCODE_with_pvalues.csv
#   HG38_SNV_motifbreakR_FactorBook_with_pvalues.csv
#   sessionInfo.txt - R version and package info
#
#
# Notes:
#   - modify snps.mb[1:12] to run the number of desired variants.
#   - P-values calculated with granularity 1e-6
#   - Motif disruption threshold: 1e-4
# ==============================================================================

# ==============================================================================
# SECTION 1: INSTALL & LOAD PACKAGES (run once if needed)
# ==============================================================================
# Uncomment and run these lines if packages are not yet installed:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(version = "3.19")
# BiocManager::install(c("motifbreakR", "BSgenome.Hsapiens.UCSC.hg38",
#                        "SNPlocs.Hsapiens.dbSNP155.GRCh38", "MotifDb"))
# install.packages("Cairo")

# ==============================================================================
# SECTION 2: LOAD MOTIFBREAKR AND DATABASES
# ==============================================================================
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(MotifDb)

available.SNPs()
cat("motifbreakR version:", as.character(packageVersion("motifbreakR")), "\n")

# ==============================================================================
# SECTION 3: READ SNV BED FILE
# ==============================================================================
setwd('/Users/tamaravujic/Desktop/BHF_1+3/Computational/ALL_NC_VARIANTS_FROM_GEL')
snps.bed.file <- "hg19_patients_SNVs.bed"

# Read SNVs from BED file (no dbSNP checking)
snps.mb <- snps.from.file(file = snps.bed.file,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38,
                          format = "bed", check.unnamed.for.rsid = FALSE)
snps.mb

# Load motif databases
data(motifbreakR_motif)
data(hocomoco)
data(homer)
data(encodemotif)
data(factorbook)

# ==============================================================================
# SECTION 4: RUN MOTIFBREAKR ACROSS ALL DATABASES
# ==============================================================================
# Note: Currently runs on first 12 variants. Modify snps.mb[1:12] to run all.

cat("Running motifbreakR analysis...\n")

# HOCOMOCO database
results.hocomoco <- motifbreakR(snpList = snps.mb[1:12], hocomoco, filterp = TRUE,
                                threshold = 1e-4, method = "ic",
                                bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                                BPPARAM = BiocParallel::SerialParam())

# HOMER database
results.homer <- motifbreakR(snpList = snps.mb[1:12], homer, filterp = TRUE,
                             threshold = 1e-4, method = "ic",
                             bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                             BPPARAM = BiocParallel::SerialParam())

# ENCODE database
results.encodemotif <- motifbreakR(snpList = snps.mb[1:12], encodemotif, filterp = TRUE,
                                   threshold = 1e-4, method = "ic",
                                   bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                                   BPPARAM = BiocParallel::SerialParam())

# FactorBook database
results.factorbook <- motifbreakR(snpList = snps.mb[1:12], factorbook, filterp = TRUE,
                                  threshold = 1e-4, method = "ic",
                                  bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                                  BPPARAM = BiocParallel::SerialParam())

# ==============================================================================
# SECTION 5: CALCULATE P-VALUES FOR MOTIF DISRUPTION
# ==============================================================================
p.hocomoco <- calculatePvalue(results.hocomoco, granularity = 1e-6)
p.encodemotif <- calculatePvalue(results.encodemotif, granularity = 1e-6)
p.homer <- calculatePvalue(results.homer, granularity = 1e-6)
p.factorbook <- calculatePvalue(results.factorbook, granularity = 1e-6)

cat("P-values calculated for all 4 databases\n")

# ==============================================================================
# SECTION 6: EXPORT RESULTS TO CSV
# ==============================================================================
library(dplyr)

# Function to export motifbreakR GRanges objects to CSV
export_motifbreakR_results <- function(gr_object, filename) {
  df <- as.data.frame(gr_object, row.names = NULL, optional = TRUE)
  df <- df %>%
    mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = ","))))
  rownames(df) <- NULL
  write.csv(df, file = filename, row.names = FALSE)
  message("Exported: ", filename)
}

# Save results from all 4 databases
results_list <- list(
  hocomoco = p.hocomoco,
  homer = p.homer,
  encodemotif = p.encodemotif,
  factorbook = p.factorbook
)

for (name in names(results_list)) {
  export_motifbreakR_results(
    results_list[[name]],
    paste0("HG38_SNV_motifbreakR_", name, "_with_pvalues.csv")
  )
}

# ==============================================================================
# SECTION 7: OPTIONAL - QUERY AND VISUALIZE SPECIFIC VARIANTS
# ==============================================================================

hocomoco <- results.hocomoco[results.hocomoco$SNP_id == "chr21:39313478:G:A"]
hocomoco

encodemotif <- results.encodemotif[results.encodemotif$SNP_id == "chr21:39313478:G:A"]
encodemotif

homer <- results.homer[results.homer$SNP_id == "chr21:39313478:G:A"]
homer

factorbook <- results.factorbook[results.factorbook$SNP_id == "chr21:39313478:G:A"]
factorbook

# Plot variant
install.packages("Cairo")
plotMB(results = results.hocomoco, rsid = "chr21:39034323:C:G", effect = "strong", altAllele = "G")

plotMB(results = results.encodemotif, rsid = "chr21:39034323:C:G", effect = "strong", altAllele = "G")

plotMB(results = results.factorbook, rsid = "chr21:39034323:C:G", effect = "strong", altAllele = "G")

plotMB(results = results.homer, rsid = "chr21:39034323:C:G", effect = "strong", altAllele = "G")

# ==============================================================================
# SECTION 8: SAVE SESSION INFO
# ==============================================================================
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
