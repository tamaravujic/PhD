# ==============================================================================
# Script: INDEL_motifbreakR_pipeline.R
#
# Purpose:
#   Complete pipeline for analyzing INDEL variants using motifbreakR.
#   Converts variant list to  anchored BED format, then predicts
#   transcription factor binding disruption across 4 motif databases.
#   
#   Indels must be "anchored" to a reference base to standardise their representation.
#   Without anchoring, the same variant can be written multiple ways (left-aligned vs
#   right-aligned), causing inconsistencies across tools. Anchoring ensures:
#   - For insertions: include the base immediately before the insertion
#   - For deletions: include the base immediately before the deleted sequence
#   This VCF-compliant format is required by motifbreakR for accurate predictions.
#
# Workflow:
#   1. Read variant list (chr_pos_ref/alt format)
#   2. Anchor variants and convert to VCF-compliant BED format
#   3. Load motifbreakR and motif databases (HOCOMOCO, HOMER, ENCODE, FactorBook)
#   4. Run motifbreakR analysis across all 4 databases
#   5. Calculate p-values for motif disruption
#   6. Export results to CSV files
#   7. Generate optional visualization plots
#
# Input:
#   variants.txt - one variant per line, format: chr_pos_ref/alt
#     Examples:
#       chr21_39313478_G/A
#       chr21_39349260-39349268_GGGCCCCCC/GCG (deletion)
#
# Output:
#   motifbreakR_anchored.bed - intermediate anchored BED file
#   HG38__indel_motifbreakR_HOCOMOCO_with_pvalues.csv
#   HG38__indel_motifbreakR_HOMER_with_pvalues.csv
#   HG38__indel_motifbreakR_ENCODE_with_pvalues.csv
#   HG38__indel_motifbreakR_FactorBook_with_pvalues.csv
#   sessionInfo.txt - R version and package info
#
# Notes:
#   - modify snps.indel[1:19] to run desired number of variants 
#   - P-values calculated with granularity 1e-6
#   - Motif disruption threshold: 1e-4
# ==============================================================================

setwd('')

# ==============================================================================
# SECTION 1: INSTALL & LOAD PACKAGES (run once if needed)
# ==============================================================================
# Uncomment and run these lines if packages are not yet installed:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(version = "3.19")
# BiocManager::install(c("motifbreakR", "BSgenome.Hsapiens.UCSC.hg38", 
#                        "SNPlocs.Hsapiens.dbSNP155.GRCh38", "MotifDb", "Biostrings"))
# install.packages("Cairo")

# ==============================================================================
# SECTION 2: VARIANT ANCHORING AND BED FILE CREATION
# ==============================================================================
# Convert variant list to motifbreakR-compatible BED file with anchoring

# Input/output files
input_variants <- "variants.txt"
anchored_bed <- "motifbreakR_anchored.bed"

# Load genome reference
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
genome <- BSgenome.Hsapiens.UCSC.hg38

# Read variant list
variants <- readLines(input_variants)
out_lines <- list()

# Process each variant and anchor to reference genome
for(ln in variants){
  ln <- trimws(ln)
  if(ln == "" || substring(ln, 1, 1) == "#") next
  
  # Parse variant string: chr_pos_ref/alt
  parts <- strsplit(ln, "_", fixed = TRUE)[[1]]
  if(length(parts) < 3) stop(paste("Cannot parse:", ln))
  
  chrom <- parts[1]
  pos <- as.integer(parts[2])
  allele_str <- paste(parts[3:length(parts)], collapse = "_")
  alleles <- strsplit(allele_str, "/", fixed = TRUE)[[1]]
  
  ref <- alleles[1]
  alts <- alleles[-1]
  
  # Process each alternate allele
  for(alt in alts){
    if(is.na(alt) || alt == ref) next
    
    # Anchor logic: add reference base to make VCF-compliant
    if(ref == "-" && alt != "-"){
      # Insertion: anchor to position, add ref base
      anchor_pos <- pos
      anchor_base <- as.character(getSeq(genome, names = chrom, start = anchor_pos, end = anchor_pos))
      REF <- anchor_base
      ALT <- paste0(anchor_base, alt)
      bed_start <- anchor_pos - 1
      bed_end <- anchor_pos
      
    } else if(ref != "-" && alt == "-"){
      # Deletion: anchor to position-1, add ref base
      anchor_pos <- pos - 1
      if(anchor_pos < 1) stop(paste("Deletion anchored before 1:", ln))
      anchor_base <- as.character(getSeq(genome, names = chrom, start = anchor_pos, end = anchor_pos))
      REF <- paste0(anchor_base, ref)
      ALT <- anchor_base
      bed_start <- anchor_pos - 1
      bed_end <- pos + nchar(ref) - 1
      
    } else {
      # Substitution/replacement
      anchor_pos <- pos - 1
      if(anchor_pos < 1) stop(paste("Variant anchored before 1:", ln))
      anchor_base <- as.character(getSeq(genome, names = chrom, start = anchor_pos, end = anchor_pos))
      REF <- paste0(anchor_base, ref)
      ALT <- paste0(anchor_base, alt)
      bed_start <- anchor_pos - 1
      bed_end <- pos + nchar(ref) - 1
    }
    
    # Create BED line
    name <- paste0(chrom, ":", anchor_pos, ":", REF, ":", ALT)
    out_lines <- append(out_lines, paste(chrom, bed_start, bed_end, name, 0, "+", sep = "\t"))
  }
}

# Write anchored BED file
writeLines(unlist(out_lines), anchored_bed)
cat("Anchored", length(out_lines), "variants to", anchored_bed, "\n")

# ==============================================================================
# SECTION 3: LOAD MOTIFBREAKR AND DATABASES
# ==============================================================================
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(MotifDb)

cat("motifbreakR version:", as.character(packageVersion("motifbreakR")), "\n")

# Load motif databases
data(motifbreakR_motif)
data(hocomoco)
data(homer)
data(encodemotif)
data(factorbook)

# Read variants from anchored BED file
snps.indel <- variants.from.file(file = anchored_bed, 
                                 search.genome = BSgenome.Hsapiens.UCSC.hg38, 
                                 format = "bed")
cat("Loaded", length(snps.indel), "variants\n")

# ==============================================================================
# SECTION 4: RUN MOTIFBREAKR ACROSS ALL DATABASES
# ==============================================================================
# Note: Currently runs on first 19 variants. Modify snps.indel[1:19] to run all.

cat("Running motifbreakR analysis...\n")

# HOCOMOCO database
results.hocomoco <- motifbreakR(snpList = snps.indel[1:19], hocomoco, filterp = TRUE,
                                threshold = 1e-4, method = "ic",
                                bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                                BPPARAM = BiocParallel::SerialParam())

# HOMER database
results.homer <- motifbreakR(snpList = snps.indel[1:19], homer, filterp = TRUE,
                             threshold = 1e-4, method = "ic",
                             bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                             BPPARAM = BiocParallel::SerialParam())

# ENCODE database
results.encodemotif <- motifbreakR(snpList = snps.indel[1:19], encodemotif, filterp = TRUE,
                                   threshold = 1e-4, method = "ic",
                                   bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                                   BPPARAM = BiocParallel::SerialParam())

# FactorBook database
results.factorbook <- motifbreakR(snpList = snps.indel[1:19], factorbook, filterp = TRUE,
                                  threshold = 1e-4, method = "ic",
                                  bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                                  BPPARAM = BiocParallel::SerialParam())

# ==============================================================================
# SECTION 5: CALCULATE P-VALUES FOR MOTIF DISRUPTION
# ==============================================================================
p.hocomoco <- calculatePvalue(results.hocomoco, granularity = 1e-6)
p.homer <- calculatePvalue(results.homer, granularity = 1e-6)
p.encodemotif <- calculatePvalue(results.encodemotif, granularity = 1e-6)
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
  HOCOMOCO = p.hocomoco,
  HOMER = p.homer,
  ENCODE = p.encodemotif,
  FactorBook = p.factorbook
)

for (name in names(results_list)) {
  export_motifbreakR_results(
    results_list[[name]],
    paste0("HG38__indel_motifbreakR_", name, "_with_pvalues.csv")
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
