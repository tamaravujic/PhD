# ==============================================================================
# Script: motifbreakR_output_sort.R
#
# Purpose:
#   Post-processing script for motifbreakR results. Summarises TF predictions
#   per variant and merges with pathogenicity scores into a single annotated file.
#
# Workflow:
#   1. Summarise motifbreakR CSV: collapse multiple TF hits per variant into
#      comma-separated lists (one row per variant)
#   2. Merge TF summary with pathogenicity scores table using variant IDs
#   3. Output single file with both TF predictions and pathogenicity scores
#
# Input:
#   variant_motifbreakr.csv - motifbreakR output (one row per TF hit)
#   priority_scores_table.csv - variant list with pathogenicity scores (CADD, DeepSEA, etc.)
#
# Output:
#   tf_summary_per_variant.csv - summarised TF predictions (one row per variant)
#   merged_variants_with_TFs.csv - final merged file with all annotations
#
# Prerequisites:
#   - Run SNV_motifbreakR_pipeline.R and/or INDEL_motifbreakR_pipeline.R first
#   - motifbreakR CSV output files available
#   - priority_scores_table.csv with variant pathogenicity data
#
# Notes:
#   - Standardizes variant ID formats (chr21:39313478:G:A → chr21_39313478_G_A)
#   - Uses left_join to preserve all variants in scores table
#   - Duplicates lines removed to prevent re-running steps
# ==============================================================================

# ==============================================================================
# PART 1: SUMMARISE MOTIFBREAKR OUTPUT
# ==============================================================================
# Collapse multiple TF hits per variant into one row per variant

library(dplyr)
library(stringr)

setwd('/Users/tamaravujic/Desktop/Variant_analysis')

# Read motifbreakR output (one row per TF hit per variant)
df <- read.csv("motifbreakr.csv")

# Summarize per variant: concatenate TF gene symbols and binding changes
summary_df <- df %>%
  group_by(SNP_id) %>%
  summarise(
    TFs_predicted = paste(geneSymbol, collapse = ", "),
    Binding_change = paste(Binding_change, collapse = ", ")
  ) %>%
  ungroup()

# Save summarized results (one row per variant)
write.csv(summary_df, "tf_summary_per_variant.csv", row.names = FALSE)
print(head(summary_df))

# ==============================================================================
# PART 2: MERGE TF SUMMARY WITH PATHOGENICITY SCORES
# ==============================================================================
# Combine TF predictions with variant pathogenicity annotations

# Load both tables
scores <- read.csv("priority_scores_table.csv")        # pathogenicity data
tfs     <- read.csv("tf_summary_per_variant.csv")      # TF motif summary

# Standardize variant ID formats to enable matching
# Convert SNP_id format (chr21:39313478:G:A) to variant_key (chr21_39313478_G_A)
tfs <- tfs %>%
  mutate(
    variant_key = SNP_id %>%
      gsub(":", "_", .) %>%
      gsub("/", "_", .)
  )

# Standardize scores table variant ID
scores <- scores %>%
  mutate(
    variant_key = gsub("/", "_", Uploaded_variation)
  )

# Merge by standardized variant key (left join keeps all score variants)
merged_df <- left_join(scores, tfs, by = "variant_key")

# Remove helper column and save final merged file
merged_df <- merged_df %>%
  select(-variant_key)

write.csv(merged_df, "merged_variants_with_TFs.csv", row.names = FALSE)

# Inspect results
head(merged_df)
