# ==============================================================================
# Script: variant_threshold_filtering.R
#
# Purpose:
#   Filters a table of genomic variants based on a set of tool-specific scoring
#   thresholds, to flag which variants show evidence of pathogenicity/regulatory
#   impact from at least one prediction tool.
#
# Workflow:
#   1. Read variant scores from multiple annotation tools
#   2. Define tool-specific thresholds for pathogenic/regulatory evidence
#   3. Evaluate each variant against all thresholds
#   4. Flag variants passing one or more thresholds
#   5. Output filtered variant table and summary statistics
#
# Input:
#   priority_scores_table.csv
#     - A CSV file containing variant-level scores from multiple annotation
#       tools (e.g. CADD, DeepSEA/Composite Priority scores, RegVar per-tissue
#       scores, and MotifBreakR transcription factor predictions).
#
# Output:
#   pathogenic_variants.csv
#     - Subset of input variants that pass at least one scoring threshold,
#       with two new columns added:
#         - thresholds_passed: count of thresholds passed per variant
#         - passed_thresholds: names of the specific thresholds passed
#   Console output: summary table showing how many variants passed 0, 1, 2, ... 
#   thresholds (distribution of filtering results)
#
# Scoring Tools Reference:
#   - CADD: Combined Annotation Dependent Depletion (variant deleteriousness)
#   - Fathmm Non_Coding_Score: Regulatory impact prediction
#   - DeepSEA Composite Priority: Machine learning variant impact scores
#   - RegVar: Tissue-specific regulatory variant prediction (6 tissues tested)
#   - MotifBreakR: Transcription factor binding disruption predictions
#
# Notes:
#   - Thresholds were pre-defined based on tool-specific literature/cutoffs
#   - RegVar scores are tissue-specific; a variant is considered "RegVar-positive"
#     if it passes the threshold in ANY one of the six tissues tested.
# ==============================================================================

# --- Setup ---
# Update this path to match your local directory structure
setwd('/Users/tamaravujic/Desktop/Variant_analysis')

# --- Read input file ---
# 'priority_scores_table.csv' contains all scores generated using various tools
# for variant analysis (CADD, Fathmm, DeepSEA, RegVar, MotifBreakR).
df <- read.csv("priority_scores_table.csv", stringsAsFactors = FALSE)

# --- Define numeric thresholds ---
# Each named value is the cutoff above which a variant is considered to "pass"
# for that particular scoring tool/category.
thresholds <- c(
  CADD_PHRED = 10,                        # CADD Phred-scaled score cutoff
  Non_Coding_Score = 0.7,                 # Fathmm regulatory impact score
  Composite_Priority_INDELS = 0.618,      # DeepSEA composite score, indels
  Composite_Priority_snv = 0.333,         # DeepSEA composite score, SNVs
  
  # RegVar thresholds (per tissue) - tissue-specific regulatory variant scores
  RegVar_Score_sigmoid_colon = 0.41,
  RegVar_Score_blood = 0.33663,
  RegVar_Score_liver = 0.4,
  RegVar_Score_lung = 0.47,
  RegVar_Score_spleen = 0.36,
  RegVar_Score_small_intestine = 0.24
)

# --- Ensure columns are numeric ---
numeric_cols <- names(thresholds)
df[numeric_cols] <- lapply(df[numeric_cols], function(x) as.numeric(as.character(x)))

# --- RegVar: TRUE if threshold passed for at least one tissue ---
regvar_cols <- grep("^RegVar_Score_", names(thresholds), value = TRUE)
regvar_pass_any <- Reduce(
  `|`,
  lapply(regvar_cols, function(col) df[[col]] > thresholds[[col]])
)

# Treat missing/NA comparisons as "did not pass" rather than NA
regvar_pass_any[is.na(regvar_pass_any)] <- FALSE

# --- Build conditions (TRUE if threshold passed) ---
# One  column per scoring category, used to evaluate how many thresholds
# each variant passes.
conditions <- data.frame(
  CADD_PHRED       = df$CADD_PHRED > thresholds["CADD_PHRED"],
  Non_Coding_Score = df$Non_Coding_Score > thresholds["Non_Coding_Score"],
  
  # Single combined RegVar condition (passes in >=1 tissue, computed above)
  RegVar_Score_any_tissue = regvar_pass_any,
  
  # DeepSEA composite priority scores, split by variant type
  Composite_Priority_INDELS = df$Composite_Priority_INDELS > thresholds["Composite_Priority_INDELS"],
  Composite_Priority_snv    = df$Composite_Priority_snv > thresholds["Composite_Priority_snv"],
  
  # MotifBreakR: TRUE if the variant is predicted to disrupt binding of any
  # TF from this list of lymphatic relevant TF families
  TF_family_match = grepl(
    "SOX|FOX|ETS|ERG|NFAT|GATA|FLI|KLF|HIF|NR2F|ELF1|ELF4|ELK2|ETV1",
    df$TFs_predicted, ignore.case = TRUE
  )
)

# --- Clean up ---
# Convert any remaining NAs to FALSE (i.e. "not passed") and ensure all
# columns are strictly logical (TRUE/FALSE) type.
conditions[is.na(conditions)] <- FALSE
conditions[] <- lapply(conditions, as.logical)

# For each variant: count how many of the 6 conditions were TRUE, and
# also record *which* thresholds were passed as a comma-separated string.
df$thresholds_passed <- rowSums(conditions)
df$passed_thresholds <- sapply(
  seq_len(nrow(conditions)),
  function(i) {
    passed <- names(conditions)[as.logical(unlist(conditions[i, , drop = TRUE]))]
    if (length(passed) == 0) "" else paste(passed, collapse = ", ")
  }
)

# --- Select variants passing at least one threshold ---
selected_variants <- df[df$thresholds_passed > 0, ]

# --- Save results ---
write.csv(selected_variants, "pathogenic_variants.csv", row.names = FALSE)

# --- Summary table ---
# Prints a frequency table: how many variants passed 0, 1, 2, ... thresholds
print(table(df$thresholds_passed))